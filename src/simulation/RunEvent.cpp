// Copyright (c) 2012-2014 Eric M. Heien, Michael K. Sachs, John B. Rundle
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include "RunEvent.h"
#include "SimError.h"

/*!
 At the end of each sweep after we have recalculated block CFF, we determine
 which blocks will have a failure due to dynamic or static stress changes.
 To avoid runaway failures we limit the number of failures for a given block in
 a sweep to 20. However, it is unclear if this ever happens in practice.
 */
void RunEvent::markBlocks2Fail(VCSimulation *sim, const FaultID &trigger_fault, VCEventSweep &current_sweep) {
    int         lid;
    BlockID     gid;
    bool        add;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        //current_sweep.setFinalStresses(gid, sim->getBlock(gid).getShearStress(), sim->getBlock(gid).getNormalStress());
        // Add this block if it has a static CFF failure or dynamic failure
        add = sim->getBlock(gid).cffFailure() || sim->getBlock(gid).dynamicFailure(trigger_fault) || sim->getBlock(gid).additionalFailure();

        //sim->console() << gid <<
        //" failure s:" << sim->getBlock(gid).cffFailure() << " d:" << sim->getBlock(gid).dynamicFailure(trigger_fault) << " a:" << sim->getBlock(gid).additionalFailure() << ", ";
        // Let each block fail at most 10 times. This is somewhat arbitrary, but it prevents runaway ruptures.
        if (add && num_failures[gid] < 10) {
            num_failures[gid] += 1;
            sim->getBlock(gid).setFailed(true);
            sim->getBlock(gid).setFailedThisSweep(true);
            blocks2fail.insert(gid);
        }
    }
}

/*!
 Process the list of blocks that failed on this node using the original friction law.
 */
void RunEvent::processBlocksOrigFrictionLaw(VCSimulation *sim, VCEventSweep &current_sweep) {
    BlockIDSet::iterator        fit;
    double                      rn, slip, stress_drop, a;
    FaultBlockMapping           failed_faults;

    sim->getFaultBlockMapping(failed_faults, blocks2fail);

    // For each block that has failed
    for (fit=blocks2fail.begin(); fit!=blocks2fail.end(); ++fit) {
        if (sim->isLocalBlockID(*fit)) {
            Block &b = sim->getBlock(*fit);

            // calculate the drop in stress from the failure
            stress_drop = b.getStressDrop() - b.getCFF();
            a = sim->getEventNoise();

            // noise in the range [1-a, 1+a)
            if (a != 0) rn = (1.0-a) + 2*a*b.state.rand.nextDouble();
            else rn = 1;

            // Slip is in m
            slip = rn * (stress_drop/b.getSelfStresses());

            // Record how much the block slipped in this sweep and initial stresses
            current_sweep.setSlipAndArea(b.getBlockID(), slip, b.get_area(), b.getMu());
            current_sweep.setInitStresses(b.getBlockID(), b.getShearStress(), b.getNormalStress());
            b.saveFStresses();

            b.state.slipDeficit += slip;
        }
    }

    //sim->console()  << std::endl;
}

/*!
 Process the list of blocks that failed on this node using the stepped friction law.
 */
void RunEvent::processBlocksNewFrictionLaw(VCSimulation *sim, VCEventSweep &current_sweep) {
    BlockIDSet::iterator        fit;
    double                      rn, slip, stress_drop, a, failure_size;
    FaultBlockMapping           failed_faults;
    FaultFailureAreaMapping     failed_faults_area;
    double                      slip_scaling_threshold;

    sim->getFaultBlockMapping(failed_faults, blocks2fail);

    sim->getFaultFailureAreaMapping(failed_faults_area, blocks2fail);

    // For each block that has failed
    for (fit=blocks2fail.begin(); fit!=blocks2fail.end(); ++fit) {
        if (sim->isLocalBlockID(*fit)) {
            Block &b = sim->getBlock(*fit);

            // calculate the drop in stress from the failure
            stress_drop = b.getStressDrop() - b.getCFF();

            a = sim->getEventNoise();

            // noise in the range [1-a, 1+a)
            if (a != 0) rn = (1.0-a) + 2*a*b.state.rand.nextDouble();
            else rn = 1;

            //rn = (1.0-a) + 2*a*b.state.rand.nextDouble();
            // Slip is in m



            failure_size = double(failed_faults[b.getFaultID()].size());
            slip_scaling_threshold = double(b.getSlipScalingThreshold());

            if (failure_size <= slip_scaling_threshold) {
                slip = rn * (stress_drop/b.getSelfStresses()) * (failure_size/slip_scaling_threshold);
            } else {
                slip = rn * (stress_drop/b.getSelfStresses());
            }

            //slip =  -(sqrt(failed_faults_area[b.getFaultID()])/b.getMu()) * stress_drop;

            //sim->console() << failed_faults_area[b.getFaultID()] << std::endl;

            //sim->console() << -(sqrt(failed_faults_area[b.getFaultID()])/b.getMu()) * stress_drop << " " << slip << std::endl;
            //sim->console() << -(sqrt(failed_faults_area[b.getFaultID()])/b.getMu()) * stress_drop << " " << -(sqrt(b.get_area())/b.getMu()) * stress_drop << " " << slip << std::endl;

            // Record how much the block slipped in this sweep and initial stresses
            current_sweep.setSlipAndArea(b.getBlockID(), slip, b.get_area(), b.getMu());
            current_sweep.setInitStresses(b.getBlockID(), b.getShearStress(), b.getNormalStress());
            b.saveFStresses();

            b.state.slipDeficit += slip;
        }
    }
}

/*!
 Given an initial failed block, propagates the failure throughout the system
 by calculating changes in stress and using static and dynamic stress
 failure functions. A single step in the failure propagation is called a sweep
 and multiple sweeps comprise an entire event.
 */
SimRequest RunEvent::run(SimFramework *_sim) {
    VCSimulation            *sim = static_cast<VCSimulation *>(_sim);
    BlockList::iterator     it;
    VCEvent::iterator       eit;
    VCEventSweep            current_sweep;
    BlockID                 triggerID, gid;
    FaultID                 trigger_fault;
    int                     lid, more_blocks_to_fail;
    EventSweeps             sweep_list;
    std::pair<BlockIDSet::const_iterator, BlockIDSet::const_iterator> nbr_start_end;
    BlockIDSet::const_iterator  nit;

    // Get the event trigger block and fault
    triggerID = sim->getCurrentEvent().getEventTrigger();
    trigger_fault = sim->getBlock(triggerID).getFaultID();

    // Clear the list of failed blocks, and add the trigger block
    blocks2fail.clear();
    num_failures.clear();

    if (sim->getCurrentEvent().getEventTriggerOnThisNode()) {
        num_failures[triggerID] += 1;
        blocks2fail.insert(triggerID);
        sim->getBlock(triggerID).setFailed(true);
        sim->getBlock(triggerID).setFailedThisSweep(true);
    }

    // Save stress information at the beginning of the event
    // This is used to determine dynamic block failure
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) sim->getBlock(sim->getGlobalBID(lid)).saveStresses();

    more_blocks_to_fail = sim->blocksToFail(!blocks2fail.empty());

    // While there are still failed blocks to handle or we're in localized failure mode
    while (more_blocks_to_fail) {
        // Clear the current sweep event
        current_sweep.clear();

        // Share the failed blocks with other processors to correctly handle
        // faults that are split among different processors
        sim->distributeFailedBlocks(blocks2fail);

        // Process the blocks that failed using the specified friction law
        switch (sim->getFrictionLaw()) {
            case FRIC_LAW_ORIG:
                processBlocksOrigFrictionLaw(sim, current_sweep);
                break;

            case FRIC_LAW_STEPPED:
                processBlocksNewFrictionLaw(sim, current_sweep);
                break;

            default:
                break;
        }

        // Set the update field for stress recalculation with
        // the total amount each block has slipped in this sweep
        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
            sim->setUpdateField(gid, current_sweep.getBlockSlip(gid));
        }

        // Set the update field for stress recalculation with
        // the total amount each block has slipped in this sweep
        // Working with the slip deficit
        //for(gid=0;gid<sim->numGlobalBlocks();++gid) {
        //  sim->setFShearStress(gid, 0.0);
        //  sim->setFNormalStress(gid, sim->getBlock(gid).getRhogd());
        //sim->setUpdateField(it->getBlockID(), it->state.slipDeficit);
        //  sim->setUpdateField(gid, sim->getBlock(gid).state.slipDeficit);
        //}

        // Distribute the update field values to other processors
        sim->distributeUpdateField();

        // Set dynamic triggering on for any blocks neighboring blocks that slipped in the last sweep
        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
            if (sim->getUpdateField(gid) > 0) {
                nbr_start_end = sim->getNeighbors(gid);

                for (nit=nbr_start_end.first; nit!=nbr_start_end.second; ++nit) {
                    sim->getBlock(*nit).dynamicOn();
                }
            }
        }

        // Calculate the new shear stresses and CFFs given the new update field values
        sim->matrixVectorMultiplyAccum(sim->getFShearStressPtr(),
                                       sim->greenShear(),
                                       sim->getUpdateFieldPtr(),
                                       false);

        if (sim->doNormalStress()) {
            sim->matrixVectorMultiplyAccum(sim->getFNormalStressPtr(),
                                           sim->greenNormal(),
                                           sim->getUpdateFieldPtr(),
                                           false);
        }

        sim->computeCFFs(true);

        // Record the final stresses of blocks that failed during this sweep
        BlockIDSet::iterator        fit;

        for (fit=blocks2fail.begin(); fit!=blocks2fail.end(); ++fit) {
            if (sim->isLocalBlockID(*fit)) {
                Block &b = sim->getBlock(*fit);
                current_sweep.setFinalStresses(*fit, b.getShearStress(), b.getNormalStress());
            }
        }

        blocks2fail.clear(); // we are done with these blocks

        // Find any blocks that fail because of the new stresses
        markBlocks2Fail(sim, trigger_fault, current_sweep);

        sim->collectEventSweep(current_sweep);

        // Add the recorded sweep to the list
        sweep_list.push_back(current_sweep);

        more_blocks_to_fail = sim->blocksToFail(!blocks2fail.empty());
    }

    // Set the completed list as the sweep list for the entire event
    sim->getCurrentEvent().addSweeps(sweep_list);

    // Record the stress in the system before and after the event
    recordEventStresses(sim);

    // Update the cumulative slip for this fault and turn off dynamic triggering
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        sim->getBlock(gid).setFailed(false);
        sim->getBlock(gid).setFailedThisSweep(false);
        sim->getBlock(gid).dynamicOff();
    }

    assertThrow(sim->getCurrentEvent().getNumSweeps() > 0, "There was a trigger but no failed blocks.");

    return SIM_STOP_OK;
}

void RunEvent::recordEventStresses(VCSimulation *sim) {
    BlockIDSet involved_blocks;
    double shear_init, shear_final, normal_init, normal_final;
    double total_shear_init, total_shear_final, total_normal_init, total_normal_final;

    sim->getCurrentEvent().getInvolvedBlocks(involved_blocks);

    sim->getInitialFinalStresses(involved_blocks, shear_init, shear_final, normal_init, normal_final);

    // If we have multiple processors, sum the values and store on the root node
#ifdef MPI_C_FOUND
    MPI_Reduce(&shear_init, &total_shear_init, 1, MPI_DOUBLE, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);
    MPI_Reduce(&shear_final, &total_shear_final, 1, MPI_DOUBLE, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);
    MPI_Reduce(&normal_init, &total_normal_init, 1, MPI_DOUBLE, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);
    MPI_Reduce(&normal_final, &total_normal_final, 1, MPI_DOUBLE, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);
#else
    total_shear_init = shear_init;
    total_shear_final = shear_final;
    total_normal_init = normal_init;
    total_normal_final = normal_final;
#endif

    sim->getCurrentEvent().setEventStresses(total_shear_init, total_shear_final, total_normal_init, total_normal_final);
}
