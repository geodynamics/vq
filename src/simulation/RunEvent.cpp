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
 */
void RunEvent::markBlocks2Fail(VCSimulation *sim, const FaultID &trigger_fault, VCEventSweep &current_sweep) {
    int         lid;
    BlockID     gid;
    bool        add;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);

        // Blocks can only fail once per event, after that they slide freely
        if (sim->getBlock(gid).getFailed()) continue;

        // Add this block if it has a static CFF failure
        add = sim->getBlock(gid).cffFailure();

        // Allow dynamic failure if the block is "loose" (next to a previously failed block)
        if (looseBlocks.count(gid) > 0) add |= sim->getBlock(gid).dynamicFailure(trigger_fault);

        if (add) {
            sim->getBlock(gid).setFailed(true);
            local_failed_blocks.insert(gid);
        }
    }
}

/*!
 Process the list of blocks that failed on this node using the original friction law.
 */
void RunEvent::processBlocksOrigFail(VCSimulation *sim, VCEventSweep &current_sweep) {
    BlockIDSet::iterator    fit;
    double                  slip, stress_drop;

    // For each block that fails in this sweep, calculate how much it slips
    for (fit=local_failed_blocks.begin(); fit!=local_failed_blocks.end(); ++fit) {
        if (sim->isLocalBlockID(*fit)) {
            Block &b = sim->getBlock(*fit);

            // calculate the drop in stress from the failure
            stress_drop = b.getCFF0() - b.getCFF();

            if (!stress_drop) stress_drop = b.getStressDrop() - b.getCFF();

            //std::cerr << stress_drop << std::endl;

            // Slip is in m
            slip = (stress_drop/b.getSelfStresses());

            if (slip < 0) slip = 0;

            if (slip > -b.state.slipDeficit) slip = -b.state.slipDeficit;

            // Record how much the block slipped in this sweep and initial stresses
            current_sweep.setSlipAndArea(b.getBlockID(), slip, b.get_area(), b.lame_mu());
            current_sweep.setInitStresses(b.getBlockID(), b.getShearStress(), b.getNormalStress());

            b.state.slipDeficit += slip;
        }
    }
}

void solve_it(int n, double *x, double *A, double *b) {
    int     i, j, k;
    double  v, f, sum;

    for (i=0; i<n; ++i) {
        v = A[i+n*i];

        for (j=i+1; j<n; ++j) {
            f = A[i+n*j]/v;

            for (k=0; k<n; ++k) {
                A[k+n*j] -= f*A[k+n*i];
            }

            b[j] -= f*b[i];
        }
    }

    for (i=n-1; i>=0; --i) {
        sum = b[i];

        for (j=i+1; j<n; ++j) {
            sum -= A[j+n*i]*x[j];
        }

        x[i] = sum/A[i+n*i];
    }
}

void RunEvent::processBlocksSecondaryFailures(VCSimulation *sim, VCEventSweep &current_sweep) {
    int             lid;
    BlockID         gid;
    unsigned int    i, n;
    BlockIDSet                      local_id_list;
    BlockIDProcMapping              global_id_list;
    BlockIDSet::const_iterator      it;
    BlockIDProcMapping::const_iterator  jt;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        Block &b = sim->getBlock(gid);

        // If the block has already failed (but not in this sweep) then adjust the slip
        if (b.getFailed() && global_failed_blocks.count(gid) == 0) {
            local_id_list.insert(gid);
        }
    }

    // Figure out how many failures there were over all processors
    sim->distributeBlocks(local_id_list, global_id_list);
    
    int num_local_failed = local_id_list.size();
    int num_global_failed = global_id_list.size();
    
    double *A = new double[num_local_failed*num_global_failed];
    double *b = new double[num_local_failed];
    double *x = new double[num_local_failed];

    for (i=0,it=local_id_list.begin(); it!=local_id_list.end();++i,++it) {
        Block &blk = sim->getBlock(*it);

        for (n=0,jt=global_id_list.begin(); jt!=global_id_list.end(); ++n,++jt) {
            A[i*num_global_failed+n] = sim->getGreenShear(*it, jt->first);

            if (sim->doNormalStress()) {
                A[i*num_global_failed+n] -= blk.friction()*sim->getGreenNormal(*it, jt->first);
            }
        }

        b[i] = blk.getCFF()+blk.friction()*blk.getRhogd();
    }

    if (sim->isRootNode()) {
        double *fullA = new double[num_global_failed*num_global_failed];
        double *fullb = new double[num_global_failed];
        double *fullx = new double[num_global_failed];
        
        // Fill in the A matrix and b vector from the various processes
        for (i=0,n=0,jt=global_id_list.begin();jt!=global_id_list.end();++jt,++i) {
            if (jt->second != sim->getNodeRank()) {
                MPI_Recv(&(fullA[i*num_global_failed]), num_global_failed, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(fullb[i]), 1, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                memcpy(&(fullA[i*num_global_failed]), &(A[n*num_global_failed]), sizeof(double)*num_global_failed);
                memcpy(&(fullb[i]), &(b[n]), sizeof(double));
                n++;
            }
        }
        
        // Solve the global system on the root node
        solve_it(num_global_failed, fullx, fullA, fullb);
        
        // Send back the resulting values from x to each process
        for (i=0,n=0,jt=global_id_list.begin();jt!=global_id_list.end();++jt,++i) {
            if (jt->second != sim->getNodeRank()) {
                MPI_Send(&(fullx[i]), 1, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD);
            } else {
                memcpy(&(x[n]), &(fullx[i]), sizeof(double));
                n++;
            }
        }
        
        // Delete the memory arrays created
        delete fullx;
        delete fullb;
        delete fullA;
    } else {
        for (i=0;i<num_local_failed;++i) {
            MPI_Send(&(A[i*num_global_failed]), num_global_failed, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(b[i]), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        
        for (i=0;i<num_local_failed;++i) {
            MPI_Recv(&(x[i]), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // Take the results of the calculation and determine how much each ruptured block slipped
    for (i=0,it=local_id_list.begin(); it!=local_id_list.end(); ++i,++it) {
        Block &block = sim->getBlock(*it);

        double slip = x[i] - block.getSlipDeficit();

        if (slip > 0) {
            // Record how much the block slipped in this sweep and initial stresses
            current_sweep.setSlipAndArea(block.getBlockID(), slip, block.get_area(), block.lame_mu());
            current_sweep.setInitStresses(block.getBlockID(), block.getShearStress(), block.getNormalStress());

            block.state.slipDeficit += slip;
        }
    }

    delete A;
    delete b;
    delete x;
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
    bool                    final_sweep = false;
    std::pair<BlockIDSet::const_iterator, BlockIDSet::const_iterator> nbr_start_end;
    BlockIDSet::const_iterator  nit;

    // Get the event trigger block and fault
    triggerID = sim->getCurrentEvent().getEventTrigger();
    trigger_fault = sim->getBlock(triggerID).getFaultID();

    // Clear the list of "loose" (can dynamically fail) blocks
    looseBlocks.clear();

    // Clear the list of failed blocks, and add the trigger block
    local_failed_blocks.clear();

    if (sim->getCurrentEvent().getEventTriggerOnThisNode()) {
        local_failed_blocks.insert(triggerID);
        looseBlocks.insert(triggerID);
        sim->getBlock(triggerID).setFailed(true);
    }

    // Save stress information at the beginning of the event
    // This is used to determine dynamic block failure
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) sim->getBlock(sim->getGlobalBID(lid)).saveStresses();

    more_blocks_to_fail = sim->blocksToFail(!local_failed_blocks.empty());

    // While there are still failed blocks to handle
    while (more_blocks_to_fail || final_sweep) {
        // Clear the current sweep event
        current_sweep.clear();

        // Share the failed blocks with other processors to correctly handle
        // faults that are split among different processors
        sim->distributeBlocks(local_failed_blocks, global_failed_blocks);

        // Process the blocks that failed
        processBlocksOrigFail(sim, current_sweep);

        // Recalculate CFF for all blocks where slipped blocks don't contribute
        for (it=sim->begin(); it!=sim->end(); ++it) {
            sim->setShearStress(it->getBlockID(), 0.0);
            sim->setNormalStress(it->getBlockID(), it->getRhogd());
            sim->setUpdateField(it->getBlockID(), (it->getFailed() ? 0 : it->state.slipDeficit));
        }

        // Distribute the update field values to other processors
        sim->distributeUpdateField();

        // Set dynamic triggering on for any blocks neighboring blocks that slipped in the last sweep
        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
            if (current_sweep.getBlockSlip(gid) > 0) {
                nbr_start_end = sim->getNeighbors(gid);

                for (nit=nbr_start_end.first; nit!=nbr_start_end.second; ++nit) {
                    looseBlocks.insert(*nit);
                }
            }
        }

        // Calculate the CFFs based on the stuck blocks
        sim->matrixVectorMultiplyAccum(sim->getShearStressPtr(),
                                       sim->greenShear(),
                                       sim->getUpdateFieldPtr(),
                                       true);

        if (sim->doNormalStress()) {
            sim->matrixVectorMultiplyAccum(sim->getNormalStressPtr(),
                                           sim->greenNormal(),
                                           sim->getUpdateFieldPtr(),
                                           true);
        }

        sim->computeCFFs();

        // For each block that has failed already, calculate the slip from the movement of blocks that just failed
        processBlocksSecondaryFailures(sim, current_sweep);

        // Set the update field to the slip of all blocks
        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
            Block &block=sim->getBlock(gid);
            sim->setShearStress(gid, 0.0);
            sim->setNormalStress(gid, block.getRhogd());
            sim->setUpdateField(gid, (block.getFailed() ? 0 : block.state.slipDeficit));
        }

        // Calculate the new shear stresses and CFFs given the new update field values
        sim->matrixVectorMultiplyAccum(sim->getShearStressPtr(),
                                       sim->greenShear(),
                                       sim->getUpdateFieldPtr(),
                                       true);

        if (sim->doNormalStress()) {
            sim->matrixVectorMultiplyAccum(sim->getNormalStressPtr(),
                                           sim->greenNormal(),
                                           sim->getUpdateFieldPtr(),
                                           true);
        }

        sim->computeCFFs();

        // Record the final stresses of blocks that failed during this sweep
        BlockIDProcMapping::iterator        fit;

        for (fit=global_failed_blocks.begin(); fit!=global_failed_blocks.end(); ++fit) {
            if (sim->isLocalBlockID(fit->first)) {
                Block &b = sim->getBlock(fit->first);
                current_sweep.setFinalStresses(fit->first, b.getShearStress(), b.getNormalStress());
            }
        }

        global_failed_blocks.clear(); // we are done with these blocks
        local_failed_blocks.clear(); // we are done with these blocks

        // Find any blocks that fail because of the new stresses
        markBlocks2Fail(sim, trigger_fault, current_sweep);

        sim->collectEventSweep(current_sweep);

        // Add the recorded sweep to the list
        sweep_list.push_back(current_sweep);

        if (final_sweep) {
            final_sweep = false;
        } else {
            more_blocks_to_fail = sim->blocksToFail(!local_failed_blocks.empty());

            if (!more_blocks_to_fail) final_sweep = true;
        }
    }

    // Set the completed list as the sweep list for the entire event
    sim->getCurrentEvent().addSweeps(sweep_list);

    // Record the stress in the system before and after the event
    recordEventStresses(sim);

    // Update the cumulative slip for this fault
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        sim->getBlock(gid).setFailed(false);
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
#ifdef HAVE_MPI
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
