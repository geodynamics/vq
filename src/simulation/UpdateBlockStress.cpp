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

#include "UpdateBlockStress.h"

/*!
 Initialize the stress calculation by setting initial block slip and stresses
 then calculating the stress in the whole system.
 */
void UpdateBlockStress::init(SimFramework *_sim) {
    BlockList::iterator nt;
    BlockID             gid;
    SectionID           sid;
    int                 lid, i;
    bool                err;
    double              stress_drop;
    double rho = 2700.0;      // density of rock in kg m^-3
    double g = 9.81;           // force of gravity in m s^-2
    double depth = 0.0;         //
    //double mean_slip_rate = 0.0;
    //double norm_velocity;
    double char_magnitude, char_slip, fault_length, fault_area, fault_width, nu, R;
    std::map<SectionID, double> section_lengths;
    std::map<SectionID, double> section_areas;
    std::map<SectionID, double>::iterator sit;
    std::map<SectionID, double> section_min_das;
    std::map<SectionID, double>::iterator ssit;
    quakelib::ModelStressSet    stress_set;
    quakelib::ModelStress       stress;
    
    sim = static_cast<Simulation *>(_sim);
    tmpBuffer = new double[sim->numGlobalBlocks()];
    
    // Read the stress input file for initial stress conditions
    std::string stress_file_type = sim->getStressInfileType();
    std::string stress_filename = sim->getStressInfile();
    std::string stress_index_filename = sim->getStressIndexInfile();

    if (stress_filename != "" && stress_file_type != "") {
        
        if (stress_file_type == "text") {
            if (stress_index_filename == "") {
                sim->errConsole() << "ERROR: Must specify stress index file " << std::endl;
                return;
            } else {
                err = stress_set.read_file_ascii(stress_index_filename, stress_filename);
            }
        } else if (stress_file_type == "hdf5") {
            err = stress_set.read_file_hdf5(stress_filename);
        } else {
            sim->errConsole() << "ERROR: unknown file type " << stress_file_type << std::endl;
            return;
        }
        
        // Schultz: Currently we just load the last event saved in the stress state file.
        stress = stress_set[stress_set.size()-1].stresses();
        // Also set the sim year to the year the stresses were saved
        sim->setYear(stress_set[stress_set.size()-1].getYear());
        sim->console() << "--- Setting initial stresses from file, starting new sim at year " << sim->getYear() << " ---" << std::endl;
    
        // If given an initial stress state, set those stresses and slip deficits.
        // Schultz: The slip deficit is really the only information used to start the sim, as we
        // recalculate stresses at the end of this init() based on slip deficits.
        for (i=0; i<stress.size(); ++i) {
            sim->setInitShearNormalStress(stress[i]._element_id, stress[i]._shear_stress, stress[i]._normal_stress);
            sim->setSlipDeficit(stress[i]._element_id, stress[i]._slip_deficit);
        }
    }

    // Determine section minimum distance along strike (required since das is defined along
    // faults and doesn't reset to 0 when entering a new section of the same fault.
    for (nt=sim->begin(); nt!=sim->end(); ++nt) {
        sid = nt->getSectionID();
        
        if (section_min_das.count(sid)) {
            sit = section_min_das.find(sid);
            // Replace the current max length with this element's distance along strike if it's smaller
            // Trying to find the lower bound for the fault
            sit->second = std::min(sit->second, nt->min_das());
            
        } else {
            // If it's not already in here, add this section
            section_min_das.insert(std::make_pair(sid, nt->min_das()));
        }
        
    }
    
    // Determine section lengths and add up the areas
    for (nt=sim->begin(); nt!=sim->end(); ++nt) {
        sid = nt->getSectionID();

        if (section_lengths.count(sid)) {
            sit = section_lengths.find(sid);
            ssit = section_min_das.find(sid);
            double min_das = ssit->second;
            // Replace the current max length with this element's distance along strike if it's larger
            sit->second = std::max(sit->second, nt->max_das()-min_das);

        } else {
            ssit = section_min_das.find(sid);
            double min_das = ssit->second;
            // If it's not already in here, add this section
            section_lengths.insert(std::make_pair(sid, nt->max_das()-min_das));
        }

        if (section_areas.count(sid)) {
            sit = section_areas.find(sid);
            // Add the current element's area to the section's total
            sit->second += nt->area();

        } else {
            // If it's not already in here, add this section
            section_areas.insert(std::make_pair(sid, nt->area()));
        }

    }
    
    // Set the simulation section areas now that we computed them
    for (sit=section_areas.begin(); sit!=section_areas.end(); ++sit) {
        sim->setSectionArea(sit->first, sit->second);
    }
    
    // Set the simulation section lengths now that we computed them
    for (sit=section_lengths.begin(); sit!=section_lengths.end(); ++sit) {
        sim->setSectionLength(sit->first, sit->second);
    }
    
    
    /*
    /////// Schultz: First we compute the mean slip rate to avoid NaN's
    for (nt=sim->begin(); nt!=sim->end(); ++nt) {
        mean_slip_rate += nt->slip_rate();
    }
    // Normalize by number of blocks
    mean_slip_rate /= sim->numGlobalBlocks();
    */

    // All processes need the friction values for all blocks, so we set rhogd here
    // and transfer stress drop values between nodes later
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        //
        // TODO: check if this negative sign is warranted and not double counted
        depth = fabs(sim->getBlock(gid).center()[2]);  // depth of block center in m

        sim->setRhogd(gid, rho*g*depth);       // kg m^-3 * m s^-2 * m = kg m^-1 * s^-2 = Pa
        //printf("**Degbug(%d): setting rhogd[%d/%dg], %f, %f, %f ==> %f\n", getpid(), lid, gid, rho, g, depth, sim->getRhogd(gid));

        sim->setDynamicVal(gid, sim->getDynamic());
        sim->setFailed(gid, false);

        // Set stresses to their specified initial values
        //printf("**Degbug(%d): shear[%d]\n", getpid(), gid);
        sim->setShearStress(gid, sim->getInitShearStress(gid));
        //printf("**Degbug(%d): normal[%d]\n", getpid(), gid);
        sim->setNormalStress(gid, sim->getInitNormalStress(gid));

        // Set the stress drop based on the Greens function calculations
        if (sim->computeStressDrops()) {
            // TODO: Instead of dividing by the local slip rate as your normalization,
            //       we need to have some minimum or mean normalization constant to
            //       handle zero slip rates that lead to stress_drop = NaN.

            /* Eric's stress drop method
            // Now compute weighted stress drops from slip rates and shear interactions
            stress_drop = 0;
            norm_velocity = sim->getBlock(gid).slip_rate();

            for (nt=sim->begin(); nt!=sim->end(); ++nt) {
                stress_drop += ((nt->slip_rate() + mean_slip_rate)/(norm_velocity + mean_slip_rate))*sim->getGreenShear(gid, nt->getBlockID());
            }

            stress_drop *= sim->getBlock(gid).max_slip();
            /////// Schultz: All stress drops must be negative
            if (stress_drop > 0) stress_drop = -1.0*fabs(stress_drop);

            /////// Schultz: Hack #2, multiply the stress drops by 5.0 to get closer to the prescribed VC stress drop values
            // TODO: Make this a parameter
            stress_drop *= 5.0;
            */

            ///// Schultz stress drop method
            // Get the section id, sid
            sid = sim->getBlock(gid).getSectionID();
            // Get fault area and length (we will compute the width)
            fault_area = sim->getSectionArea(sid);
            fault_length = sim->getSectionLength(sid);
            fault_width = fault_area/fault_length; // This way we get the mean width

            // Use Wells & Coppersmith scaling to find the characteristic magnitude and slip given the fault geometry
            char_magnitude = 4.0+log10(fault_area*1e-6) + sim->stressDropFactor();
            char_slip = pow(10, (3.0/2.0)*(char_magnitude+10.7))/(1e7*sim->getBlock(gid).lame_mu()*fault_area);

            // Compute stress drop from geometry and expected slip/mag
            nu = 0.5*sim->getBlock(gid).lame_lambda()/(sim->getBlock(gid).lame_mu() + sim->getBlock(gid).lame_lambda());
            R  = sqrt(fault_width*fault_width + fault_length*fault_length);

            stress_drop = -2*sim->getBlock(gid).lame_mu()*char_slip*( (1-nu)*fault_length/fault_width + fault_width/fault_length )/( (1-nu)*M_PI*R ) ;

            sim->setStressDrop(gid, stress_drop);
            sim->setMaxStressDrop(gid, stress_drop);

        } else {
            sim->setStressDrop(gid, sim->getBlock(gid).stress_drop());
            sim->setMaxStressDrop(gid, sim->getBlock(gid).stress_drop());
        }

        // Initialize element slips to equilibrium position, slip=0
        // Unless we are reading in a stress input file, then we have already set the slip deficit higher up in this method
        if (sim->getStressInfile() == "" && sim->getStressIndexInfile() == "") sim->setSlipDeficit(gid, 0);

        if (sim->isLocalBlockID(gid)) {
            sim->decompressNormalRow(gid);
            //sim->setGreenNormal(bid, bid, 0.0);  // Erase diagonal for normal Greens matrix
            sim->compressNormalRow(gid, 0.7);
        }
    }

#ifdef MPI_C_FOUND

    // Transfer stress drop values between nodes
    // also broadcast the rhogd value.
    // note that the array rhogd[] is a protected array, so we cannot access it directly (aka, using MPI_Reduce with a custom mpi funciton to
    // gather/add/whatever, but handle nan), but since this only runs once at the sim initialization, it's not a huge problem.
    for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
        // initialize these to avoid memcheck errors:
        double       stress_drop=std::numeric_limits<float>::quiet_NaN();
        double       max_stress_drop=std::numeric_limits<float>::quiet_NaN();
        double       tmp_rhogd = std::numeric_limits<float>::quiet_NaN();
        double       slip_deficit = std::numeric_limits<float>::quiet_NaN();

        if (sim->isLocalBlockID(gid)) {
            stress_drop = sim->getStressDrop(gid);
            max_stress_drop = sim->getMaxStressDrop(gid);
            //
            tmp_rhogd = sim->getRhogd(gid);
            slip_deficit = sim->getSlipDeficit(gid);
        }

        //
        MPI_Bcast(&stress_drop, 1, MPI_DOUBLE, sim->getBlockNode(gid), MPI_COMM_WORLD);
        MPI_Bcast(&max_stress_drop, 1, MPI_DOUBLE, sim->getBlockNode(gid), MPI_COMM_WORLD);
        MPI_Bcast(&slip_deficit, 1, MPI_DOUBLE, sim->getBlockNode(gid), MPI_COMM_WORLD);
        //
        MPI_Bcast(&tmp_rhogd, 1, MPI_DOUBLE, sim->getBlockNode(gid), MPI_COMM_WORLD);

        //
        if (!sim->isLocalBlockID(gid)) {
            sim->setStressDrop(gid, stress_drop);
            sim->setMaxStressDrop(gid, max_stress_drop);
            //
            sim->setRhogd(gid, tmp_rhogd);
            sim->setSlipDeficit(gid, slip_deficit);
        }
    }

#endif

    // Compute initial stress on all blocks
    stressRecompute();
    
//    Debug output
//    if (sim->isRootNode()) {
//        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
//            std::cout << gid << "  " << sim->getShearStress(gid) << "  " << sim->getNormalStress(gid) << "  " << sim->getSlipDeficit(gid) <<std::endl;
//        }
//    }

}

/*!
 Calculate the stress on all blocks and determine which block will be
 the next to fail and when it will fail. Use this to determine the
 slipDeficit on all blocks at the time of failure. Finally, use the
 new slipDeficit to recalculate the stress on all blocks at the time of failure.
 */
SimRequest UpdateBlockStress::run(SimFramework *_sim) {
    // Put a stress load on all blocks and determine which block will fail first
    int                     lid;
    BlockVal                next_static_fail, next_aftershock, next_event, next_event_global;
    quakelib::Conversion    convert;
    quakelib::ModelEvent    new_event;

    // Calculate the current rates of stress change on all blocks
    stressRecompute();

    // Given the rates of change, determine which block will fail next
    nextStaticFailure(next_static_fail);

    // Get the next aftershock event time
    nextAftershock(next_aftershock);

    // Take whichever is sooner, with ties going in favor of aftershocks
    if (next_static_fail.val < next_aftershock.val) {
        next_event.val = next_static_fail.val;
        next_event.block_id = next_static_fail.block_id;
    } else {
        next_event.val = next_aftershock.val;
        next_event.block_id = next_aftershock.block_id;
    }

    // Each node now has the time before the first failure among its blocks
    // Determine the time to first failure over all nodes
    // (MPI_Reduce operation executed if MPI is present).
    sim->allReduceBlockVal(next_event, next_event_global, BLOCK_VAL_MIN);

    // If we didn't find any static failures or aftershocks, abort the simulation
    if (sim->isRootNode()) assertThrow(next_event_global.val < DBL_MAX, "System stuck, no blocks to move.");

    // Increment the simulation year to the next failure time and
    // update the slip on all other blocks
    sim->incrementYear(next_event_global.val);

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        BlockID gid = sim->getGlobalBID(lid);
        Block &local_block = sim->getBlock(gid);
        double cur_slip_deficit = sim->getSlipDeficit(gid);
        sim->setSlipDeficit(gid, cur_slip_deficit-local_block.slip_rate()*convert.year2sec(next_event_global.val)*(1.0-local_block.aseismic()));
    }

    // Recalculate the stress on all blocks with the new slip deficits
    stressRecompute();

    // Record the current event
    new_event.setEventTriggerOnThisNode(next_event_global.block_id==next_static_fail.block_id);
    new_event.setEventTrigger(next_event_global.block_id);
    new_event.setEventYear(sim->getYear());
    new_event.setEventNumber(sim->getEventCount());
    sim->addEvent(new_event);

    if (sim->getYear() > sim->getSimDuration()) return SIM_STOP_REQUIRED;
    else return SIM_CONTINUE;
}

/*!
 Calculate the time to the next aftershock. Aftershocks are only held
 on the root node so other nodes ignore this.
 */
void UpdateBlockStress::nextAftershock(BlockVal &next_aftershock) {
    if (sim->isRootNode() && sim->numAftershocksToProcess() > 0) {
        next_aftershock.val = sim->nextAftershockTime() - sim->getYear();
    } else {
        next_aftershock.val = DBL_MAX;
    }

    next_aftershock.block_id = UNDEFINED_ELEMENT_ID;
}

/*!
 Determine the next time step in the simulation when a failure occurs.
 Return the block ID of the block responsible for the failure and the timestep until the failure.
 */
void UpdateBlockStress::nextStaticFailure(BlockVal &next_static_fail) {
    BlockList::iterator     it;
    double                  ts;
    BlockID                 gid;
    int                     lid;
    quakelib::Conversion    convert;

    // Set up the temporary buffer and update field
    for (it=sim->begin(); it!=sim->end(); ++it) {
        tmpBuffer[it->getBlockID()] = 0.0;

        // Set the update field to be the slip rate of each block (note: these are local blockID values.)
        sim->setUpdateField(it->getBlockID(), it->slip_rate());
    }

    // update the temporary buffer with the Greens function applied to the block slip rates
    sim->matrixVectorMultiplyAccum(tmpBuffer,
                                   sim->greenShear(),
                                   sim->getUpdateFieldPtr(),
                                   true);

    //
    if (sim->doNormalStress()) {
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();
            sim->setUpdateField(gid, -sim->getFriction(gid)*it->slip_rate());
        }

        sim->matrixVectorMultiplyAccum(tmpBuffer,
                                       sim->greenNormal(),
                                       sim->getUpdateFieldPtr(),
                                       true);
    }

    //
    // Go through the blocks and find which one will fail first
    next_static_fail.val = DBL_MAX;
    next_static_fail.block_id = UNDEFINED_ELEMENT_ID;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        Block &block = sim->getBlock(gid);

        // Calculate the time until this block will fail
        // If the block has aseismic slip, the calculation is the exact solution of the
        // differential equation d(cff)/dt = rate_of_stress_change + cff*aseismic_frac*self_shear/recurrence
        if (block.aseismic() > 0) {
            double      A, B, K;
            A = -tmpBuffer[gid];
            B = -block.aseismic()*sim->getSelfStresses(gid)/sim->getRecurrence(gid);
            K = -log(A+B*sim->getCFF(gid))/B;
            ts = K + log(A)/B;
        } else {
            ts = convert.sec2year(sim->getCFF(gid)/tmpBuffer[gid]);
        }

        // Blocks with negative timesteps are skipped. These effectively mean the block
        // should have failed already but didn't. This can happen in the starting phase
        // of a simulation due to initial stresses. These blocks will fail anyway in the
        // event propagation section, so we can safely ignore them here. However, negative
        // timesteps should not happen after the simulation has progressed for a while.
        if (ts <= 0) continue;

        // If the time to slip is less than the current shortest time, record the block
        // To ensure reproducibility with multiple processes, if multiple blocks fail
        // at the same time then we choose the block with the lowest ID over all the processes
        // yoder: regarding a tie for failure time:: this is probably fine, but for a given instantiation of geometry, it
        // favors certain fault segments for failure. we should probably determine these randomly... unless we very specifically want
        // this part of the simulation to remain deterministic.
        if (ts < next_static_fail.val) {
            next_static_fail.block_id = gid;
            next_static_fail.val = ts;
        } else if (ts == next_static_fail.val) {
            next_static_fail.block_id = (gid < next_static_fail.block_id ? gid : next_static_fail.block_id);
        }
    }
}

/*!
 Recompute the stress on each block based on the slip deficits of all the other blocks.
 */
void UpdateBlockStress::stressRecompute(void) {
    BlockList::iterator it;

    // Reset the shear and normal stress, and set the update field to be the slip deficit
    for (it=sim->begin(); it!=sim->end(); ++it) {
        BlockID gid = it->getBlockID();
        sim->setShearStress(gid, 0.0);
        sim->setNormalStress(gid, sim->getRhogd(gid));
        sim->setUpdateField(gid, sim->getSlipDeficit(gid));
    }

    // Distribute the new update field over all nodes
    // (MPI_ calls when MPI enabled)
    sim->distributeUpdateField();

    // Multiply the Greens shear function by the slipDeficit vector
    // to get shear on local blocks
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

    // Recompute the CFF on blocks based on the new shear/normal stresses
    sim->computeCFFs();

}

void UpdateBlockStress::finish(SimFramework *_sim) {
    delete [] tmpBuffer;
}
