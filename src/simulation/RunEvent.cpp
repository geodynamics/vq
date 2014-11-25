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

/*!
 At the end of each sweep after we have recalculated block CFF, we determine
 which blocks will have a failure due to dynamic or static stress changes.
 */
void RunEvent::markBlocks2Fail(Simulation *sim, const FaultID &trigger_fault) {
    int         lid;
    BlockID     gid;
    bool        add;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);

        // Blocks can only fail once per event, after that they slide freely
        if (sim->getBlock(gid).getFailed()) continue;

        // Add this block if it has a static CFF failure
        add = sim->cffFailure(gid);

        // Allow dynamic failure if the block is "loose" (next to a previously failed block)
        if (loose_elements.count(gid) > 0) add |= sim->dynamicFailure(gid, trigger_fault);

        if (add) {
            sim->getBlock(gid).setFailed(true);
            local_failed_elements.insert(gid);
        }
    }
}

/*!
 Process the list of blocks that failed on this node using the original friction law.
 */
void RunEvent::processBlocksOrigFail(Simulation *sim, quakelib::ModelSweeps &sweeps) {
    quakelib::ElementIDSet::iterator    fit;
    double                              slip, stress_drop;

    // For each block that fails in this sweep, calculate how much it slips
    for (fit=local_failed_elements.begin(); fit!=local_failed_elements.end(); ++fit) {
        if (sim->isLocalBlockID(*fit)) {
            BlockID gid = *fit;
            Block &b = sim->getBlock(*fit);

            // calculate the drop in stress from the failure
            stress_drop = sim->getCFF0(gid) - sim->getCFF(gid);

            if (!stress_drop) stress_drop = sim->getStressDrop(gid) - sim->getCFF(gid);

            // Slip is in m
            slip = (stress_drop/sim->getSelfStresses(gid));

            if (slip < 0) slip = 0;

            // Record how much the block slipped in this sweep and initial stresses
            sweeps.setSlipAndArea(sweep_num,
                                  b.getBlockID(),
                                  slip,
                                  b.area(),
                                  b.lame_mu());
            sweeps.setInitStresses(sweep_num,
                                   b.getBlockID(),
                                   sim->getShearStress(gid),
                                   sim->getNormalStress(gid));

            sim->setSlipDeficit(gid, sim->getSlipDeficit(gid)+slip);
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

void RunEvent::processBlocksSecondaryFailures(Simulation *sim, quakelib::ModelSweeps &sweeps) {
    int             lid;
    BlockID         gid;
    unsigned int    i, n;
    quakelib::ElementIDSet          local_id_list;
    BlockIDProcMapping              global_id_list;
    quakelib::ElementIDSet::const_iterator      it;
    BlockIDProcMapping::const_iterator  jt;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        Block &b = sim->getBlock(gid);

        // If the block has already failed (but not in this sweep) then adjust the slip
        if (b.getFailed() && global_failed_elements.count(gid) == 0) {
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

    for (i=0,it=local_id_list.begin(); it!=local_id_list.end(); ++i,++it) {
        for (n=0,jt=global_id_list.begin(); jt!=global_id_list.end(); ++n,++jt) {
            A[i*num_global_failed+n] = sim->getGreenShear(*it, jt->first);

            if (sim->doNormalStress()) {
                A[i*num_global_failed+n] -= sim->getFriction(*it)*sim->getGreenNormal(*it, jt->first);
            }
        }

        b[i] = sim->getCFF(*it)+sim->getFriction(*it)*sim->getRhogd(*it);
    }

    if (sim->isRootNode()) {
        double *fullA = new double[num_global_failed*num_global_failed];
        double *fullb = new double[num_global_failed];
        double *fullx = new double[num_global_failed];

        // Fill in the A matrix and b vector from the various processes
        for (i=0,n=0,jt=global_id_list.begin(); jt!=global_id_list.end(); ++jt,++i) {
            if (jt->second != sim->getNodeRank()) {
#ifdef MPI_C_FOUND
                MPI_Recv(&(fullA[i*num_global_failed]), num_global_failed, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(fullb[i]), 1, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
                assertThrow(false, "Single processor version of code, but faults mapped to multiple processors.");
#endif
            } else {
                memcpy(&(fullA[i*num_global_failed]), &(A[n*num_global_failed]), sizeof(double)*num_global_failed);
                memcpy(&(fullb[i]), &(b[n]), sizeof(double));
                n++;
            }
        }

        // Solve the global system on the root node
        solve_it(num_global_failed, fullx, fullA, fullb);

        // Send back the resulting values from x to each process
        for (i=0,n=0,jt=global_id_list.begin(); jt!=global_id_list.end(); ++jt,++i) {
            if (jt->second != sim->getNodeRank()) {
#ifdef MPI_C_FOUND
                MPI_Send(&(fullx[i]), 1, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD);
#else
                assertThrow(false, "Single processor version of code, but faults mapped to multiple processors.");
#endif
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
#ifdef MPI_C_FOUND

        for (i=0; i<num_local_failed; ++i) {
            MPI_Send(&(A[i*num_global_failed]), num_global_failed, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(b[i]), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        for (i=0; i<num_local_failed; ++i) {
            MPI_Recv(&(x[i]), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

#else
        assertThrow(false, "Single processor version of code, but processor MPI rank is non-zero.");
#endif
    }

    // Take the results of the calculation and determine how much each ruptured block slipped
    for (i=0,it=local_id_list.begin(); it!=local_id_list.end(); ++i,++it) {
        Block &block = sim->getBlock(*it);

        double slip = x[i] - sim->getSlipDeficit(*it);;

        if (slip > 0) {
            // Record how much the block slipped in this sweep and initial stresses
            sweeps.setSlipAndArea(sweep_num,
                                  *it,
                                  slip,
                                  block.area(),
                                  block.lame_mu());
            sweeps.setInitStresses(sweep_num,
                                   *it,
                                   sim->getShearStress(*it),
                                   sim->getNormalStress(*it));

            sim->setSlipDeficit(*it, sim->getSlipDeficit(*it)+slip);
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
void RunEvent::processStaticFailure(Simulation *sim) {
    BlockList::iterator     it;
    quakelib::ModelSweeps   event_sweeps;
    BlockID                 triggerID, gid;
    FaultID                 trigger_fault;
    int                     more_blocks_to_fail;
    bool                    final_sweep = false;
    std::pair<quakelib::ElementIDSet::const_iterator, quakelib::ElementIDSet::const_iterator> nbr_start_end;
    quakelib::ElementIDSet::const_iterator  nit;

    // Get the event trigger block and fault
    triggerID = sim->getCurrentEvent().getEventTrigger();
    trigger_fault = sim->getBlock(triggerID).getFaultID();
    sweep_num = 0;

    // Clear the list of "loose" (can dynamically fail) blocks
    loose_elements.clear();

    // Clear the list of failed blocks, and add the trigger block
    local_failed_elements.clear();

    if (sim->getCurrentEvent().getEventTriggerOnThisNode()) {
        local_failed_elements.insert(triggerID);
        loose_elements.insert(triggerID);
        sim->getBlock(triggerID).setFailed(true);
    }

    more_blocks_to_fail = sim->blocksToFail(!local_failed_elements.empty());

    // While there are still failed blocks to handle
    while (more_blocks_to_fail || final_sweep) {
        // Share the failed blocks with other processors to correctly handle
        // faults that are split among different processors
        sim->distributeBlocks(local_failed_elements, global_failed_elements);

        // Process the blocks that failed
        processBlocksOrigFail(sim, event_sweeps);

        // Recalculate CFF for all blocks where slipped blocks don't contribute
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();
            sim->setShearStress(gid, 0.0);
            sim->setNormalStress(gid, sim->getRhogd(gid));
            sim->setUpdateField(gid, (it->getFailed() ? 0 : sim->getSlipDeficit(gid)));
        }

        // Distribute the update field values to other processors
        sim->distributeUpdateField();

        // Set dynamic triggering on for any blocks neighboring blocks that slipped in the last sweep
        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
            // Add block neighbors if the block has slipped
            if (sim->getUpdateField(gid) > 0) {
                nbr_start_end = sim->getNeighbors(gid);

                for (nit=nbr_start_end.first; nit!=nbr_start_end.second; ++nit) {
                    loose_elements.insert(*nit);
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
        processBlocksSecondaryFailures(sim, event_sweeps);

        // Set the update field to the slip of all blocks
        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
            Block &block=sim->getBlock(gid);
            sim->setShearStress(gid, 0.0);
            sim->setNormalStress(gid, sim->getRhogd(gid));
            sim->setUpdateField(gid, (block.getFailed() ? 0 : sim->getSlipDeficit(gid)));
        }

        sim->distributeUpdateField();

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

        for (fit=global_failed_elements.begin(); fit!=global_failed_elements.end(); ++fit) {
            if (sim->isLocalBlockID(fit->first)) {
                event_sweeps.setFinalStresses(sweep_num,
                                              fit->first,
                                              sim->getShearStress(fit->first),
                                              sim->getNormalStress(fit->first));
            }
        }

        global_failed_elements.clear(); // we are done with these blocks
        local_failed_elements.clear(); // we are done with these blocks

        // Find any blocks that fail because of the new stresses
        markBlocks2Fail(sim, trigger_fault);

        if (final_sweep) {
            final_sweep = false;
        } else {
            more_blocks_to_fail = sim->blocksToFail(!local_failed_elements.empty());

            if (!more_blocks_to_fail) final_sweep = true;
        }

        sweep_num++;
    }

    // Set the completed list as the sweep list for the entire event
    sim->collectEventSweep(event_sweeps);
    sim->getCurrentEvent().setSweeps(event_sweeps);
}

/*!
 Process the next aftershock. This involves determining a suitable rupture area from an empirical relationship, finding the nearest elements,
 */
void RunEvent::processAftershock(Simulation *sim) {
    std::map<double, BlockID>                   as_elem_dists;
    std::map<double, BlockID>::const_iterator   it;
    std::map<BlockID, double>                   elem_slips;
    EventAftershock                           as;
    BlockID                                     gid;
    quakelib::ElementIDSet                      id_set;
    quakelib::ElementIDSet::const_iterator      bit;
    quakelib::Conversion                        convert;
    quakelib::ModelSweeps                       event_sweeps;

    // Pop the next aftershock off the list
    as = sim->popAftershock();

    // Calculate the distance from the aftershock to all elements
    for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
        double as_to_elem_dist = sim->getBlock(gid).center().dist(as.loc());
        as_elem_dists.insert(std::make_pair(as_to_elem_dist, gid));
    }

    // Determine the target rupture area given the aftershock magnitude
    // TODO: make this equation user specifiable?
    double rupture_area = pow(10, -3.49+0.91*as.mag);
    double selected_rupture_area = 0;
    double selected_rupture_area_mu = 0;

    // Go through the elements, closest first, until we find enough to match the rupture area
    for (it=as_elem_dists.begin(); it!=as_elem_dists.end(); ++it) {
        Block &b=sim->getBlock(it->second);
        selected_rupture_area += convert.sqm2sqkm(b.area());
        selected_rupture_area_mu += b.area()*b.lame_mu();
        id_set.insert(it->second);

        if (selected_rupture_area > rupture_area) break;
    }

    // Determine the amount of slip needed to match the aftershock magnitude
    // The contribution of each block to moment is based on its fraction of total area*mu
    double total_moment = pow(10, (as.mag + 10.7)*(3.0/2.0))/1e7;

    for (bit=id_set.begin(); bit!=id_set.end(); ++bit) {
        Block &b=sim->getBlock(*bit);

        // Calculate the slip based on the earthquake moment
        double element_moment = total_moment*(b.area()*b.lame_mu()/selected_rupture_area_mu);
        double element_slip = element_moment/(b.lame_mu()*b.area());

        // Adjust the slip deficit appropriately
        sim->setSlipDeficit(*bit, sim->getSlipDeficit(*bit)+element_slip);

        // Create the sweep describing this aftershock
        // Since we don't distinguish sweeps, every slip occurs in sweep 0
        event_sweeps.setSlipAndArea(0, *bit, element_slip, b.area(), b.lame_mu());
        event_sweeps.setInitStresses(0, *bit, sim->getShearStress(*bit), sim->getNormalStress(*bit));
    }

    // Set the update field to the slip of all blocks
    for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
        sim->setShearStress(gid, 0.0);
        sim->setNormalStress(gid, sim->getRhogd(gid));
        sim->setUpdateField(gid, sim->getSlipDeficit(gid));
    }

    sim->distributeUpdateField();

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

    // Record final stresses on each block involved in the aftershock
    for (bit=id_set.begin(); bit!=id_set.end(); ++bit) {
        event_sweeps.setFinalStresses(0, *bit, sim->getShearStress(*bit), sim->getNormalStress(*bit));
    }

    // Add sweeps to list
    sim->getCurrentEvent().setSweeps(event_sweeps);
}

SimRequest RunEvent::run(SimFramework *_sim) {
    Simulation            *sim = static_cast<Simulation *>(_sim);
    int                     lid;

    // Save stress information at the beginning of the event
    // This is used to determine dynamic block failure
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) sim->saveStresses(sim->getGlobalBID(lid));

    // If there's a specific block that triggered the event, it's a static stress failure type event
    if (sim->getCurrentEvent().getEventTrigger() != UNDEFINED_ELEMENT_ID) {
        processStaticFailure(sim);
    } else {
        // Otherwise it's an aftershock
        processAftershock(sim);
    }

    // Record the stress in the system before and after the event
    recordEventStresses(sim);

    // Update the cumulative slip for this fault
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        BlockID gid = sim->getGlobalBID(lid);
        sim->getBlock(gid).setFailed(false);
    }

    // TODO: reinstate this check
    // TODO: currently fails because processors may locally have no failures
    //assertThrow(sim->getCurrentEvent().size() > 0, "There was a trigger but no failed blocks.");

    return SIM_STOP_OK;
}

void RunEvent::recordEventStresses(Simulation *sim) {
    quakelib::ElementIDSet involved_blocks;
    double shear_init, shear_final, normal_init, normal_final;
    double total_shear_init, total_shear_final, total_normal_init, total_normal_final;

    involved_blocks = sim->getCurrentEvent().getInvolvedElements();

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
