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

        // Blocks can only fail once per event, after that they slide freely (in secondary ruptures).
        // Also, we don't want to add elements that have over-slipped and have very negative CFFs (i.e. negative shear stress).
        // These elements won't statically fail, but they may meet the dynamic failure criteria, so skip them.
        if (sim->getFailed(gid) || sim->getCFF(gid) < sim->getMaxStressDrop(gid)) continue;

        // Add this block if it has a static or dynamic CFF failure
        add = sim->cffFailure(gid) ||  sim->dynamicFailure(gid, trigger_fault);

        if (add) {
            sim->setFailed(gid, true);
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

            ///// Schultz: This has been moved to the function called getEffectiveStressDrop in Simulation.h
            //stress_drop = sim->getStressDrop(gid) - sim->getCFF(gid);

            // Slip is in m
            slip = sim->getEffectiveStressDrop(gid)/sim->getSelfStresses(gid);
            /// Schultz: The hack below implements Eric's version. It's non-ideal, but it works. To be improved later.
            //slip = -1.0*(sim->getSlipDeficit(gid));


            ////// Schultz:
            // Do not allow negative slips, it means you are allowing an element to gain stress during an event.
            if (slip > 0) {

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
            } else {
              // Schultz; If slip <= 0, then CFF <= max_stress_drop and we must have over-slipped.
              //   Therefore we should consider it as not failed, so it won't contribute any later in the rupture.
                if (sim->getCFF(gid) <= sim->getMaxStressDrop(gid)) {
                    sim->setFailed(gid, false);
                }
            }
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
    // yoder:  This bit of code is a likely candidate for the heisenbug/heisen_hang problem. basically, i think the is_root(),send/receive
    // logic loop has a tendency to get hung up for complex operations. revise that code block, nominally into two "isRoot()" blocks.
    // 1) first, distribute the A,B arrays (an array and a vector) between the nodes.
    // 2) then, multiply, etc.
    // 3) then redistribute the result back to the various nodes.
    // basically move the second part of the isRoot() (don't recall how the not isRoot() block looks) outside the send/recv block.
    //
    int             lid;
    BlockID         gid;
    unsigned int    i, n;
    quakelib::ElementIDSet          local_secondary_id_list;  // lists of local/global secondary failures.
    BlockIDProcMapping              global_secondary_id_list;
    //
    quakelib::ElementIDSet::const_iterator      it;
    BlockIDProcMapping::const_iterator  jt;

    //
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        //for (quakelib::ModelSweeps::iterator s_it=sweeps.begin(); s_it!=sweeps.end(); ++s_it) {
        // would a faster way to do this step be to look through the current event_sweeps list?
        //yes, but we're assuming that "original failure" has been processed,
        // which i think is a pretty safe bet. BUT, let's leave the original looping code in comment, to facilitate an easy recovery if this is a mistake.
        // another possible concern is keepting track of local/global blocks. for now, let's leave this alone. it is a (relatively) small matter of optimization.
        //lid = s_it->_element_id;
        gid = sim->getGlobalBID(lid);

        //
        // If the block has already failed (but not in this sweep) then adjust the slip
        if (sim->getFailed(gid) && global_failed_elements.count(gid) == 0) {
            //local_id_list.insert(gid);
            local_secondary_id_list.insert(gid);
        }
    }

    //
    // use: global/local_secondary_id_list;
    // Figure out how many failures there were over all processors
    //sim->distributeBlocks(local_id_list, global_id_list);
    // can this somehow distribute a block to global_failed_elements twice? (multiple copies of same value?)
    // yoder (note): after we distributeBlocks(), we can check to see that all items in local_ exist in global_ exactly once.
    // if not, throw an exception... and then we'll figure out how this is happening. remember, local_ is like [gid, gig, gid...]
    // global_ is like [(gid, p_rank), (gid, p_rank)...], and each pair item is accessed like global_[rw_num]->first /->second
    sim->distributeBlocks(local_secondary_id_list, global_secondary_id_list);

    // ==== DYNAMIC STRESS DROPS ========== 
    // Schultz: now that we know how many elements are involved, assign dynamic stress drops
    if (sim->doDynamicStressDrops()) {
        double dynamicStressDrop;
        quakelib::ElementIDSet::const_iterator cit;
        BlockIDProcMapping::const_iterator  bit;

        // Compute the current event area
        // Add in global_failed_elements
        for (bit=global_failed_elements.begin(); bit!=global_failed_elements.end(); ++bit) {
            // Avoid double counting
            if (!all_event_blocks.count(bit->first)) {
                current_event_area += sim->getBlock(bit->first).area();
                all_event_blocks.insert(bit->first);
            }
        }

        // Also add in the area from the secondary failed elements
        for (bit=global_secondary_id_list.begin(); bit!=global_secondary_id_list.end(); ++bit) {
            // Avoid double counting
            if (!all_event_blocks.count(bit->first)) {
                current_event_area += sim->getBlock(bit->first).area();
                all_event_blocks.insert(bit->first);
            }
        }

        for (cit=all_event_blocks.begin(); cit!=all_event_blocks.end(); ++cit) {
            if (current_event_area < sim->getFaultArea(sim->getBlock(*cit).getFaultID())) {
                // If the current area is smaller than the section area, scale the stress drop
                dynamicStressDrop = sim->computeDynamicStressDrop(*cit, current_event_area);
                
                // Try to make the triggering element slip enough to trigger others.
                if (*cit == sim->getCurrentEvent().getEventTrigger() && sweep_num == 0) {
                    sim->setStressDrop(*cit, sim->getMaxStressDrop(*cit), false);
                } else {
                    sim->setStressDrop(*cit, dynamicStressDrop, false);
                }
                
            } else {
                sim->setStressDrop(*cit, sim->getMaxStressDrop(*cit), false);
            }
        }
    }
    // Note: For multiprocessing, we do not need to distribute/communicate these changes to the stress drops.
    // We have computed new stress drops for our local elements, and when we communicate the b-vector around,
    // the appropriate stress drops will be communicated to other processors. 


    //int num_local_failed = local_id_list.size();
    //int num_global_failed = global_id_list.size();
    int num_local_failed = local_secondary_id_list.size();
    int num_global_failed = global_secondary_id_list.size();

    double *A = new double[num_local_failed*num_global_failed];
    double *b = new double[num_local_failed];
    double *x = new double[num_local_failed];

    //
    // stress transfer (greens functions) between each local element and all global elements.
    for (i=0,it=local_secondary_id_list.begin(); it!=local_secondary_id_list.end(); ++i,++it) {
        for (n=0,jt=global_secondary_id_list.begin(); jt!=global_secondary_id_list.end(); ++n,++jt) {

            A[i*num_global_failed+n] = sim->getGreenShear(*it, jt->first);

            if (sim->doNormalStress()) {
                A[i*num_global_failed+n] -= sim->getFriction(*it)*sim->getGreenNormal(*it, jt->first);
            }
        }

        ///// Schultz:
        // Even if we are doing dynamic stress drops, they've already been set. Check processStaticFailure() and
        // the beginning of this method
        //b[i] = sim->getStressDrop(*it) - sim->getCFF(*it);
        b[i] = sim->getEffectiveStressDrop(*it);
    }

    //
    // so A,b are calculated for each local node (with dimension n_local x n_global and now they'll be consolidated on the root node. note that
    // the corresponding mpi_send comes after this block. the root node blocks until child nodes have sent A,b and root_node has received A,b.

    /////////////
    //
    if (sim->isRootNode()) {
        double *fullA = new double[num_global_failed*num_global_failed];
        double *fullb = new double[num_global_failed];
        double *fullx = new double[num_global_failed];

        // Fill in the A matrix and b vector from the various processes
        //for (i=0,n=0,jt=global_id_list.begin(); jt!=global_id_list.end(); ++jt,++i) {
        for (i=0,n=0,jt=global_secondary_id_list.begin(); jt!=global_secondary_id_list.end(); ++jt,++i) {
            if (jt->second != sim->getNodeRank()) {
#ifdef MPI_C_FOUND
                //
                // element NOT local to this (root) node, so receive this data element from the remote process with mpi_rank jt->second:
                // note: jt-> first: global_id, jt->second: node_id/rank
                // MPI_Recv(in_buff{out}, in_len, mpi_dtype, src_node, tag, comm, status{out})
                // note that in_buff and status are technically "output" parameters for the MPI_Recv function; we read data into in_buff by
                // calling MPI_Recv()
                MPI_Recv(&(fullA[i*num_global_failed]), num_global_failed, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(fullb[i]), 1, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
                assertThrow(false, "Single processor version of code, but faults mapped to multiple processors.");
#endif
            } else {
                // YES, element is local to this node; copy from the "local" array to the "global" array:
                memcpy(&(fullA[i*num_global_failed]), &(A[n*num_global_failed]), sizeof(double)*num_global_failed);
                memcpy(&(fullb[i]), &(b[n]), sizeof(double));
                n++;
            }
        }

        //
        // Solve the global system on the root node (we're inside an if (sim->isRootNode()) {} block )
        solve_it(num_global_failed, fullx, fullA, fullb);

        // Send back the resulting values from x to each process
        //for (i=0,n=0,jt=global_id_list.begin(); jt!=global_id_list.end(); ++jt,++i) {
        for (i=0,n=0,jt=global_secondary_id_list.begin(); jt!=global_secondary_id_list.end(); ++jt,++i) {
            if (jt->second != sim->getNodeRank()) {
#ifdef MPI_C_FOUND
                // send these values to node-rank jt->second:
                // yoder: try using synchronous send, MPI_Ssend()
                //MPI_Send(&(fullx[i]), 1, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD);
                /// THIS WAS THE FIX FOR HEISENBUG ////////////////////////
                MPI_Ssend(&(fullx[i]), 1, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD);
                //////////////////////////////////////////////////////////
#else
                assertThrow(false, "Single processor version of code, but faults mapped to multiple   processors.");
#endif
            } else {
                // yoder: so we're copying the revised xfull[j] -- x[n] ('global' to 'local'); x[n] are the local
                // node's values. aka, node sends these values back to itself.
                memcpy(&(x[n]), &(fullx[i]), sizeof(double));
                n++;
            }
        }

        //
        // Delete the memory arrays created (use delete [] for arrays)
        delete [] fullx;
        delete [] fullb;
        delete [] fullA;
    } else {
        // NOT root_node:
#ifdef MPI_C_FOUND
        // send these values to root (rank 0) node:
        // Child Nodes do their sending bit:
        for (i=0; i<num_local_failed; ++i) {
            // yoder: We must use synchronous MPI_Ssend() (this waits for all processors to report back before proceeding).
            //MPI_Send(&(A[i*num_global_failed]), num_global_failed, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            //MPI_Send(&(b[i]), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            /// THIS WAS THE FIX FOR HEISENBUG ////////////////////////
            MPI_Ssend(&(A[i*num_global_failed]), num_global_failed, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Ssend(&(b[i]), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }


        //
        for (i=0; i<num_local_failed; ++i) {
            MPI_Recv(&(x[i]), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

#else
        assertThrow(false, "Single processor version of code, but processor MPI rank is non-zero.");
#endif
    }

    // Take the results of the calculation and determine how much each ruptured block slipped
    //for (i=0,it=local_id_list.begin(); it!=local_id_list.end(); ++i,++it) {
    for (i=0,it=local_secondary_id_list.begin(); it!=local_secondary_id_list.end(); ++i,++it) {
        Block &block = sim->getBlock(*it);
        ///////////////
        // Schultz:: The matrix solution solves for the slip, not the final slip deficit.
        double slip = x[i];
        ///////////////
        
        ////// Schultz:
        // Must not allow negative slips. It prevents periodicity in single fault sims.
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
            //
            sim->setSlipDeficit(*it, sim->getSlipDeficit(*it)+slip);
        } else {
              // Schultz; If slip <= 0, then CFF <= stress_drop and we must have over-slipped.
              //   Therefore we should consider it as not failed, so it won't contribute any later in the rupture.
                if (sim->getCFF(gid) <= sim->getMaxStressDrop(gid)) {
                    sim->setFailed(gid, false);
                }
        }
    }

    //
    // delete/de-allocate arrays (use "delete []" for arrays, as opposed to "delete" for single objects)
    delete [] A;
    delete [] b;
    delete [] x;
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
    BlockID                 triggerID;          // limit variable gid to local loop scopes.
    FaultID                 trigger_fault;
    int                     more_blocks_to_fail;
    bool                    final_sweep = false;
    std::pair<quakelib::ElementIDSet::const_iterator, quakelib::ElementIDSet::const_iterator> nbr_start_end;
    quakelib::ElementIDSet::const_iterator  nit;
    // Get the event trigger block and fault
    triggerID = sim->getCurrentEvent().getEventTrigger();
    trigger_fault = sim->getBlock(triggerID).getFaultID();
    sweep_num = 0;

    // Clear the list of failed blocks, and add the trigger block
    local_failed_elements.clear();
    
    if (sim->getCurrentEvent().getEventTriggerOnThisNode()) {
        local_failed_elements.insert(triggerID);
        sim->setFailed(triggerID, true);
    }
    
    // Keep track of all the elements that have slipped in the event.
    // Use this to determine current event area.
    double current_event_area = sim->getBlock(triggerID).area();
    all_event_blocks.clear();
    all_event_blocks.insert(triggerID);
    
    
    // yoder (note): Comm::blocksToFail() executes a single MPI_Allreduce() (when MPI is present)
    more_blocks_to_fail = sim->blocksToFail(!local_failed_elements.empty());


    // While there are still failed blocks to handle
    while (more_blocks_to_fail || final_sweep) {
        // Share the failed blocks with other processors to correctly handle
        // faults that are split among different processors
        sim->distributeBlocks(local_failed_elements, global_failed_elements);


        ///////////////////////////////////////////////////////////////////
        // Schultz:: Uncomment the following to write out simulation variables during a simulation
        //     that are other wise unobservable.
        // -----  OUTPUT HACK, only works on 1 proc --------------------
        //        for (unsigned int gid=0; gid<sim->numGlobalBlocks(); ++gid) {
        //            sim->console() << sweep_num << "  " << gid << "  " << sim->getShearStress(gid) << "  " << sim->getNormalStress(gid) << "  " << sim->getCFF(gid) << "  " << sim->getStressDrop(gid) << std::endl;
        //        }
        ///////////////////////////////////////////////////////////////////


        // ==== DYNAMIC STRESS DROPS ========== 
        // Schultz: now that we know how many elements are involved, assign dynamic stress drops
        if (sim->doDynamicStressDrops()) {
            double dynamicStressDrop;
            quakelib::ElementIDSet::const_iterator cit;
            BlockIDProcMapping::const_iterator  bit;

            // Compute the current event area
            // Add in global_failed_elements
            for (bit=global_failed_elements.begin(); bit!=global_failed_elements.end(); ++bit) {
                // Avoid double counting
                if (!all_event_blocks.count(bit->first)) {
                    current_event_area += sim->getBlock(bit->first).area();
                    all_event_blocks.insert(bit->first);
                }
            }

            for (cit=all_event_blocks.begin(); cit!=all_event_blocks.end(); ++cit) {
                if (current_event_area < sim->getFaultArea(sim->getBlock(*cit).getFaultID())) {
                    // If the current area is smaller than the section area, scale the stress drop
                    dynamicStressDrop = sim->computeDynamicStressDrop(*cit, current_event_area);
                    // Try to make the triggering element slip enough to trigger others.
                    if (*cit == triggerID && sweep_num == 0) {
                        sim->setStressDrop(*cit, sim->getMaxStressDrop(*cit), false);
                    } else {
                        sim->setStressDrop(*cit, dynamicStressDrop, false);
                    }
                } else {
                    sim->setStressDrop(*cit, sim->getMaxStressDrop(*cit), false);
                }
            }
        }

        // Process the slips for the blocks that have failed.
        // note: setInitStresses() called in processBlocksOrigFail().
        // note: processBlocksOrigFail() is entirely local (no MPI).
        processBlocksOrigFail(sim, event_sweeps);

        // Now, each process has updated slips for its elements but only that process has the update.
        // We must communicate these updates among all processes.
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();
            sim->setShearStress(gid, 0.0);
            sim->setNormalStress(gid, sim->getRhogd(gid));
            //sim->setUpdateField(gid, (sim->getFailed(gid) ? 0 : sim->getSlipDeficit(gid)));
            ///////// Schultz:
            // Update the stresses using the current slip of all elements, or else we throw away the slips 
            //   computed in processBlocksOrigFail().
            sim->setUpdateField(gid, sim->getSlipDeficit(gid));
            // Although we are adding slip deficits for all elements not just the local ones, when we execute the 
            //   distributeUpdateField() command below only the local elements are selected.
        }

        // Distribute the update field values to other processors
        sim->distributeUpdateField();

        
        
        // Calculate the new CFFs based on the slips computed in processBlocksOrigFail()
        // multiply greenSchear() x getUpdateFieldPtr() --> getShearStressPtr() ... right?
        // assign stress values (shear stresses at this stage are all set to 0; normal stresses are set to sim->getRhogd(gid) -- see code a couple paragraphs above.
        //
        // The following matrix operation adds in normal/shear changes due to the current slip deficits and the Greens function interaction matrix.
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


        //
        // Create the matrix equation, including interactions, and solve the system for slips
        processBlocksSecondaryFailures(sim, event_sweeps);

        // Set the update field to the slip of all blocks
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();
            sim->setShearStress(gid, 0.0);
            sim->setNormalStress(gid, sim->getRhogd(gid));

            //////// Schultz: try setting updatefield 0 for newly failed elements this sweep
            //sim->setUpdateField(gid, ((sim->getFailed(gid) && global_failed_elements.count(gid) == 0) ? 0 : sim->getSlipDeficit(gid)));
            // We need to ensure our slip economics books are balanced. I suspect we need here
            // instead: sim->setUpdateField(gid, sim->getSlipDeficit(gid) ). Update the stresses using
            // the current slip of all elements, or else we throw away the slip information from failed elements.
            sim->setUpdateField(gid, sim->getSlipDeficit(gid));
        }


        // Communicate the slip deficits between processors.
        sim->distributeUpdateField();
        
        // Calculate the new shear stresses and CFFs given the new update field values
        sim->matrixVectorMultiplyAccum(sim->getShearStressPtr(),
                                       sim->greenShear(),
                                       sim->getUpdateFieldPtr(),
                                       true);

        //
        if (sim->doNormalStress()) {
            sim->matrixVectorMultiplyAccum(sim->getNormalStressPtr(),
                                           sim->greenNormal(),
                                           sim->getUpdateFieldPtr(),
                                           true);
        }

        //
        sim->computeCFFs();

        // ------------------------------------------------------------------------------------------------------------
        
        
        
        // and this loop updates final_stress values by looping directly over the current sweeps list.
        //  note that event_sweeps is of type quakelib::ModelSweeps, which contains a vector<SweepData> _sweeps.
        for (quakelib::ModelSweeps::iterator s_it=event_sweeps.begin(); s_it!=event_sweeps.end(); ++s_it) {
            //
            // yoder: as per request by KS, change std::isnan() --> std::isnan(); std::isnan() appears to throw an error on some platforms.
            // Eric: Probably don't need this if check
            if (isnan(s_it->_shear_final) and isnan(s_it->_normal_final)) {
                // note: the stress entries are initialized with nan values, but if there are cases where non nan values need to be updated,
                // this logic should be revisited.
                event_sweeps.setFinalStresses(sweep_num,
                                              s_it->_element_id,
                                              sim->getShearStress(s_it->_element_id),
                                              sim->getNormalStress(s_it->_element_id));

            }
        }

        //
        global_failed_elements.clear(); // we are done with these blocks
        local_failed_elements.clear();  // we are done with these blocks

        
        // Find any blocks that fail because of the new stresses (all local; no MPI).
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

    // Schultz: Dynamic stress drops
    // Now that the event is over, reset the stress drops to the inter-event values
    if (sim->doDynamicStressDrops()) {
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();
            sim->setStressDrop(gid, sim->getMaxStressDrop(gid));
        }
    }


    // Write stress state to the stress output file if we're
    // at a multiple of N_events = sim->getStressOutInterval().
    unsigned int evnum = sim->getCurrentEvent().getEventNumber();

    if (evnum >= sim->getStressOutInterval() && evnum%sim->getStressOutInterval() == 0) {
        sim->output_stress(sim->getCurrentEvent().getEventNumber());
        sim->console() << std::endl << "--- Writing sim stress state to file after event " << sim->getCurrentEvent().getEventNumber() << " ---" << std::endl << std::flush;
    }

}

/*!
 Process the next aftershock. This involves determining a suitable rupture area from an empirical
 relationship, finding the nearest elements, choosing enough elements to match the empirical
 relationship, calculating the slip needed to generate the aftershock, and updating the stress
 field appropriately.
 */
void RunEvent::processAftershock(Simulation *sim) {
    std::map<double, BlockID>                   as_elem_dists;
    std::map<double, BlockID>::const_iterator   it;
    std::map<BlockID, double>                   elem_slips;
    EventAftershock                             as;
    BlockID                                     gid;
    quakelib::ElementIDSet                      id_set;
    quakelib::ElementIDSet::const_iterator      bit;
    quakelib::Conversion                        convert;
    quakelib::ModelSweeps                       event_sweeps;

    // Set the update field to the slip on each block
    for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
        sim->setShearStress(gid, 0.0);
        sim->setNormalStress(gid, sim->getRhogd(gid));
        sim->setUpdateField(gid, sim->getSlipDeficit(gid));
    }

    // And distribute this around
    sim->distributeUpdateField();

    // Only process the aftershock stress effects on the root node
    if (sim->isRootNode()) {
        // Pop the next aftershock off the list
        as = sim->popAftershock();

        // Calculate the distance from the aftershock to all elements
        ////// Schultz: This is very inefficient. Makes aftershock simulations slow.
        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
            double as_to_elem_dist = sim->getBlock(gid).center().dist(as.loc());
            as_elem_dists.insert(std::make_pair(as_to_elem_dist, gid));
        }

        // Determine the target rupture area given the aftershock magnitude
        // TODO:
        // user_defined_constants (flag this for later revisions in which we move these contant definitions to a parameters file).
        double rupture_area = convert.sqkm2sqm(pow(10, as.mag-4.0));
        // Schultz: This scaling relation returns rupture area in km^2. Lets keep everything in M.K.S. units.
        // Scaling relations come from Leonard 2010 paper:
        // "Earthquake Fault Scaling: Self-Consistent Relating of Rupture Length, Width, Average Displacement, and Moment Release"
        double selected_rupture_area = 0;
        double selected_rupture_area_mu = 0;

        // Go through the elements, closest first, until we find enough to match the rupture area
        for (it=as_elem_dists.begin(); it!=as_elem_dists.end(); ++it) {
            Block &b=sim->getBlock(it->second);
            selected_rupture_area += b.area();
            selected_rupture_area_mu += b.area()*b.lame_mu();
            id_set.insert(it->second);

            // If this is the first aftershock element, assign it as the event trigger
            if (it==as_elem_dists.begin()) sim->getCurrentEvent().setEventTrigger(it->second);

            if (selected_rupture_area > rupture_area) break;
        }

        // Determine the amount of slip needed to match the aftershock magnitude
        // The contribution of each block to moment is based on its fraction of total area*mu
        double total_moment = pow(10, (as.mag + 6.0)*(3.0/2.0));

        for (bit=id_set.begin(); bit!=id_set.end(); ++bit) {
            Block &b=sim->getBlock(*bit);

            // Calculate the slip based on the earthquake moment
            double element_moment = total_moment*(b.area()*b.lame_mu()/selected_rupture_area_mu);
            double element_slip = element_moment/(b.lame_mu()*b.area());

            // Adjust the slip deficit appropriately
            sim->setUpdateField(*bit, sim->getUpdateField(*bit)+element_slip);

            // Create the sweep describing this aftershock
            // Since we don't distinguish sweeps, every slip occurs in sweep 0
            event_sweeps.setSlipAndArea(0, *bit, element_slip, b.area(), b.lame_mu());
            event_sweeps.setInitStresses(0, *bit, sim->getShearStress(*bit), sim->getNormalStress(*bit));
        }

    }

    // Broadcast the new slip deficit from the root node to all processes
    sim->broadcastUpdateField();

    // And update the slip deficit on each process to take this into account
    for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
        sim->setSlipDeficit(gid, sim->getUpdateField(gid));
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

    if (sim->getCurrentEvent().getEventTrigger() != UNDEFINED_ELEMENT_ID) {
        processStaticFailure(sim);
    } else {
        // Otherwise it's an aftershock
        processAftershock(sim);
    }

    // Record the stress in the system before and after the event.
    recordEventStresses(sim);

    // Reset the failed status for each local block
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        BlockID gid = sim->getGlobalBID(lid);
        sim->setFailed(gid, false);
    }

    // TODO: reinstate this check
    // TODO: currently fails because processors may locally have no failures
    // TODO: notes to self(s) then: this single line works in SPP mode, or more specifically for a single node, so we can use an MPI_reduce() call to get the max or sum
    // ... and could this be causing heisen_hang? would this create a scenario where a child node would send/receive block data to root
    // but the root node would not send back anything (block not in list of failed blocks), so that child node would just hang there?
    // maybe, but i think that by this time, we'd be long since past that point.
    //       of all local sim->getCurrentEvent().size() values.
    //
    //assertThrow(sim->getCurrentEvent().size() > 0, "There was a trigger but no failed blocks.");
    //
#ifdef MPI_C_FOUND
    // so let's get the total event size by doing an MPI sum. note this can be an MPI_Reduce() to the root node only, or we can count on all processors.
    // this should be an equivalent operation (though slower in the latter case); i guess then, if we gather only to the root noode, we do the assert/throw
    // only on root node.:
    int local_event_size = sim->getCurrentEvent().size();
    int global_event_size = 0;
    //
    // aggregate and assert on root node:
    MPI_Reduce(&local_event_size, &global_event_size, 1, MPI_INT, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);

    if (sim->isRootNode()) {
        assertThrow(global_event_size > 0, "There was a trigger but no failed blocks.");
    };

    //
    //// aggregate and assert on all nodes:
    //MPI_Allreduce(&local_event_size, &global_event_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //assertThrow(sim->getCurrentEvent().size() > 0, "There was a trigger but no failed blocks.");
#else
    //global_event_size=local_event_size;
    int global_event_size = sim->getCurrentEvent().size();

    assertThrow(global_event_size > 0, "There was a trigger but no failed blocks. (" << getpid() << "/" << sim->getNodeRank() << ")");

#endif
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
