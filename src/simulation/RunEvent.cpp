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
        if (sim->getFailed(gid)) continue;

        // Add this block if it has a static CFF failure
        add = sim->cffFailure(gid);

        // Allow dynamic failure if the block is "loose" (next to a previously failed block)
        if (loose_elements.count(gid) > 0) add |= sim->dynamicFailure(gid, trigger_fault);

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
            //

            ///// Schultz:
            // Avoid the scenario that CFF is slightly changed due to other failed blocks,
            //  so that we never satisfy the condition if (!stress_drop). Must check this
            //  in the future.
            
            // calculate the drop in stress from the failure
            //stress_drop = sim->getCFF0(gid) - sim->getCFF(gid);
            //if (!stress_drop) stress_drop = sim->getStressDrop(gid) - sim->getCFF(gid);
            
            ////// Schultz
            stress_drop = sim->getStressDrop(gid) - sim->getCFF(gid);
            
            // Slip is in m
            slip = (stress_drop/sim->getSelfStresses(gid));

            ////// Schultz:
            // Relaxing the restriction on negative slip for now. We will investigate the cause/effect
            // of negative slips in the future.
            //if (slip < 0) slip = 0;

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
        // do we have the correct list of global failed elements?
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
    //sim->barrier();
    // yoder (note): after we distributeBlocks(), we can check to see that all items in local_ exist in global_ exactly once.
    // if not, throw an exception... and then we'll figure out how this is happening. remember, local_ is like [gid, gig, gid...]
    // global_ is like [(gid, p_rank), (gid, p_rank)...], and each pair item is accessed like global_[rw_num]->first /->second
    sim->distributeBlocks(local_secondary_id_list, global_secondary_id_list);
    //sim->barrier();   //yoder: (debug)

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

        //b[i] = sim->getCFF(*it)+sim->getFriction(*it)*sim->getRhogd(*it);
        ////// Schultz:
        // Need to incorporate stress drop here somehow, we need secondary
        // rupture physics to be same as primary. This incorporates the physics that
        // for static failure (CFF=0) the amount of stress to release is the stress drop.
        // For dynamic failure this amount of stress is reduced by however much stress is still
        // needed to reach static failure. Or if the element has accumulated a ton of stress (CFF>0)
        // then this adds that extra amount to the stress needed to be released.
        b[i] = sim->getStressDrop(*it) - sim->getCFF(*it);
    }

    //
    // heisenbug/heisen_hang likely candidate...
    //sim->barrier();     // yoder: (debug)
    //
    // so A,b are calculated for each local node (with dimension n_local x n_global and now they'll be consolidated on the root node. note that
    // the corresponding mpi_send comes after this block. the root node blocks until child nodes have sent A,b and root_node has received A,b.
    //
    /*
    #ifdef MPI_C_FOUND
        // yoder (debug):
        // (but this does not appear to be the problem; this exception is not being thrown, and we're still getting heisen_hang).
        // before we start our loops, let's see that we're all on the same event (and?) sweep.
        // use MPI_Allgather(*send, send_count, mpi_type_send, *rec, rec_count, mpi_type_rec, mpi_comm)
        // for now, check on all processes that each process thinks all other processes are on the same event, sweep.
        // use sim->getWorldSize() to get worldsize.
        int n_procs = sim->getWorldSize();
        int process_event_ids[n_procs];
        int process_sweep_ids[n_procs];
        int this_event_id = sim->getCurrentEvent().getEventNumber();
        int this_sweep_id = sweep_num;      // declared on the RunEvent class level. we don't really need to copy it.
        // to test on all processes:
        //MPI_Allgather(&this_event_id, 1, MPI_INT, &process_event_ids, 1, MPI_INT, MPI_COMM_WORLD);
        //MPI_Allgather(&this_sweep_id, 1, MPI_INT, &process_sweep_ids, 1, MPI_INT, MPI_COMM_WORLD);
        // to test on only root_node:
        MPI_Gather(&this_event_id, 1, MPI_INT, &process_event_ids, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&this_sweep_id, 1, MPI_INT, &process_sweep_ids, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //
        if (sim->getNodeRank()==0) {
            // and nominally, this test could also be if (sim->isRootNode() ), but really, any node can be "root" for the purposes of "gather", so we should be specific.
            // it might be slick to assign events and sweeps to different nodes, but hopefully we'll be removing this bit eventually.
            // now, see that all event,sweep values are the same:
            //printf("**Debug(%d/%d/%d):: event_ids: ", sim->getNodeRank(), getpid(), n_procs);
            for (int j=0; j<n_procs;++j) {
                //
                //printf("(%d/%d), ", process_event_ids[0], process_event_ids[j]);
                assertThrow(process_event_ids[j]==process_event_ids[0], "Processes do not have same EventID: (0: " << process_event_ids[0] << "), (" << j << ": " << process_event_ids[j] << "\n");
                };
                //printf("\n");
            for (int j=0; j<n_procs;++j) assertThrow(process_sweep_ids[j]==process_sweep_ids[0], "Processes do not have same SweepID: (0: " << process_sweep_ids[0] << "), (" << j << ": " << process_sweep_ids[j] << "\n");
        };
        // and everybody waits until we've evaluated this.
        //
    #endif
    */

    /////////////
    //
    if (sim->isRootNode()) {
        double *fullA = new double[num_global_failed*num_global_failed];
        double *fullb = new double[num_global_failed];
        double *fullx = new double[num_global_failed];

        // Fill in the A matrix and b vector from the various processes
        // note: for an empty global_id_list, this list will do nothing, so MPI_Recv() will not execute, and we'll probably end up with a hanging MPI_Send().
        //       can that ever happen? global_secondary_id_list is empty on one node but not on another? maybe between iterations?
        //sim->barrier();       //yoder: (debug)        note: this (and the other one too... sim->barrier() is paried with an ifRoot==False set.
        //
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
                //printf("**Debug[%d/%d]: A[],b MPI_Recv root-node blocking...\n", sim->getNodeRank(), getpid());
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
        //
        // root node regroups (barrier() ), then does its sending bit:
        //sim->barrier();       //yoder: (debug)

        /*
                // Debug:
        #ifdef MPI_C_FOUND
                // how many elements do we expect to send back to processes?
                //int id_counts[sim->getWorldSize()];
                int id_counts[n_procs];
                for (int j=0; j<n_procs; ++j) id_counts[j]=0;   // noting that we initilize our own count,j=0, which we won't use.
                for (jt=global_secondary_id_list.begin(); jt!=global_secondary_id_list.end(); ++jt) ++id_counts[jt->second];
                for (int w=1; w<sim->getWorldSize();++w) {
                    MPI_Send(&(id_counts[w]), 1, MPI_INT, w, 0, MPI_COMM_WORLD);
                    }
        #endif
                // END DEBUG
        */
        // Send back the resulting values from x to each process
        //for (i=0,n=0,jt=global_id_list.begin(); jt!=global_id_list.end(); ++jt,++i) {
        for (i=0,n=0,jt=global_secondary_id_list.begin(); jt!=global_secondary_id_list.end(); ++jt,++i) {
            if (jt->second != sim->getNodeRank()) {
#ifdef MPI_C_FOUND
                // send these values to node-rank jt->second:
                // yoder: try using synchronous send, MPI_Ssend()
                //MPI_Send(&(fullx[i]), 1, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD);
                MPI_Ssend(&(fullx[i]), 1, MPI_DOUBLE, jt->second, 0, MPI_COMM_WORLD);
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
        // ... and note that if there are no local failures, these MPI_Send() commands will not execute. do we end up with a hanging send/receive pair?
        // Child Nodes do their sending bit:
        //sim->barrier();       //yoder: (debug)        note: this (and the other one too... sim->barrier() is paried with an ifRoot==True set.
        for (i=0; i<num_local_failed; ++i) {
            //printf("**Debug[%d/%d]: A[],b MPI_Send child-node blocking...\n", sim->getNodeRank(), getpid());
            // yoder: try using synchronous MPI_Ssend()
            //MPI_Send(&(A[i*num_global_failed]), num_global_failed, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            //MPI_Send(&(b[i]), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Ssend(&(A[i*num_global_failed]), num_global_failed, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Ssend(&(b[i]), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        //
        // regroup (barrier()), then child nodes do their receiving bit:
        //sim->barrier();       //yoder: (debug)
        // fetch x[i] from root (rank 0) node:
        // yoder (debugging):
        // check to see that the local ids exist in the global ids:

        /*
        // Debugging:
        //bool loc_glob_ok = true;
        // start with the count. how many entries does the root_node expect to send?
        int expected_recv_count = 0;
        MPI_Recv(&expected_recv_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (expected_recv_count!=num_local_failed) {
            printf("**ERROR(%d/%d): expected recv_count does not match root node: %d/%d\n", sim->getNodeRank(), getpid(), expected_recv_count, num_local_failed);
            assertThrow(0, "send/receive count does not match in secondary blocks.\n");
            };
        ////
        */

        //
        for (i=0; i<num_local_failed; ++i) {
            // yoder: (debug) this is, apparently, where heisen_hang happens.
            //    has apparently moved on; it looks like having finished the secondary loop.
            //    so it seems likely that the problem is that one or more secondary blocks ends up on multiple processors, but is listed only once in the globals (??)
            //    so let's try to verify that our list of local_failed blocks is compatible with the root list of gloabal_failed. we might, yet, end up just
            //    rewriting this section with MPI_Gatherv()...
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
        //
        double slip = x[i] - sim->getSlipDeficit(*it);

        //
        ////// Schultz:
        // We may be destabilizing the system here if the solution includes negative slipped elements.
        // We cannot solve the whole system then throw out a few elements.
        
        ////// Schultz:
        // Relaxing the restriction on negative slip for now. We will investigate the cause/effect
        // of negative slips in the future.
        //if (slip > 0) {
        // Record how much the block slipped in this sweep and initial stresse
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
        //}
    }

    //
    // delete/de-allocate arrays (use "delete []" for arrays, as opposed to "delete" for single objects)
    delete [] A;
    delete [] b;
    delete [] x;
    //sim->barrier();
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
    //BlockID                 triggerID, gid;
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
    //
    // Clear the list of "loose" (can dynamically fail) blocks
    loose_elements.clear();
    //
    // Clear the list of failed blocks, and add the trigger block
    local_failed_elements.clear();

    //
    if (sim->getCurrentEvent().getEventTriggerOnThisNode()) {
        local_failed_elements.insert(triggerID);
        loose_elements.insert(triggerID);
        sim->setFailed(triggerID, true);
    }

    //
    // yoder (note): Comm::blocksToFail() executes a single MPI_Allreduce() (when MPI is present),
    // ... but since it's only the one "all-hands" call, it is not a likely candidate for heisen_hang.
    // also note that blocksToFail() involves some bool/int conversions that might warrant further inspection and clarification.
    more_blocks_to_fail = sim->blocksToFail(!local_failed_elements.empty());
    //
    // use this iterator/counter to efficiently walk through the event_sweeps list when we update stresses:
    unsigned int event_sweeps_pos = 0;

    //
    // While there are still failed blocks to handle
    //sim->barrier();
    while (more_blocks_to_fail || final_sweep) {
        // write stress, slip, etc. to events and sweeps output (text or hdf5).
        sim->output_stress(sim->getCurrentEvent().getEventNumber(), sweep_num);

        // Share the failed blocks with other processors to correctly handle
        // faults that are split among different processors
        //sim->barrier();    // yoder: (debug)   (we're probably safe without this barrier() )... but at some point, i was able to generate a hang during distributeBlocks()
        // so let's try it with this in place...
        sim->distributeBlocks(local_failed_elements, global_failed_elements);
        //sim->barrier(); // yoder: (debug)
        //
        // Process the blocks that failed.
        // note: setInitStresses() called in processBlocksOrigFail().
        // note: processBlocksOrigFail() is entirely local (no MPI).
        processBlocksOrigFail(sim, event_sweeps);

        // Recalculate CFF for all blocks where slipped blocks don't contribute
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();
            sim->setShearStress(gid, 0.0);
            sim->setNormalStress(gid, sim->getRhogd(gid));
            //sim->setUpdateField(gid, (sim->getFailed(gid) ? 0 : std::isnan(sim->getSlipDeficit(gid)) ? 0 :sim->getSlipDeficit(gid) )); // ... also check for nan values
            sim->setUpdateField(gid, (sim->getFailed(gid) ? 0 : sim->getSlipDeficit(gid)));

            ////////// Schultz:
            // We need to ensure our slip economics books are balanced. I suspect we need here
            // instead: sim->setUpdateField(gid, sim->getSlipDeficit(gid) ). Update the stresses using
            // the current slip of all elements, or else we throw away the slip information from primary ruptures.
            //sim->setUpdateField(gid, sim->getSlipDeficit(gid));
        }

        // Distribute the update field values to other processors
        sim->distributeUpdateField();

        // Set dynamic triggering on for any blocks neighboring blocks that slipped in the last sweep
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();

            // Add block neighbors if the block has slipped
            if (sim->getFailed(gid)) {
                nbr_start_end = sim->getNeighbors(gid);

                for (nit=nbr_start_end.first; nit!=nbr_start_end.second; ++nit) {
                    loose_elements.insert(*nit);
                }
            }
        }

        //
        // Calculate the CFFs based on the stuck blocks
        // multiply greenSchear() x getUpdateFieldPtr() --> getShearStressPtr() ... right?
        // assign stress values (shear stresses at this stage are all set to 0; normal stresses are set to (i think) sim->getRhogd(gid) -- see code a couple paragraphs above.
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
        // For each block that has failed already, calculate the slip from the movement of blocks that just failed
        //
        processBlocksSecondaryFailures(sim, event_sweeps);
        //

        // Set the update field to the slip of all blocks
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();
            sim->setShearStress(gid, 0.0);
            sim->setNormalStress(gid, sim->getRhogd(gid));
            //
            // if this is a current/original failure, then 0 else...
            sim->setUpdateField(gid, (global_failed_elements.count(gid)>0 ? 0 : sim->getSlipDeficit(gid)));
            //sim->setUpdateField(gid, (global_failed_elements.count(gid)>0 ? 0 : std::isnan(sim->getSlipDeficit(gid)) ? 0 : sim->getSlipDeficit(gid)));
            
            //////// Schultz: try setting updatefield 0 for newly failed elements this sweep
            //sim->setUpdateField(gid, ((sim->getFailed(gid) && global_failed_elements.count(gid) == 0) ? 0 : sim->getSlipDeficit(gid)));            
            // We need to ensure our slip economics books are balanced. I suspect we need here
            // instead: sim->setUpdateField(gid, sim->getSlipDeficit(gid) ). Update the stresses using
            // the current slip of all elements, or else we throw away the slip information from failed elements.
            //sim->setUpdateField(gid, sim->getSlipDeficit(gid));
        }

        //
        //
        // could we be getting MPI conflicts here? if one process is trying to distribute original failures (using update_field), and another process is already
        // trying to distribute secondary failures, this could (???) cause a conflict and hang???
        // let's try an MPI_Barrier here (i think we've got this wrapped up in an ifmpi block somewhere...):
        //
        //sim->barrier();    // yoder: (debug)
        sim->distributeUpdateField();
        //sim->barrier();    // yoder: (debug)
        //
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
        //
        // Record the final stresses of blocks that failed during this sweep
        BlockIDProcMapping::iterator        fit;

        //
        // yoder: there are a  couple ways to see that we loop over all the necessary elements. here's an approach that starts with the original code.
        // i recommend leaving these comments in place for a revision or two, just until we make a final decision about how we're going to do this
        // (process these two populations of block failures).
        //
        // this is the original code, to loop over all new failures. nominally, we want to do this plus loop over secondary failures as well.
        // this loop ignores all secondary slip events. do we loop over a different list, or do we need to maintain a collective global_failed_element list (aka,
        // not re-initialize it every sweep)? ... and of course, we must ask: is this by design? can we loop over event_sweeps (or at least new entries therein) directly?
        /*
        for (fit=global_failed_elements.begin(); fit!=global_failed_elements.end(); ++fit) {
            if (sim->isLocalBlockID(fit->first)) {
                printf("**Debug: setFinalStresses(): sweep: %d, gid: %d, ss: %f, ns: %f\n", sweep_num, fit->first, sim->getShearStress(fit->first), sim->getNormalStress(fit->first));
                event_sweeps.setFinalStresses(sweep_num,
                                              fit->first,
                                              sim->getShearStress(fit->first),
                                              sim->getNormalStress(fit->first));
            }
        }
        // yoder: now, this loops over a global list of secondary failures (which has/d been moved to a class-wide declaration, if this is how we want to handle this).
        for (fit=global_secondary_id_list.begin(); fit!=global_secondary_id_list.end(); ++fit) {
            if (sim->isLocalBlockID(fit->first)) {
                printf("**Debug: setFinalStresses_secondary(): sweep: %d, gid: %d, ss: %f, ns: %f\n", sweep_num, fit->first, sim->getShearStress(fit->first), sim->getNormalStress(fit->first));
                //if (not fit->_slip>0) {continue;};        // note: i think nan will always throw a false, so this should work for both nan and <=0.
                event_sweeps.setFinalStresses(sweep_num,
                                              fit->first,
                                              sim->getShearStress(fit->first),
                                              sim->getNormalStress(fit->first));
            }
        }
        */
        //
        //
        // and this loop updates final_stress values by looping directly over the current sweeps list.
        //  note that event_sweeps is of type quakelib::ModelSweeps, which contains a vector<SweepData> _sweeps.
        for (quakelib::ModelSweeps::iterator s_it=event_sweeps.begin(); s_it!=event_sweeps.end(); ++s_it, ++event_sweeps_pos) {
            //
            // yoder: as per request by KS, change isnan() --> std::isnan(); isnan() appears to throw an error on some platforms.
            if (std::isnan(s_it->_shear_final) and std::isnan(s_it->_normal_final)) {
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
        //
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

    //
    // output_stress() for final item in list.
    sim->output_stress(sim->getCurrentEvent().getEventNumber(), sweep_num);

    // Set the completed list as the sweep list for the entire event
    sim->collectEventSweep(event_sweeps);
    sim->getCurrentEvent().setSweeps(event_sweeps);
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
    EventAftershock                           as;
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
        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
            double as_to_elem_dist = sim->getBlock(gid).center().dist(as.loc());
            as_elem_dists.insert(std::make_pair(as_to_elem_dist, gid));
        }

        // Determine the target rupture area given the aftershock magnitude
        // TODO:
        // user_defined_constants (flag this for later revisions in which we move these contant definitions to a parameters file).
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

    // If there's a specific block that triggered the event, it's a static stress failure type event
    //sim->barrier();
    if (sim->getCurrentEvent().getEventTrigger() != UNDEFINED_ELEMENT_ID) {
        processStaticFailure(sim);
    } else {
        // Otherwise it's an aftershock
        processAftershock(sim);
    }

    //sim->barrier();

    // Record the stress in the system before and after the event.
    // yoder: note that recordEventStresses() emloyes a bit of MPI action, so it might be advisable to
    // add some barrier() blocking (which might have been added above to barrier()-wrap the process_{earthquake type}() call).
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
