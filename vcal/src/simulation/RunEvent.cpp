// Copyright (c) 2010, John B. Rundle <rundle@cse.ucdavis.edu>, 
// All rights reserved.
// 
// Redistribution and use of this code or any derivative works are
// permitted provided that the following conditions are met:
// 
// * Redistributions may not be sold, nor may they be used in a
// commercial product or activity.
// 
// * Redistributions that are modified from the original source must
// include the complete source code, including the source code for all
// components used by a binary built from the modified
// sources. However, as a special exception, the source code
// distributed need not include anything that is normally distributed
// (in either source or binary form) with the major components
// (compiler, kernel, and so on) of the operating system on which the
// executable runs, unless that component itself accompanies the
// executable.
// 
// * Redistributions must reproduce the above copyright notice, this list
// of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "RunEvent.h"
#include "SimError.h"

/*!
 At the end of each sweep after we have recalculated block CFF, we determine
 which blocks will have a failure due to dynamic or static stress changes.
 To avoid runaway failures we limit the number of failures for a given block in
 a sweep to 20. However, it is unclear if this ever happens in practice.
 */
void RunEvent::markBlocks2Fail(VCSimulation *sim, const FaultID &trigger_fault, VCEventSweep &current_sweep) {
	int			lid;
	BlockID		gid;
	bool		add;
	
    for(lid=0;lid<sim->numLocalBlocks();++lid) {
		gid = sim->getGlobalBID(lid);
        //current_sweep.setFinalStresses(gid, sim->getBlock(gid).getShearStress(), sim->getBlock(gid).getNormalStress());
		// Add this block if it has a static CFF failure or dynamic failure
		add = sim->getBlock(gid).cffFailure() || sim->getBlock(gid).dynamicFailure(trigger_fault) || sim->getBlock(gid).additionalFailure();
		//sim->console() << gid << 
        //" failure s:" << sim->getBlock(gid).cffFailure() << " d:" << sim->getBlock(gid).dynamicFailure(trigger_fault) << " a:" << sim->getBlock(gid).additionalFailure() << ", ";
		// Let each block fail at most 10 times. This is somewhat arbitrary, but it prevents runaway ruptures.
        if(add && num_failures[gid] < 10) {
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
	BlockIDSet::iterator		fit;
	double						rn, slip, stress_drop, a;
	FaultBlockMapping			failed_faults;
    
    sim->getFaultBlockMapping(failed_faults, blocks2fail);
    
	// For each block that has failed
	for(fit=blocks2fail.begin();fit!=blocks2fail.end();++fit) {
		if (sim->isLocalBlockID(*fit)) {
			Block& b = sim->getBlock(*fit);
			
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
	BlockIDSet::iterator		fit;
	double						rn, slip, stress_drop, a, failure_size;
	FaultBlockMapping			failed_faults;
	FaultFailureAreaMapping		failed_faults_area;
	double						slip_scaling_threshold;
    
    sim->getFaultBlockMapping(failed_faults, blocks2fail);
	
	sim->getFaultFailureAreaMapping(failed_faults_area, blocks2fail);
    
	// For each block that has failed
	for(fit=blocks2fail.begin();fit!=blocks2fail.end();++fit) {
		if (sim->isLocalBlockID(*fit)) {
			Block& b = sim->getBlock(*fit);
			
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
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	BlockList::iterator		it;
	VCEvent::iterator		eit;
	VCEventSweep			current_sweep;
	BlockID					triggerID, gid;
	FaultID					trigger_fault;
	unsigned int			i;
	int						lid, more_blocks_to_fail, num_spec_exec, start_sweep_size;
	SpecExecStage			stage;
	bool					local_blocks_failed;
	BlockVal				num_local_sweeps, max_sweeps;
	EventSweeps				sweep_list;
	std::pair<BlockIDSet::const_iterator, BlockIDSet::const_iterator> nbr_start_end;
	BlockIDSet::const_iterator	nit;
    
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
    for(lid=0;lid<sim->numLocalBlocks();++lid) sim->getBlock(sim->getGlobalBID(lid)).saveStresses();
	
	// Initially we assume the failure can affect all nodes
	stage = NORMAL_OPERATION;
	sim->startEvent();
	
	more_blocks_to_fail = sim->blocksToFail(!blocks2fail.empty());
	
	// While there are still failed blocks to handle or we're in localized failure mode
    while(stage == LOCALIZED_FAILURE || more_blocks_to_fail) {
		// Given the blocks that we need to process, predict whether this will be a localized failure
		switch (stage) {
			case NORMAL_OPERATION:
				// If this is the only node with failing blocks and we predict the rupture
				// will remain local, change to LOCALIZED_FAILURE mode
				if (more_blocks_to_fail == 1 && sim->isLocalizedFailure(blocks2fail)) {
					stage = LOCALIZED_FAILURE;
					start_sweep_size = sweep_list.size();
				}
				break;
			case LOCALIZED_FAILURE:
				// If we are in LOCALIZED_FAILURE mode and we continue to predict it will remain that way, do nothing
				// Otherwise, we move back to NORMAL_OPERATION by checking the results
				if (!sim->isLocalizedFailure(blocks2fail)) {
					stage = CHECK_SELF_IGNORE;
				}
				// If we are in LOCALIZED_FAILURE mode and we have lots of localized sweeps, recalculate stress to avoid
				// stress buildups which cause rewinding
				if (sweep_list.size() - start_sweep_size > 50) {
					stage = CHECK_SELF_IGNORE;
				}
				break;
			default:
				assertThrow(false, "Error occurred in speculative execution.");
				break;
		}

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
		// Original method
		for(gid=0;gid<sim->numGlobalBlocks();++gid) {
			sim->setUpdateField(gid, current_sweep.getBlockSlip(gid));
		}
		
		// Set the update field for stress recalculation with
		// the total amount each block has slipped in this sweep
		// Working with the slip deficit
		//for(gid=0;gid<sim->numGlobalBlocks();++gid) {
		//	sim->setFShearStress(gid, 0.0);
		//	sim->setFNormalStress(gid, sim->getBlock(gid).getRhogd());
			//sim->setUpdateField(it->getBlockID(), it->state.slipDeficit);
		//	sim->setUpdateField(gid, sim->getBlock(gid).state.slipDeficit);
		//}

		
		// Share the update field with all other nodes (needed for stress recalculations)
		if (stage != LOCALIZED_FAILURE) {
			// If we were previously in localized mode, sum the slips that occurred in this mode
			// and add them to the update field for transmission to other nodes
			if (stage == CHECK_SELF_IGNORE) {
				for (i=start_sweep_size;i<sweep_list.size();++i) {
					for(lid=0;lid<sim->numLocalBlocks();++lid) {
						gid = sim->getGlobalBID(lid);
						sim->setUpdateField(gid, sim->getUpdateField(gid)+sweep_list[i].getBlockSlip(gid));
					}
				}
			}
			
			// Distribute the update field values to other processors
			num_spec_exec = sim->distributeUpdateField(stage == CHECK_SELF_IGNORE);
			
			// If more than one node did speculative execution, treat it as a failure (though this should never happen)
			assertThrow(num_spec_exec <= 1, "Only one node should go into localized mode at a time.");
			
			// If only one node did speculative execution and we're not that node, go into check mode
			if (num_spec_exec == 1 && stage == NORMAL_OPERATION) {
				stage = CHECK_IF_SELF_FAILED;
			}
			
			// After transmission, restore the previous values of the update field for the node that ran localized
			if (stage == CHECK_SELF_IGNORE) {
				for(gid=0;gid<sim->numGlobalBlocks();++gid) {
					sim->setUpdateField(gid, current_sweep.getBlockSlip(gid));
				}
			}
		}
		
		// Set dynamic triggering on for any blocks neighboring blocks that slipped in the last sweep
		for(gid=0;gid<sim->numGlobalBlocks();++gid) {
			if (sim->getUpdateField(gid) > 0) {
				nbr_start_end = sim->getNeighbors(gid);
				
				for (nit=nbr_start_end.first;nit!=nbr_start_end.second;++nit) {
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
		BlockIDSet::iterator		fit;
		for(fit=blocks2fail.begin();fit!=blocks2fail.end();++fit) {
			if (sim->isLocalBlockID(*fit)) {
				Block& b = sim->getBlock(*fit);
				current_sweep.setFinalStresses(*fit, b.getShearStress(), b.getNormalStress());
			}
		}
		blocks2fail.clear(); // we are done with these blocks
		
		// Find any blocks that fail because of the new stresses
        markBlocks2Fail(sim, trigger_fault, current_sweep);
		
		// If we were doing localized failure before, make sure that no block elsewhere failed
		max_sweeps.val = 0;
		if (stage == CHECK_SELF_IGNORE) {
			if (sim->blocksToFail(false)) {
				stage = REWIND_ALL;
				fail_sweep_sizes.push_back(sweep_list.size()-start_sweep_size);
			} else {
				succ_sweep_sizes.push_back(sweep_list.size()-start_sweep_size);
				sim->speculationSuccess(current_sweep);
				stage = NORMAL_OPERATION;
				num_local_sweeps.val = sweep_list.size() - start_sweep_size;
				num_local_sweeps.block_id = UNDEFINED_BLOCK_ID;
				sim->allReduceBlockVal(num_local_sweeps, max_sweeps, BLOCK_VAL_MAX);
			}
		} else if (stage == CHECK_IF_SELF_FAILED) {
			local_blocks_failed = !blocks2fail.empty();
			if (sim->blocksToFail(local_blocks_failed)) {
				stage = REWIND_ALL;
			} else {
				sim->speculationSuccess(current_sweep);
				stage = NORMAL_OPERATION;
				num_local_sweeps.val = 1;
				num_local_sweeps.block_id = UNDEFINED_BLOCK_ID;
				sim->allReduceBlockVal(num_local_sweeps, max_sweeps, BLOCK_VAL_MAX);
				for (i=0;i<max_sweeps.val;++i) sweep_list.push_back(VCEventSweep());
			}
		}
		
		// Collect the sweep history from all nodes onto the root
		if (stage != LOCALIZED_FAILURE && stage != REWIND_ALL) {
			for (i=max_sweeps.val;i>0;--i) sim->collectEventSweep(*(sweep_list.end()-i));
			
			sim->collectEventSweep(current_sweep);
		}
		
		// Add the recorded sweep to the list
		sweep_list.push_back(current_sweep);
		
		// If we need to rewind, then reset everything to the values at the event start
		if (stage == REWIND_ALL) {
			// Note that the current speculation failed
			sim->speculationFailed(current_sweep);
			
			// Clear the sweeps and block fail lists
			sweep_list.clear();
			blocks2fail.clear();
			
			// Return to the initial state
			if (sim->getCurrentEvent().getEventTriggerOnThisNode()) {
                blocks2fail.insert(triggerID);
            }
			for(lid=0;lid<sim->numLocalBlocks();++lid) sim->getBlock(sim->getGlobalBID(lid)).restoreStresses();
			stage = NORMAL_OPERATION;
		}
		
		if (stage != LOCALIZED_FAILURE) more_blocks_to_fail = sim->blocksToFail(!blocks2fail.empty());
    }
	
	// Set the completed list as the sweep list for the entire event
	sim->getCurrentEvent().addSweeps(sweep_list);
	
    // Record the stress in the system before and after the event
    recordEventStresses(sim);
    
	// Update the cumulative slip for this fault and turn off dynamic triggering
	for(lid=0;lid<sim->numLocalBlocks();++lid) {
		gid = sim->getGlobalBID(lid);
        sim->getBlock(gid).state.slipCumulative += sim->getCurrentEvent().getEventSlip(gid);
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

void RunEvent::finish(SimFramework *_sim) {
	/*std::vector<int>::iterator		it;
	int		i;
	
	for (i=0;i<_sim->getWorldSize();++i) {
		if (i == _sim->getNodeRank()) {
			std::cerr << "Node " << i << std::endl;
			std::cerr << "Success: ";
			for (it=succ_sweep_sizes.begin();it!=succ_sweep_sizes.end();++it) {
				std::cerr << *it << " ";
			}
			std::cerr << std::endl << "Failure: ";
			for (it=fail_sweep_sizes.begin();it!=fail_sweep_sizes.end();++it) {
				std::cerr << *it << " ";
			}
			std::cerr << std::endl;
		}
		_sim->barrier();
	}*/
}
