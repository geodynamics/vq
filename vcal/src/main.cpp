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

#include "VCSimulation.h"

// Plugin #includes
#include "HDF5DataShare.h"

// File parsing and output related
#include "EqSimFileOutput.h"
#include "EqSimFileParse.h"
#include "GreensFileOutput.h"
#include "CheckpointFileOutput.h"
#include "CheckpointFileParse.h"
#include "SystemFileOutput.h"

// Simulation related
#include "AddAsperities.h"
#include "AddNoise.h"
#include "BadFaultKill.h"
#include "BASSAftershocks.h"
#include "BGPoissonEvents.h"
#include "BlockValCompute.h"
#include "DepthDepVelocity.h"
#include "EventRecorder.h"
#include "GracefulQuit.h"
#include "GreensInit.h"
#include "GreensKillInteraction.h"
#include "ProgressMonitor.h"
#include "RunEvent.h"
#include "SanityCheck.h"
#include "UpdateBlockStress.h"
#include "VCInitBlocks.h"

int main (int argc, char **argv) {
	PluginID		read_eqsim_file, init_blocks, block_val_compute;
	PluginID		greens_init, greens_outfile, system_output_file, add_noise, add_asperities, bad_fault_kill;
	PluginID		depth_dependent_velocity, greens_kill, update_block_stress, run_event;
	PluginID		sanity_checking, bass_model_aftershocks, bg_poisson_events, display_progress, state_output_file;
	PluginID		eqsim_output_file, h5_data_share, vc_events_output_file, graceful_quit;
	VCSimulation	*vc_sim;
	
	vc_sim = new VCSimulation(argc, argv);
	
	// ************************************************************
	// ** Define plugins and whether they are active
	// ************************************************************
	// EqSim files are parsed if a geometry file name is specified
	read_eqsim_file = vc_sim->registerPlugin(new EqSimFileParse, !vc_sim->getEqSimGeometryFile().empty());
	
	// Calculate block values if we aren't using EqSim files
	block_val_compute = vc_sim->registerPlugin(new BlockValCompute, vc_sim->getEqSimGeometryFile().empty());
	
	// Calculate Greens function if a calculation method is specified
	greens_init = vc_sim->registerPlugin(new GreensInit, vc_sim->getGreensCalcMethod() != GREENS_CALC_NONE);
	
	// Write the Greens values to a file if a file name is specified
	greens_outfile = vc_sim->registerPlugin(new GreensFileOutput, !vc_sim->getGreensOutfile().empty());
	
	// Write the system information to a file if a file name is specified
	system_output_file = vc_sim->registerPlugin(new SystemFileOutput, !vc_sim->getSystemOutfile().empty());
	
	// Add stress noise if the stress noise level is above 0
	add_noise = vc_sim->registerPlugin(new AddNoise, vc_sim->getStressNoise() > 0);
	
	// Add asperities if there are more than 0 asperities requested
	add_asperities = vc_sim->registerPlugin(new AddAsperities, vc_sim->getNumAsperities() > 0);
	
	// Kill faults that drop below a certain CFF value
	bad_fault_kill = vc_sim->registerPlugin(new BadFaultKill, vc_sim->getFaultKillCFF() < 0);
	
	// Use depth dependent velocity if requested
	depth_dependent_velocity = vc_sim->registerPlugin(new DepthDepVelocity, vc_sim->doDepthDependentVelocity());
	
	// Implement Greens matrix killing if the kill distance is greater than 0
	greens_kill = vc_sim->registerPlugin(new GreensKillInteraction, vc_sim->getGreensKillDistance() > 0);
	
	// Implement sanity checking if requested
	sanity_checking = vc_sim->registerPlugin(new SanityCheck, vc_sim->doSanityCheck());
	
	// Implement BASS aftershocks if the number of BASS generations is more than 0
	bass_model_aftershocks = vc_sim->registerPlugin(new BASSAftershocks, vc_sim->getBASSMaxGenerations() != 0);
	
	// Implement background Poisson process events if the interevent time is more than 0
	bg_poisson_events = vc_sim->registerPlugin(new BGPoissonEvents, vc_sim->getBGEventMeanInterevent() > 0);
	
	// Display progress if the progress display period is more than 0
	display_progress = vc_sim->registerPlugin(new ProgressMonitor, vc_sim->getProgressPeriod() > 0);
	
	// Write the simulation state if the period is more than zero
	state_output_file = vc_sim->registerPlugin(new CheckpointFileOutput, vc_sim->getCheckpointPeriod() > 0);
	
	// Write an EqSim event file if the output file is specified
	eqsim_output_file = vc_sim->registerPlugin(new EqSimFileOutput, !vc_sim->getEqSimOutputFile().empty());
	
	// Write events in the old output format if the file name is specified
	vc_events_output_file = vc_sim->registerPlugin(new EventRecorder, !vc_sim->getEventsFile().empty());
	
	// These plugins are always active in a simulation
	init_blocks = vc_sim->registerPlugin(new VCInitBlocks, true);
	update_block_stress = vc_sim->registerPlugin(new UpdateBlockStress, true);
	run_event = vc_sim->registerPlugin(new RunEvent, true);
	h5_data_share = vc_sim->registerPlugin(new HDF5DataShare, true);
	graceful_quit = vc_sim->registerPlugin(new GracefulQuit, true);
	
	// ************************************************************
	// ** Specify plugin dependencies
	// ************************************************************
	// Block initialization must occur after the blocks are read in
	vc_sim->registerDependence(init_blocks, read_eqsim_file, DEP_OPTIONAL);

	// Check for quit file before doing any intensive calculation
	vc_sim->registerDependence(graceful_quit, init_blocks, DEP_OPTIONAL);
	
	// Greens calculation occurs first, then initial Greens/system value output
	vc_sim->registerDependence(greens_init, graceful_quit, DEP_OPTIONAL);
	
	// Noise/asperities/other block modifications occur after block initialization
	vc_sim->registerDependence(block_val_compute, greens_init, DEP_OPTIONAL);
	vc_sim->registerDependence(greens_outfile, greens_init, DEP_REQUIRE);
	vc_sim->registerDependence(greens_kill, greens_init, DEP_REQUIRE);
	vc_sim->registerDependence(greens_outfile, greens_kill, DEP_OPTIONAL);
	vc_sim->registerDependence(depth_dependent_velocity, block_val_compute, DEP_OPTIONAL);
	vc_sim->registerDependence(add_noise, depth_dependent_velocity, DEP_OPTIONAL);
	vc_sim->registerDependence(add_asperities, depth_dependent_velocity, DEP_OPTIONAL);
	
	// Output the system file after we have applied noise/asperities
	vc_sim->registerDependence(system_output_file, add_asperities, DEP_OPTIONAL);
	vc_sim->registerDependence(system_output_file, add_noise, DEP_OPTIONAL);
	
	// Just to be safe we allow either noise or asperities but not both
	vc_sim->registerDependence(add_asperities, add_noise, DEP_EXCLUSIVE);
	
	// The core of the simulation, which must occur after Greens function is calculated and after asperities/noise
	vc_sim->registerDependence(update_block_stress, system_output_file, DEP_OPTIONAL);
	vc_sim->registerDependence(update_block_stress, block_val_compute, DEP_OPTIONAL);
	vc_sim->registerDependence(update_block_stress, greens_kill, DEP_OPTIONAL);
	vc_sim->registerDependence(run_event, update_block_stress, DEP_REQUIRE);
	
	// Background events and BASS aftershocks are added after the initial event is processed
	vc_sim->registerDependence(bg_poisson_events, run_event, DEP_REQUIRE);
	vc_sim->registerDependence(bass_model_aftershocks, run_event, DEP_REQUIRE);
	
	// BASS aftershocks after background events
	vc_sim->registerDependence(bass_model_aftershocks, bg_poisson_events, DEP_OPTIONAL);
	
	// Sanity checking occurs after the events are run
	vc_sim->registerDependence(sanity_checking, run_event, DEP_OPTIONAL);
	vc_sim->registerDependence(bad_fault_kill, run_event, DEP_OPTIONAL);
	
	// Output occurs after events (including background) are finished
	vc_sim->registerDependence(display_progress, run_event, DEP_OPTIONAL);
	vc_sim->registerDependence(state_output_file, run_event, DEP_OPTIONAL);
	vc_sim->registerDependence(eqsim_output_file, bg_poisson_events, DEP_OPTIONAL);
	vc_sim->registerDependence(h5_data_share, bg_poisson_events, DEP_OPTIONAL);
	vc_sim->registerDependence(vc_events_output_file, bg_poisson_events, DEP_OPTIONAL);
	
	//vc_sim->writeDOT(vc_sim->errConsole());
	
	// ************************************************************
	// ** Initialize and run the simulation
	// ************************************************************
	vc_sim->init();
	vc_sim->run();
	vc_sim->finish();
	
	vc_sim->printTimers();
	
	delete vc_sim;
	
    return 0;
}
