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

#include "CheckpointFileOutput.h"
#include <sstream>

void CheckpointFileOutput::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# Saving checkpoint every " << sim->getCheckpointPeriod() << " events." << std::endl;
}

void CheckpointFileOutput::writeCheckpoint(const std::string &ckpt_file_name, const VCSimulation *sim) {
	CheckpointSet			checkpoint_set;
	int						i;
	BlockID					bid;
	
	// Get current state for local blocks
	for (i=0;i<sim->numLocalBlocks();++i) {
		bid = sim->getGlobalBID(i);
		checkpoint_set[bid] = sim->getBlock(bid).state.readCheckpointData();
	}
	HDF5CheckpointWriter checkpoint_file(ckpt_file_name,
										 sim->numGlobalBlocks(),
										 sim->getYear(),
										 sim->getEventCount(),
										 checkpoint_set);
}

SimRequest CheckpointFileOutput::run(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	std::stringstream		ss;
	
	// Periodically checkpoint simulation state to a file
	if (sim->getCheckpointPeriod() > 0 && sim->getEventCount() % sim->getCheckpointPeriod() == 0) {
		ss << sim->getCheckpointPrefix() << checkpoint_num << ".h5";
		writeCheckpoint(ss.str(), sim);
		checkpoint_num++;
	}
	
	return SIM_STOP_OK;
}

void CheckpointFileOutput::finish(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	std::stringstream		ss;
	
	ss << sim->getCheckpointPrefix() << "final" << ".h5";
	writeCheckpoint(ss.str(), sim);
}
