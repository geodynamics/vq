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

#include "HDF5DataShare.h"
#include "HDF5Data.h"

bool HDF5DataShare::pauseFileExists(void) {
	FILE			*fp;
	if ((fp = fopen(HDF5_PAUSE_FILE_NAME, "r"))) {
		fclose(fp);
		return true;
	}
	return false;
}

void HDF5DataShare::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# To access the HDF5 file during the simulation, pause it " << std::endl;
	sim->console() << "# by creating the file " << HDF5_PAUSE_FILE_NAME ".  Delete the file to resume." << std::endl;
}

void HDF5DataShare::init(SimFramework *_sim) {
	VCSimulation				*sim = static_cast<VCSimulation*>(_sim);
	BlockList::const_iterator	it;
	
	h5_data = NULL;
	next_pause_check = sim->itersPerSecond();
	
	// Only the root node writes to the output file
	if (sim->isRootNode()) {
		h5_data = new HDF5DataWriter(sim->getHDF5File(), sim->numGlobalBlocks());
		
		// Set the start and end years of the simulation
		h5_data->setStartEndYears(sim->getYear(), sim->getSimDuration());
		
		// Set the base latitude/longitude of the simulation
		h5_data->setLatLon0(sim->getBaseLatLon());
		
		// Record block information
		for (it=sim->begin();it!=sim->end();++it) h5_data->setBlockInfo(*it);
	}
}

/*!
 During the simulation, this framework just writes events to the HDF5 file.
 */
SimRequest HDF5DataShare::run(SimFramework *_sim) {
	VCSimulation		*sim = static_cast<VCSimulation*>(_sim);
	
	if (h5_data) h5_data->writeEvent(sim->getCurrentEvent(), sim->getBGEvents());
	
	// Check if the pause file exists each second
	if (sim->getEventCount() >= next_pause_check) {
		next_pause_check += sim->itersPerSecond();
		if (sim->isRootNode() && pauseFileExists()) {
			// Flush out the HDF5 data
			h5_data->flush();
			sim->console() << "# Pausing simulation due to presence of file " << HDF5_PAUSE_FILE_NAME << std::endl;
			while (pauseFileExists()) sleep(1);
		}
		// Other processes wait until the pause file is gone
		sim->barrier();
	}
	
	return SIM_STOP_OK;
}

/*!
 Finish by unattaching from the shared memory structure.
 */
void HDF5DataShare::finish(SimFramework *_sim) {
	if (h5_data) delete h5_data;
}
