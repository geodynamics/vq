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

#include "EventRecorder.h"
#include "SimError.h"
#include <limits.h>
#include <iomanip>

void EventRecorder::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# Events output file: " << sim->getEventsFile() << std::endl;
}

/*!
 Set up the event recorder by defining the fault boundaries (in terms of block IDs).
 */
void EventRecorder::init(SimFramework *_sim) {
	VCSimulation				*sim = static_cast<VCSimulation*>(_sim);
	std::string					file_name = sim->getEventsFile().c_str();
	BlockList::const_iterator	it;
	
	if (!sim->isRootNode()) return;
	
    eFile.open(file_name.c_str());
    assertThrow(eFile, "Failed to open events file: " + file_name);
    eFile << std::fixed << std::setprecision(2);
	
    // SETTING FID BOUNDS, ID OF THE FIRST AND THE LAST BLOCK
    FaultID prev_fid = UNDEFINED_FAULT_ID;
	BlockID start_block_id = UNDEFINED_BLOCK_ID;
    for(it=sim->begin();it!=sim->end();++it) {
        if(it->getFaultID() != prev_fid) {
			if (start_block_id != UNDEFINED_BLOCK_ID) {
				faultBlockIDs.push_back(std::make_pair(start_block_id, it->getBlockID()-1));
			}
			start_block_id = it->getBlockID();
            prev_fid = it->getFaultID();
        }
    }
	faultBlockIDs.push_back(std::make_pair(start_block_id, sim->numGlobalBlocks()-1));
	assertThrow(start_block_id!=UNDEFINED_BLOCK_ID, "Could not find any blocks in simulation.");
}

/*!
 Record the year, trigger block and slip amount of each ruptured block in the event.
 */
SimRequest EventRecorder::run(SimFramework *_sim) {
	VCSimulation		*sim = static_cast<VCSimulation*>(_sim);
	VCEvent::iterator	it;

	if (!sim->isRootNode()) return SIM_STOP_OK;
	
	// Write the event year, number, and which block triggered the event
    eFile << "YEAR: " << sim->getCurrentEvent().getEventYear() << '\t' << sim->getEventCount() << '\n';
    eFile << "TRIG: " << sim->getCurrentEvent().getEventTrigger() << '\n';
	
	// For each ruptured block, write the total slip of the block
    eFile << "SLIP: ";
    for(it=sim->getCurrentEvent().begin();it!=sim->getCurrentEvent().end();++it) {
        eFile << it->first << ':' << it->second.slip << '\t';
    }
    eFile << "\n";
	
	// For each fault in the event, write a character graphic showing which elements ruptured
    for(unsigned f=0; f<faultBlockIDs.size(); ++f) {
        bool is_event = false;
        for(unsigned i=faultBlockIDs[f].first; i<=faultBlockIDs[f].second; ++i) {
            bool is_slip = (sim->getCurrentEvent().getEventSlip(i) != 0);
			
            if(is_event) {
                eFile << (is_slip ? "*" : " ");
            } else if(is_slip) {
				is_event = true;
				eFile << sim->getBlock(i).getFaultID() << "\t[";
				for(unsigned j=0; j<i-faultBlockIDs[f].first; j++) eFile << " ";
				eFile << "*";
            }
        }
        if(is_event) eFile << "]\n";
    }
    eFile << "\n";
    eFile << "\n";
	
	return SIM_STOP_OK;
}

void EventRecorder::finish(SimFramework *_sim) {
    eFile.close();
}
