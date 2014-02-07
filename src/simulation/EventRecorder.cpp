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
