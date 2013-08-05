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

#include "VCBlock.h"
#include "VCEvent.h"

#ifndef _VCSIM_DATA_EVENTS_H_
#define _VCSIM_DATA_EVENTS_H_

class VCSimDataEvents {
private:
	//! Current event in the simulation (older events are discarded)
	VCEvent						cur_event;
	
	//! Set of background events that occurred between current event and previous event
	VCGeneralEventSet			bg_events;
	
	//! Current count of events
    unsigned int				event_cnt;

public:
	VCSimDataEvents(void) : event_cnt(0) {};
	
	int getEventCount(void) const { return event_cnt; };
	VCEvent &getCurrentEvent(void) { return cur_event; };
	void addEvent(const VCEvent &new_event) { cur_event = new_event; event_cnt++; };
	
	void clearBackgroundEvents(void) { bg_events.clear(); };
	void addBackgroundEvent(const VCGeneralEvent &new_bg_event) { bg_events.push_back(new_bg_event); };
	VCGeneralEventSet &getBGEvents(void) { return bg_events; };
};

#endif
