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

#include "VCBlock.h"
#include "VCEvent.h"

#ifndef _VCSIM_DATA_EVENTS_H_
#define _VCSIM_DATA_EVENTS_H_

class VCSimDataEvents {
    private:
        //! Current event in the simulation (older events are discarded)
        VCEvent                     cur_event;

        //! Set of background events that occurred between current event and previous event
        VCGeneralEventSet           bg_events;

        //! Current count of events
        unsigned int                event_cnt;

    public:
        VCSimDataEvents(void) : event_cnt(0) {};

        int getEventCount(void) const {
            return event_cnt;
        };
        VCEvent &getCurrentEvent(void) {
            return cur_event;
        };
        void addEvent(const VCEvent &new_event) {
            cur_event = new_event;
            event_cnt++;
        };

        void clearBackgroundEvents(void) {
            bg_events.clear();
        };
        void addBackgroundEvent(const VCGeneralEvent &new_bg_event) {
            bg_events.push_back(new_bg_event);
        };
        VCGeneralEventSet &getBGEvents(void) {
            return bg_events;
        };
};

#endif
