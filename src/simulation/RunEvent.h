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

#include "VCSimulation.h"

#ifndef _RUN_EVENT_H_
#define _RUN_EVENT_H_

/*!
 Starts with an initial failed block and propagates the failure
 throughout the system using static and dynamic failure functions.
 */
class RunEvent : public SimPlugin {
    private:
        quakelib::ElementIDSet  local_failed_blocks;
        BlockIDProcMapping      global_failed_blocks;
        quakelib::ElementIDSet  looseBlocks;
        unsigned int            sweep_num;

        void processBlocksOrigFail(VCSimulation *sim, VCEventSweep &current_sweep);
        void processBlocksSecondaryFailures(VCSimulation *sim, VCEventSweep &current_sweep);
        virtual void markBlocks2Fail(VCSimulation *sim, const FaultID &trigger_fault);
        void recordEventStresses(VCSimulation *sim);

    public:
        virtual std::string name() const {
            return "Propagate event ruptures";
        };
        virtual void initDesc(const SimFramework *_sim) const {};
        virtual bool needsTimer(void) const {
            return true;
        };
        virtual SimRequest run(SimFramework *_sim);
};

#endif
