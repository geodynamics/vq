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

#ifndef _UPDATE_BLOCK_STRESS_H_
#define _UPDATE_BLOCK_STRESS_H_

/*!
 Calculates the rate of stress change for all local blocks and
 uses this to determine which block will be the next to fail.
 */
class UpdateBlockStress : public SimPlugin {
    public:
        virtual std::string name() const {
            return "Block stress calculation";
        };
        virtual void initDesc(const SimFramework *_sim) const {};
        virtual bool needsTimer(void) const {
            return true;
        }

        virtual void init(SimFramework *_sim);
        virtual SimRequest run(SimFramework *_sim);
        virtual void finish(SimFramework *_sim);

    private:
        void nextTimeStep(BlockVal &fail_time);
        void stressRecompute(void);
        //void stressRecompute(double backslip);
        //void stressRecomputeB(void);
        //void stressRecompute(BlockVal &fail_time);

        double          *tmpBuffer;
        VCSimulation    *sim;
};

#endif
