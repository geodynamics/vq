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

#ifndef _BLOCK_VAL_COMPUTE_H_
#define _BLOCK_VAL_COMPUTE_H_

/*!
 This plugin initializes internal block values before simulation starts.
 Currently it is used to calculate the stress drop of a block given
 a specified recurrence interval. This cannot be done when the block is
 read because the stress drop will depend on the self shear value of the
 Greens function.
 */
class BlockValCompute : public SimPlugin {
    public:
        virtual std::string name() const {
            return "Initialize block internal values";
        };
        virtual void initDesc(const SimFramework *_sim) const {};

        virtual void init(SimFramework *_sim);
};

#endif
