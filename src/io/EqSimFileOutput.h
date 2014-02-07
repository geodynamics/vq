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
#include <fstream>
#include "QuakeLib.h"

#ifndef _EQSIM_FILE_OUTPUT_H_
#define _EQSIM_FILE_OUTPUT_H_

/*!
 A plugin to write VC events in the EqSim event format.
 */
class EqSimFileOutput : public SimPlugin {
private:
	quakelib::EQSimEventWriter	eqsim_event_file;
	
public:
    virtual std::string name(void) const { return "EQSim file output"; }
	virtual void initDesc(const SimFramework *_sim) const;
	
	virtual void init(SimFramework *_sim);
    virtual SimRequest run(SimFramework *_sim);
	virtual void finish(SimFramework *_sim);
};

#endif
