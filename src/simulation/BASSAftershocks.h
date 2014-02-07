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

#ifndef _BASS_AFTERSHOCKS_H_
#define _BASS_AFTERSHOCKS_H_

/*!
 This plugin creates VCAftershockEvent objects associated with the
 current event using the BASS model and the user supplied parameters.
 The first generation aftershocks can occur anywhere near the ruptured
 blocks on the fault, later generation locations are based on the previous
 generation locations.
 */
class BASSAftershocks : public SimPlugin {
private:
	// BASS parameters
	float _Mm;		//! Minimum magnitude
	float _dM;		//! Strength of aftershock sequence (intensity of aftershocks)
	float _b;		//! Scaling of frequency magnitude
	float _c;		//! Start of aftershocks (days)
	float _p;		//! Decay rate of aftershocks (time)
	float _d;		//! Distance of aftershocks (meters)
	float _q;		//! Decay rate of aftershocks (distance)
	
	// Member data
	AftershockSet	*events_to_process;
	BlockIDSet		event_blocks;
	bool			first;
	
	// Work horse function
	unsigned int generateAftershocks(VCSimulation *sim, VCEventAftershock seed);
	
public:
    virtual std::string name(void) const { return "BASS model aftershocks"; }
	virtual void initDesc(const SimFramework *_sim) const;
	virtual bool needsTimer(void) const { return true; };
	
	virtual SimRequest run(SimFramework *_sim);
};

#endif
