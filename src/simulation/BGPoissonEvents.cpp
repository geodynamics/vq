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

#include "BGPoissonEvents.h"
#include <algorithm>

void BGPoissonEvents::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# Creating background Poisson events (mean interevent time "
	<< sim->getBGEventMeanInterevent() << " years)." << std::endl;
}

void BGPoissonEvents::init(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	
	last_year = sim->getSimStart() + (-log(_sim->randDouble())/(1.0/sim->getBGEventMeanInterevent()));
}

SimRequest BGPoissonEvents::run(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	double					cur_year, mag, r, bg_x, bg_y, theta;
	int						selected_ind;
	VCGeneralEvent			bg_event;
	quakelib::Conversion	convert;
	
	if (!_sim->isRootNode()) return SIM_STOP_OK;
	
	// Get background event parameters
	_l = 1.0/sim->getBGEventMeanInterevent();
	_mag_min = sim->getBGEventMinMagnitude();
	_mag_max = sim->getBGEventMaxMagnitude();
	_d = sim->getBGEventDistance();
	_q = sim->getBGEventDistanceDecay();
	
	// Get the current simulation year and clear old events
	cur_year = sim->getCurrentEvent().getEventYear();
	sim->clearBackgroundEvents();
	
	// Generate background events until we catch up
	while (last_year < cur_year) {
		// Determine a magnitude within the Gutenberg-Richter range
		mag = _mag_max+1;
		while (mag > _mag_max) {
			mag = _mag_min - log10(_sim->randDouble());
		}
		
		// Randomly select a block to be the basis for determining location
		selected_ind = sim->randInt(sim->numGlobalBlocks());
		
		// Choose an angle and distance from the selected block using the BASS model
		theta = 2.0 * M_PI * sim->randFloat();
		r = _d*powf(sim->randFloat(), -1.0/(_q-1.0)) - _d;
		
		// Set the x, y position of the new event to be an offset from the randomly chosen block
		bg_x = sim->getBlock(selected_ind).center()[0] + convert.m2km(r*cos(theta));
		bg_y = sim->getBlock(selected_ind).center()[1] + convert.m2km(r*sin(theta));
		
		bg_event = VCGeneralEvent(mag, last_year, bg_x, bg_y);
		sim->addBackgroundEvent(bg_event);
		last_year += (-log(_sim->randDouble())/_l);
	}
	
	return SIM_STOP_OK;
}
