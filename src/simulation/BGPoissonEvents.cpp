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
