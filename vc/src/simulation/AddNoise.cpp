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

#include "AddNoise.h"

void AddNoise::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# Adding noise (" << sim->getStressNoiseResolution()
		<< "km resolution) to blocks." << std::endl;
}

void AddNoise::init(SimFramework *_sim) {
	VCSimulation						*sim = static_cast<VCSimulation*>(_sim);
	double								stress_noise, stress_res_dist, cur_noise, old_stress_drop;
	std::map<quakelib::Vec<3>, double>	centers_noise;
	quakelib::Vec<3>					mid_pt;
	std::map<quakelib::Vec<3>, double>::iterator	cit;
	BlockList::iterator					it;
	
	stress_noise = sim->getStressNoise();
	stress_res_dist = sim->getStressNoiseResolution();
	
	//! For each block in the simulation, get the midpoint and record the noise added to the block.
	//! Each block within a specified radius of the selected block will have the same level of noise added.
	//! The noise is defined as a modification of the stress drop from x to x +/- stress_noise.
	for (it=sim->begin();it!=sim->end();++it) {
		cur_noise = 1.0e10;
		mid_pt = it->center();
		for (cit=centers_noise.begin();cit!=centers_noise.end();++cit) {
			if (mid_pt.dist(cit->first) <= stress_res_dist) {
				cur_noise = cit->second;
				break;
			}
		}
		if (cur_noise > 1.0e9) {
			cur_noise = (1.0-stress_noise) + 2*stress_noise*sim->randDouble();
			centers_noise.insert(std::make_pair(mid_pt, cur_noise));
		}
		old_stress_drop = it->getStressDrop();
		it->setStressDrop(old_stress_drop*cur_noise);
	}
}
