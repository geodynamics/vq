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

#include "AddAsperities.h"

void AddAsperities::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# Randomly adding " << sim->getNumAsperities() << " asperities (width "
		<< sim->getAsperityWidth() << " blocks, magnitude: " << sim->getAsperityMagnitude() << ")." << std::endl;
}

void AddAsperities::init(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	int						num_asperities, i, n, selection, asperity_width;
	BlockIDList				avail_faults;
	BlockID					bid;
	double					asperity_mag;
	BlockIDSet				selected_faults;
	BlockIDSet::iterator	it;
	bool					first;
	
	num_asperities = sim->getNumAsperities();
	asperity_width = sim->getAsperityWidth();
	asperity_mag = sim->getAsperityMagnitude();
	
	if(sim->numGlobalBlocks() < num_asperities*asperity_width) throw "Asperity size mismatch.";
	
	//! First, make a list of potential surface blocks to put asperities on.
	//! For simplicity, this just allows asperities on block IDs which are multiples of the width.
	for (i=0;i<sim->numSurfaceBlocks()/asperity_width-1;++i) avail_faults.push_back(i*asperity_width);

	//! Next, select the blocks to add asperities to and set the asperity
	//! (stress drop multiplication) for all blocks in the layer.
	for (i=0;i<num_asperities;++i) {
		selection = sim->randInt(avail_faults.size());
		bid = avail_faults[selection];
		selected_faults.insert(bid);
		for (n=0;n<asperity_width;++n) {
			//for (int p=0;p<sim->numLayers();++p) {
			//	BlockID selected_bid = (bid+n)+(p*sim->numSurfaceBlocks());
			//	double old_stress_drop = sim->getBlock(selected_bid).getStressDrop();
			//	sim->getBlock(selected_bid).setStressDrop(old_stress_drop*asperity_mag);
			//}
		}
		
		//! Remove any selected faults from the list to avoid double selection. 
		avail_faults.erase(avail_faults.begin()+selection);
	}
	
	//! Finally, notify the user which blocks were modified.
	if (sim->isRootNode()) {
		first = true;
		sim->console() << "# Asperities added to blocks: ";
		for (it=selected_faults.begin();it!=selected_faults.end();++it) {
			sim->console() << (first?"":", ") << *it;
			if (asperity_width > 1) sim->console() << "-" << *it+asperity_width-1;
			first = false;
		}
		sim->console() << std::endl;
	}
}
