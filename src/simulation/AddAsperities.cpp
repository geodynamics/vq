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

#include "AddAsperities.h"

void AddAsperities::initDesc(const SimFramework *_sim) const {
    const VCSimulation          *sim = static_cast<const VCSimulation *>(_sim);

    sim->console() << "# Randomly adding " << sim->getNumAsperities() << " asperities (width "
                   << sim->getAsperityWidth() << " blocks, magnitude: " << sim->getAsperityMagnitude() << ")." << std::endl;
}

void AddAsperities::init(SimFramework *_sim) {
    VCSimulation            *sim = static_cast<VCSimulation *>(_sim);
    int                     num_asperities, i, n, selection, asperity_width;
    BlockIDList             avail_faults;
    BlockID                 bid;
    double                  asperity_mag;
    BlockIDSet              selected_faults;
    BlockIDSet::iterator    it;
    bool                    first;

    num_asperities = sim->getNumAsperities();
    asperity_width = sim->getAsperityWidth();
    asperity_mag = sim->getAsperityMagnitude();

    if (sim->numGlobalBlocks() < num_asperities*asperity_width) throw "Asperity size mismatch.";

    //! First, make a list of potential surface blocks to put asperities on.
    //! For simplicity, this just allows asperities on block IDs which are multiples of the width.
    for (i=0; i<sim->numSurfaceBlocks()/asperity_width-1; ++i) avail_faults.push_back(i*asperity_width);

    //! Next, select the blocks to add asperities to and set the asperity
    //! (stress drop multiplication) for all blocks in the layer.
    for (i=0; i<num_asperities; ++i) {
        selection = sim->randInt(avail_faults.size());
        bid = avail_faults[selection];
        selected_faults.insert(bid);

        for (n=0; n<asperity_width; ++n) {
            //for (int p=0;p<sim->numLayers();++p) {
            //  BlockID selected_bid = (bid+n)+(p*sim->numSurfaceBlocks());
            //  double old_stress_drop = sim->getBlock(selected_bid).getStressDrop();
            //  sim->getBlock(selected_bid).setStressDrop(old_stress_drop*asperity_mag);
            //}
        }

        //! Remove any selected faults from the list to avoid double selection.
        avail_faults.erase(avail_faults.begin()+selection);
    }

    //! Finally, notify the user which blocks were modified.
    if (sim->isRootNode()) {
        first = true;
        sim->console() << "# Asperities added to blocks: ";

        for (it=selected_faults.begin(); it!=selected_faults.end(); ++it) {
            sim->console() << (first?"":", ") << *it;

            if (asperity_width > 1) sim->console() << "-" << *it+asperity_width-1;

            first = false;
        }

        sim->console() << std::endl;
    }
}
