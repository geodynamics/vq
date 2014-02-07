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

#include "GreensKillInteraction.h"

void GreensKillInteraction::initDesc(const SimFramework *_sim) const {
    const VCSimulation          *sim = static_cast<const VCSimulation *>(_sim);

    sim->console() << "# Greens kill distance: " << sim->getGreensKillDistance() << " km." << std::endl;
}

/*!
 Test how many Greens values will be killed for the dry run.
 */
void GreensKillInteraction::dryRun(SimFramework *_sim) {
    VCSimulation                *sim = static_cast<VCSimulation *>(_sim);

    sim->console() << "# Total number of Greens values: " << sim->numLocalBlocks()*double(sim->numGlobalBlocks()) << std::endl;
    sim->console() << "# Number of Greens values killed: " << killInteraction(_sim, true) << std::endl;
}

/*!
 Actively remove Greens function values during the simulation.
 In future simulations, we might reallocate the sparse matrix to improve memory usage.
 */
void GreensKillInteraction::init(SimFramework *_sim) {
    killInteraction(_sim, false);
}

/*!
 Set a Greens function interaction value to 0 for blocks
 over a specified distance from each other.
 */
double GreensKillInteraction::killInteraction(SimFramework *_sim, const bool &dry_run) {
    VCSimulation                *sim = static_cast<VCSimulation *>(_sim);
    double                      num_kill;
    BlockList::const_iterator   it, jit;
    int                         lid;

    num_kill = 0;

    //! For each combination of fault segments, check if the 3D midpoint
    //! distances are over the kill distance.
    //! If so, set the Greens function value to 0.
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        Block &local_block = sim->getBlock(sim->getGlobalBID(lid));

        for (jit=sim->begin(); jit!=sim->end(); ++jit) {
            if (local_block.center_distance(*jit) > sim->getGreensKillDistance()) {
                num_kill++;

                if (!dry_run) {
                    sim->setGreens(local_block.getBlockID(), jit->getBlockID(), 0.0, 0.0);
                }
            }
        }
    }

    return num_kill;
}
