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

#include "VCInitBlocks.h"
#include "SimError.h"

void VCInitBlocks::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# Initializing blocks." << std::endl;
}

/*!
 Set up the dry run for parallel systems by partitioning
 the blocks and creating MPI related data types.
 */
void VCInitBlocks::dryRun(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	sim->partitionBlocks();
}

/*!
 Initialize the simulation by partitioning the model, allocating the necessary memory,
 setting value pointers and creating the needed MPI data types.
 */
void VCInitBlocks::init(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	int						i;
	BlockList::iterator		bit;
	
	assertThrow(sim->numGlobalBlocks() > 0, "Simulation must include at least 1 block.");
	
	sim->partitionBlocks();
	
	// Initialize simulation arrays and classes
	sim->setupArrays(sim->numGlobalBlocks(),
							  sim->numLocalBlocks(),
							  // use compressed array for Barnes Hut style Greens function calculations
							  sim->getGreensCalcMethod()==GREENS_CALC_BARNES_HUT,
							  // transposed array for faster sweep calculations
							  sim->useTransposedMatrix());
	
	// Set the block pointers to the storage arrays
	for (i=0,bit=sim->begin();bit!=sim->end();++bit,++i)
		bit->setStatePtrs(sim->getShearStressPtr(i),
						  sim->getNormalStressPtr(i),
                          sim->getFShearStressPtr(i),
						  sim->getFNormalStressPtr(i),
						  sim->getUpdateFieldPtr(i));
	
	// Set the starting year of the simulation
	sim->setYear(sim->getSimStart());

	// Determine neighbors of blocks
    sim->determineBlockNeighbors();
}
