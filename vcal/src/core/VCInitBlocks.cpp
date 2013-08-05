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
