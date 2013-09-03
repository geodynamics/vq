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

#include "ProgressMonitor.h"
#include <iomanip>
#include <numeric>

#ifdef HAVE_FLOAT_H
#include <float.h>
#endif

#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif

void ProgressMonitor::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# Displaying simulation progress every ";
	if (sim->getProgressPeriod() == 1) sim->console() << "second." << std::endl;
	else sim->console() << sim->getProgressPeriod() << " seconds." << std::endl;
}

/*!
 Print initial information about the simulation environment and model.
 */
void ProgressMonitor::init(SimFramework *_sim) {
	VCSimulation	*sim = static_cast<VCSimulation*>(_sim);
	BlockVal		min, avr, max;
	int				width = 30;
	
	getStats(sim, min, avr, max);
	
	next_event_count = sim->itersPerSecond();
	
	if (!sim->isRootNode()) return;
	
	sim->console() << "#" << std::endl;
	sim->console() << "# **********************************************" << std::endl;
	sim->console() << "# ***" << std::endl;
	sim->console() << std::setw(width) << std::left << "# *** Blocks" << ": " << sim->numGlobalBlocks() << std::endl;
	//sim->console() << std::setw(width) << std::left << "# *** Layers" << ": " << sim->numLayers() << std::endl;
	sim->console() << std::setw(width) << std::left << "# *** Faults" << ": " << sim->numFaults() << std::endl;
	sim->console() << std::setw(width) << std::left << "# ***" << std::endl;
	sim->console() << std::setw(width) << std::left << "# *** Present Time (years)" << ": " << sim->getYear() << std::endl;
	sim->console() << std::setw(width) << std::left << "# *** Min  cff" << ": " << min.val << std::endl;
	sim->console() << std::setw(width) << std::left << "# *** Mean cff" << ": " << avr.val << std::endl;
	sim->console() << std::setw(width) << std::left << "# *** Max  cff" << ": " << max.val << std::endl;
	sim->console() << std::setw(width) << std::left << "# ***" << std::endl;
	sim->console() << std::right;
}

/*!
 Gather statistics over all computing nodes about the CFF values.
 */
void ProgressMonitor::getStats(VCSimulation *sim, BlockVal &min, BlockVal &avg, BlockVal &max) {
	int			lid;
	BlockID		gid;
	BlockVal	min_val, max_val, sum_val;
	double		cur_cff;
	
	min_val.val = DBL_MAX;
	max_val.val = -DBL_MAX;
	sum_val.val = 0;
	min_val.block_id = max_val.block_id = sum_val.block_id = UNDEFINED_BLOCK_ID;
	
	// Get the minimum, maximum and sum CFF values on this node
	for (lid=0;lid<sim->numLocalBlocks();++lid) {
		gid = sim->getGlobalBID(lid);
		cur_cff = sim->getBlock(gid).getCFF();
		if (cur_cff < min_val.val) {
			min_val.val = cur_cff;
			min_val.block_id = gid;
		}
		sum_val.val += cur_cff;
		if (cur_cff > max_val.val) {
			max_val.val = cur_cff;
			max_val.block_id = gid;
		}
	}
	// Reduce to the min/avg/max over all nodes
	sim->allReduceBlockVal(min_val, min, BLOCK_VAL_MIN);
	sim->allReduceBlockVal(max_val, max, BLOCK_VAL_MAX);
	sim->allReduceBlockVal(sum_val, avg, BLOCK_VAL_SUM);
	avg.val /= sim->numGlobalBlocks();
}

/*!
 Periodically (roughly every getProgressPeriod() seconds) print simulation information.
 This includes the number of events and current stress values for the system.
 */
SimRequest ProgressMonitor::run(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	BlockVal				min, avg, max;
	int						max_bid_width;
	std::ios_base::fmtflags	fmt_flags;
	
	// Print the header on the first time
	if (first) {
		sim->console() << "# events        year      minCff[index]      avrCff      maxCff[index]" << std::endl;
		first = false;
	}
	
	if(sim->getEventCount() >= next_event_count) {
		// Gather simulation-wide statistics
		getStats(sim, min, avg, max);
		
		// Estimate how many events until next progress report and notify other processes
		next_event_count += (sim->getProgressPeriod()*sim->itersPerSecond())+1;
		
		// Write the current progress
		max_bid_width = (int)(log10(sim->numGlobalBlocks())+1);
		
		fmt_flags = sim->console().flags();
		sim->console()	<< std::fixed << std::setw(8) << sim->getEventCount()		// number of events
						<< std::setw(12) << std::setprecision(1) << sim->getYear()	// current simulation year
						<< std::scientific << std::setw(17-max_bid_width) << std::setprecision(3) << min.val
						<< "[" << std::setw(max_bid_width) << min.block_id << "]"	// min cff
						<< std::setw(12) << std::setprecision(3) << avg.val			// average cff
						<< std::setw(17-max_bid_width) << std::setprecision(3) << max.val
						<< "[" << std::setw(max_bid_width) << max.block_id << "]"	// max cff
						<< std::endl;
		sim->console().flags(fmt_flags);
	}
	
	return SIM_STOP_OK;
}
