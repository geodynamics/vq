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

#include "GreensInit.h"
#include <sstream>

/*!
 Calculate how much memory the Greens function matrices will require
 given the number of blocks and number of nodes.
 */
void GreensInit::dryRun(SimFramework *_sim) {
	VCSimulation		*sim = static_cast<VCSimulation*>(_sim);
	double				init_memory, running_memory, init_memory_per_node, running_memory_per_node;
	std::stringstream	rm_ss, rmpn_ss, im_ss, impn_ss;
	int					rm_ind, rmpn_ind, im_ind, impn_ind, nblocks, num_nodes;
	std::string			space_vals[] = {"bytes", "kilobytes", "megabytes", "gigabytes", "terabytes", "petabytes"};
	
	num_nodes = sim->getWorldSize();
	nblocks = sim->numGlobalBlocks();
	
	// Calculate how much memory we will use during the run
	running_memory_per_node = 2.0*sizeof(GREEN_VAL)*nblocks*nblocks/num_nodes;
	init_memory_per_node = running_memory_per_node*(2-(1.0/float(num_nodes)));
	running_memory = running_memory_per_node*num_nodes;
	init_memory = init_memory_per_node*num_nodes;
	
	// Convert bytes to KB/MB/etc format
	rm_ind = rmpn_ind = im_ind = impn_ind = 0;
	while (running_memory > 1024) { running_memory /= 1024; rm_ind++; }
	while (init_memory > 1024) { init_memory /= 1024; im_ind++; }
	while (running_memory_per_node > 1024) { running_memory_per_node /= 1024; rmpn_ind++; }
	while (init_memory_per_node > 1024) { init_memory_per_node /= 1024; impn_ind++; }
	
	rm_ss << running_memory << " " << space_vals[rm_ind];
	rmpn_ss << running_memory_per_node << " " << space_vals[rmpn_ind];
	im_ss << init_memory << " " << space_vals[im_ind];
	impn_ss << init_memory_per_node << " " << space_vals[impn_ind];

	sim->console() << "# Greens function memory during initialization: " << im_ss.str() << " (" << impn_ss.str() << " per CPU)" << std::endl;
	sim->console() << "# Greens function memory during simulation: " << rm_ss.str() << " (" << rmpn_ss.str() << " per CPU)" << std::endl;
}

/*!
 Calculates Greens functions for given parameters and type of calculation.
 */
void GreensInit::init(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	GreensFuncCalcBarnesHut	bh_calc;
    GreensFuncCalc2011      g2011_calc;
	GreensFuncFileParse		file_parse;
	double					start_time;
	double					shear_bytes, normal_bytes;
	int						shear_ind, norm_ind;
	std::string				space_vals[] = {"bytes", "kilobytes", "megabytes", "gigabytes", "terabytes", "petabytes"};
	
	start_time = sim->curTime();
	
	switch (sim->getGreensCalcMethod()) {
		case GREENS_FILE_PARSE:
			sim->console() << "# Reading Greens function data from file " << sim->getGreensInputfile() << std::flush;
			file_parse.CalculateGreens(sim);
			break;
		case GREENS_CALC_BARNES_HUT:
			sim->console() << "# Calculating Greens function with Barnes Hut technique" << std::flush;
			bh_calc.CalculateGreens(sim);
			break;
        case GREENS_CALC_2011:
			sim->console() << "# Calculating Greens function with the 2011 Okada class" << std::flush;
			g2011_calc.CalculateGreens(sim);
			break;
		default:
			exit(-1);
			break;
	}
	
	// Write out the number of seconds it took for the Greens function calculation
	sim->console() << std::endl << "# Greens function took " << sim->curTime() - start_time << " seconds." << std::endl;
	
	// Determine Green's function matrix memory usage
	// TODO: global reduction for parallel simulations
	shear_bytes = sim->greenShear()->mem_bytes();
	normal_bytes = sim->greenNormal()->mem_bytes();
	shear_ind = norm_ind = 0;
	while (shear_bytes > 1024.0) { shear_bytes /= 1024.0; shear_ind++; }
	while (normal_bytes > 1024.0) { normal_bytes /= 1024.0; norm_ind++; }
	sim->console() << "# Greens shear matrix takes " << shear_bytes << " " << space_vals[shear_ind] << "." << std::endl;
	sim->console() << "# Greens normal matrix takes " << normal_bytes << " " << space_vals[norm_ind] << "." << std::endl;
}
