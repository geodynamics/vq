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

#include "GreensInit.h"
#include <sstream>

/*!
 Calculate how much memory the Greens function matrices will require
 given the number of blocks and number of nodes.
 */
void GreensInit::dryRun(SimFramework *_sim) {
    VCSimulation        *sim = static_cast<VCSimulation *>(_sim);
    double              init_memory, running_memory, init_memory_per_node, running_memory_per_node;
    std::stringstream   rm_ss, rmpn_ss, im_ss, impn_ss;
    int                 rm_ind, rmpn_ind, im_ind, impn_ind, nblocks, num_nodes;
    std::string         space_vals[] = {"bytes", "kilobytes", "megabytes", "gigabytes", "terabytes", "petabytes"};

    num_nodes = sim->getWorldSize();
    nblocks = sim->numGlobalBlocks();

    // Calculate how much memory we will use during the run
    running_memory_per_node = 2.0*sizeof(GREEN_VAL)*nblocks*nblocks/num_nodes;
    init_memory_per_node = running_memory_per_node*(2-(1.0/float(num_nodes)));
    running_memory = running_memory_per_node*num_nodes;
    init_memory = init_memory_per_node*num_nodes;

    // Convert bytes to KB/MB/etc format
    rm_ind = rmpn_ind = im_ind = impn_ind = 0;

    while (running_memory > 1024) {
        running_memory /= 1024;
        rm_ind++;
    }

    while (init_memory > 1024) {
        init_memory /= 1024;
        im_ind++;
    }

    while (running_memory_per_node > 1024) {
        running_memory_per_node /= 1024;
        rmpn_ind++;
    }

    while (init_memory_per_node > 1024) {
        init_memory_per_node /= 1024;
        impn_ind++;
    }

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
    VCSimulation            *sim = static_cast<VCSimulation *>(_sim);
    GreensFuncCalcBarnesHut bh_calc;
    GreensFuncCalc2011      g2011_calc;
    GreensFuncFileParse     file_parse;
    double                  start_time;
    double                  shear_bytes, normal_bytes;
    int                     shear_ind, norm_ind;
    std::string             space_vals[] = {"bytes", "kilobytes", "megabytes", "gigabytes", "terabytes", "petabytes"};

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

    while (shear_bytes > 1024.0) {
        shear_bytes /= 1024.0;
        shear_ind++;
    }

    while (normal_bytes > 1024.0) {
        normal_bytes /= 1024.0;
        norm_ind++;
    }

    sim->console() << "# Greens shear matrix takes " << shear_bytes << " " << space_vals[shear_ind] << "." << std::endl;
    sim->console() << "# Greens normal matrix takes " << normal_bytes << " " << space_vals[norm_ind] << "." << std::endl;
}
