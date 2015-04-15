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
    Simulation        *sim = static_cast<Simulation *>(_sim);
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
    Simulation            *sim = static_cast<Simulation *>(_sim);
    GreensFuncCalcBarnesHut bh_calc;
    GreensFuncCalcStandard  gstandard_calc;
    GreensFuncFileParse     file_parse;
    std::string             space_vals[] = {"bytes", "kilobytes", "megabytes", "gigabytes", "terabytes", "petabytes"};

    double start_time = sim->curTime();

    switch (sim->getGreensCalcMethod()) {
        case GREENS_FILE_PARSE:
            sim->console() << "# Reading Greens function data from file " << sim->getGreensInputfile() << std::flush;
            file_parse.CalculateGreens(sim);
            break;

        case GREENS_CALC_BARNES_HUT:
            sim->console() << "# Calculating Greens function with Barnes Hut technique" << std::flush;
            bh_calc.CalculateGreens(sim);
            break;

        case GREENS_CALC_STANDARD:
            sim->console() << "# Calculating Greens function with the standard Okada class" << std::flush;
            gstandard_calc.CalculateGreens(sim);
            break;

        default:
            exit(-1);
            break;
    }

    // Write out the number of seconds it took for the Greens function calculation
    sim->console() << std::endl << "# Greens function took " << sim->curTime() - start_time << " seconds." << std::endl;

    // Determine Green's function matrix memory usage
    double shear_bytes = sim->greenShear()->mem_bytes();
    double normal_bytes = sim->greenNormal()->mem_bytes();
    int shear_ind = (log(shear_bytes)/log(2))/10;
    int norm_ind = (log(shear_bytes)/log(2))/10;
    double abbr_shear_bytes = shear_bytes/pow(2,shear_ind*10);
    double abbr_normal_bytes = normal_bytes/pow(2,norm_ind*10);

    sim->console() << "# Greens shear matrix takes " << abbr_shear_bytes << " " << space_vals[shear_ind] << std::endl;
    sim->console() << "# Greens normal matrix takes " << abbr_normal_bytes << " " << space_vals[norm_ind] << std::endl;
	//
<<<<<<< HEAD
	// yoder: and now, if GF limits have been set, let's spin through the GF values and truncate non-physical entries:
	// (note we can also do this during the various **.CalculateGreens(). 
=======
	// yoder: print some greens max/min/mean stats:
    double shear_min, shear_max, shear_mean, normal_min, normal_max, normal_mean;
    getGreensStats(sim, shear_min, shear_max, shear_mean, normal_min, normal_max, normal_mean);
    sim->console() << "# Greens Shear:\n max: " << shear_max << "\n min: " << shear_min << "\n mean: " << shear_mean << std::endl << std::endl;
    sim->console() << "# Greens Normal:\n max: " << normal_max << "\n min: " << normal_min << "\n mean: " << normal_mean << std::endl << std::endl;
>>>>>>> d4184e4a58897c1316ec8ec2a906ca0035ba2810

#ifdef MPI_C_FOUND
    //
    if (sim->getWorldSize() > 1) {
        //
        // initialize these to suppress memcheck errors:
        double global_shear_bytes = std::numeric_limits<float>::quiet_NaN();
        double global_normal_bytes = std::numeric_limits<float>::quiet_NaN();
        //double global_shear_bytes, global_normal_bytes;
        //
        MPI_Reduce(&shear_bytes, &global_shear_bytes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&normal_bytes, &global_normal_bytes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        //
        if (sim->isRootNode()) {
            int global_shear_ind = (log(global_shear_bytes)/log(2))/10;
            int global_norm_ind = (log(global_normal_bytes)/log(2))/10;

            double abbr_global_shear_bytes = global_shear_bytes/pow(2,global_shear_ind*10);
            double abbr_global_normal_bytes = global_normal_bytes/pow(2,global_norm_ind*10);
            //
            sim->console() << "# Global Greens shear matrix takes " << abbr_global_shear_bytes << " " << space_vals[global_shear_ind] << "." << std::endl;
            sim->console() << "# Global Greens normal matrix takes " << abbr_global_normal_bytes << " " << space_vals[global_norm_ind] << "." << std::endl;
        };

    }

#endif
    

}

// yoder:
void GreensInit::getGreensStats(Simulation *sim, double &shear_min, double &shear_max, double &shear_mean, double &normal_min, double &normal_max, double &normal_mean) {
	// gather max/min values for greens functions (which we may have defined in the parameter file).
	// note: we (see ProgressMonitor.cpp/h; we could use BlockVal declarations (is there a GreensVal object? we could create one) to
	// combine the block_id, or in this case the block_id pair, with the gr_values.
	//
	int         lid;
    BlockID     gid;
    double      shear_min_l, shear_max_l, normal_min_l, normal_max_l;	// local-node copies of min/max. we'll mpi_reduce --> shear_min, shear_max, etc.
    double      shear_sum, normal_sum;
    double      cur_shear, cur_normal;
    BlockList::iterator nt;

    shear_min  = DBL_MAX;
    normal_min = DBL_MAX;
    shear_max  = -DBL_MAX;
    normal_max = -DBL_MAX;
    //
    shear_sum = normal_sum = 0;
    
    //min_val.block_id = max_val.block_id = sum_val.block_id = UNDEFINED_ELEMENT_ID;

    // Get the minimum, maximum and sum CFF values on this node
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        //
        // now, loop over the sim events to get greens-pair values:
        for (nt=sim->begin(); nt!=sim->end(); ++nt) {
            //stress_drop += (nt->slip_rate()/norm_velocity)*sim->getGreenShear(gid, nt->getBlockID());
            cur_shear  = sim->getGreenShear(gid, nt->getBlockID());
            cur_normal = sim->getGreenNormal(gid, nt->getBlockID());
            //
            if (cur_shear<shear_min_l) shear_min_l = cur_shear;
            if (cur_shear>shear_max_l) shear_max_l = cur_shear;
            shear_sum += cur_shear;
            //
            if (cur_normal<normal_min_l) normal_min_l = cur_normal;
            if (cur_normal>normal_max_l) normal_max_l = cur_normal;
            normal_sum += cur_normal;
        };
    };
    

    // Reduce to the min/avg/max over all nodes
    sim->allReduceBlockVal(shear_min_l, shear_min, BLOCK_VAL_MIN);
    sim->allReduceBlockVal(shear_max_l, shear_max, BLOCK_VAL_MAX);
    //
    sim->allReduceBlockVal(normal_min_l, normal_min, BLOCK_VAL_MIN);
    sim->allReduceBlockVal(normal_max_l, normal_max, BLOCK_VAL_MAX);
    
    sim->allReduceBlockVal(shear_sum, shear_mean, BLOCK_VAL_SUM);
    sim->allReduceBlockVal(normal_sum, normal_mean, BLOCK_VAL_SUM);
    shear_mean  /= sim->numGlobalBlocks();
    normal_mean /= sim->numGlobalBlocks();
};












