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
	// yoder: print some greens max/min/mean stats:
    double shear_min, shear_max, shear_mean, normal_min, normal_max, normal_mean;
    getGreensStats(sim, shear_min, shear_max, shear_mean, normal_min, normal_max, normal_mean);
    sim->console() << "# Greens Shear:\n max: " << shear_max << "\n min: " << shear_min << "\n mean: " << shear_mean << std::endl << std::endl;
    sim->console() << "# Greens Normal:\n max: " << normal_max << "\n min: " << normal_min << "\n mean: " << normal_mean << std::endl << std::endl;
    //
    // yoder: and now, get Greens Stats for (off)diagonal elements separately.
    double shear_diag_min, shear_diag_max, shear_diag_mean, normal_diag_min, normal_diag_max, normal_diag_mean;
    double shear_offdiag_min, shear_offdiag_max, shear_offdiag_mean, normal_offdiag_min, normal_offdiag_max, normal_offdiag_mean;
    getGreensDiagStats(sim, shear_diag_min, shear_diag_max, shear_diag_mean, normal_diag_min, normal_diag_max, normal_diag_mean, shear_offdiag_min, shear_offdiag_max, shear_offdiag_mean, normal_offdiag_min, normal_offdiag_max, normal_offdiag_mean);
    //
    sim->console() << "# Greens DiagShear:: " << shear_diag_min << " -- " << shear_diag_max << " (" << shear_diag_mean << ")\n"; // << std::endl << std::endl;
    sim->console() << "# Greens DiagNormal:: " << normal_diag_min << " -- " << normal_diag_max << " (" << normal_diag_mean << ")\n"; // std::endl << std::endl;
    sim->console() << "# Greens offDiagShear:: " << shear_offdiag_min << " -- " << shear_offdiag_max << " (" << shear_offdiag_mean << ")\n"; // std::endl << std::endl;
    sim->console() << "# Greens offDiagNormal:: " << normal_offdiag_min << " -- " << normal_offdiag_max << " (" << normal_offdiag_mean << ")\n\n"; //  std::endl << std::endl;

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
    //
    shear_min_l = normal_min_l =  DBL_MAX;
    shear_max_l = normal_max_l = -DBL_MAX;
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
    #ifdef MPI_C_FOUND
    int mpi_return=0;	// holder variable; mpi_reduce returns an int. maybe we'll have on efor each?
    // construct MPI_reduce calls here:
    mpi_return = MPI_Allreduce(&shear_min_l, &shear_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&shear_max_l, &shear_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&shear_sum, &shear_mean,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    mpi_return = MPI_Allreduce(&normal_min_l, &normal_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&normal_max_l, &normal_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&normal_sum, &normal_mean,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    #else
    shear_min = shear_min_l;
    shear_max = shear_max_l;
    shear_mean = shear_sum;
    
    normal_min = normal_min_l;
    normal_max = normal_max_l;
    normal_mean = normal_sum;
    #endif
    shear_mean  /= sim->numGlobalBlocks();
    normal_mean /= sim->numGlobalBlocks();
};

void GreensInit::getGreensDiagStats(Simulation *sim, double &shear_diag_min, double &shear_diag_max, double &shear_diag_mean, double &normal_diag_min, double &normal_diag_max, double &normal_diag_mean, double &shear_offdiag_min, double &shear_offdiag_max, double &shear_offdiag_mean, double &normal_offdiag_min, double &normal_offdiag_max, double &normal_offdiag_mean) {
	// gather max/min values for greens functions (which we may have defined in the parameter file).
	// note: we (see ProgressMonitor.cpp/h; we could use BlockVal declarations (is there a GreensVal object? we could create one) to
	// combine the block_id, or in this case the block_id pair, with the gr_values.
	// for (off) diagonal elements
	//
	int         lid;
    BlockID     gid;
    double      shear_diag_min_l, shear_diag_max_l, normal_diag_min_l, normal_diag_max_l;   // local-node copies...
    double      shear_diag_sum, normal_diag_sum;
    //
    double      shear_offdiag_min_l, shear_offdiag_max_l, normal_offdiag_min_l, normal_offdiag_max_l;	// local-node copies
    double      shear_offdiag_sum, normal_offdiag_sum;
    //
    double      cur_shear, cur_normal;
    //
    BlockList::iterator nt;
    //
    shear_diag_min_l = shear_offdiag_min_l = normal_diag_min_l = normal_offdiag_min_l = DBL_MAX;
    shear_diag_max_l = shear_offdiag_max_l = normal_diag_max_l = normal_offdiag_max_l = -DBL_MAX;
    //
    shear_diag_sum = normal_diag_sum = 0;
    shear_offdiag_sum = normal_offdiag_sum = 0;
    
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
            if (gid==nt->getBlockID()) {
            	// diagonal elements:
            	//printf("diag: %d, %d:: %f/%f\n", gid, nt->getBlockID(), cur_shear, cur_normal);
	            shear_diag_min_l = std::min(shear_diag_min_l, cur_shear);
	            shear_diag_max_l = std::max(shear_diag_max_l, cur_shear);
	            shear_diag_sum += cur_shear;
	            //
	            normal_diag_min_l = std::min(normal_diag_min_l, cur_normal);
	            normal_diag_max_l = std::max(normal_diag_max_l, cur_normal);
	            normal_diag_sum += cur_normal;
            }
            //
            else { 
            	//if (gid!=nt->getBlockID()) {
            	// off-diagonal elements:
    	        //printf("off-diag: %d, %d:: %f/%f\n", gid, nt->getBlockID(), cur_shear, cur_normal);
    	        shear_offdiag_min_l = std::min(shear_offdiag_min_l, cur_shear);
    	        shear_offdiag_max_l = std::max(shear_offdiag_max_l, cur_shear);
    	        shear_offdiag_sum += cur_shear;
    	        //
    	        normal_offdiag_min_l = std::min(normal_offdiag_min_l, cur_normal);
    	        normal_offdiag_max_l = std::max(normal_offdiag_max_l, cur_normal);
    	        normal_offdiag_sum += cur_normal;
    	    };
        };
    };
    

    // Reduce to the min/avg/max over all nodes
    #ifdef MPI_C_FOUND
    int mpi_return=0;	// holder variable; mpi_reduce returns an int. maybe we'll have on efor each?
    // construct MPI_reduce calls here:
    mpi_return = MPI_Allreduce(&shear_diag_min_l, &shear_diag_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&shear_diag_max_l, &shear_diag_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&shear_diag_sum, &shear_diag_mean,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    mpi_return = MPI_Allreduce(&shear_offdiag_min_l, &shear_offdiag_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&shear_offdiag_max_l, &shear_offdiag_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&shear_offdiag_sum, &shear_offdiag_mean,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    mpi_return = MPI_Allreduce(&normal_diag_min_l, &normal_diag_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&normal_diag_max_l, &normal_diag_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&normal_diag_sum, &normal_diag_mean,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    mpi_return = MPI_Allreduce(&normal_offdiag_min_l, &normal_offdiag_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&normal_offdiag_max_l, &normal_offdiag_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    mpi_return = MPI_Allreduce(&normal_offdiag_sum, &normal_offdiag_mean,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
    #else
    shear_diag_min = shear_diag_min_l;
    shear_diag_max = shear_diag_max_l;
    shear_diag_mean = shear_diag_sum;
    
    shear_offdiag_min = shear_offdiag_min_l;
    shear_offdiag_max = shear_offdiag_max_l;
    shear_offdiag_mean = shear_offdiag_sum;

    normal_diag_min = normal_diag_min_l;
    normal_diag_max = normal_diag_max_l;
    normal_diag_mean = normal_diag_sum;
    
    normal_offdiag_min = normal_offdiag_min_l;
    normal_offdiag_max = normal_offdiag_max_l;
    normal_offdiag_mean = normal_offdiag_sum;
    #endif
    shear_diag_mean  /= sim->numGlobalBlocks();
    normal_diag_mean /= sim->numGlobalBlocks();
    
    shear_offdiag_mean  /= sim->numGlobalBlocks();
    normal_offdiag_mean /= sim->numGlobalBlocks();
};











