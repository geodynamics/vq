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
    
    // Debug:
    sim->console() << "\n\n** Degug: Begin debugging GreensInit::init() MPI bit\n";
    sim->console() << "**Debug: RootNode: " << getpid() << "\n";
    // debugging comment: basically, anything that uses the << operator seems to crash big-time.

#ifdef MPI_C_FOUND
    //
    if (sim->getWorldSize() > 1) {
    	//printf("**Debug(%d): begin MPI debug\n", getpid());
    	//std::cout << "** Debug("<<getpid()<<") Begin GreensInit::init() MPI bit\n";		// even this appears to fail on child nodes, so
    	                                                                              // i'm guessing that a seg-fault long since made a mess.
        // yoder: try initializing these (just to suppress memcheck errors):
        double global_shear_bytes = std::numeric_limits<float>::quiet_NaN();
        double global_normal_bytes = std::numeric_limits<float>::quiet_NaN();
        //double global_shear_bytes, global_normal_bytes;
        //
        MPI_Reduce(&shear_bytes, &global_shear_bytes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&normal_bytes, &global_normal_bytes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //
        int global_shear_ind = (log(global_shear_bytes)/log(2))/10;
        int global_norm_ind = (log(global_normal_bytes)/log(2))/10;
        double abbr_global_shear_bytes = global_shear_bytes/pow(2,global_shear_ind*10);
        double abbr_global_normal_bytes = global_normal_bytes/pow(2,global_norm_ind*10);
        
        printf("\n** GreensInit Debugging(%d): %f, %f ## %f, %f \n", getpid(), global_shear_bytes, global_normal_bytes, abbr_global_shear_bytes, abbr_global_normal_bytes); 
        //
        /*
        printf("**Debug(%d): now printf globals (%d)\n", getpid(), getpid());
        printf("**Debug(%d): ... and one more write...\n", getpid());
        printf("**#Debug(%d): \n", getpid());
        
        std::cout <<"**Debug(" << getpid() << "): \n";
        std::cout << "       is_root: " << sim->isRootNode() << "\n";
        sim->console() << "**DebugConsole(" << getpid() << "): now ->console() globals: " << getpid() << "\n";
        
        // Debug:
        
        std::cout << "** Debug("<<getpid()<<") cout:: try one more write to console():\n";
        sim->console() << "** Debug("<<getpid()<<") trying one more write to console()\n";
        std::cout << "** Debug("<<getpid()<<") cout:: wrote to console again...\n";
        
        std::cout << "**Debug " << getpid() << "# Global Greens shear matrix takes " << abbr_global_shear_bytes << " " << "appropri-bytes" << "." << std::endl;
        std::cout << "**Debug " << getpid() << "# Global Greens normal matrix takes " << abbr_global_normal_bytes << " " << "appropri-bytes" << "." << std::endl;
        
        sim->console() << "# Global Greens shear matrix takes " << abbr_global_shear_bytes << " " << "appropri-bytes" << "." << std::endl;
        sim->console() << "# Global Greens normal matrix takes " << abbr_global_normal_bytes << " " << "appropri-bytes" << "." << std::endl;
        */
        sim->console() << "# Global Greens shear matrix takes " << abbr_global_shear_bytes << " " << space_vals[global_shear_ind] << "." << std::endl;
        sim->console() << "# Global Greens normal matrix takes " << abbr_global_normal_bytes << " " << space_vals[global_norm_ind] << "." << std::endl;
        
    }

#endif
    sim->console() << "** Degug: (" << getpid() << ") END debugging GreensInit::init() MPI bit\n";
}
