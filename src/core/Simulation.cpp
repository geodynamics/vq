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

#include "Simulation.h"
#include "SimFramework.h"

#ifdef VQ_HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef VQ_HAVE_STRING_H
#include <string.h>
#endif

#include <list>
#include <vector>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 Initialize the simulation by reading the parameter file and checking the validity of parameters.
 */
Simulation::Simulation(int argc, char **argv) : SimFramework(argc, argv) {
    srand(time(0));

    // Ensure we are given the parameter file name
    assertThrow(argc == 2, "usage: vc param_file");

    read_params(argv[argc-1]);

    // Check validity of parameters
    // TODO: add more checks here
    assertThrow(!getVersion().compare("2.0"),
                "sim.version: Parameter file version must be 2.0");
    assertThrow(getSimStart()>=0,
                "sim.start_year: Start year must be at least 0.");
    assertThrow(getSimStart() < getSimDuration(),
                "sim.start_year: Start year must be before end year.");
    assertThrow(getGreensCalcMethod() != GREENS_CALC_UNDEFINED,
                "Greens calculation method must be either standard, Barnes Hut or file based.");

    // Now that we have the parameters, write them out to a file
    // on the root node for record keeping purposes
    if (isRootNode()) {
        std::stringstream    param_file_name;
        param_file_name << "vq_params_" << getPID() << ".prm";
        this->write_params(param_file_name.str());
    }
}

void Simulation::output_stress(quakelib::UIndex event_num, quakelib::UIndex sweep_num) {
    quakelib::ModelStress       stress;
    quakelib::ModelStressState  stress_state;

    if (getStressOutfileType() == "") return;

    stress_state.setYear(getYear());
    stress_state.setEventNum(event_num);
    stress_state.setSweepNum(sweep_num);
    stress_state.setStartEndRecNums(num_stress_recs, num_stress_recs+numGlobalBlocks());

    // Store stress values
    for (unsigned int i=0; i<numGlobalBlocks(); ++i) {
        stress.add_stress_entry(i, shear_stress[i], normal_stress[i]);
    }

    num_stress_recs += numGlobalBlocks();

    if (getStressOutfileType() == "text") {
        // Write the stress state details
        stress_state.write_ascii(stress_index_outfile);
        stress_index_outfile.flush();

        // Write the stress details
        stress.write_ascii(stress_outfile);
        stress_outfile.flush();
    } else if (getStressOutfileType() == "hdf5") {
        // Write the stress state details
        stress_state.append_stress_state_hdf5(stress_data_file);

        // Write the stress details
        stress.append_stress_hdf5(stress_data_file);
    }
}

/*!
 Finish the simulation by deallocating memory and freeing MPI related structures.
 */
Simulation::~Simulation(void) {
#ifdef DEBUG
    console() << "Number matrix multiplies: " << num_mults << std::endl;
#endif

    if (mult_buffer) free( mult_buffer );

    if (decompress_buf) free( decompress_buf );

    //
    deallocateArrays();
}

/*!
 Initialize the VC specific part of the simulation by creating communication timers.
 */
void Simulation::init(void) {
#ifdef DEBUG
    reduce_comm_timer = initTimer("Reduce Comm", false, false);
    fail_comm_timer = initTimer("Fail Synch Comm", false, false);
    dist_comm_timer = initTimer("Distribute Comm", false, false);
    sweep_comm_timer = initTimer("Sweep Comm", false, false);
    mult_timer = initTimer("Vec-Mat Mults", false, false);
    num_mults = 0;
#endif
    SimFramework::init();

    if (getStressOutfileType() == "text") {
        if (getStressOutfile() == "" || getStressIndexOutfile() == "") {
            errConsole() << "ERROR: Stress file names cannot be blank." << std::endl;
            exit(-1);
        }

        stress_index_outfile.open(getStressIndexOutfile().c_str());
        stress_outfile.open(getStressOutfile().c_str());

        if (!stress_index_outfile.good()) {
            errConsole() << "ERROR: Could not open output file " << getStressIndexOutfile() << std::endl;
            exit(-1);
        }

        if (!stress_outfile.good()) {
            errConsole() << "ERROR: Could not open output file " << getStressOutfile() << std::endl;
            exit(-1);
        }

        quakelib::ModelStressState::write_ascii_header(stress_index_outfile);
        quakelib::ModelStress::write_ascii_header(stress_outfile);
    } else if (getStressOutfileType() == "hdf5") {
#ifdef HDF5_FOUND
        open_stress_hdf5_file(getStressOutfile());
#else
        sim->errConsole() << "ERROR: HDF5 library not linked, cannot use HDF5 output files." << std::endl;
        exit(-1);
#endif
    } else if (!(getStressOutfileType() == "")) {
        errConsole() << "ERROR: Unknown stress output file type " << getStressOutfileType() << std::endl;
        exit(-1);
    }

    num_stress_recs = 0;
}

#ifdef HDF5_FOUND
void Simulation::open_stress_hdf5_file(const std::string &hdf5_file_name) {
    hid_t   plist_id;

    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    if (plist_id < 0) exit(-1);

#ifdef HDF5_IS_PARALLEL
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
    // Create the data file, overwriting any old files
    stress_data_file = H5Fcreate(hdf5_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    if (stress_data_file < 0) exit(-1);

    // Create the stress index table
    quakelib::ModelStressState::setup_stress_state_hdf5(stress_data_file);

    // Create the stress table
    quakelib::ModelStress::setup_stress_hdf5(stress_data_file);

    H5Pclose(plist_id);
}
#endif

/*!
 Calculate the number of faults in the simulation based on
 the number of unique block fault IDs.
 */
int Simulation::numFaults(void) const {
    BlockList::const_iterator   it;
    FaultIDSet                  fault_set;

    for (it=begin(); it!=end(); ++it) fault_set.insert(it->getFaultID());

    return fault_set.size();
}

/*!
 Calculate the stress before and after an event of the blocks that failed during the event.
 */
void Simulation::getInitialFinalStresses(const quakelib::ElementIDSet &block_set, double &shear_init, double &shear_final, double &normal_init, double &normal_final) const {
    quakelib::ElementIDSet::const_iterator      it;

    shear_init = shear_final = normal_init = normal_final = 0.0;

    // For each specified block
    for (it=block_set.begin(); it!=block_set.end(); ++it) {
        // Add the before/after stresses if it is on this node
        // Non-local blocks will have incorrect stress data
        if (isLocalToNode(*it)) {
            shear_init += shear_stress0[*it];
            shear_final += getShearStress(*it);
            normal_init += normal_stress0[*it];
            normal_final += getNormalStress(*it);
        }
    }
}

void Simulation::sumStresses(const quakelib::ElementIDSet &block_set,
                             double &shear_stress_sum,
                             double &shear_stress0_sum,
                             double &normal_stress_sum,
                             double &normal_stress0_sum) const {
    quakelib::ElementIDSet::const_iterator      it;

    shear_stress_sum = shear_stress0_sum = normal_stress_sum = normal_stress0_sum = 0;

    for (it=block_set.begin(); it!=block_set.end(); ++it) {
        shear_stress_sum += getShearStress(*it);
        shear_stress0_sum += shear_stress0[*it];
        normal_stress_sum += getNormalStress(*it);
        normal_stress0_sum += normal_stress0[*it];
    }
}

void Simulation::printTimers(void) {
    if (!dry_run) printAllTimers(console(), world_size, node_rank, ROOT_NODE_RANK);
}

void Simulation::determineBlockNeighbors(void) {
    BlockList::iterator     bit, iit;
    quakelib::ElementIDSet  all_sweeps;
    double                  block_size;

    for (bit=begin(); bit!=end(); ++bit) {
        for (iit=begin(); iit!=end(); ++iit) {
            block_size  = floor(bit->largest_dimension() * 1e-3) * 1e3;

            if (bit->getBlockID() != iit->getBlockID() &&           // ensure blocks are not the same
                bit->getFaultID() == iit->getFaultID() &&           // ensure blocks are on the same fault
                bit->center().dist(iit->center()) < block_size * 2.0) { // ensure blocks are "nearby" each other
                neighbor_map[bit->getBlockID()].insert(iit->getBlockID());
                neighbor_map[iit->getBlockID()].insert(bit->getBlockID());
            }
        }
    }
}

std::pair<quakelib::ElementIDSet::const_iterator, quakelib::ElementIDSet::const_iterator> Simulation::getNeighbors(const BlockID &bid) const {
    std::map<BlockID, quakelib::ElementIDSet>::const_iterator       it;
    quakelib::ElementIDSet      empty_set;

    it = neighbor_map.find(bid);

    if (it != neighbor_map.end()) {
        return std::make_pair(it->second.begin(), it->second.end());
    } else {
        return std::make_pair(empty_set.end(), empty_set.end());
    }
}

/*!
 Computes CFF values for all blocks on this node.
 */
void Simulation::computeCFFs(void) {
    int         i;

    for (i=0; i<numLocalBlocks(); ++i) {
        BlockID gid = getGlobalBID(i);
        calcCFF(gid);
    }
}

//! Calculates and stores the CFF of this block.
void Simulation::calcCFF(const BlockID gid) {
    cff[gid] = fabs(shear_stress[gid]) - fabs(friction[gid]*normal_stress[gid]);
}

/*!
 Performs a matrix-vector multiply (C += A * B) where C is Nx1, A is MxN and B is Nx1.
 This directly accesses the matrix A and assumes it is stored in transpose format.
 It assumes c should be referenced by the local-global map.
 dense specifies whether the vector is likely mostly non-zero and just used for accounting purposes.
 The SSE version of this function is still under testing and should not yet be used since it assumes
 certain things about memory alignment and matrix padding that are not true for all models.
 */

// Whether to perform sparse multiplications - these have no effect on the simulation
// and will only slow things down, generally used for testing purposes
//#define PERFORM_SPARSE_MULTIPLIES

void Simulation::matrixVectorMultiplyAccum(double *c, const quakelib::DenseMatrix<GREEN_VAL> *a, const double *b, const bool dense) {
    int         x, y, l, width, height, array_dim;
    double      val;
#ifdef DEBUG

    if (dense) {
        num_mults++;
        startTimer(mult_timer);
    }

#endif

    height = numLocalBlocks();
    width = numGlobalBlocks();
    array_dim = localSize();

    //
    if (!decompress_buf) decompress_buf = (GREEN_VAL *)valloc(sizeof(GREEN_VAL)*array_dim);

    //decompress_buf = (GREEN_VAL *)valloc(sizeof(GREEN_VAL)*array_dim);

    /*if (!mult_buffer) mult_buffer = (double *)valloc(sizeof(double)*array_dim);
    for (x=0;x<height;++x) mult_buffer[x] = 0;
    for (y=0;y<width;++y) {
        if (b[y]) {
            for (x=0;x<height;++x) mult_buffer[x] += a->val(x,y)*b[y];
        }
    }
    for (x=0;x<height;++x) {
        l = getGlobalBID(x);
        c[l] += mult_buffer[x];
    }*/

    if (a->transpose()) {
        // This works by calculating the contribution of each input vector (b) element
        // to the final answer rather than calculating each element of the answer
        // individually. In the case where most of the vector elements are zero this
        // will be much faster than normal multiplication techniques. For VC about
        // 80% or more of the matrix-vector multiplications have sparse vectors.
        //
        // TODO:
        if (!mult_buffer) mult_buffer = (double *)valloc(sizeof(double)*array_dim);

        //mult_buffer = (double *)valloc(sizeof(double)*array_dim);    //??

        // Reset the temporary buffer
        for (x=0; x<height; ++x) mult_buffer[x] = 0;

        // Perform the multiplication
        if (dense) {
            for (y=0; y<width; ++y) {
                multiplyRow(mult_buffer, &(b[y]), a->getCol(decompress_buf, y), array_dim);
            }
        } else {
            for (y=0; y<width; ++y) {
                val = b[y];
#ifndef PERFORM_SPARSE_MULTIPLIES

                if (!val) continue;

#endif
                multiplyRow(mult_buffer, &val, a->getCol(decompress_buf, y), array_dim);
            }
        }

        // Add the temporary buffer values into the result array
        for (x=0; x<height; ++x) {
            l = getGlobalBID(x);
            c[l] += mult_buffer[x];
        }
    } else {
        for (x=0; x<height; ++x) {
            val = 0;
            multiplySumRow(&val, b, a->getRow(decompress_buf, x), width, dense);
            c[x] += val;
        }
    }

#ifdef DEBUG

    if (dense) stopTimer(mult_timer);

#endif
}

// SSE old execution notes
// Original speed (no OpenMP)
//  41.5 events per second
//  3662 MB/sec bandwidth
//  914 Mflop/sec
// SSE
//  41.3 events per second
//  3829 MB/sec bandwidth
//  955 Mflop/sec
//  4.5% improvement over original
// SSE loop unroll without prefetch
//  43.4 events per second
//  4059 MB/sec bandwidth
//  1013 Mflop/sec
//  10.8% improvement over original
// SSE loop unroll with a* prefetch
//  46.4 events per second
//  4316 MB/sec bandwidth
//  1077 Mflop/sec
//  17.8% improvement over original

// New SSE execution notes (desktop iMac)
// Benchmark max mem speed: 12300 MB/sec
// Benchmark max comp speed: 5058 MFlop/sec
// 9350 mults, 906 elems
// BW = 906*907*8/(9350/t)
// flops = 906*906*2/(9350/t)
// Original speed (no OpenMP, no SSE, one CPU)
//  10.007 secs
//  5857 MB/sec bandwidth
//  1462 Mflop/sec
//
// SSE loop unroll without prefetch
//  5.28896 secs
//  11083 MB/sec bandwidth
//  2767 Mflop/sec
// SSE loop unroll with a[] prefetch
//  5.14001 secs
//  11404 MB/sec bandwidth
//  2847 Mflop/sec


// NOTES: the SSE version is mostly memory bound. Performance could be improved by
// changing to floats.

//#define USE_SSE
//#define PREFETCH_A_VEC
//#define PREFETCH_B_VEC
//#define PREFETCH_C_VEC
#define SSE_LOOP_UNROLL     8

#ifdef USE_SSE  // SSE version
#include <xmmintrin.h>

void Simulation::multiplySumRow(double *c, const double *b, const GREEN_VAL *a, const int n, const bool dense) {
    __m128d     aval, bval, cval, tmpval;
    double      tmp[2];
#if SSE_LOOP_UNROLL != 16 && SSE_LOOP_UNROLL != 8 && SSE_LOOP_UNROLL != 6 && SSE_LOOP_UNROLL != 4 && SSE_LOOP_UNROLL != 2
#error "Invalid value of SSE_LOOP_UNROLL"
#endif
    cval = _mm_setzero_pd();

    for (int x=0; x<n; x+=SSE_LOOP_UNROLL) {
#ifdef PREFETCH_A_VEC
        _mm_prefetch(&a[x+SSE_LOOP_UNROLL], _MM_HINT_T0);
#endif
#ifdef PREFETCH_B_VEC
        _mm_prefetch(&b[x+SSE_LOOP_UNROLL], _MM_HINT_T0);
#endif
        aval = _mm_load_pd(&a[x]);
        bval = _mm_load_pd(&b[x]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);

#if SSE_LOOP_UNROLL > 2
        aval = _mm_load_pd(&a[x+2]);
        bval = _mm_load_pd(&b[x+2]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
#endif

#if SSE_LOOP_UNROLL > 4
        aval = _mm_load_pd(&a[x+4]);
        bval = _mm_load_pd(&b[x+4]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
#endif

#if SSE_LOOP_UNROLL > 6
        aval = _mm_load_pd(&a[x+6]);
        bval = _mm_load_pd(&b[x+6]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
#endif
    }

    _mm_store_pd(tmp, cval);
    c[0] += tmp[0] + tmp[1];
}

#else   // Non-SSE version

void Simulation::multiplySumRow(double *c, const double *b, const GREEN_VAL *a, const int n, const bool dense) {
    double val = 0;

    if (dense) {
        for (int x=0; x<n; ++x) {
            val += a[x]*b[x];
        }
    } else {
        for (int x=0; x<n; ++x) {
#ifndef PERFORM_SPARSE_MULTIPLIES

            if (!b[x]) continue;

#endif
            val += a[x]*b[x];
        }
    }

    c[0] += val;
}

#endif

/*!
 Multiplies each value in a by b and adds the result to c.
 Using this function is faster than leaving the code in matrixVectorMultiplyAccum.
 */
#ifdef USE_SSE
void Simulation::multiplyRow(double *c, const double *b, const GREEN_VAL *a, const int n) {
    __m128d     aval, bval, cval, tmpval;
#if SSE_LOOP_UNROLL != 16 && SSE_LOOP_UNROLL != 8 && SSE_LOOP_UNROLL != 6 && SSE_LOOP_UNROLL != 4 && SSE_LOOP_UNROLL != 2
#error "Invalid value of SSE_LOOP_UNROLL"
#endif
    bval = _mm_load1_pd(b);

    for (int x=0; x<n; x+=SSE_LOOP_UNROLL) {
#ifdef PREFETCH_A_VEC
        _mm_prefetch(&a[x+SSE_LOOP_UNROLL], _MM_HINT_T0);
#endif
#ifdef PREFETCH_C_VEC
        _mm_prefetch(&c[x+SSE_LOOP_UNROLL], _MM_HINT_T0);
#endif
        aval = _mm_load_pd(&a[x]);
        cval = _mm_load_pd(&c[x]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
        _mm_store_pd(&c[x], cval);

#if SSE_LOOP_UNROLL > 2
        aval = _mm_load_pd(&a[x+2]);
        cval = _mm_load_pd(&c[x+2]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
        _mm_store_pd(&c[x+2], cval);
#endif

#if SSE_LOOP_UNROLL > 4
        aval = _mm_load_pd(&a[x+4]);
        cval = _mm_load_pd(&c[x+4]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
        _mm_store_pd(&c[x+4], cval);
#endif

#if SSE_LOOP_UNROLL > 6
        aval = _mm_load_pd(&a[x+6]);
        cval = _mm_load_pd(&c[x+6]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
        _mm_store_pd(&c[x+6], cval);
#endif

#if SSE_LOOP_UNROLL > 8
        aval = _mm_load_pd(&a[x+8]);
        cval = _mm_load_pd(&c[x+8]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
        _mm_store_pd(&c[x+8], cval);
#endif

#if SSE_LOOP_UNROLL > 10
        aval = _mm_load_pd(&a[x+10]);
        cval = _mm_load_pd(&c[x+10]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
        _mm_store_pd(&c[x+10], cval);
#endif

#if SSE_LOOP_UNROLL > 12
        aval = _mm_load_pd(&a[x+12]);
        cval = _mm_load_pd(&c[x+12]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
        _mm_store_pd(&c[x+12], cval);
#endif

#if SSE_LOOP_UNROLL > 14
        aval = _mm_load_pd(&a[x+14]);
        cval = _mm_load_pd(&c[x+14]);
        tmpval = _mm_mul_pd(aval, bval);
        cval = _mm_add_pd(tmpval, cval);
        _mm_store_pd(&c[x+14], cval);
#endif
    }
}
#else
void Simulation::multiplyRow(double *c, const double *b, const GREEN_VAL *a, const int n) {
    for (int x=0; x<n; ++x) c[x] += b[0]*a[x];
}
#endif

/*!
 Distributes the local part of the update field to other nodes and
 receives their local update fields.
 */
void Simulation::distributeUpdateField(void) {
#ifdef MPI_C_FOUND
#ifdef DEBUG
    startTimer(dist_comm_timer);
#endif
    int     i;
    BlockID bid;

    // Copy the local update field values to the send buffer
    for (i=0; i<numLocalBlocks(); ++i) {
        bid = updateFieldSendIDs[i];
        updateFieldSendBuf[i] = getUpdateFieldPtr()[bid];       // getUpdateField{Send/Recv}Buff[] declared in core/Comm.h as double * . note that it is "new"
        // allocated as type GREEN_VAL, which is macro-defined as #define GREEN_VAL       double in core/Block.h
    }

    // check the buffer allocations for correct size. what about numLocalBlocks() what if this is 0?
    MPI_Allgatherv(updateFieldSendBuf,
                   numLocalBlocks(),
                   MPI_DOUBLE,
                   updateFieldRecvBuf,
                   updateFieldCounts,
                   updateFieldDisps,
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);

    // Copy the received values from the buffer to the update field
    for (i=0; i<numGlobalBlocks(); ++i) {
        bid = updateFieldRecvIDs[i];
        getUpdateFieldPtr()[bid] = updateFieldRecvBuf[i];
    }

#ifdef DEBUG
    stopTimer(dist_comm_timer);
#endif
#endif
}

/*!
 Broadcast the update field from the root node to other nodes.
 This is used for the aftershock slip adjustment calculations.
 */
void Simulation::broadcastUpdateField(void) {
#ifdef MPI_C_FOUND
#ifdef DEBUG
    startTimer(dist_comm_timer);
#endif

    MPI_Bcast(update_field, numGlobalBlocks(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef DEBUG
    stopTimer(dist_comm_timer);
#endif
#endif
}

/*!
 Distributes a list of blocks among all processors. Used for determining failed blocks in a sweep.
 */
void Simulation::distributeBlocks(const quakelib::ElementIDSet &local_id_list, BlockIDProcMapping &global_id_list) {
    quakelib::ElementIDSet::const_iterator  it;
#ifdef MPI_C_FOUND
    int                                     i, n, p;
    int                         *proc_block_count = new int[world_size];
    int                         *proc_block_disps = new int[world_size];
    BlockID                     *local_block_ids = new BlockID[local_id_list.size()];

#ifdef DEBUG
    startTimer(dist_comm_timer);
#endif

    // Let all procs know how many blocks there were on other processors
    int local_block_count = local_id_list.size();
    MPI_Allgather(&local_block_count, 1, MPI_INT, proc_block_count, 1, MPI_INT, MPI_COMM_WORLD);

#ifdef DEBUG
    stopTimer(dist_comm_timer);
#endif

    // Create list of local block IDs
    // why do we do this? aren't we just copying local_id_list[] to local_block_ids[]? can't we just use local_id_list[] ? same size and everything...
    for (i=0,it=local_id_list.begin(); it!=local_id_list.end(); ++i,++it) local_block_ids[i] = *it;

    // Count total, displacement of block IDs
    int total_blocks = 0;

    for (i=0; i<world_size; ++i) {
        proc_block_disps[i] = total_blocks;
        total_blocks += proc_block_count[i];
    }

    BlockID                     *block_ids = new BlockID[total_blocks];

#ifdef DEBUG
    startTimer(dist_comm_timer);
#endif
    MPI_Allgatherv(local_block_ids,
                   local_block_count,
                   MPI_INT,
                   block_ids,
                   proc_block_count,
                   proc_block_disps,
                   MPI_INT,
                   MPI_COMM_WORLD);
#ifdef DEBUG
    stopTimer(dist_comm_timer);
#endif
    //
    n = 0;

    for (p=0; p<world_size; ++p) {
        for (i=0; i<proc_block_count[p]; ++i) {
            // global_id_list is a collection of pairs like <block_id, node_rank_id> .
            global_id_list.insert(std::make_pair(block_ids[n], p));
            ++n;
        }
    }

    //
    // use delete [] for arrays...
    delete [] block_ids;
    delete [] proc_block_count;
    delete [] proc_block_disps;
    delete [] local_block_ids;
#else   // MPI_C_FOUND

    // Copy the local IDs into the global list just for the single processor
    for (it=local_id_list.begin(); it!=local_id_list.end(); ++it) {
        global_id_list.insert(std::make_pair(*it, 0));
    }

#endif
}

/*!
 Collect the individual event sweeps spread through all nodes
 on to the root node in a single sweep.
 */
void Simulation::collectEventSweep(quakelib::ModelSweeps &sweeps) {
#ifdef MPI_C_FOUND
    // note: this does nothing if not MPI_C_FOUND
    int                             *sweep_counts, *sweep_offsets;
    int                             num_local_sweeps, i, total_sweep_count;
    BlockSweepVals                  *sweep_vals, *all_sweeps;
    quakelib::ModelSweeps::const_iterator    it;

#ifdef DEBUG
    startTimer(sweep_comm_timer);
#endif

    // Gather the number of sweeps per node at the root
    num_local_sweeps = sweeps.size();

    if (isRootNode()) {
        sweep_counts = new int[world_size];
        sweep_offsets = new int[world_size];
    } else {
        sweep_counts = sweep_offsets = NULL;
    }

    MPI_Gather(&num_local_sweeps, 1, MPI_INT, sweep_counts, 1, MPI_INT, ROOT_NODE_RANK, MPI_COMM_WORLD);

    // Record the number of blocks the root will receive from each node
    if (isRootNode()) {
        total_sweep_count = 0;

        for (i=0; i<world_size; ++i) {
            // note: this gives the starting position of sweeps from each node (first node starts at position i=0,
            // second at i=n_0, third at i=n_1, etc.)
            sweep_offsets[i] = total_sweep_count;
            total_sweep_count += sweep_counts[i];
        }
    }

    // TODO: Change BlockSweepVals to use the Quakelib data structure
    // Record the values of each block in this sweep
    sweep_vals = new BlockSweepVals[num_local_sweeps];

    if (isRootNode()) all_sweeps = new BlockSweepVals[total_sweep_count];

    //
    for (i=0,it=sweeps.begin(); it!=sweeps.end(); ++it,++i) {
        sweep_vals[i].element_id = it->_element_id;
        sweep_vals[i].sweep_num = it->_sweep_number;
        sweep_vals[i].slip = it->_slip;
        sweep_vals[i].init_shear = it->_shear_init;
        sweep_vals[i].init_normal = it->_normal_init;
        sweep_vals[i].final_shear = it->_shear_final;
        sweep_vals[i].final_normal = it->_normal_final;
    }

    // Gather the sweep info at the root node
    // (aka, consolidate sweep_vals into all_sweeps)
    MPI_Gatherv(sweep_vals, num_local_sweeps, element_sweep_type,
                all_sweeps, sweep_counts, sweep_offsets, element_sweep_type,
                ROOT_NODE_RANK, MPI_COMM_WORLD);

    // Record the received blocks into the current sweep on the root node
    if (isRootNode()) {
        for (i=0; i<total_sweep_count; ++i) {
            BlockID         bid = all_sweeps[i].element_id;
            unsigned int    sweep_num = all_sweeps[i].sweep_num;
            sweeps.setSlipAndArea(sweep_num,
                                  bid,
                                  all_sweeps[i].slip,
                                  getBlock(bid).area(),
                                  getBlock(bid).lame_mu());
            sweeps.setInitStresses(sweep_num,
                                   bid,
                                   all_sweeps[i].init_shear,
                                   all_sweeps[i].init_normal);
            sweeps.setFinalStresses(sweep_num,
                                    bid,
                                    all_sweeps[i].final_shear,
                                    all_sweeps[i].final_normal);
        }

        // use delete [] for arrays:
        if (isRootNode()) {
            delete [] sweep_counts;
            delete [] sweep_offsets;
            delete [] all_sweeps;
        } else {
            // note that in this case, sweep_counts=NULL; sweep_offsets=NULL; and so far, all_sweeps is uninitialized or maybe populated by an MPI call?
            delete sweep_counts;
            delete sweep_offsets;
            delete all_sweeps;
        }
    }

    delete [] sweep_vals;

#ifdef DEBUG
    stopTimer(sweep_comm_timer);
#endif
#endif
}

/*!
 Partition the blocks over different nodes using a simple blocked partition scheme.
 */
void Simulation::partitionBlocks(void) {
    int                     i;
#ifdef MPI_C_FOUND
    PartitionMethod                 part_method = PARTITION_DISTANCE;
    int                             world_size, num_global_blocks, num_local_blocks, local_rank, j, n;
    std::multimap<int, BlockID>::iterator   it, it_start, it_end;
    std::set<BlockID>               cur_assigns;
    std::set<BlockID>::iterator     bit;
    std::set<BlockID>               avail_ids;
    //BlockID                         *assign_array;        // declare in line?
    int                             num_assign;
    BlockList::iterator             git;
    bool                            more_to_assign;
    FaultID                         cur_fault;
    std::multimap<double, BlockID>  dist_map;
    BlockID                         base_id;

    world_size = getWorldSize();
    num_global_blocks = numGlobalBlocks();
    num_local_blocks = num_global_blocks/world_size;
    local_rank = getNodeRank();

    //
    // Make a set of available BlockIDs
    for (git=begin(); git!=end(); ++git) avail_ids.insert(git->getBlockID());

    // Segment->Node assignment is made at the root node, then transmitted to other nodes
    // Since it is relatively simple, it won't take much time
    for (i=0; i<world_size; ++i) {
        cur_assigns.clear();

        if (isRootNode()) {
            switch (part_method) {
                case PARTITION_BLOCK:
                    if (i == world_size-1) num_local_blocks = avail_ids.size();

                    for (n=0; n<num_local_blocks; ++n) {
                        cur_assigns.insert(*(avail_ids.begin()));
                        avail_ids.erase(avail_ids.begin());
                    }

                    break;

                case PARTITION_DISTANCE:
                    base_id = UNDEFINED_ELEMENT_ID;

                    if (i == world_size-1) num_local_blocks = avail_ids.size();

                    more_to_assign = true;

                    while (more_to_assign) {
                        // Find a starting segment
                        if (base_id == UNDEFINED_ELEMENT_ID) {
                            base_id = *(avail_ids.begin());
                            cur_fault = getBlock(base_id).getFaultID();
                            dist_map.clear();

                            // Calculate the distances to all blocks on the same fault from the starting block
                            for (bit=avail_ids.begin(); bit!=avail_ids.end(); ++bit) {
                                if (getBlock(*bit).getFaultID() == cur_fault) {
                                    double      dist;
                                    dist = getBlock(*bit).center().dist(getBlock(base_id).center());
                                    dist_map.insert(std::make_pair(dist, *bit));
                                }
                            }
                        }

                        // Assign the closest element and remove it from the distance map
                        cur_assigns.insert(dist_map.begin()->second);
                        avail_ids.erase(dist_map.begin()->second);
                        dist_map.erase(dist_map.begin());

                        // If we have enough assignments for this node, quit
                        if (cur_assigns.size() >= num_local_blocks) more_to_assign = false;

                        // If we're out of blocks to assign, start another fault
                        if (dist_map.size() == 0) base_id = UNDEFINED_ELEMENT_ID;
                    }

                    break;

                default:
                    throw std::logic_error("Unknown partitioning method.");
                    break;
            }
        }

        // Send the BlockIDs to all nodes
        num_assign = cur_assigns.size();
        MPI_Bcast(&num_assign, 1, MPI_INT, ROOT_NODE_RANK, MPI_COMM_WORLD);
        //
        // try delclaring/deleting inline...
        //assign_array = new BlockID[num_assign];
        BlockID *assign_array = new BlockID[num_assign];

        //
        // Put the assigned IDs into the array on the root node
        if (isRootNode()) {
            //
            for (n=0,bit=cur_assigns.begin(); bit!=cur_assigns.end(); ++n,++bit) {
                assign_array[n] = *bit;
            }
        }


        // Broadcast the assignments to all nodes
        // what is the correct MPI data type to use here? int has size 4; MPI_INT 8.
        //MPI_Bcast(assign_array, num_assign, MPI_INT, ROOT_NODE_RANK, MPI_COMM_WORLD);
        MPI_Bcast(assign_array, num_assign, MPI_UNSIGNED, ROOT_NODE_RANK, MPI_COMM_WORLD);

        //
        for (n=0; n<num_assign; ++n) {
            if (local_rank == i) {
                local_block_ids.push_back(assign_array[n]);
                global_block_ids.insert(std::make_pair(assign_array[n], n));
            }

            block_node_map.insert(std::make_pair(assign_array[n], i));
            node_block_map.insert(std::make_pair(i, assign_array[n]));
        }

        delete [] assign_array;
    }

    if (isRootNode()) assertThrow(avail_ids.size()==0, "Did not assign all blocks in partitioning.");

    //
    // these are declared and deallocated in core/Comm.h.
    updateFieldRecvIDs = new BlockID[num_global_blocks];
    updateFieldRecvBuf = new GREEN_VAL[num_global_blocks];

    updateFieldSendIDs = new BlockID[numLocalBlocks()];
    updateFieldSendBuf = new GREEN_VAL[numLocalBlocks()];

    updateFieldCounts = new int[world_size];
    updateFieldDisps = new int[world_size];

    failBlockRecvBuf = new int[num_global_blocks];
    failBlockSendBuf = new int[numLocalBlocks()];
    failBlockCounts = new int[world_size];
    failBlockDisps = new int[world_size];

    // Get the counts of elements from each node and order them in the receive ID list
    for (j=0,i=0; i<world_size; ++i) {
        updateFieldCounts[i] = node_block_map.count(i);
        failBlockCounts[i] = node_block_map.count(i);

        it_start = node_block_map.equal_range(i).first;
        it_end = node_block_map.equal_range(i).second;

        // Figure out what IDs we will get in the receive buffer
        for (it=it_start; it!=it_end; ++it,++j) updateFieldRecvIDs[j] = it->second;

        // If we're on the local node, also record the order of IDs for the send buffer
        if (i == local_rank) {
            for (n=0,it=it_start; it!=it_end; ++it,++n) updateFieldSendIDs[n] = it->second;
        }
    }

    // Create the displacement map for receiving update field and block failure values
    updateFieldDisps[0] = failBlockDisps[0] = 0;

    for (i=1; i<world_size; ++i) {
        updateFieldDisps[i] = updateFieldDisps[i-1] + updateFieldCounts[i-1];
        failBlockDisps[i] = failBlockDisps[i-1] + failBlockCounts[i-1];
    }

#else

    for (i=0; i<numGlobalBlocks(); ++i) {
        local_block_ids.push_back(i);
        global_block_ids.insert(std::make_pair(i, i));
        block_node_map.insert(std::make_pair(i, 0));
    }

#endif
    //
    mult_buffer = NULL;
    decompress_buf = NULL;
}

// yoder:
void Simulation::setGreens(const BlockID &r, const BlockID &c, const double &new_green_shear, const double &new_green_normal) {
    // TODO:
    // yoder: add extreme greens-value checking here. see Params.h/.cpp, double getGreensShearMax(void), etc.
    // looks like we place a condition on the values of {new_green_shear, new_green_normal}
    // it may also be desirable to distinguish between self-stress (diagonal elements) and off-diagonals. however, if it comes
    // to that, we'll probably want to use Python/quakelib tools to modify a Greens output file directly and just load those values.
    //
    // update code modified to use max/min thresholds for Greens values:
    //greenShear()->setVal(getLocalInd(r), c, std::max(getGreenShearMin(), std::min(new_green_shear, getGreenShearMax())));
    //greenNormal()->setVal(getLocalInd(r), c, std::max(getGreenNormalMin(), std::min(new_green_normal, getGreenNormalMax())));
    //
    // yoder: ... and now, specify blockID values in max/min() to distinguish (off)diag greens elements.
    //
    greenShear()->setVal(getLocalInd(r), c, std::max(getGreenShearMin(r,c), std::min(new_green_shear, getGreenShearMax(r,c))));
    greenNormal()->setVal(getLocalInd(r), c, std::max(getGreenNormalMin(r,c), std::min(new_green_normal, getGreenNormalMax(r,c))));
    
    // original update code:
    /*
    //greenShear()->setVal(getLocalInd(r), c, new_green_shear);
    //greenNormal()->setVal(getLocalInd(r), c, new_green_normal);

    //
    //if (r == c) setSelfStresses(r, new_green_shear, new_green_normal);
    */
    //if (r == c) setSelfStresses(r, std::max(getGreenShearMin(), std::min(new_green_shear, getGreenShearMax())), std::max(getGreenNormalMin(), std::min(new_green_normal, getGreenNormalMax())));
    if (r == c) setSelfStresses(r, std::max(getGreenShearMin(r,c), std::min(new_green_shear, getGreenShearMax(r,c))), std::max(getGreenNormalMin(r,c), std::min(new_green_normal, getGreenNormalMax(r,c))));
};

// yoder:
void Simulation::debug_out(std::string str_in) {
    // simple debug output code; print the inout string plus the process_id, node_rank.
    printf("**Debug(%d/%d)[ev: %d]: %s", getNodeRank(), getpid(), getCurrentEvent().getEventNumber(), str_in.c_str());
};
