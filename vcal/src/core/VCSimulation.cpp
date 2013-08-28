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

#include "VCSimulation.h"
#include "SimFramework.h"
#include "SimError.h"
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <vector>

#undef HAVE_PAPI_H

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 Initialize the simulation by reading the parameter file and checking the validity of parameters.
 */
VCSimulation::VCSimulation(int argc, char **argv) : SimFramework(argc, argv) {
	srand(time(0));
	
	// TODO: more thorough argument checking
	if (argc == 1) read_params("params.d");
	else read_params(argv[argc-1]);
	
	// Check validity of parameters
	// TODO: add more checks here
	assertThrow(!getVersion().compare("2.0"),
				"sim.version: Parameter file version must be 2.0");
	assertThrow(getSimStart()>=0,
				"sim.start_year: Start year must be at least 0.");
	assertThrow(getSimStart() < getSimDuration(),
				"sim.start_year: Start year must be before end year.");
	assertThrow(getEventNoise()>=0 && getEventNoise()<=1,
				"sim.noise.event: Stress noise value must be between 0 and 1.");
	assertThrow(getSlipDeficitNoise()>=0 && getSlipDeficitNoise()<=1,
				"sim.noise.slip_deficit: Initial slip deficit noise value must be between 0 and 1.");
	assertThrow(getStressNoise()>=0 && getStressNoise()<=1,
				"sim.noise.stress: Stress drop noise value must be between 0 and 1.");
	assertThrow(getGreensCalcMethod() != GREENS_CALC_UNDEFINED,
				"Greens calculation method must be either 2011, Barnes Hut or file based.");
	assertThrow(getSpecExecMethod() != SPEC_EXEC_UNDEFINED,
				"Speculative execution method must be either none, fixed or adaptive.");
	
#ifdef HAVE_PAPI_H
	int		retval;
	
	assertThrow(PAPI_query_event (PAPI_TOT_INS) == PAPI_OK, "No instruction counter");
	multiply_event_set = PAPI_NULL;
	assertThrow(PAPI_create_eventset(&multiply_event_set) == PAPI_OK, PAPI_strerror(retval));
	assertThrow(PAPI_add_event(multiply_event_set, PAPI_TOT_INS) == PAPI_OK, PAPI_strerror(retval));
#endif
}

/*!
 Finish the simulation by deallocating memory and freeing MPI related structures.
 */
VCSimulation::~VCSimulation(void) {
#ifdef DEBUG
	console() << "Number matrix multiplies: " << num_mults << std::endl;
#endif
	if (getSpecExecMethod() != SPEC_EXEC_NONE) {
		console() << "Number of predictions: " << num_predictions << std::endl;
		console() << "Number predicted local: " << num_predicted_local << std::endl;
		console() << "Number predictions failed: " << num_predictions_failed << std::endl;
		console() << "Number predictions succeed: " << num_predictions_success << std::endl;
	}
	
	// Deallocate PAPI event set
#ifdef HAVE_PAPI_H
	console() << "Number instructions: " << total_inst << std::endl;
	PAPI_cleanup_eventset(multiply_event_set);
	PAPI_destroy_eventset(&multiply_event_set);
#endif
	
	if (mult_buffer) delete mult_buffer;
	if (decompress_buf) delete decompress_buf;
	deallocateArrays();
	
	//CFF_out_file.close();
}

/*!
 Initialize the VC specific part of the simulation by creating communication timers.
 */
void VCSimulation::init(void) {
#ifdef DEBUG
	reduce_comm_timer = initTimer("Reduce Comm", false, false);
	fail_comm_timer = initTimer("Fail Synch Comm", false, false);
	dist_comm_timer = initTimer("Distribute Comm", false, false);
	sweep_comm_timer = initTimer("Sweep Comm", false, false);
	mult_timer = initTimer("Vec-Mat Mults", false, false);
	num_mults = 0;
#endif
	SimFramework::init();
	
	//CFF output hack.
	std::string file_prepend;
	size_t pos;
	
	pos = getHDF5File().find(".");
	file_prepend = getHDF5File().substr(0,pos);
	//sprintf(CFF_out_filename, "%s_CFF.dat", file_prepend.c_str());
	//CFF_out_filename_str = CFF_out_filename;
	//CFF_out_file.open(CFF_out_filename_str.c_str());
	
	//BlockList::iterator	it;
    //CFF_out_file << "-1 -1 ";
    //for(it=begin();it!=end();++it) CFF_out_file << it->getBlockID() << " ";
	//CFF_out_file << std::endl;
}

/*!
 Calculate the number of layers in the simulation based on
 the number of unique block top depths.
 
 This is deprecated.

int VCSimulation::numLayers(void) const {
	BlockList::const_iterator	it;
	std::set<double>			depth_set;
	
	//for(it=begin();it!=end();++it) depth_set.insert(it->getTop());
	
	//return depth_set.size();
    return;
}
  */

/*!
 Calculate the number of faults in the simulation based on
 the number of unique block fault IDs.
 */
int VCSimulation::numFaults(void) const {
	BlockList::const_iterator	it;
	FaultIDSet					fault_set;
	
	for(it=begin();it!=end();++it) fault_set.insert(it->getFaultID());
	
	return fault_set.size();
}

/*!
 Calculate the total shear stress in the simulation.
 */
double VCSimulation::shearStress(void) {
	BlockList::const_iterator	it;
    double shear_stress = 0.0;
	
	for(it=begin();it!=end();++it) shear_stress += it->getShearStress();
	
	return shear_stress;
}

/*!
 Calculate the total normal stress in the simulation.
 */
double VCSimulation::normalStress(void) {
	BlockList::const_iterator	it;
    double normal_stress = 0.0;
	
	for(it=begin();it!=end();++it) normal_stress += it->getNormalStress();
	
	return normal_stress;
}

/*!
 Calculate the stress before and after an event of the blocks that failed during the event.
 */
void VCSimulation::getInitialFinalStresses(const BlockIDSet &block_set, double &shear_init, double &shear_final, double &normal_init, double &normal_final) const {
	BlockIDSet::const_iterator		it;
	
	shear_init = shear_final = normal_init = normal_final = 0.0;
	
	// For each specified block
	for (it=block_set.begin();it!=block_set.end();++it) {
		// Add the before/after stresses if it is on this node
		// Non-local blocks will have incorrect stress data
		if (isLocalToNode(*it)) {
			shear_init += getBlock(*it).getStressS0();
			shear_final += getBlock(*it).getShearStress();
			normal_init += getBlock(*it).getStressN0();
			normal_final += getBlock(*it).getNormalStress();
		}
	}
}

/*!
 Blocks must have trace flags defined for this to work
*/
int VCSimulation::numSurfaceBlocks(void) const {
	//return numGlobalBlocks()/numLayers();
    BlockList::const_iterator	it;
	BlockIDSet					block_set;
    
    for(it=begin();it!=end();++it) {
        if (it->on_trace()) {
            block_set.insert(it->getFaultID());
        }
    }
    
    return block_set.size();
}

/*!
 Get the set of FaultIDs associated with the given set of blocks.
 */
void VCSimulation::getBlockFaultIDs(FaultIDSet &fault_ids, const BlockIDSet &block_ids) const {
	BlockIDSet::const_iterator		it;
	
	for (it=block_ids.begin();it!=block_ids.end();++it) {
		fault_ids.insert(getBlock(*it).getFaultID());
	}
}

/*!
 Get a mapping of exactly which faults are associated with a specified set of blocks.
 */
void VCSimulation::getFaultBlockMapping(FaultBlockMapping &fault_block_mapping, const BlockIDSet &event_blocks) const {
	BlockIDSet::const_iterator		it;
	FaultID							fault_id;
	
	for (it=event_blocks.begin();it!=event_blocks.end();++it) {
		fault_id = getBlock(*it).getFaultID();
		fault_block_mapping[fault_id].insert(*it);
	}
}

/*!
 Get a mapping of faults to current failure area.
 */
void VCSimulation::getFaultFailureAreaMapping(FaultFailureAreaMapping &fault_failure_area_mapping, const BlockIDSet &event_blocks) const {
	BlockIDSet::const_iterator		it;
	FaultID							fault_id;
	double							block_area;
	
	for (it=event_blocks.begin();it!=event_blocks.end();++it) {
		fault_id = getBlock(*it).getFaultID();
		block_area = getBlock(*it).get_area();
		if (fault_failure_area_mapping.find(fault_id) == fault_failure_area_mapping.end() )
		{
			fault_failure_area_mapping[fault_id] = block_area;
		} else {
			fault_failure_area_mapping[fault_id] += block_area;
		}
	}
}

/*!
 Get a mapping of exactly which sections are associated with a specified set of blocks.
 */
void VCSimulation::getSectionBlockMapping(SectionBlockMapping &section_block_mapping, const BlockIDSet &event_blocks) const {
	BlockIDSet::const_iterator		it;
	quakelib::SectionID						section_id;
	
	for (it=event_blocks.begin();it!=event_blocks.end();++it) {
		section_id = getBlock(*it).getSectionID();
		section_block_mapping[section_id].insert(*it);
	}
}

/*!
 Determines the boundaries of depth and DAS in the specified set of blocks.
 Since depth is negative below ground, the lowest depth is the lowest (most negative) value.
 */
void VCSimulation::depthDASBounds(const BlockIDSet &block_set, double &low_depth, double &high_depth, double &low_das, double &high_das) const {
	BlockIDSet::const_iterator		it;
	high_das = high_depth = -DBL_MAX;
	low_das = low_depth = DBL_MAX;
	
	for (it=block_set.begin();it!=block_set.end();++it) {
		low_depth = fmin(low_depth, getBlock(*it).max_depth());
		high_depth = fmax(high_depth, getBlock(*it).min_depth());
		low_das = fmin(low_das, getBlock(*it).min_das());
		high_das = fmax(high_das, getBlock(*it).max_das());
	}
}

void VCSimulation::sumStresses(const BlockIDSet &block_set,
							   double &shear_stress,
							   double &shear_stress0,
							   double &normal_stress,
							   double &normal_stress0) const {
	BlockIDSet::const_iterator		it;
	
	shear_stress = shear_stress0 = normal_stress = normal_stress0 = 0;
	
	for (it=block_set.begin();it!=block_set.end();++it) {
		shear_stress += getShearStress(*it);
		shear_stress0 += getBlock(*it).getStressS0();
		normal_stress += getNormalStress(*it);
		normal_stress0 += getBlock(*it).getStressN0();
	}
}

/*!
 Determines whether a specified block is on top of all blocks it overlaps.
 Used for writing slip maps in EqSim output files.
 */
bool VCSimulation::isTopOfSlipRectangle(const BlockID &bid, const BlockIDSet &block_set) {
	BlockIDSet::const_iterator		it;
	Block &main_block = getBlock(bid);
	
	for (it=block_set.begin();it!=block_set.end();++it) {
		Block &target_block = getBlock(*it);
		if (!main_block.is_above(target_block) && main_block.overlaps(target_block)) {
			return false;
		}
	}
	return true;
}

void VCSimulation::getSlipRectangleBlocks(BlockIDSet &slip_rect_blocks, const BlockID &bid, const BlockIDSet &event_blocks) {
	BlockIDSet::const_iterator		it;
	Block &main_block = getBlock(bid);
	
	slip_rect_blocks.insert(bid);
	for (it=event_blocks.begin();it!=event_blocks.end();++it) {
		Block &target_block = getBlock(*it);
		if (main_block.overlaps(target_block)) slip_rect_blocks.insert(*it);
	}
}

void VCSimulation::printTimers(void) {
	if (!dry_run) printAllTimers(console(), world_size, node_rank, ROOT_NODE_RANK);
}

void VCSimulation::determineBlockNeighbors(void) {
	BlockList::iterator	bit, iit;
    BlockIDSet			all_blocks;
	double				block_size;
    
	for (bit=begin();bit!=end();++bit) {
		for (iit=begin();iit!=end();++iit) {
			block_size  = floor(bit->largest_dimension() * 1e-3) * 1e3;
			if (bit->getBlockID() != iit->getBlockID() &&			// ensure blocks are not the same
				bit->getFaultID() == iit->getFaultID() &&			// ensure blocks are on the same fault
				bit->center().dist(iit->center()) < block_size * 2.0) {	// ensure blocks are "nearby" each other
					addNeighbor(bit->getBlockID(), iit->getBlockID());
			}
		}
	}
}

void VCSimulation::addNeighbor(const BlockID &b1, const BlockID &b2) {
	neighbor_map[b1].insert(b2);
	neighbor_map[b2].insert(b1);
}

std::pair<BlockIDSet::const_iterator, BlockIDSet::const_iterator> VCSimulation::getNeighbors(const BlockID &bid) const {
	std::map<BlockID, BlockIDSet>::const_iterator		it;
	BlockIDSet		empty_set;
	
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
void VCSimulation::computeCFFs(bool in_event) {
	int			i;
	
	
	for (i=0;i<numLocalBlocks();++i) {
		getBlock(getGlobalBID(i)).calcCFF(in_event);
	}
	
	/* uncomment to have slip deficit and cff dumped to a text file.*/
	static bool inited = false;
	
	if (!inited) {
		block_dat_out_file.open("block_info.dat");
		//console() << "in_event" << " ";
		printHeaders();
		inited = true;
	}

    printAll();
	//end text dump. comment to here to stop text dumping
	
		
}

void VCSimulation::finish(void) {
	//block_dat_out_file.close();
}

void VCSimulation::printAll(void) {
	BlockList::iterator	it;
    block_dat_out_file << std::setprecision(6) << getYear() << " " << std::flush;
    for(it=begin();it!=end();++it) {
		block_dat_out_file 
			<< std::setprecision(6) << it->getSlipDeficit() << " "
			<< std::setprecision(6) << it->getCFF() << " " << std::flush;
    }
    block_dat_out_file << std::endl << std::flush;
}

void VCSimulation::printHeaders(void) {
	BlockList::iterator	it;
    block_dat_out_file << "year" << " " << std::flush;
    for(it=begin();it!=end();++it) {
		block_dat_out_file
			<< it->getBlockID() << "_slip_deficit" << " "
			<< it->getBlockID() << "_cff" << " " << std::flush;
    }
    block_dat_out_file << std::endl << std::flush;
}

void VCSimulation::printStresses(void) {
	BlockList::iterator	it;
    console() << getYear() << " " << getEventCount() << " ";
    for(it=begin();it!=end();++it) {
        if (it->getSectionID() == 13) {
        
            console() 
                << it->getShearStress() << " " << it->getNormalStress() << " " << it->getCFF() << " " << it->getFCFF() << " ";
        }
    }
    console() << std::endl;
}


void VCSimulation::printSlipDeficits(void) {
	BlockList::iterator	it;
    console() << getYear() << " ";
    for(it=begin();it!=end();++it) console() << it->getSlipDeficit() << " ";
	console() << std::endl;
}

void VCSimulation::printSlipCumulative(void) {
	BlockList::iterator	it;
    console() << getYear() << " ";
    for(it=begin();it!=end();++it) console() << it->getSlipCumulative() << " ";
	console() << std::endl;
}


void VCSimulation::printCFFs(void) {
	BlockList::iterator	it;
    //CFF_out_file << getYear() << " ";
    //for(it=begin();it!=end();++it) CFF_out_file << it->getCFF() << " ";
	//CFF_out_file << std::endl;
}

void VCSimulation::printShearStress(void) {
	BlockList::iterator	it;
    console() << getYear() << " ";
    //for(it=begin();it!=end();++it) console() << it->getShearStress() << " " << it->getFShearStress() << " ";
    for(it=begin();it!=end();++it) console() << it->getShearStress() << " ";
	console() << std::endl;
}

void VCSimulation::printNormalStress(void) {
	BlockList::iterator	it;
    console() << getYear() << " ";
    for(it=begin();it!=end();++it) console() << it->getNormalStress() << " ";
	console() << std::endl;
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

void VCSimulation::matrixVectorMultiplyAccum(double *c, const quakelib::DenseMatrix<GREEN_VAL> *a, const double *b, const bool dense) {
	int			x, y, l, width, height, array_dim;
	double		val;
#ifdef DEBUG
	if (dense) {
		num_mults++;
		startTimer(mult_timer);
	}
#endif
#ifdef HAVE_PAPI_H
	int		retval;
	retval = PAPI_start(multiply_event_set);
	assertThrow(retval == PAPI_OK, PAPI_strerror(retval));
#endif
	
	height = numLocalBlocks();
	width = numGlobalBlocks();
	array_dim = localSize();
	if (!decompress_buf) decompress_buf = (GREEN_VAL*)valloc(sizeof(GREEN_VAL)*array_dim);
	
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
		
		if (!mult_buffer) mult_buffer = (double *)valloc(sizeof(double)*array_dim);
		
		// Reset the temporary buffer
		for (x=0;x<height;++x) mult_buffer[x] = 0;
		
		// Perform the multiplication
		if (dense) {
			for (y=0;y<width;++y) {
				multiplyRow(mult_buffer, &(b[y]), a->getCol(decompress_buf, y), array_dim);
			}
		} else {
			for (y=0;y<width;++y) {
				val = b[y];
#ifndef PERFORM_SPARSE_MULTIPLIES
				if (!val) continue;
#endif
				multiplyRow(mult_buffer, &val, a->getCol(decompress_buf, y), array_dim);
			}
		}
		
		// Add the temporary buffer values into the result array
		for (x=0;x<height;++x) {
			l = getGlobalBID(x);
			c[l] += mult_buffer[x];
		}
	} else {
		for (x=0;x<height;++x) {
			val = 0;
			multiplySumRow(&val, b, a->getRow(decompress_buf, x), width, dense);
			c[x] += val;
		}
	}
	
#ifdef HAVE_PAPI_H
	long_long	num_inst;
	retval = PAPI_stop(multiply_event_set, &num_inst);
	//assertMsg(retval == PAPI_OK, PAPI_strerror(retval));
	total_inst += num_inst;
#endif
	
#ifdef DEBUG
	if (dense) stopTimer(mult_timer);
#endif
}

// SSE old execution notes
// Original speed (no OpenMP)
//	41.5 events per second
//  3662 MB/sec bandwidth
//  914 Mflop/sec
// SSE
//	41.3 events per second
//  3829 MB/sec bandwidth
//  955 Mflop/sec
//	4.5% improvement over original
// SSE loop unroll without prefetch
//  43.4 events per second
//  4059 MB/sec bandwidth
//  1013 Mflop/sec
//	10.8% improvement over original
// SSE loop unroll with a* prefetch
//  46.4 events per second
//  4316 MB/sec bandwidth
//  1077 Mflop/sec
//	17.8% improvement over original

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
#define SSE_LOOP_UNROLL		8

#ifdef USE_SSE	// SSE version
#include <xmmintrin.h>

void VCSimulation::multiplySumRow(double *c, const double *b, const GREEN_VAL *a, const int n, const bool dense) {
	__m128d		aval, bval, cval, tmpval;
	double		tmp[2];
#if SSE_LOOP_UNROLL != 16 && SSE_LOOP_UNROLL != 8 && SSE_LOOP_UNROLL != 6 && SSE_LOOP_UNROLL != 4 && SSE_LOOP_UNROLL != 2
#error "Invalid value of SSE_LOOP_UNROLL"
#endif
	cval = _mm_setzero_pd();
	for (int x=0;x<n;x+=SSE_LOOP_UNROLL) {
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

#else	// Non-SSE version

void VCSimulation::multiplySumRow(double *c, const double *b, const GREEN_VAL *a, const int n, const bool dense) {
	double val = 0;
	if (dense) {
		for (int x=0;x<n;++x) {
			val += a[x]*b[x];
		}
	} else {
		for (int x=0;x<n;++x) {
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
void VCSimulation::multiplyRow(double *c, const double *b, const GREEN_VAL *a, const int n) {
	__m128d		aval, bval, cval, tmpval;
#if SSE_LOOP_UNROLL != 16 && SSE_LOOP_UNROLL != 8 && SSE_LOOP_UNROLL != 6 && SSE_LOOP_UNROLL != 4 && SSE_LOOP_UNROLL != 2
#error "Invalid value of SSE_LOOP_UNROLL"
#endif
	bval = _mm_load1_pd(b);
	for (int x=0;x<n;x+=SSE_LOOP_UNROLL) {
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
void VCSimulation::multiplyRow(double *c, const double *b, const GREEN_VAL *a, const int n) {
	for (int x=0;x<n;++x) c[x] += b[0]*a[x];
}
#endif

/*!
 Distributes the local part of the update field to other nodes and
 receives their local update fields.
 */
int VCSimulation::distributeUpdateField(const bool &did_spec_exec) {
	double		num_spec_exec = 0;
	
#ifdef MPI_C_FOUND
#ifdef DEBUG
	startTimer(dist_comm_timer);
#endif
	int		i;
	BlockID	bid;
	// Copy the local update field values to the send buffer
	for (i=0;i<numLocalBlocks();++i) {
		bid = updateFieldSendIDs[i];
		updateFieldSendBuf[i] = getUpdateFieldPtr()[bid];
	}
	updateFieldSendBuf[i] = (did_spec_exec ? 1.0 : 0.0);
	MPI_Allgatherv(updateFieldSendBuf,
				   numLocalBlocks()+1,
				   MPI_DOUBLE,
				   updateFieldRecvBuf,
				   updateFieldCounts,
				   updateFieldDisps,
				   MPI_DOUBLE,
				   MPI_COMM_WORLD);
	// Copy the received values from the buffer to the update field
	for (i=0;i<numGlobalBlocks()+getWorldSize();++i) {
		bid = updateFieldRecvIDs[i];
		if (bid != UINT_MAX) {
			getUpdateFieldPtr()[bid] = updateFieldRecvBuf[i];
		} else {
			num_spec_exec += updateFieldRecvBuf[i];
		}
	}
#ifdef DEBUG
	stopTimer(dist_comm_timer);
#endif
#endif
	return num_spec_exec;
}

/*!
 Distributes a bit field of which blocks failed in the current sweep.
 Needed for slip calculation if faults are spread over multiple processors.
 */
void VCSimulation::distributeFailedBlocks(BlockIDSet &failed_blocks) {
#ifdef MPI_C_FOUND
	int							i;
	BlockIDSet::iterator		it;
	
	for (i=0;i<numGlobalBlocks();++i) failBlockRecvBuf[i] = 0;
	for (i=0;i<numLocalBlocks();++i) failBlockSendBuf[i] = -1;
	for (i=0,it=failed_blocks.begin();it!=failed_blocks.end();++i,++it) {
		if (isLocalBlockID(*it)) {
			failBlockSendBuf[i] = *it;
		}
	}
	
#ifdef DEBUG
	startTimer(dist_comm_timer);
#endif
	MPI_Allgatherv(failBlockSendBuf,
				   numLocalBlocks(),
				   MPI_INT,
				   failBlockRecvBuf,
				   failBlockCounts,
				   failBlockDisps,
				   MPI_INT,
				   MPI_COMM_WORLD);
#ifdef DEBUG
	stopTimer(dist_comm_timer);
#endif
	for (i=0;i<numGlobalBlocks();++i) {
		if (failBlockRecvBuf[i] >= 0) {
			failed_blocks.insert(failBlockRecvBuf[i]);
		}
	}
#endif
}

/*!
 Try to predict whether a rupture on the given faults will propagate
 to another node during this sweep.
 */
bool VCSimulation::isLocalizedFailure(const BlockIDSet &fail_set) {
	BlockIDSet::const_iterator		it;
	BlockList::iterator				bit;
	double							min_dist, fixed_dist;
	int								i;
	bool							predict_local;
	
	// If we're using the fixed distance method, calculate the boundary distances
	if ((getSpecExecMethod() == SPEC_EXEC_FIXED_DIST || getSpecExecMethod() == SPEC_EXEC_ADAPTIVE) && boundary_dist_map.empty()) {
		// For each local block
		for (i=0;i<numLocalBlocks();++i) {
			Block		&local_block = getBlock(getGlobalBID(i));
			min_dist = DBL_MAX;
			// Compute the distance to the boundary over all global blocks
			for (bit=begin();bit!=end();++bit) {
				// For non-local blocks, compute the minimum distance
				if (!isLocalBlockID(bit->getBlockID())) {
					min_dist = fmin(min_dist, bit->center().dist(local_block.center()));
				}
			}
			boundary_dist_map.insert(std::make_pair(getGlobalBID(i), min_dist));
		}
	}
	
	// If our last prediction failed, just do the normal
	// execution method for the rest of the event
	if (last_prediction_failed) return false;
	
	if (fail_set.empty()) return false;
	
	num_predictions++;
	
	switch (getSpecExecMethod()) {
		case SPEC_EXEC_NONE:
			predict_local = false;
			break;
		case SPEC_EXEC_FIXED_DIST:
		case SPEC_EXEC_ADAPTIVE:
			fixed_dist = getSpecExecDistance();
			predict_local = true;
			for (it=fail_set.begin();it!=fail_set.end();++it) {
				if (boundary_dist_map.at(*it) < fixed_dist) {
					predict_local = false;
					break;
				}
			}
			break;
		default:
			predict_local = false;
			break;
	}
	
	if (predict_local) num_predicted_local++;
	
	return predict_local;
}

/*!
 Collect the individual event sweeps spread through all nodes
 on to the root node in a single sweep.
 NOTE: This does not transfer stresses, these will be zeroed in the output file
 */
void VCSimulation::collectEventSweep(VCEventSweep &cur_sweep) {
#ifdef MPI_C_FOUND
	int								*block_counts, *block_offsets;
	int								num_blocks, i, total_block_count;
	VCEventSweep::const_iterator	it;
	BlockSweepVals					*sweep_blocks, *all_blocks;
	
#ifdef DEBUG
	startTimer(sweep_comm_timer);
#endif
	
	// Gather the number of blocks per node at the root
	num_blocks = cur_sweep.size();
	if (isRootNode()) {
		block_counts = new int[world_size];
		block_offsets = new int[world_size];
	}
	MPI_Gather(&num_blocks, 1, MPI_INT, block_counts, 1, MPI_INT, ROOT_NODE_RANK, MPI_COMM_WORLD);
	
	// Record the number of blocks the root will receive from each node
	if (isRootNode()) {
		total_block_count = 0;
		for (i=0;i<world_size;++i) {
			block_offsets[i] = total_block_count;
			total_block_count += block_counts[i];
		}
	}
	
	// Record the values of each block in this sweep
	sweep_blocks = new BlockSweepVals[num_blocks];
	if (isRootNode()) all_blocks = new BlockSweepVals[total_block_count];
	for (i=0,it=cur_sweep.begin();it!=cur_sweep.end();++it,++i) {
		sweep_blocks[i].block_id = it->first;
		sweep_blocks[i].slip = it->second.slip;
		sweep_blocks[i].init_shear = it->second.shear_init;
		sweep_blocks[i].init_normal = it->second.normal_init;
		sweep_blocks[i].final_shear = it->second.shear_final;
		sweep_blocks[i].final_normal = it->second.normal_final;
	}
	
	// Gather the sweep info at the root node
	MPI_Gatherv(sweep_blocks, num_blocks, block_sweep_type,
				all_blocks, block_counts, block_offsets, block_sweep_type,
				ROOT_NODE_RANK, MPI_COMM_WORLD);
	
	// Record the received blocks into the current sweep on the root node
	if (isRootNode()) {
		for (i=0;i<total_block_count;++i) {
			BlockID		bid;
			bid = all_blocks[i].block_id;
			cur_sweep.setSlipAndArea(bid,
									 all_blocks[i].slip,
									 getBlock(bid).get_area(),
									 getBlock(bid).getMu());
			cur_sweep.setInitStresses(bid,
									all_blocks[i].init_shear,
									all_blocks[i].init_normal);
			cur_sweep.setFinalStresses(bid,
									 all_blocks[i].final_shear,
									 all_blocks[i].final_normal);
		}
		delete block_counts;
		delete block_offsets;
		delete all_blocks;
	}
	
	delete sweep_blocks;
	
#ifdef DEBUG
	stopTimer(sweep_comm_timer);
#endif
#endif
}

/*!
 Partition the blocks over different nodes using a simple blocked partition scheme.
 */
void VCSimulation::partitionBlocks(void) {
	int						i;
#ifdef MPI_C_FOUND
	PartitionMethod					part_method = PARTITION_DISTANCE;
	int								world_size, num_global_blocks, num_local_blocks, local_rank, j, n;
	std::multimap<int, BlockID>::iterator	it, it_start, it_end;
	std::set<BlockID>				cur_assigns;
	std::set<BlockID>::iterator		bit;
	std::set<BlockID>				avail_ids;
	BlockID							*assign_array;
	int								num_assign;
	BlockList::iterator				git;
	bool							more_to_assign;
	FaultID							cur_fault;
	std::multimap<double, BlockID>	dist_map;
	BlockID							base_id;
	
	world_size = getWorldSize();
	num_global_blocks = numGlobalBlocks();
	num_local_blocks = num_global_blocks/world_size;
	local_rank = getNodeRank();
	
	// Make a set of available BlockIDs
	for (git=begin();git!=end();++git) avail_ids.insert(git->getBlockID());
	
	// Segment->Node assignment is made at the root node, then transmitted to other nodes
	// Since it is relatively simple, it won't take much time
	for (i=0;i<world_size;++i) {
		cur_assigns.clear();
		if (isRootNode()) {
			switch (part_method) {
				case PARTITION_BLOCK:
					if (i == world_size-1) num_local_blocks = avail_ids.size();
					for (n=0;n<num_local_blocks;++n) {
						cur_assigns.insert(*(avail_ids.begin()));
						avail_ids.erase(avail_ids.begin());
					}
					break;
				case PARTITION_DISTANCE:
					base_id = UNDEFINED_BLOCK_ID;
					if (i == world_size-1) num_local_blocks = avail_ids.size();
					more_to_assign = true;
					while (more_to_assign) {
						// Find a starting segment
						if (base_id == UNDEFINED_BLOCK_ID) {
							base_id = *(avail_ids.begin());
							cur_fault = getBlock(base_id).getFaultID();
							dist_map.clear();
							
							// Calculate the distances to all blocks on the same fault from the starting block
							for (bit=avail_ids.begin();bit!=avail_ids.end();++bit) {
								if (getBlock(*bit).getFaultID() == cur_fault) {
									double		dist;
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
						if (dist_map.size() == 0) base_id = UNDEFINED_BLOCK_ID;
					}
					break;
				default:
					throw VCException("Unknown partitioning method.");
					break;
			}
		}
		
		// Send the BlockIDs to all nodes
		num_assign = cur_assigns.size();
		MPI_Bcast(&num_assign, 1, MPI_INT, ROOT_NODE_RANK, MPI_COMM_WORLD);
		assign_array = new BlockID[num_assign];
		
		// Put the assigned IDs into the array on the root node
		if (isRootNode()) {
			for (n=0,bit=cur_assigns.begin();bit!=cur_assigns.end();++n,++bit) {
				assign_array[n] = *bit;
			}
		}
		
		// Broadcast the assignments to all nodes
		MPI_Bcast(assign_array, num_assign, MPI_INT, ROOT_NODE_RANK, MPI_COMM_WORLD);
		for (n=0;n<num_assign;++n) {
			if (local_rank == i) {
				local_block_ids.push_back(assign_array[n]);
				global_block_ids.insert(std::make_pair(assign_array[n], n));
			}
			block_node_map.insert(std::make_pair(assign_array[n], i));
			node_block_map.insert(std::make_pair(i, assign_array[n]));
		}
		delete assign_array;
	}
	if (isRootNode()) assertThrow(avail_ids.size()==0, "Did not assign all blocks in partitioning.");
	
	// One extra field for speculative execution information
	updateFieldRecvIDs = new BlockID[num_global_blocks+world_size];
	updateFieldRecvBuf = new GREEN_VAL[num_global_blocks+world_size];
	
	updateFieldSendIDs = new BlockID[numLocalBlocks()+1];
	updateFieldSendBuf = new GREEN_VAL[numLocalBlocks()+1];
	
	updateFieldCounts = new int[world_size];
	updateFieldDisps = new int[world_size];
	
	failBlockRecvBuf = new int[num_global_blocks];
	failBlockSendBuf = new int[numLocalBlocks()];
	failBlockCounts = new int[world_size];
	failBlockDisps = new int[world_size];
	
	// Get the counts of elements from each node and order them in the receive ID list
	for (j=0,i=0;i<world_size;++i) {
		updateFieldCounts[i] = node_block_map.count(i)+1;
		failBlockCounts[i] = node_block_map.count(i);
		
		it_start = node_block_map.equal_range(i).first;
		it_end = node_block_map.equal_range(i).second;
		// Figure out what IDs we will get in the receive buffer
		for (it=it_start;it!=it_end;++it,++j) updateFieldRecvIDs[j] = it->second;
		updateFieldRecvIDs[j++] = UINT_MAX;
		// If we're on the local node, also record the order of IDs for the send buffer
		if (i == local_rank) {
			for (n=0,it=it_start;it!=it_end;++it,++n) updateFieldSendIDs[n] = it->second;
			updateFieldSendIDs[n] = UINT_MAX;
		}
	}
	
	// Create the displacement map for receiving update field and block failure values
	updateFieldDisps[0] = failBlockDisps[0] = 0;
	for (i=1;i<world_size;++i) {
		updateFieldDisps[i] = updateFieldDisps[i-1] + updateFieldCounts[i-1];
		failBlockDisps[i] = failBlockDisps[i-1] + failBlockCounts[i-1];
	}
#else
	for (i=0;i<numGlobalBlocks();++i) {
		local_block_ids.push_back(i);
		global_block_ids.insert(std::make_pair(i, i));
		block_node_map.insert(std::make_pair(i, 0));
	}
#endif
	
	mult_buffer = NULL;
	decompress_buf = NULL;
}
