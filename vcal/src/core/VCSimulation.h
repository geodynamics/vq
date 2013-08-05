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

#include "SimFramework.h"
#include "VCParams.h"
#include "VCBlock.h"
#include "VCEvent.h"
#include "VCSimData.h"
#include "VCComm.h"
#include "VCCommPartition.h"
#include "VCCommSpecExec.h"
#include "HDF5Data.h"

#ifndef _VCSIMULATION_H_
#define _VCSIMULATION_H_

#include <string>

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iomanip>

enum PartitionMethod {
	PARTITION_UNDEFINED,
	PARTITION_BLOCK,
	PARTITION_DISTANCE
};

/*!
 VCSimulation is an instantiation of a Virtual California simulation using the SimFramework.
 It contains functions related to the simulation, parameter retrieval functions and basic
 functions for earthquake analysis.
 */
class VCSimulation : public SimFramework, public VCParams, public VCSimData, public VCCommPartition, public VCCommSpecExec, public VCComm {
public:
	VCSimulation(int argc, char **argv);
	~VCSimulation(void);
	
	void init(void);
	void finish(void);

	double getYear(void) const { return year; };
	void setYear(const double &new_year) { year = new_year; };
	void incrementYear(const double &year_step) { year += year_step; };

	int numFaults(void) const;
    double shearStress(void);
    double normalStress(void);
    void getInitialFinalStresses(const BlockIDSet &block_set, double &shear_init, double &shear_final, double &normal_init, double &normal_final) const;
	//int numLayers(void) const;
	int numSurfaceBlocks(void) const;
	void getFaultNames(std::map<FaultID, std::string> &fault_names) const;
	void getBlockFaultIDs(FaultIDSet &fault_ids, const BlockIDSet &block_ids) const;
	void getFaultBlockMapping(FaultBlockMapping &fault_block_mapping, const BlockIDSet &event_blocks) const;
	void getFaultFailureAreaMapping(FaultFailureAreaMapping &fault_failure_area_mapping, const BlockIDSet &event_blocks) const;
    void getSectionBlockMapping(SectionBlockMapping &section_block_mapping, const BlockIDSet &event_blocks) const;
	
	void depthDASBounds(const BlockIDSet &block_set, double &low_depth, double &high_depth, double &low_das, double &high_das) const;
	void sumStresses(const BlockIDSet &block_set, double &shear_stress, double &shear_stress0, double &normal_stress, double &normal_stress0) const;
	bool isTopOfSlipRectangle(const BlockID &bid, const BlockIDSet &block_set);
	void getSlipRectangleBlocks(BlockIDSet &slip_rect_blocks, const BlockID &bid, const BlockIDSet &event_blocks);
	
	double sumGreenShear(const BlockID &r) const { double sum = 0; for (int i=0;i<numGlobalBlocks();++i) sum += greenShear()->val(i, r); return sum; };
	double getGreenShear(const BlockID &r, const BlockID &c) const { return greenShear()->val(getLocalInd(r), c); };
	void setGreens(const BlockID &r, const BlockID &c, const double &new_green_shear, const double &new_green_normal) {
		greenShear()->setVal(getLocalInd(r), c, new_green_shear);
		greenNormal()->setVal(getLocalInd(r), c, new_green_normal);
		if (r == c) getBlock(r).setSelfStresses(new_green_shear, new_green_normal);
	};
	double getGreenNormal(const BlockID &r, const BlockID &c) const { return greenNormal()->val(getLocalInd(r), c); };
	bool compressShearRow(const unsigned int &row, const float &ratio) { return greenShear()->compressRow(row, ratio); };
	bool compressNormalRow(const unsigned int &row, const float &ratio) { return greenNormal()->compressRow(row, ratio); };
	bool decompressShearRow(const unsigned int &row) { return greenShear()->decompressRow(row); };
	bool decompressNormalRow(const unsigned int &row) { return greenNormal()->decompressRow(row); };
	void allocateShearNormalRows(const unsigned int &row) { greenNormal()->allocateRow(row); greenShear()->allocateRow(row); };
	
    void determineBlockNeighbors(void);
    
	void computeCFFs(bool in_event);
	void printHeaders(void);
    void printAll(void);
	void printStresses(void);
	void printCFFs(void);
    void printSlipDeficits(void);
	void printShearStress(void);
    void printNormalStress(void);
    void printSlipCumulative(void);
	void matrixVectorMultiplyAccum(double *c, const quakelib::DenseMatrix<GREEN_VAL> *a, const double *b, const bool dense);
	void multiplySumRow(double *c, const double *b, const GREEN_VAL *a, const int n, const bool dense);
	void multiplyRow(double *c, const double *b, const GREEN_VAL *a, const int n);
	int distributeUpdateField(const bool &did_spec_exec);
	void distributeFailedBlocks(BlockIDSet &failed_blocks);
	void collectEventSweep(VCEventSweep &cur_sweep);
	
	void addNeighbor(const BlockID &b1, const BlockID &b2);
	std::pair<BlockIDSet::const_iterator, BlockIDSet::const_iterator> getNeighbors(const BlockID &bid) const;
	void printTimers(void);
	
	bool isLocalizedFailure(const BlockIDSet &fail_set);
	
	bool isLocalBlockID(const BlockID &block_id) const { return (block_node_map.at(block_id) == node_rank); };
	
	void partitionBlocks(void);
    
    static bool distanceCompare(BlockVal first, BlockVal second) {
        return (first.val < second.val);
    };
	
	// This stuff is so I can dump CFF data to a file. This is just a hack.
	//char block_dat_out_filename[100];
    //std::string block_dat_out_filename_str;
	std::ofstream block_dat_out_file;
	
private:
#ifdef DEBUG
	int							mult_timer;
	
	//! Counter for number of matrix multiplies performed
	int							num_mults;
#endif
	//! Current simulation year
    double						year;

	//! Temporary buffer used to speed up calculations
	double						*mult_buffer;
	GREEN_VAL					*decompress_buf;
	
	//! Map of which blocks have which neighbors
	std::map<BlockID, BlockIDSet>	neighbor_map;
	
#ifdef HAVE_PAPI_H
	int							multiply_event_set;
	long_long					total_inst;
#endif
};

#endif
