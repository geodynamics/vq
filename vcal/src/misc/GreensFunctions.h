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

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

// TODO: use specified Lame parameters (e.g. EqSim) in Greens function calculation

#ifndef _GREENS_FUNCTIONS_H_
#define _GREENS_FUNCTIONS_H_

/*!
 Class to encapsulate strain values used in Greens function calculation.
 */
class GreensValsSparseRow {
private:
	std::vector<GREEN_VAL> row_data;
public:
	GreensValsSparseRow(void) {};
	
	GreensValsSparseRow(const unsigned int &width) {
		row_data.resize(width);
	};
	
	GREEN_VAL &operator[](const unsigned int &col) throw(std::out_of_range) {
		if(col >= row_data.size()) throw std::out_of_range("GreensValsSparseRow[]");
		return row_data.at(col);
	};
	//friend std::ostream& operator<<(std::ostream& os, const GreensValsSparseRow& m);
};

// As used in the Green's function calculations, each sparse matrix will take:
// 8*N^2*(2-1/M)/M bytes, where N is the number of blocks and M is the number of CPUs
// For example, with 100,000 blocks on 64 machines, each matrix will take ~2.3 GB
// The GreensValsSparseMatrix is needed because to symmetrize the matrix we need to
// calculate (j,i) for each (i,j) the CPU will use.
// Example of sparse matrix (X=computed and used by this machine, O=computed but not used):
// |-|-|-|-|
// | |O| | |
// |X|X|X|X|
// | |O| | |
// | |O| | |
// |-|-|-|-|

class GreensValsSparseMatrix {
private:
	std::vector<GreensValsSparseRow>	matrix;
public:
	//GreensValsSparseMatrix(const GreensValsSparseMatrix &x) { assert(false); };
	
	GreensValsSparseMatrix(const std::vector<int> &row_widths) {
		std::vector<int>::const_iterator		it;
		
		// Allocate matrix rows with specified sizes
		matrix.clear();
		for (it=row_widths.begin();it!=row_widths.end();++it) {
			matrix.push_back(GreensValsSparseRow(*it));
		}
	};
	
	~GreensValsSparseMatrix(void) {
		matrix.clear();
	};
	
	GreensValsSparseRow &operator[](const unsigned int &row) throw(std::out_of_range) {
		if(row >= matrix.size()) throw std::out_of_range("GreensValsSparseMatrix[]");
		return matrix.at(row);
	};
	friend std::ostream& operator<<(std::ostream& os, const GreensValsSparseMatrix& m);
};

/*!
 Abstract class representing a generic method of calculating Greens function.
 Particular methods are instantiated as descendent classes.
 */
class GreensFuncCalc {
private:
    double			last_update;
	int				outcnt;
	
public:
	GreensFuncCalc(void) : last_update(0.0), outcnt(0) {};
	
	void progressBar(VCSimulation *sim, const int &thread_num, const int &num_done_blocks);
	
	virtual void CalculateGreens(VCSimulation *sim) = 0;
	
	void symmetrizeMatrix(VCSimulation *sim,
						  GreensValsSparseMatrix &ssh);
};

class GreensFuncFileParse : public GreensFuncCalc {
public:
	void CalculateGreens(VCSimulation *sim);
};

class GreensFuncCalcBarnesHut : public GreensFuncCalc {
public:
	void CalculateGreens(VCSimulation *sim);
	void bhInnerCalc(VCSimulation *sim, quakelib::Octree<3> *tree, const BlockID &bid);
};

class GreensFuncCalc2011 : public GreensFuncCalc {
public:
    void CalculateGreens(VCSimulation *sim);
    void InnerCalc2011(VCSimulation *sim,
                       const BlockID &bnum,
                       GreensValsSparseMatrix &ssh,
                       GreensValsSparseMatrix &snorm);
};

#endif
