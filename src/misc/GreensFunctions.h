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

#include "VCSimulation.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef VC_HAVE_STDLIB_H
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
            if (col >= row_data.size()) throw std::out_of_range("GreensValsSparseRow[]");

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
        std::vector<GreensValsSparseRow>    matrix;
    public:
        //GreensValsSparseMatrix(const GreensValsSparseMatrix &x) { assert(false); };

        GreensValsSparseMatrix(const std::vector<int> &row_widths) {
            std::vector<int>::const_iterator        it;

            // Allocate matrix rows with specified sizes
            matrix.clear();

            for (it=row_widths.begin(); it!=row_widths.end(); ++it) {
                matrix.push_back(GreensValsSparseRow(*it));
            }
        };

        ~GreensValsSparseMatrix(void) {
            matrix.clear();
        };

        GreensValsSparseRow &operator[](const unsigned int &row) throw(std::out_of_range) {
            if (row >= matrix.size()) throw std::out_of_range("GreensValsSparseMatrix[]");

            return matrix.at(row);
        };
        friend std::ostream &operator<<(std::ostream &os, const GreensValsSparseMatrix &m);
};

/*!
 Abstract class representing a generic method of calculating Greens function.
 Particular methods are instantiated as descendent classes.
 */
class GreensFuncCalc {
    private:
        double          last_update;
        int             outcnt;

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

class GreensFuncCalcStandard : public GreensFuncCalc {
    public:
        void CalculateGreens(VCSimulation *sim);
        void InnerCalcStandard(VCSimulation *sim,
                               const BlockID &bnum,
                               GreensValsSparseMatrix &ssh,
                               GreensValsSparseMatrix &snorm);
};

#endif
