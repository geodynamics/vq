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

#include "VCSimDataBlocks.h"
#include "VCSimDataEvents.h"

#ifndef _VCSIM_DATA_H_
#define _VCSIM_DATA_H_

/*!
 VCSimData stores shared simulation data in a commonly accessible structure.
 The Greens matrices are stored in a transpose format where successive entries
 are successive matrix rows rather than columns. This greatly increases speed
 for the matrix-vector multiplication. However, it also means you must be
 careful when accessing these matrices directly since reading element (x, y)
 in the array will actually be reading element (y, x) in the Greens matrix.
 For parallel computations only a segment of the matrix is stored which makes
 access even more complicated. Unless speed is required, it is recommended to
 use the VCSimData accessor functions to read/write the Greens matrices.
 */
class VCSimData : public VCSimDataBlocks, public VCSimDataEvents {
    private:
        unsigned int            global_size, local_size;
        quakelib::DenseMatrix<GREEN_VAL>    *green_shear, *green_normal;

        double                  *shear_stress;
        double                  *normal_stress;
        double                  *update_field;

    public:
        VCSimData(void) : global_size(0), local_size(0), green_shear(NULL), green_normal(NULL),
            shear_stress(NULL), normal_stress(NULL), update_field(NULL) {};

        void setupArrays(const unsigned int &global_sys_size,
                         const unsigned int &local_sys_size,
                         const bool &compressed,
                         const bool &transposed);
        void deallocateArrays(void);

        unsigned int localSize(void) const {
            return local_size;
        };
        unsigned int globalSize(void) const {
            return global_size;
        };

        quakelib::DenseMatrix<GREEN_VAL> *greenShear(void) const {
            return green_shear;
        };
        quakelib::DenseMatrix<GREEN_VAL> *greenNormal(void) const {
            return green_normal;
        };
        double *getUpdateFieldPtr(const unsigned int &elem=0) const {
            return &(update_field[elem]);
        };
        double *getShearStressPtr(const unsigned int &elem=0) const {
            return &(shear_stress[elem]);
        };
        double *getNormalStressPtr(const unsigned int &elem=0) const {
            return &(normal_stress[elem]);
        };

        double getUpdateField(const BlockID &b) const {
            return update_field[b];
        };
        void setUpdateField(const BlockID &b, const double &new_val) {
            update_field[b] = new_val;
        };

        double getNormalStress(const BlockID &b) const {
            return normal_stress[b];
        };
        void setNormalStress(const BlockID &b, const double &new_val) {
            normal_stress[b] = new_val;
        };
        void addToNormalStress(const BlockID &b, const double &add_val) {
            normal_stress[b] += add_val;
        };

        double getShearStress(const BlockID &b) const {
            return shear_stress[b];
        };
        void setShearStress(const BlockID &b, const double &new_val) {
            shear_stress[b] = new_val;
        };
        void addToShearStress(const BlockID &b, const double &add_val) {
            shear_stress[b] += add_val;
        };
};

#endif
