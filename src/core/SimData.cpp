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

#include "SimData.h"
#include <stdlib.h>

/*!
 Allocate and initialize the arrays needed for a VC simulation.
 These include the shear and normal Greens function matrices
 and stress value arrays.
 */
void VCSimData::setupArrays(const unsigned int &global_sys_size,
                            const unsigned int &local_sys_size,
                            const bool &compressed,
                            const bool &transposed) {
    deallocateArrays();

    global_size = global_sys_size;
    // Make the local size a factor of 16 to allow SSE loop unrolling
    local_size = local_sys_size+(16-local_sys_size%16);

    if (compressed) {
        if (transposed) {
            green_shear = new quakelib::CompressedRowMatrixTranspose<GREEN_VAL>(local_size, global_size);
            green_normal = new quakelib::CompressedRowMatrixTranspose<GREEN_VAL>(local_size, global_size);
        } else {
            green_shear = new quakelib::CompressedRowMatrixStraight<GREEN_VAL>(local_size, global_size);
            green_normal = new quakelib::CompressedRowMatrixStraight<GREEN_VAL>(local_size, global_size);
        }
    } else {
        if (transposed) {
            green_shear = new quakelib::DenseStdTranspose<GREEN_VAL>(local_size, global_size);
            green_normal = new quakelib::DenseStdTranspose<GREEN_VAL>(local_size, global_size);
        } else {
            green_shear = new quakelib::DenseStdStraight<GREEN_VAL>(local_size, global_size);
            green_normal = new quakelib::DenseStdStraight<GREEN_VAL>(local_size, global_size);
        }
    }

    shear_stress = (double *)malloc(sizeof(double)*global_size);
    assertThrow(shear_stress, "Not enough memory to allocate shear stress array.");
    normal_stress = (double *)malloc(sizeof(double)*global_size);
    assertThrow(normal_stress, "Not enough memory to allocate normal stress array.");
    update_field = (double *)malloc(sizeof(double)*global_size);
    assertThrow(update_field, "Not enough memory to allocate update field array.");
    slip_deficit = (double *)malloc(sizeof(double)*global_size);
    assertThrow(slip_deficit, "Not enough memory to allocate slip deficit array.");
    rhogd = (double *)malloc(sizeof(double)*global_size);
    assertThrow(rhogd, "Not enough memory to allocate rhogd array.");
    stress_drop = (double *)malloc(sizeof(double)*global_size);
    assertThrow(stress_drop, "Not enough memory to allocate stress drop array.");
    max_stress_drop = (double *)malloc(sizeof(double)*global_size);
    assertThrow(max_stress_drop, "Not enough memory to allocate max stress drop array.");
    cff = (double *)malloc(sizeof(double)*global_size);
    assertThrow(cff, "Not enough memory to allocate cff array.");
    friction = (double *)malloc(sizeof(double)*global_size);
    assertThrow(friction, "Not enough memory to allocate friction array.");
    cff0 = (double *)malloc(sizeof(double)*global_size);
    assertThrow(cff0, "Not enough memory to allocate cff0 array.");
    self_shear = (double *)malloc(sizeof(double)*global_size);
    assertThrow(self_shear, "Not enough memory to allocate self_shear array.");
    self_normal = (double *)malloc(sizeof(double)*global_size);
    assertThrow(self_normal, "Not enough memory to allocate self_normal array.");
    shear_stress0 = (double *)malloc(sizeof(double)*global_size);
    assertThrow(shear_stress0, "Not enough memory to allocate shear_stress0 array.");
    normal_stress0 = (double *)malloc(sizeof(double)*global_size);
    assertThrow(normal_stress0, "Not enough memory to allocate normal_stress0 array.");
    dynamic_val = (double *)malloc(sizeof(double)*global_size);
    assertThrow(dynamic_val, "Not enough memory to allocate dynamic_val array.");
    failed = (bool *)malloc(sizeof(bool)*global_size);
    assertThrow(failed, "Not enough memory to allocate failed array.");

    for (BlockID i=0; i<global_size; ++i) {
        shear_stress[i] = normal_stress[i] = std::numeric_limits<float>::quiet_NaN();
        update_field[i] = slip_deficit[i] = std::numeric_limits<float>::quiet_NaN();
        rhogd[i] = stress_drop[i] = max_stress_drop[i] = cff[i] = std::numeric_limits<float>::quiet_NaN();
        friction[i] = cff0[i] = self_shear[i] = self_normal[i] = std::numeric_limits<float>::quiet_NaN();
        shear_stress0[i] = normal_stress0[i] = dynamic_val[i] = std::numeric_limits<float>::quiet_NaN();
        failed[i] = false;
    }
}

// bh_theta = 0.1

// Results for 6 layers 1 km
// No compression: 112MB in Greens
// Compression: 56MB in Greens

// Results for 12 layers 1 km
// No compression: 433MB in Greens, 306MB after
// Compression: 142MB in Greens

// Results for 24 layers 1 km
// No compression: 1689MB in Greens
// Compression: 465MB in Greens

/*!
 Deallocate previously created arrays.
 */
void VCSimData::deallocateArrays(void) {
    // note: green_shear and green_normal are container objects, not arrays, so we use delete, not delete[]
    if (green_shear) delete green_shear;

    if (green_normal) delete green_normal;

    if (shear_stress) free(shear_stress);

    if (normal_stress) free(normal_stress);

    if (update_field) free(update_field);

    if (slip_deficit) free(slip_deficit);

    if (rhogd) free(rhogd);

    if (stress_drop) free(stress_drop);

    if (max_stress_drop) free(max_stress_drop);

    if (cff) free(cff);

    if (friction) free(friction);

    if (cff0) free(cff0);

    if (self_shear) free(self_shear);

    if (self_normal) free(self_normal);

    if (shear_stress0) free(shear_stress0);

    if (normal_stress0) free(normal_stress0);

    if (dynamic_val) free(dynamic_val);

    if (failed) free(failed);

    green_shear = green_normal = NULL;
    shear_stress = normal_stress = update_field = NULL;
    slip_deficit = rhogd = stress_drop = max_stress_drop = cff = friction = NULL;
    cff0 = self_shear = self_normal = shear_stress0 = normal_stress0 = dynamic_val = NULL;
    failed = NULL;
}
