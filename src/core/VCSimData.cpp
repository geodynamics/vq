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

#include "VCSimData.h"
#include "SimError.h"
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
    int             i;

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
    f_shear_stress = (double *)malloc(sizeof(double)*global_size);
    assertThrow(f_shear_stress, "Not enough memory to allocate shear stress array.");
    f_normal_stress = (double *)malloc(sizeof(double)*global_size);
    assertThrow(f_normal_stress, "Not enough memory to allocate normal stress array.");
    update_field = (double *)malloc(sizeof(double)*global_size);
    assertThrow(update_field, "Not enough memory to allocate update field array.");

    for (i=0; i<global_size; ++i) {
        shear_stress[i] = normal_stress[i] = 0;
        f_shear_stress[i] = f_normal_stress[i] = 0;
        update_field[i] = 0;
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
    if (green_shear) delete green_shear;

    if (green_normal) delete green_normal;

    if (shear_stress) free(shear_stress);

    if (normal_stress) free(normal_stress);

    if (f_shear_stress) free(f_shear_stress);

    if (f_normal_stress) free(f_normal_stress);

    if (update_field) free(update_field);

    green_shear = green_normal = NULL;
    shear_stress = normal_stress = f_shear_stress = f_normal_stress = update_field = NULL;
}
