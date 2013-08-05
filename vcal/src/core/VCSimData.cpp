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
	int				i;
	
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
	
	shear_stress = (double*)malloc(sizeof(double)*global_size);
	assertThrow(shear_stress, "Not enough memory to allocate shear stress array.");
	normal_stress = (double*)malloc(sizeof(double)*global_size);
	assertThrow(normal_stress, "Not enough memory to allocate normal stress array.");
    f_shear_stress = (double*)malloc(sizeof(double)*global_size);
	assertThrow(f_shear_stress, "Not enough memory to allocate shear stress array.");
	f_normal_stress = (double*)malloc(sizeof(double)*global_size);
	assertThrow(f_normal_stress, "Not enough memory to allocate normal stress array.");
	update_field = (double*)malloc(sizeof(double)*global_size);
	assertThrow(update_field, "Not enough memory to allocate update field array.");
	
	for (i=0;i<global_size;++i) {
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
