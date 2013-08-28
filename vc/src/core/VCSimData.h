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
	unsigned int			global_size, local_size;
	quakelib::DenseMatrix<GREEN_VAL>	*green_shear, *green_normal;
	
	double					*shear_stress;
	double					*normal_stress;
    double					*f_shear_stress;
	double					*f_normal_stress;
	double					*update_field;
	
public:
	VCSimData(void) : global_size(0), local_size(0), green_shear(NULL), green_normal(NULL),
		shear_stress(NULL), normal_stress(NULL), f_shear_stress(NULL), f_normal_stress(NULL), update_field(NULL) {};
	
	void setupArrays(const unsigned int &global_sys_size,
					 const unsigned int &local_sys_size,
					 const bool &compressed,
					 const bool &transposed);
	void deallocateArrays(void);
	
	unsigned int localSize(void) const { return local_size; };
	unsigned int globalSize(void) const { return global_size; };
	
	quakelib::DenseMatrix<GREEN_VAL> *greenShear(void) const { return green_shear; };
	quakelib::DenseMatrix<GREEN_VAL> *greenNormal(void) const { return green_normal; };
	double *getUpdateFieldPtr(const unsigned int &elem=0) const { return &(update_field[elem]); };
	double *getShearStressPtr(const unsigned int &elem=0) const { return &(shear_stress[elem]); };
	double *getNormalStressPtr(const unsigned int &elem=0) const { return &(normal_stress[elem]); };
    double *getFShearStressPtr(const unsigned int &elem=0) const { return &(f_shear_stress[elem]); };
	double *getFNormalStressPtr(const unsigned int &elem=0) const { return &(f_normal_stress[elem]); };
	
	double getUpdateField(const BlockID &b) const { return update_field[b]; };
	void setUpdateField(const BlockID &b, const double &new_val) { update_field[b] = new_val; };
	
	double getNormalStress(const BlockID &b) const { return normal_stress[b]; };
	void setNormalStress(const BlockID &b, const double &new_val) { normal_stress[b] = new_val; };
	void addToNormalStress(const BlockID &b, const double &add_val) { normal_stress[b] += add_val; };
	
	double getShearStress(const BlockID &b) const { return shear_stress[b]; };
	void setShearStress(const BlockID &b, const double &new_val) { shear_stress[b] = new_val; };
	void addToShearStress(const BlockID &b, const double &add_val) { shear_stress[b] += add_val; };
    
    double getFNormalStress(const BlockID &b) const { return f_normal_stress[b]; };
	void setFNormalStress(const BlockID &b, const double &new_val) { f_normal_stress[b] = new_val; };
	void addToFNormalStress(const BlockID &b, const double &add_val) { f_normal_stress[b] += add_val; };
	
	double getFShearStress(const BlockID &b) const { return f_shear_stress[b]; };
	void setFShearStress(const BlockID &b, const double &new_val) { f_shear_stress[b] = new_val; };
	void addToFShearStress(const BlockID &b, const double &add_val) { f_shear_stress[b] += add_val; };
};

#endif
