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

#include "VCBlock.h"
#include "QuakeLibOkada.h"
#include <math.h>
#include <sstream>
#include <string.h>
#include <iomanip>
#include <list>

#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif

#ifdef HAVE_FLOAT_H
#include <float.h>
#endif

#define IM1 2147483563
#define IM2 2147483399

#define IA1 40014
#define IA2 40692

#define IQ1 53668
#define IQ2 52774

#define IR1 12211
#define IR2 3791

#define IMM1 (IM1-1)
#define NDIV (1+IMM1/NTAB)
#define NORM 1./IM1

/*!
 New methods defined for the updated okada functions
 */

#define TRIG_TOLERANCE 0.0001

void Block::get_rake_and_normal_stress_due_to_block(double stresses[2], const Block &source_block) const {
    // test block vectors. This block is the test block
	quakelib::Vec<3> rake_vec, normal_vec, center_vec;
    // source block vectors
    quakelib::Vec<3> source_normal_vec, source_center_vec;
    // modified test vectors
    quakelib::Vec<3> mrake_vec, mnormal_vec, mcenter_vec;
    // other stuff
    quakelib::Vec<3> rot_axis, stress_vec, shift_vec, xy_projected_source_normal;
	quakelib::Tensor<3,3> stress_tensor;
    double theta;
    
    normal_vec = normal();
    rake_vec = rake_vector();
    center_vec = center();
    
    source_normal_vec = source_block.normal();
    source_center_vec = source_block.center();
    
    // first shift all points
    shift_vec = quakelib::Vec<3>(source_block._vert[1][0],source_block._vert[1][1],0.0);
    mcenter_vec = center_vec - shift_vec;
    
    // now we need to perform a 2d rotation in the x-y plane so that the new x axis
    // aligns with the pt3-pt2 vector of the source block. this is okada's coord sys
    xy_projected_source_normal[0] = source_normal_vec[0];
    xy_projected_source_normal[1] = source_normal_vec[1];
    xy_projected_source_normal[2] = 0.0;
    
    theta = xy_projected_source_normal.vector_angle(quakelib::Vec<3>(0.0, -1.0, 0.0));
    
    if (normal_vec[0] >= 0.0) {
        rot_axis = quakelib::Vec<3>(0.0, 0.0, 1.0);
    } else {
        rot_axis = quakelib::Vec<3>(0.0, 0.0, -1.0);
    }
    
    mnormal_vec = normal_vec. rotate_around_axis(rot_axis, theta);
    mrake_vec   = rake_vec.   rotate_around_axis(rot_axis, theta);
    mcenter_vec = mcenter_vec.rotate_around_axis(rot_axis, theta);
    
    stress_tensor = source_block.calc_stress_tensor(mcenter_vec, unit_slip, getLambda(), getMu());
    
	stress_vec = stress_tensor*mnormal_vec;
    
    stresses[0] = stress_vec.dot_product(mrake_vec);
    stresses[1] = stress_vec.dot_product(mnormal_vec);
}

/*!*/

VCRand::VCRand(const VCRand::State& state) {
    idum1 = state.idum1;
    idum2 = state.idum2;
    iy    = state.iy;
    for(int i=0; i<NTAB; i++) iv[i] = state.iv[i];
}

void VCRand::save(VCRand::State& state) {
    state.idum1 = idum1;
    state.idum2 = idum2;
    state.iy    = iy;
    for(int i=0; i<NTAB; i++) state.iv[i] = iv[i];
}

void VCRand::init(int seed) {
    idum2 = idum1 = (seed < 1 ? 1 : seed);
    for (int j=NTAB+7; j>=0; j--) {
        int k = idum1/IQ1;
        idum1  = IA1*(idum1 - k*IQ1) - k*IR1;
        if(idum1 < 0) idum1 += IM1;
        if(j < NTAB) iv[j] = idum1;
    }
    iy=iv[0];
}

void VCRand::save(int* state) {
    state[0] = idum1;
    state[1] = idum2;
    state[2] = iy;
	
    int shift = 3;
    for(int i=0; i<VCRand::NTAB; i++) state[shift+i] = iv[i];
}

void VCRand::init(int* state) {
    idum1 = state[0];
    idum2 = state[1];
    iy    = state[2];
    
    int shift = 3;
    for(int i=0; i<VCRand::NTAB; i++) iv[i] = state[shift+i];
}

double VCRand::nextDouble() {
    int k = idum1/IQ1;
    idum1  = IA1*(idum1 - k*IQ1) - k*IR1;
    if(idum1 < 0) idum1 += IM1;
	
    k = idum2/IQ2;
    idum2  = IA2*(idum2 - k*IQ2) - k*IR2;
    if(idum2 < 0) idum2 += IM2;
	
    int tmp = iy/NDIV;
    iy = iv[tmp] - idum2;
    iv[tmp] = idum1;
	
    if(iy < 1) iy += IMM1;
	
    if(iy > IMM1) return 0.99999999;
    else          return iy*NORM;
}

std::ostream& operator<<(std::ostream& os, const BlockVal& bv) {
	if (bv.block_id == UNDEFINED_BLOCK_ID) os << "(Undef, ";
	else os << "(B" << bv.block_id << ", ";
	os << bv.val << ")";
	return os;
}

std::ostream& operator<<(std::ostream& os, const Block& b) {
	quakelib::Conversion		convert;
	quakelib::Vec<3>		center = b.center();
	
	os  << std::setw(7) << b.getBlockID()
		<< std::setw(5) << b.getFaultID()
        << std::setw(8) << std::setprecision(4) << center[0]
		<< std::setw(8) << std::setprecision(4) << center[1]
        << std::setw(8) << std::setprecision(4) << center[2]
		<< std::setw(12) << std::setprecision(8) << b.getRecurrence()
		<< std::setw(11) << std::setprecision(6) << b.slip_rate()
		<< std::setw(13) << std::setprecision(9) << b.getStressDrop()
		<< std::setw(9) << std::setprecision(5) << b.friction()
		<< std::setw(12) << std::setprecision(5) << b.getSelfStresses()
		<< std::setw(9) << std::setprecision(5) << b.getDynamicVal()
		<< std::setw(7) << std::setprecision(4) << b.aseismic()
		<< std::setw(8) << std::setprecision(4) << convert.rad2deg(b.rake())
		<< std::setw(8) << std::setprecision(4) << convert.rad2deg(b.dip())
        << std::setw(8) << std::setprecision(4) << b.getLambda()
        << std::setw(8) << std::setprecision(4) << b.getMu()
		//<< std::setw(9) << std::setprecision(4) << convert.rad2deg(b.getDipDir())
		//<< std::setw(8) << std::setprecision(4) << b.getMaxDepth()
		//<< std::setw(8) << std::setprecision(4) << b.getMinDepth()
		//<< std::setw(11) << std::setprecision(7) << b.getWestVector().X()
		//<< std::setw(11) << std::setprecision(7) << b.getWestVector().Y()
		//<< std::setw(11) << std::setprecision(7) << b.getEastVector().X()
		//<< std::setw(11) << std::setprecision(7) << b.getEastVector().Y()
		<< " " << b.getFaultName();
	return os;
}

std::string Block::header(void) const {
	std::stringstream		ss;
	ss	<< std::setw(6) << "id"
		<< std::setw(5) << "fid"
        << std::setw(5) << "center x"
        << std::setw(5) << "center y"
        << std::setw(5) << "center z"
		<< std::setw(12) << "recurrence"
		<< std::setw(11) << "slip_rate"
		<< std::setw(13) << "stress_drop"
		<< std::setw(9) << "friction"
		<< std::setw(12) << "self_stress"
		<< std::setw(9) << "dynamic"
		<< std::setw(7) << "aseis"
		<< std::setw(8) << "rake"
		<< std::setw(8) << "dip"
		//<< std::setw(9) << "dip_dir"
		//<< std::setw(8) << "bottom"
		//<< std::setw(8) << "top"
		//<< std::setw(11) << "west_x"
		//<< std::setw(11) << "west_y"
		//<< std::setw(11) << "east_x"
		//<< std::setw(11) << "east_y"
		<< " name";
	return ss.str();
}

std::string Block::headerUnits(void) const {
	std::stringstream		ss;
	ss	<< std::setw(6) << "" 
		<< std::setw(5) << ""
        << std::setw(5) << "(meters)"
        << std::setw(5) << "(meters)"
        << std::setw(5) << "(meters)"
		<< std::setw(12) << "(seconds)"
		<< std::setw(11) << "(meters/second)"
		<< std::setw(13) << "(Pascals)"
		<< std::setw(9) << ""
		<< std::setw(12) << "(Pascals)"
		<< std::setw(9) << ""
		<< std::setw(7) << ""
		<< std::setw(8) << "(deg)"
		<< std::setw(8) << "(deg)";
		//<< std::setw(9) << "(deg)"
		//<< std::setw(8) << "(km)"
		//<< std::setw(8) << "(km)"
		//<< std::setw(11) << "(km)"
		//<< std::setw(11) << "(km)"
		//<< std::setw(11) << "(km)"
		//<< std::setw(11) << "(km)";
	return ss.str();
}

void Block::clear(void) {
	quakelib::Element<4>::clear();
	is_valid = false;
	id = UNDEFINED_BLOCK_ID;
	fid = UNDEFINED_FAULT_ID;
    sid = UNDEFINED_SECTION_ID;
	self_shear = self_normal = dynamic_val = stress_drop = init_shear_stress = 0;
    init_normal_stress = rhogd = friction_val = section_layers = lambda = mu = unit_slip = 0;
	fault_name.clear();
    active_dynamic_val = 1.0;
}

void Block::calcFriction(void) {
    friction_val = (fabs(getStressDrop())/getRhogd());
}

bool Block::cffFailure(void) const {
	return (state.cff >= 0);
}

bool Block::dynamicFailure(const FaultID &event_fault) const {
        return (fid == event_fault &&
            state.cff > state.cff0 &&
			state.cff > getStressDrop() &&
			(state.cff0-state.cff)/state.cff0 > active_dynamic_val);
}

bool Block::additionalFailure(void) const {
	return (state.stressS[0] > 0.0 && (state.cff-state.Fcff)/state.cff > dynamic_val);
}

void Block::calcCFF(bool in_event) {
    if (in_event) {
        if (failed) {
            if (failed_this_sweep) {
				if (state.FstressS[0] > state.stressS[0] ) {
					state.stressS[0] =  2.0 * state.stressS[0] - state.FstressS[0];
				} else {
					state.stressS[0] = state.FstressS[0];
				}
                //double new_stress = state.FstressS[0];//state.stressS[0] - (state.FstressS0 - state.FstressS[0]);
                if (state.stressS[0] < 0) {
                    state.stressS[0] = 0;
                }
				
				

                ////state.stressS[0] -= state.FstressS0 - state.FstressS[0];
                //state.stressN[0] -= state.FstressN0 - state.FstressN[0];
                //state.FstressS[0] = state.stressS[0];
                //state.FstressN[0] = state.stressN[0];
				state.stressN[0] = state.FstressN[0];
				state.FstressS[0] = state.stressS[0];
                state.cff = state.stressS[0] - friction()*state.stressN[0];
                failed_this_sweep = false;
            }
			
			state.Fcff = state.FstressS[0] - friction()*state.FstressN[0];
			
			//if (state.Fcff > 0 ) {
			//	state.Fcff = 0;
			//}
			

        }
        else {
            state.stressS[0] = state.FstressS[0];
            state.stressN[0] = state.FstressN[0];
            state.cff = state.stressS[0] - friction()*state.stressN[0];
            state.Fcff = state.cff;
        }        
    }
    else {
        state.FstressS[0] = state.stressS[0];
        state.FstressN[0] = state.stressN[0];
        state.cff = state.stressS[0] - friction()*state.stressN[0];
        state.Fcff = state.FstressS[0] - friction()*state.FstressN[0];
    }
	
    //state.cff = state.stressS[0] - stress_drop;
}

void Block::setStatePtrs(double *stressS, double *stressN, double *FstressS, double *FstressN, double *update_field) {
	state.FstressS = FstressS;
	state.FstressN = FstressN;
    state.stressS = stressS;
	state.stressN = stressN;
	state.updateField = update_field;
}

/*!
 Returns a BlockInfo structure for the info of this block.
 */
BlockInfo Block::getBlockInfo(void) const {
	BlockInfo		block_info;
	unsigned int	i;
	
	block_info.bid = getBlockID();
	block_info.fid = getFaultID();
    block_info.sid = getSectionID();
	
	for (i=0;i<4;++i) {
		block_info.x_pt[i] = _vert[i][0];
		block_info.y_pt[i] = _vert[i][1];
		block_info.z_pt[i] = _vert[i][2];
		block_info.das_pt[i] = _das[i];
		block_info.trace_flag_pt[i] = _trace_flag[i];
	}
	
	block_info.slip_velocity = _slip_rate;
	block_info.aseismicity = _aseis_factor;
    block_info.rake = _rake;
    block_info.dip = dip();
	block_info.dynamic_strength = dynamic_strength;
	block_info.static_strength = static_strength;
	block_info.lame_mu = mu;
	block_info.lame_lambda = lambda;
	
	strncpy(block_info.fault_name, getFaultName().c_str(), FAULT_NAME_MAX_LEN);
	
	return block_info;
}

StateCheckpointData State::readCheckpointData(void) const {
	StateCheckpointData		ckpt_data;
	
	ckpt_data.slipDeficit = slipDeficit;
	ckpt_data.slipCumulative = slipCumulative;
	ckpt_data.cff = cff;
	ckpt_data.stressS = *stressS;
	ckpt_data.stressN = *stressN;
	ckpt_data.updateField = *updateField;
	
	return ckpt_data;
}

void State::storeCheckpointData(const StateCheckpointData &ckpt_data) {
	slipDeficit = ckpt_data.slipDeficit;
	slipCumulative = ckpt_data.slipCumulative;
	cff = ckpt_data.cff;
	*stressS = ckpt_data.stressS;
	*stressN = ckpt_data.stressN;
	*updateField = ckpt_data.updateField;
}

std::ostream& operator<<(std::ostream& os, const State& s) {
    os << std::setw(10) << s.cff
		<< std::setw(12) << s.updateField[0]
		<< std::setw(10) << s.stressS[0]
		<< std::setw(10) << s.stressN[0];
	
	return os;
}
