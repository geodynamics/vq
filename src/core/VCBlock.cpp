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

VCRand::VCRand(const VCRand::State &state) {
    idum1 = state.idum1;
    idum2 = state.idum2;
    iy    = state.iy;

    for (int i=0; i<NTAB; i++) iv[i] = state.iv[i];
}

void VCRand::save(VCRand::State &state) {
    state.idum1 = idum1;
    state.idum2 = idum2;
    state.iy    = iy;

    for (int i=0; i<NTAB; i++) state.iv[i] = iv[i];
}

void VCRand::init(int seed) {
    idum2 = idum1 = (seed < 1 ? 1 : seed);

    for (int j=NTAB+7; j>=0; j--) {
        int k = idum1/IQ1;
        idum1  = IA1*(idum1 - k*IQ1) - k*IR1;

        if (idum1 < 0) idum1 += IM1;

        if (j < NTAB) iv[j] = idum1;
    }

    iy=iv[0];
}

void VCRand::save(int *state) {
    state[0] = idum1;
    state[1] = idum2;
    state[2] = iy;

    int shift = 3;

    for (int i=0; i<VCRand::NTAB; i++) state[shift+i] = iv[i];
}

void VCRand::init(int *state) {
    idum1 = state[0];
    idum2 = state[1];
    iy    = state[2];

    int shift = 3;

    for (int i=0; i<VCRand::NTAB; i++) iv[i] = state[shift+i];
}

double VCRand::nextDouble() {
    int k = idum1/IQ1;
    idum1  = IA1*(idum1 - k*IQ1) - k*IR1;

    if (idum1 < 0) idum1 += IM1;

    k = idum2/IQ2;
    idum2  = IA2*(idum2 - k*IQ2) - k*IR2;

    if (idum2 < 0) idum2 += IM2;

    int tmp = iy/NDIV;
    iy = iv[tmp] - idum2;
    iv[tmp] = idum1;

    if (iy < 1) iy += IMM1;

    if (iy > IMM1) return 0.99999999;
    else          return iy*NORM;
}

std::ostream &operator<<(std::ostream &os, const BlockVal &bv) {
    if (bv.block_id == UNDEFINED_BLOCK_ID) os << "(Undef, ";
    else os << "(B" << bv.block_id << ", ";

    os << bv.val << ")";
    return os;
}

void Block::clear(void) {
    quakelib::SimElement::clear();
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
            //  state.Fcff = 0;
            //}


        } else {
            state.stressS[0] = state.FstressS[0];
            state.stressN[0] = state.FstressN[0];
            state.cff = state.stressS[0] - friction()*state.stressN[0];
            state.Fcff = state.cff;
        }
    } else {
        state.FstressS[0] = state.stressS[0];
        state.FstressN[0] = state.stressN[0];
        state.cff = state.stressS[0] - friction()*state.stressN[0];
        state.Fcff = state.FstressS[0] - friction()*state.FstressN[0];
    }

    //state.cff = state.stressS[0] - stress_drop;
}

void Block::setStatePtrs(double *stressS, double *stressN, double *update_field) {
    state.stressS = stressS;
    state.stressN = stressN;
    state.updateField = update_field;
}

StateCheckpointData State::readCheckpointData(void) const {
    StateCheckpointData     ckpt_data;

    ckpt_data.slipDeficit = slipDeficit;
    ckpt_data.cff = cff;
    ckpt_data.stressS = *stressS;
    ckpt_data.stressN = *stressN;
    ckpt_data.updateField = *updateField;

    return ckpt_data;
}

void State::storeCheckpointData(const StateCheckpointData &ckpt_data) {
    slipDeficit = ckpt_data.slipDeficit;
    cff = ckpt_data.cff;
    *stressS = ckpt_data.stressS;
    *stressN = ckpt_data.stressN;
    *updateField = ckpt_data.updateField;
}

std::ostream &operator<<(std::ostream &os, const State &s) {
    os << std::setw(10) << s.cff
       << std::setw(12) << s.updateField[0]
       << std::setw(10) << s.stressS[0]
       << std::setw(10) << s.stressN[0];

    return os;
}
