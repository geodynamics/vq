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

#include "Block.h"
#include "QuakeLibOkada.h"
#include "QuakeLibIO.h"
#include <sstream>
#include <string.h>
#include <iomanip>
#include <list>

/*!
 New methods defined for the updated okada functions
 */

void Block::get_rake_and_normal_stress_due_to_block(double stresses[2], const double &sample_dist, const Block &source_block) const {
    // test block vectors. This block is the test block
    quakelib::Vec<3> rake_vec, normal_vec, center_vec;
    // source block vectors
    quakelib::Vec<3> source_normal_vec, source_center_vec;
    // modified test vectors
    quakelib::Vec<3> mrake_vec, mnormal_vec, mcenter_vec;
    // other stuff
    quakelib::Vec<3> rot_axis, stress_vec, shift_vec, xy_projected_source_normal;
    quakelib::Tensor<3,3> stress_tensor;
    // Points to sample the Greens function at
    std::vector<quakelib::Vec<3> >  sample_points;
    unsigned int                    i, n, n_horiz_samples, n_vert_samples;
    double                          horiz_step, vert_step, horiz_pos, vert_pos;
    double                          theta;

    normal_vec = normal();
    rake_vec = rake_vector();
    center_vec = center();

    source_normal_vec = source_block.normal();
    source_center_vec = source_block.center();

    // Take samples of the Greens function at minimum distances and average the results
    // This allows better convergence between models with few large blocks and models with many small blocks
    n_horiz_samples = fmax((_vert[2] - _vert[0]).mag()/sample_dist, 1);
    n_vert_samples = fmax((_vert[1] - _vert[0]).mag()/sample_dist, 1);

    // Select a grid of N x M evenly spaced sample points on the block
    horiz_step = 1.0/n_horiz_samples;
    vert_step = 1.0/n_vert_samples;
    horiz_pos = horiz_step/2;

    for (i=0; i<n_horiz_samples; ++i) {
        vert_pos = vert_step/2;

        for (n=0; n<n_vert_samples; ++n) {
            sample_points.push_back(interpolate_point(horiz_pos, vert_pos));
            vert_pos += vert_step;
        }

        horiz_pos += horiz_step;
    }

    stress_vec = quakelib::Vec<3>();

    for (i=0; i<sample_points.size(); ++i) {
        // first shift all points
        shift_vec = quakelib::Vec<3>(source_block._vert[1][0],source_block._vert[1][1],0.0);
        //mcenter_vec = center_vec - shift_vec;
        mcenter_vec = sample_points[i] - shift_vec;

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

        // Assume unit slip of 1.0
        stress_tensor = source_block.calc_stress_tensor(mcenter_vec, 1.0, lame_lambda(), lame_mu());

        stress_vec += stress_tensor*mnormal_vec;
    }

    stress_vec *= 1.0/sample_points.size();
    stresses[0] = stress_vec.dot_product(mrake_vec);
    stresses[1] = stress_vec.dot_product(mnormal_vec);
}

/*!*/

std::ostream &operator<<(std::ostream &os, const BlockVal &bv) {
    if (bv.block_id == UNDEFINED_ELEMENT_ID) os << "(Undef, ";
    else os << "(B" << bv.block_id << ", ";

    os << bv.val << ")";
    return os;
}

void Block::clear(void) {
    quakelib::SimElement::clear();
    id = UNDEFINED_ELEMENT_ID;
    fid = UNDEFINED_FAULT_ID;
    sid = UNDEFINED_SECTION_ID;
}
