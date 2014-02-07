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

#include "AddNoise.h"

void AddNoise::initDesc(const SimFramework *_sim) const {
    const VCSimulation          *sim = static_cast<const VCSimulation *>(_sim);

    sim->console() << "# Adding noise (" << sim->getStressNoiseResolution()
                   << "km resolution) to blocks." << std::endl;
}

void AddNoise::init(SimFramework *_sim) {
    VCSimulation                        *sim = static_cast<VCSimulation *>(_sim);
    double                              stress_noise, stress_res_dist, cur_noise, old_stress_drop;
    std::map<quakelib::Vec<3>, double>  centers_noise;
    quakelib::Vec<3>                    mid_pt;
    std::map<quakelib::Vec<3>, double>::iterator    cit;
    BlockList::iterator                 it;

    stress_noise = sim->getStressNoise();
    stress_res_dist = sim->getStressNoiseResolution();

    //! For each block in the simulation, get the midpoint and record the noise added to the block.
    //! Each block within a specified radius of the selected block will have the same level of noise added.
    //! The noise is defined as a modification of the stress drop from x to x +/- stress_noise.
    for (it=sim->begin(); it!=sim->end(); ++it) {
        cur_noise = 1.0e10;
        mid_pt = it->center();

        for (cit=centers_noise.begin(); cit!=centers_noise.end(); ++cit) {
            if (mid_pt.dist(cit->first) <= stress_res_dist) {
                cur_noise = cit->second;
                break;
            }
        }

        if (cur_noise > 1.0e9) {
            cur_noise = (1.0-stress_noise) + 2*stress_noise*sim->randDouble();
            centers_noise.insert(std::make_pair(mid_pt, cur_noise));
        }

        old_stress_drop = it->getStressDrop();
        it->setStressDrop(old_stress_drop*cur_noise);
    }
}
