// Copyright (c) 2012-2013 Eric M. Heien
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

#include "QuakeLib.h"

quakelib::Vec<3> quakelib::SimElement::calc_displacement_vector(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
    Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;

    if (loc_lambda <= 0 || loc_mu <= 0) {
        throw std::invalid_argument("Lambda and mu must be greater than zero.");
    }

    if (!_is_quad) {
        throw std::invalid_argument("Displacement calculation currently only supported for rectangular elements.");
    }

    cos_result = cos(rake());
    sin_result = sin(rake());

    if (fabs(cos_result) < TRIG_TOLERANCE) {
        cos_result = 0.0;
    }

    if (fabs(sin_result) < TRIG_TOLERANCE) {
        sin_result = 0.0;
    }

    US = unit_slip * cos_result;
    UD = unit_slip * sin_result;
    UT = 0.0;

    L = (_vert[2] - _vert[0]).mag();
    W = (_vert[1] - _vert[0]).mag();
    c = fabs(max_depth());

    return block_okada.calc_displacement_vector(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}

quakelib::Vec<3> quakelib::SimElement::calc_dudx(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
    Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;

    if (loc_lambda <= 0 || loc_mu <= 0) {
        throw std::invalid_argument("Lambda and mu must be greater than zero.");
    }

    if (!_is_quad) {
        throw std::invalid_argument("Stress tensor calculation currently only supported for rectangular elements.");
    }

    cos_result = cos(rake());
    sin_result = sin(rake());

    if (fabs(cos_result) < TRIG_TOLERANCE) {
        cos_result = 0.0;
    }

    if (fabs(sin_result) < TRIG_TOLERANCE) {
        sin_result = 0.0;
    }

    US = unit_slip * cos_result;
    UD = unit_slip * sin_result;
    UT = 0.0;

    L = (_vert[2] - _vert[0]).mag();
    W = (_vert[1] - _vert[0]).mag();
    c = fabs(max_depth());

    return block_okada.calc_dudx(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}

quakelib::Vec<3> quakelib::SimElement::calc_dudy(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
    Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;

    if (loc_lambda <= 0 || loc_mu <= 0) {
        throw std::invalid_argument("Lambda and mu must be greater than zero.");
    }

    if (!_is_quad) {
        throw std::invalid_argument("Stress tensor calculation currently only supported for rectangular elements.");
    }

    cos_result = cos(rake());
    sin_result = sin(rake());

    if (fabs(cos_result) < TRIG_TOLERANCE) {
        cos_result = 0.0;
    }

    if (fabs(sin_result) < TRIG_TOLERANCE) {
        sin_result = 0.0;
    }

    US = unit_slip * cos_result;
    UD = unit_slip * sin_result;
    UT = 0.0;

    L = (_vert[2] - _vert[0]).mag();
    W = (_vert[1] - _vert[0]).mag();
    c = fabs(max_depth());

    return block_okada.calc_dudy(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}

quakelib::Vec<3> quakelib::SimElement::calc_dudz(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
    Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;

    if (loc_lambda <= 0 || loc_mu <= 0) {
        throw std::invalid_argument("Lambda and mu must be greater than zero.");
    }

    if (!_is_quad) {
        throw std::invalid_argument("Stress tensor calculation currently only supported for rectangular elements.");
    }

    cos_result = cos(rake());
    sin_result = sin(rake());

    if (fabs(cos_result) < TRIG_TOLERANCE) {
        cos_result = 0.0;
    }

    if (fabs(sin_result) < TRIG_TOLERANCE) {
        sin_result = 0.0;
    }

    US = unit_slip * cos_result;
    UD = unit_slip * sin_result;
    UT = 0.0;

    L = (_vert[2] - _vert[0]).mag();
    W = (_vert[1] - _vert[0]).mag();
    c = fabs(max_depth());

    return block_okada.calc_dudz(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}

quakelib::Tensor<3,3> quakelib::SimElement::calc_stress_tensor(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
    Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;

    if (loc_lambda <= 0 || loc_mu <= 0) {
        throw std::invalid_argument("Lambda and mu must be greater than zero.");
    }

    if (!_is_quad) {
        throw std::invalid_argument("Stress tensor calculation currently only supported for rectangular elements.");
    }

    cos_result = cos(rake());
    sin_result = sin(rake());

    if (fabs(cos_result) < TRIG_TOLERANCE) {
        cos_result = 0.0;
    }

    if (fabs(sin_result) < TRIG_TOLERANCE) {
        sin_result = 0.0;
    }

    US = unit_slip * cos_result;
    UD = unit_slip * sin_result;
    UT = 0.0;

    L = (_vert[2] - _vert[0]).mag();
    W = (_vert[1] - _vert[0]).mag();
    c = fabs(max_depth());

    return block_okada.calc_stress_tensor(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}
