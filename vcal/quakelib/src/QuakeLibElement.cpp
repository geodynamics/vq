// Copyright (c) 2012 Eric Heien <emheien@ucdavis.edu>
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

// Method of separating axes
// http://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
template <unsigned int nverts>
bool quakelib::Element<nverts>::check_separation(const Vec<3> &sep_vec, const Element<nverts> &block) const {
	double			high_b, high_a, low_b, low_a;
	unsigned int	i;
	
	high_b = high_a = -DBL_MAX;
	low_b = low_a = DBL_MAX;
	
	for (i=0;i<nverts;++i) {
		high_a = fmax(high_a, block._vert[i].dot_product(sep_vec));
		low_a = fmin(low_a, block._vert[i].dot_product(sep_vec));
		high_b = fmax(high_b, _vert[i].dot_product(sep_vec));
		low_b = fmin(low_b, _vert[i].dot_product(sep_vec));
	}
	
	return (low_a > high_b || low_b > high_a);
};

//! Checks if this block overlaps another. Coincident edges do not count as overlap.
// Method of separating axes
// http://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
template <unsigned int nverts>
bool quakelib::Element<nverts>::overlaps(const Element<nverts> &block) const {
	bool					normal_a_sep, normal_b_sep;
	std::vector<Vec<3> >	normals;
	unsigned int			i, j;
	
	normals.push_back(block.normal());
	normals.push_back(normal());
	
	normal_a_sep = check_separation(normals[0], block);
	normal_b_sep = check_separation(normals[1], block);
	
	if (normal_a_sep || normal_b_sep) {
		//there is no intersection
		return false;
	} else {
		//there could be an intersection
		//test if they are co-planar
		Vec<3> cross = normal().cross(block.normal());
		std::vector<Vec<3> > edges_a;
		std::vector<Vec<3> > edges_b;
		std::vector<Vec<3> > sep_vecs;
		
		edges_a.push_back(block._vert[1] - block._vert[0]);
		edges_a.push_back(block._vert[2] - block._vert[1]);
		edges_a.push_back(block._vert[3] - block._vert[2]);
		edges_a.push_back(block._vert[0] - block._vert[3]);
		
		edges_b.push_back(_vert[1] - _vert[0]);
		edges_b.push_back(_vert[2] - _vert[1]);
		edges_b.push_back(_vert[3] - _vert[2]);
		edges_b.push_back(_vert[0] - _vert[3]);
		
		if (cross.mag() > 0) {
			//they are not co-planar
			for (i=0;i<edges_a.size();i++) {
				for (j=0;j<edges_b.size();j++) {
					sep_vecs.push_back(edges_a[i].cross(edges_b[j]).unit_vector());
				}
			}
		} else {
			//they are coplanar. the case where they could be parallel is tested with the normal separations.
			//this is now the 2d case.
			for (i=0;i<edges_a.size();i++) {
				sep_vecs.push_back(normals[0].cross(edges_a[i]).unit_vector());
			}
			for (i=0;i<edges_b.size();i++) {
				sep_vecs.push_back(normals[0].cross(edges_b[i]).unit_vector());
			}
		}
		
		for (i=0; i < sep_vecs.size(); i++) {
			// if any of the results are seperated the elements don't intersect
			if (check_separation(sep_vecs[0], block)) {
				return false;
			}
		}
		
		return true;
	}
}

template <unsigned int nverts>
quakelib::Vec<3> quakelib::Element<nverts>::calc_displacement_vector(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
	Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;
    
	if (loc_lambda <= 0 || loc_mu <= 0) {
		throw std::invalid_argument("Lambda and mu must be greater than zero.");
	}
    
	if (nverts != 4) {
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
    
    L = (_vert[3] - _vert[0]).mag();
    W = (_vert[3] - _vert[2]).mag();
    c = fabs(_vert[1][2]);
    
    return block_okada.calc_displacement_vector(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}

template <unsigned int nverts>
quakelib::Vec<3> quakelib::Element<nverts>::calc_dudx(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
	Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;
    
	if (loc_lambda <= 0 || loc_mu <= 0) {
		throw std::invalid_argument("Lambda and mu must be greater than zero.");
	}
    
	if (nverts != 4) {
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
    
    L = (_vert[3] - _vert[0]).mag();
    W = (_vert[3] - _vert[2]).mag();
    c = fabs(_vert[1][2]);
    
    return block_okada.calc_dudx(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}

template <unsigned int nverts>
quakelib::Vec<3> quakelib::Element<nverts>::calc_dudy(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
	Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;
    
	if (loc_lambda <= 0 || loc_mu <= 0) {
		throw std::invalid_argument("Lambda and mu must be greater than zero.");
	}
    
	if (nverts != 4) {
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
    
    L = (_vert[3] - _vert[0]).mag();
    W = (_vert[3] - _vert[2]).mag();
    c = fabs(_vert[1][2]);
    
    return block_okada.calc_dudy(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}

template <unsigned int nverts>
quakelib::Vec<3> quakelib::Element<nverts>::calc_dudz(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
	Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;
    
	if (loc_lambda <= 0 || loc_mu <= 0) {
		throw std::invalid_argument("Lambda and mu must be greater than zero.");
	}
    
	if (nverts != 4) {
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
    
    L = (_vert[3] - _vert[0]).mag();
    W = (_vert[3] - _vert[2]).mag();
    c = fabs(_vert[1][2]);
    
    return block_okada.calc_dudz(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}

template <unsigned int nverts>
quakelib::Tensor<3,3> quakelib::Element<nverts>::calc_stress_tensor(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument) {
    Okada block_okada;
    double US, UD, UT, L, W, c, cos_result, sin_result;
    
	if (loc_lambda <= 0 || loc_mu <= 0) {
		throw std::invalid_argument("Lambda and mu must be greater than zero.");
	}
    
	if (nverts != 4) {
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
    
    L = (_vert[3] - _vert[0]).mag();
    W = (_vert[3] - _vert[2]).mag();
    c = fabs(_vert[1][2]);
    
    return block_okada.calc_stress_tensor(location, c, dip(), L, W, US, UD, UT, loc_lambda, loc_mu);
}

// Instantiate empty classes to ensure templates are filled
namespace quakelib {
	template class quakelib::Element<3>;
	template class quakelib::Element<4>;
}
