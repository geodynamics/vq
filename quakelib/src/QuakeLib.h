// Copyright (c) 2012-2014 Eric M. Heien
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

#include "quakelib_config.h"

#ifdef QUAKELIB_HAVE_FLOAT_H
#include <float.h>
#endif

#ifdef QUAKELIB_HAVE_STRING_H
#include <string.h>
#endif

#ifdef QUAKELIB_HAVE_MATH_H
#include <math.h>
#endif

#include "QuakeLibOkada.h"

#ifndef _QUAKELIB_H_
#define _QUAKELIB_H_

// This will result in a smallest calculated displacement of ~0.0001 m
#define DIST_SQRT_AREA_RATIO_CUTOFF_DISPLACEMENTS                46.5

// This will result in a smallest calculated gravity change of ~4e-9 
#define DIST_SQRT_AREA_RATIO_CUTOFF_GRAVITY               8.0

namespace quakelib {
    // Function to obtain the version/git information for the compiled quakelib
    static std::string quakelib_info(void) {
        std::stringstream ss;
        ss << "QuakeLib " << QUAKELIB_VERSION_STR;
#ifdef QUAKELIB_GIT_SHA1
        ss << " Git revision " << QUAKELIB_GIT_SHA1;
#endif
        return ss.str();
    };

    //! Represents a complete triangular or rectangular element for use in a simulation.
    class SimElement {
        protected:
            //! Coordinates of element vertices in meters
            Vec<3>      _vert[3];
            //! Whether the vertices describe a triangle (false) or parallelogram (quadrilateral) (true)
            bool        _is_quad;
            //! Rake angle (radians, 0.0 = left lateral, PI/2 = positive side moves up)
            double      _rake;
            //! Element slip rate (meters/second)
            double      _slip_rate;
            //! Element aseismicity factor (in [0, 1])
            double      _aseis_factor;
            //! Lame mu parameter (Pascals)
            double      _lame_mu;
            //! Lame lambda parameter (Pascals)
            double      _lame_lambda;
            //! Static yield strength (Pascals)
            double      _static_strength;
            //! Dynamic sliding strength (Pascals)
            double      _dynamic_strength;
            //! Maximum slip distance of this element (meters)
            double      _max_slip;

        public:
            SimElement(void) {
                clear();
            }
            //! Calculate the stress tensor at a location with Lame parameters lambda and mu
            //! caused by this element moving unit_slip meters.
            Tensor<3,3> calc_stress_tensor(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);
            //! Calculate the derivatives at a location with Lame parameters lambda and mu
            //! caused by this element moving unit_slip meters.
            Vec<3> calc_dudx(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);
            Vec<3> calc_dudy(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);
            Vec<3> calc_dudz(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);
            //! Calculate the displacement vector at a location with Lame parameters lambda and mu
            //! caused by this element moving unit_slip meters.
            Vec<3> calc_displacement_vector(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);

            Vec<3> interpolate_point(const double &i, const double &j) const throw(std::invalid_argument) {
                if (!_is_quad) throw std::invalid_argument("Only supports quads currently.");

                return _vert[0] + (_vert[1]-_vert[0])*i + (_vert[2]-_vert[0])*j;
            }

            //! Set the location of a vertex.
            void set_vert(const unsigned int &vert, const Vec<3> &new_vert) throw(std::out_of_range) {
                if (vert>=3) throw std::out_of_range("quakelib::Element::set_vert");

                _vert[vert] = new_vert;
            };
            //! Get the location of a vertex.
            Vec<3> vert(const unsigned int &vert) const throw(std::out_of_range) {
                if (vert>=3) throw std::out_of_range("quakelib::Element::vert");

                return _vert[vert];
            };
            // For quadrilateral elements, calculate the implicit 4th point
            Vec<3> implicit_vert(void) const throw(std::out_of_range) {
                if (!_is_quad) throw std::out_of_range("quakelib::Element::implicit_vert");

                return _vert[2]+(_vert[1]-_vert[0]);
            };

            bool is_quad(void) const {
                return _is_quad;
            };
            void set_is_quad(const bool &is_quad) {
                _is_quad = is_quad;
            };

            //! Set the slip rate in m/s for this block.
            void set_slip_rate(const double &new_slip_rate) throw(std::invalid_argument) {
                if (isnan(new_slip_rate)) throw std::invalid_argument("quakelib::Element::set_slip_rate");

                _slip_rate = new_slip_rate;
            };
            //! Get the slip rate in cm/year for this block.
            double slip_rate(void) const {
                return _slip_rate;
            };

            //! Set the rake angle of this block in radians.
            void set_rake(const double &new_rake) throw(std::invalid_argument) {
                if (isnan(new_rake)) throw std::invalid_argument("quakelib::Element::set_rake");

                _rake = new_rake;
            };
            //! Get the rake angle of this block in radians.
            double rake(void) const {
                return _rake;
            };

            //! Get the fraction of slip which is aseismic for this element.
            double aseismic(void) const {
                return _aseis_factor;
            };
            //! Set the fraction of slip which is aseismic for this element.
            void set_aseismic(const double &new_aseismic) throw(std::invalid_argument) {
                if (new_aseismic < 0 || new_aseismic > 1 || isnan(new_aseismic)) throw std::invalid_argument("quakelib::Element::set_aseismic");

                _aseis_factor = new_aseismic;
            };

            //! Get the Lame mu parameter for this element
            double lame_mu(void) const {
                return _lame_mu;
            };
            //! Set the Lame mu parameter for this element
            void set_lame_mu(const double &new_lame_mu) throw(std::invalid_argument) {
                if (new_lame_mu < 0 || isnan(new_lame_mu)) throw std::invalid_argument("quakelib::Element::set_lame_mu");

                _lame_mu = new_lame_mu;
            };

            //! Get the Lame lambda parameter for this element
            double lame_lambda(void) const {
                return _lame_mu;
            };
            //! Set the Lame lambda parameter for this element
            void set_lame_lambda(const double &new_lame_lambda) throw(std::invalid_argument) {
                if (new_lame_lambda < 0 || isnan(new_lame_lambda)) throw std::invalid_argument("quakelib::Element::set_lame_lambda");

                _lame_lambda = new_lame_lambda;
            };

            //! Get the maximum slip distance for this element
            double max_slip(void) const {
                return _max_slip;
            }
            //! Set the maximum slip distance for this element
            void set_max_slip(const double &new_max_slip) throw(std::invalid_argument) {
                if (new_max_slip < 0 || isnan(new_max_slip)) throw std::invalid_argument("quakelib::Element::set_max_slip");

                _max_slip = new_max_slip;
            }

            //! Clear all variables for this element.
            void clear(void) {
                _vert[0] = Vec<3>::nan_vec();
                _vert[1] = Vec<3>::nan_vec();
                _vert[2] = Vec<3>::nan_vec();
                _is_quad = false;
                _rake = _slip_rate = _aseis_factor = std::numeric_limits<double>::quiet_NaN();
                _lame_mu = _lame_lambda = std::numeric_limits<double>::quiet_NaN();
                _static_strength = _dynamic_strength = _max_slip = std::numeric_limits<double>::quiet_NaN();
            };

            //! Returns a unit vector along the direction of fault dip.
            double dip(void) const {
                return normal().vector_angle(Vec<3>(0,0,1));
            };

            //! Returns unit vector along the direction of fault rake.
            Vec<3> rake_vector(void) const {
                Vec<3> vec;
                vec = (_vert[2]-_vert[0]).unit_vector();
                return vec.rotate_around_axis(normal()*(-1.0), _rake);
            };

            //! Get the height of this block (in m).
            double height(void) const {
                return max_depth()-min_depth();
            };

            //! Get center point of this block
            Vec<3> center() const {
                Vec<3> c;

                for (unsigned int i=0; i<3; ++i) c += _vert[i];

                if (_is_quad) {
                    c += implicit_vert();
                    return c / 4.0;
                } else {
                    return c / 3.0;
                }
            };

            //! Calculates the area of this block based on vertices and whether it is a triangle or quad.
            double area(void) const {
                Vec<3> a,b;
                a=_vert[1]-_vert[0];
                b=_vert[2]-_vert[0];
                return (_is_quad ? 1.0 : 0.5) * a.cross(b).mag();
            };

            //! Returns length along largest dimension
            double largest_dimension(void) const {
                return fmax((_vert[2]-_vert[0]).mag(),(_vert[1]-_vert[0]).mag());
            };

            //! Determine minimum depth of this block.
            //! Note that depth is negative, so max depth returns the most negative point and min depth the least negative
            double min_depth(void) const {
                double d = -DBL_MAX;

                for (unsigned int i=0; i<3; i++) d = fmax(_vert[i][2], d);

                return d;
            };
            //! Determine maximum depth of this block.
            double max_depth(void) const {
                double d = DBL_MAX;

                for (unsigned int i=0; i<3; i++) d = fmin(_vert[i][2], d);

                return d;
            };

            //! Calculates the Euclidean distance between the 3D midpoint of this block and another block.
            double center_distance(const SimElement &other) const {
                return (other.center() - this->center()).mag();
            };

            //! Whether this block is above another block (defined by the top of the blocks).
            bool is_above(const SimElement &block) const {
                return (block.min_depth() <= min_depth());
            };

            //! Get the normal unit vector to the plane of this element.
            Vec<3> normal(void) const {
                Vec<3> a,b;
                a=_vert[1]-_vert[0];
                b=_vert[2]-_vert[0];
                return a.cross(b).unit_vector();
            };
            
            //! Returns the angle of the element relative to north. Positive rotating west, negitive rotating east.
            double strike(void) const {
                Vec<3> v1, v2, norm;
                norm = normal();
                v1[0] = norm[0];
                v1[1] = norm[1];
                v2[1] = 1.0;
                return (v1.vector_angle(v2) - M_PI/2.0);
            };
    };
	
    //Extend the SimElement class to include slip, for computing fields (gravity, displacement)
	class SlippedElement: public SimElement {
	protected:
		double _slip;
        unsigned int _id;
        unsigned int _section_id; 
	public:
		double slip(void) const { return _slip; };
		void set_slip(const double &new_slip) { _slip = new_slip; };
        double id(void) const { return _id; };
		void set_id(const double &new_id) { _id = new_id; };
        double section_id(void) const { return _section_id; };
		void set_section_id(const double &new_id) { _section_id = new_id; };
	};
    
	typedef std::vector< SlippedElement > SlippedElementList;
    
    class SlipMap {
    private:
        SlippedElementList involved_elements;
        //LatLonDepth _base;
    public:
        SlipMap() { };
        void add_element(const SlippedElement &element) {involved_elements.push_back(element);};
		void add_elements(const SlippedElementList involved_elements) {
			for (unsigned int i=0; i < involved_elements.size(); i++)
				add_element(involved_elements[i]);
		};
        //void set_base(LatLonDepth &base) { _base = base; };
        VectorList displacements(const VectorList &points, const float &lambda, const float &mu, const float &cutoff=DIST_SQRT_AREA_RATIO_CUTOFF_DISPLACEMENTS);
		FloatList gravity_changes(const VectorList &points, const float &lambda, const float &mu, const float &cutoff=DIST_SQRT_AREA_RATIO_CUTOFF_GRAVITY);
        FloatList dilat_gravity_changes(const VectorList &points, const float &lambda, const float &mu, const float &cutoff=DIST_SQRT_AREA_RATIO_CUTOFF_GRAVITY);
    };
    
}

#endif
