// Copyright (c) 2015 Eric Heien <emheien@ucdavis.edu>, Kasey Schultz <kwschultz@ucdavis.edu>
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

#ifdef QUAKELIB_HAVE_MATH_H
#include <math.h>
#endif

// Gravity change equations taken from Okubo 1992
quakelib::FloatList quakelib::SlipMap::gravity_changes(const VectorList &points, const float &lambda, const float &mu, const float &cutoff) {
	quakelib::FloatList gravity_changes;
	Okada block_okada;
    double gravity_change, slip, US, UD, UT, L, W, c, rake_cos, rake_sin, strike_cos, strike_sin, dip, strike, xp0, yp0, xp3, yp3, x, y, xp, yp;
	
	if (lambda <= 0 || mu <= 0) {
		throw std::invalid_argument("Lambda and mu must be greater than zero.");
	}
			
	for(VectorList::size_type point_id = 0; point_id != points.size(); point_id++) {
		gravity_changes.push_back(0.0);
	}
    
	for(SlippedElementList::size_type ele_id = 0; ele_id != involved_elements.size(); ele_id++) {
		slip = involved_elements[ele_id].slip();
		c = involved_elements[ele_id].max_depth();
		rake_cos = cos(involved_elements[ele_id].rake());
		rake_sin = sin(involved_elements[ele_id].rake());
		if (fabs(rake_cos) < TRIG_TOLERANCE) {
			rake_cos = 0.0;
		}
		if (fabs(rake_sin) < TRIG_TOLERANCE) {
			rake_sin = 0.0;
		}
		US = slip * rake_cos;
		UD = slip * rake_sin;
		UT = 0.0;
		
		L = (involved_elements[ele_id].implicit_vert() - involved_elements[ele_id].vert(0)).mag();
		W = (involved_elements[ele_id].implicit_vert() - involved_elements[ele_id].vert(2)).mag();
		c = fabs(involved_elements[ele_id].vert(1)[2]);
		
		dip = involved_elements[ele_id].dip();
		strike = involved_elements[ele_id].strike();
		
		strike_cos = cos(strike);
		strike_sin = sin(strike);
		if (fabs(strike_cos) < TRIG_TOLERANCE) {
			strike_cos = 0.0;
		}
		if (fabs(strike_sin) < TRIG_TOLERANCE) {
			strike_sin = 0.0;
		}
		
		xp0 = involved_elements[ele_id].vert(0)[0];
		yp0 = involved_elements[ele_id].vert(0)[1];
		
		xp3 = involved_elements[ele_id].implicit_vert()[0];
		yp3 = involved_elements[ele_id].implicit_vert()[1];
		
		for(VectorList::size_type point_id = 0; point_id != points.size(); point_id++) {
			x = points[point_id][0];
			y = points[point_id][1];
			
			if (sqrt(pow((x-(xp0 + xp3)/2.0),2)+pow((y-(yp0+yp3)/2.0),2))/sqrt(L*W) > (cutoff + slip - 1.0) ) {
				gravity_change = 0.0;
			} else {
			
				xp = (x-xp0) * strike_sin - (y-yp0) * strike_cos;
				yp = (x-xp0) * strike_cos + (y-yp0) * strike_sin;
				
				gravity_change = block_okada.calc_dg(quakelib::Vec<2>(xp,yp), c, dip, L, W, US, UD, UT, lambda, mu);

			}
			
			//std::cout << sqrt(pow((x-(xp0 + xp3)/2.0),2)+pow((y-(yp0 + yp3)/2.0),2))/sqrt(L*W) << " " << gravity_change << std::endl;
        
			gravity_changes[point_id] += gravity_change;
		}
	}
	
	return gravity_changes;

};

//The dilatational gravity changes arise from only compression/dilatation
quakelib::FloatList quakelib::SlipMap::dilat_gravity_changes(const VectorList &points, const float &lambda, const float &mu, const float &cutoff) {
	quakelib::FloatList gravity_changes;
	Okada block_okada;
    double gravity_change, slip, US, UD, UT, L, W, c, rake_cos, rake_sin, strike_cos, strike_sin, dip, strike, xp0, yp0, xp3, yp3, x, y, xp, yp;
	
	if (lambda <= 0 || mu <= 0) {
		throw std::invalid_argument("Lambda and mu must be greater than zero.");
	}
			
	for(VectorList::size_type point_id = 0; point_id != points.size(); point_id++) {
		gravity_changes.push_back(0.0);
	}
    
	for(SlippedElementList::size_type ele_id = 0; ele_id != involved_elements.size(); ele_id++) {
		slip = involved_elements[ele_id].slip();
		c = involved_elements[ele_id].max_depth();
		rake_cos = cos(involved_elements[ele_id].rake());
		rake_sin = sin(involved_elements[ele_id].rake());
		if (fabs(rake_cos) < TRIG_TOLERANCE) {
			rake_cos = 0.0;
		}
		if (fabs(rake_sin) < TRIG_TOLERANCE) {
			rake_sin = 0.0;
		}
		US = slip * rake_cos;
		UD = slip * rake_sin;
		UT = 0.0;
		
		L = (involved_elements[ele_id].implicit_vert() - involved_elements[ele_id].vert(0)).mag();
		W = (involved_elements[ele_id].implicit_vert() - involved_elements[ele_id].vert(2)).mag();
		c = fabs(involved_elements[ele_id].vert(1)[2]);
		
		dip = involved_elements[ele_id].dip();
		strike = involved_elements[ele_id].strike();
		
		strike_cos = cos(strike);
		strike_sin = sin(strike);
		if (fabs(strike_cos) < TRIG_TOLERANCE) {
			strike_cos = 0.0;
		}
		if (fabs(strike_sin) < TRIG_TOLERANCE) {
			strike_sin = 0.0;
		}
		
		xp0 = involved_elements[ele_id].vert(0)[0];
		yp0 = involved_elements[ele_id].vert(0)[1];
		
		xp3 = involved_elements[ele_id].implicit_vert()[0];
		yp3 = involved_elements[ele_id].implicit_vert()[1];
		
		for(VectorList::size_type point_id = 0; point_id != points.size(); point_id++) {
			x = points[point_id][0];
			y = points[point_id][1];
			
			//if (pow(x-event_center()[0], 2) + pow(y-event_center()[1], 2) > pow( event_radius() * cutoff ,2) ) {
			// Gotta figure the cutoff for gravity changes out
			
			if (sqrt(pow((x-(xp0 + xp3)/2.0),2)+pow((y-(yp0+yp3)/2.0),2))/sqrt(L*W) > (cutoff + slip - 1.0) ) {
				gravity_change = 0.0;
			} else {
			
				xp = (x-xp0) * strike_sin - (y-yp0) * strike_cos;
				yp = (x-xp0) * strike_cos + (y-yp0) * strike_sin;
				
				gravity_change = block_okada.calc_dg_dilat(quakelib::Vec<2>(xp,yp), c, dip, L, W, US, UD, UT, lambda, mu);
			}
			
			//std::cout << sqrt(pow((x-(xp0 + xp3)/2.0),2)+pow((y-(yp0 + yp3)/2.0),2))/sqrt(L*W) << " " << gravity_change << std::endl;
			
			gravity_changes[point_id] += gravity_change;
		}
	}
	return gravity_changes;
}

// Displacements computed from equations given in Okada 1995
quakelib::VectorList quakelib::SlipMap::displacements(const VectorList &points, const float &lambda, const float &mu, const float &cutoff)  {
	quakelib::VectorList displacements;
	Okada block_okada;
	quakelib::Vec<3> displacement;
    double slip, US, UD, UT, L, W, c, rake_cos, rake_sin, strike_cos, strike_sin, dip, strike, xp0, yp0, x, y, xp, yp, dx, dy, dz;
	
	if (lambda <= 0 || mu <= 0) {
		throw std::invalid_argument("Lambda and mu must be greater than zero.");
	}
			
	for(VectorList::size_type point_id = 0; point_id != points.size(); point_id++) {
		displacements.push_back(quakelib::Vec<3>(0.0,0.0,0.0));
	}
    
	for(SlippedElementList::size_type ele_id = 0; ele_id != involved_elements.size(); ele_id++) {
		slip = involved_elements[ele_id].slip();
		c = involved_elements[ele_id].max_depth();
		rake_cos = cos(involved_elements[ele_id].rake());
		rake_sin = sin(involved_elements[ele_id].rake());
		if (fabs(rake_cos) < TRIG_TOLERANCE) {
			rake_cos = 0.0;
		}
		if (fabs(rake_sin) < TRIG_TOLERANCE) {
			rake_sin = 0.0;
		}
		US = slip * rake_cos;
		UD = slip * rake_sin;
		UT = 0.0;
		
		L = (involved_elements[ele_id].implicit_vert() - involved_elements[ele_id].vert(0)).mag();
		W = (involved_elements[ele_id].implicit_vert() - involved_elements[ele_id].vert(2)).mag();
		c = fabs(involved_elements[ele_id].vert(1)[2]);
		
		dip = involved_elements[ele_id].dip();
		strike = involved_elements[ele_id].strike();
		
		strike_cos = cos(strike);
		strike_sin = sin(strike);
		if (fabs(strike_cos) < TRIG_TOLERANCE) {
			strike_cos = 0.0;
		}
		if (fabs(strike_sin) < TRIG_TOLERANCE) {
			strike_sin = 0.0;
		}
		
		xp0 = involved_elements[ele_id].vert(0)[0];
		yp0 = involved_elements[ele_id].vert(0)[1];
		
		/*
		std::cout << "v1: <" << involved_elements[ele_id].vert(0)[0] << ", " << involved_elements[ele_id].vert(0)[1] << ", " << involved_elements[ele_id].vert(0)[2] << ">" << std::endl;
		std::cout << "v2: <" << involved_elements[ele_id].vert(1)[0] << ", " << involved_elements[ele_id].vert(1)[1] << ", " << involved_elements[ele_id].vert(1)[2] << ">" << std::endl;
		std::cout << "v3: <" << involved_elements[ele_id].vert(2)[0] << ", " << involved_elements[ele_id].vert(2)[1] << ", " << involved_elements[ele_id].vert(2)[2] << ">" << std::endl;
		std::cout << "v4: <" << involved_elements[ele_id].implicit_vert()[0] << ", " << involved_elements[ele_id].implicit_vert()[1] << ", " << involved_elements[ele_id].implicit_vert()[2] << ">" << std::endl;
		std::cout << "L: " << L << std::endl;
		std::cout << "W: " << W << std::endl;
		std::cout << "c: " << c << std::endl;
		std::cout << "US: " << US << std::endl;
		std::cout << "UD: " << UD << std::endl;
		std::cout << "UT: " << UT << std::endl;
		std::cout << "dip: " << dip << std::endl;
		std::cout << "strike: " << strike << std::endl;
		std::cout << "xp0: " << xp0 << std::endl;
		std::cout << "yp0: " << yp0 << std::endl;
		*/
		
		for(VectorList::size_type point_id = 0; point_id != points.size(); point_id++) {
			x = points[point_id][0];
			y = points[point_id][1];
			
			//if (pow(x-event_center()[0], 2) + pow(y-event_center()[1], 2) > pow( event_radius() * cutoff ,2) ) {
			if (sqrt(pow((x-xp0),2)+pow((y-yp0),2))/sqrt(L*W) > (cutoff + slip - 1.0) ) {
				dx = 0.0;
				dy = 0.0;
				dz = 0.0;
			} else {
				xp = (x-xp0) * strike_sin - (y-yp0) * strike_cos;
				yp = (x-xp0) * strike_cos + (y-yp0) * strike_sin;
				
				displacement = block_okada.calc_displacement_vector(quakelib::Vec<3>(xp,yp,0.0), c, dip, L, W, US, UD, UT, lambda, mu);
				
				dx =  displacement[0] * strike_sin + displacement[1] * strike_cos;
				dy = -displacement[0] * strike_cos + displacement[1] * strike_sin;
				dz =  displacement[2];

			}
			
			//std::cout << pow(x-event_center[0], 2) + pow(y-event_center[1], 2) << " " << pow(event_radius * cutoff,2) << std::endl;
			
			displacements[point_id][0] += dx;
			displacements[point_id][1] += dy;
			displacements[point_id][2] += dz;
		}
		
	}
	
	return displacements;
}





