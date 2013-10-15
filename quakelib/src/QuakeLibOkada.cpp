// Copyright (c) 2012 Eric Heien <emheien@ucdavis.edu>,
//                    Michael K. Sachs <mksachs@ucdavis.edu>
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

#include "QuakeLibOkada.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace quakelib {
	OpCount op;
};

std::ostream& quakelib::operator<<(std::ostream& os, const OpCount& opc) {
#ifdef COUNT_FLOPS
	os	<< "Add: " << opc._add << "; Sub: " << opc._sub
		<< "; Mult: " << opc._mult << "; Div: " << opc._div
		<< "; Sqrt: " << opc._sqrt << "; Cmp: " << opc._cmp
		<< "; And: " << opc._and << "; Abs: " << opc._abs
		<< "; Log_e: " << opc._log << "; Total: "
		<< opc._add+opc._sub+opc._mult+opc._div+opc._sqrt+opc._cmp+opc._and+opc._abs+opc._log;
#else
	os << "Flops not counted.";
#endif
	return os;
}

// [sxx sxy sxz]
// [syx syy syz]
// [szx szy szz]
quakelib::Tensor<3,3> quakelib::Okada::calc_stress_tensor(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument) {
    double				exx, eyy, ezz, exy, exz, eyz;
	double				x, y, z;
	quakelib::Tensor<3,3>	tensor;
    
	if (mu <= 0) throw std::invalid_argument("Mu must be greater than zero.");
	
	x = location[0];
	y = location[1];
	z = location[2];
	
	precalc(dip, lambda, mu);
    
    // strain tensor components
    exx = get_exx(x, y, z, c, L, W, US, UD, UT);
    eyy = get_eyy(x, y, z, c, L, W, US, UD, UT);
    ezz = get_ezz(x, y, z, c, L, W, US, UD, UT);
    exy = get_exy(x, y, z, c, L, W, US, UD, UT);
    exz = get_exz(x, y, z, c, L, W, US, UD, UT);
    eyz = get_eyz(x, y, z, c, L, W, US, UD, UT);
    
    tensor[0][0] = lambda * (eyy + ezz) + (2.0 * mu + lambda) * exx;
    tensor[1][1] = lambda * (exx + ezz) + (2.0 * mu + lambda) * eyy;
    tensor[2][2] = lambda * (exx + eyy) + (2.0 * mu + lambda) * ezz;
    tensor[0][1] = tensor[1][0] = 2.0 * mu * exy;
    tensor[0][2] = tensor[2][0] = 2.0 * mu * exz;
    tensor[1][2] = tensor[2][1] = 2.0 * mu * eyz;
	
	return tensor;
}

// [duxdx,duydx,duzdx]
quakelib::Vec<3> quakelib::Okada::calc_dudx(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument) {
	double	x, y, z;
	Vec<3>	dudx;
    
	if (mu <= 0) throw std::invalid_argument("Mu must be greater than zero.");
    
	x = location[0];
	y = location[1];
	z = location[2];
	
	precalc(dip, lambda, mu);
    
    // strain tensor components
    dudx[0] = duxdx(x, y, z, c, L, W, US, UD, UT);
    dudx[1] = duydx(x, y, z, c, L, W, US, UD, UT);
    dudx[2] = duzdx(x, y, z, c, L, W, US, UD, UT);
	
	return dudx;
}

// [duxdy,duydy,duzdy]
quakelib::Vec<3> quakelib::Okada::calc_dudy(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument) {
	double	x, y, z;
	Vec<3>	dudy;
    
	if (mu <= 0) throw std::invalid_argument("Mu must be greater than zero.");
    
	x = location[0];
	y = location[1];
	z = location[2];
	
	precalc(dip, lambda, mu);
    
    // strain tensor components
    dudy[0] = duxdy(x, y, z, c, L, W, US, UD, UT);
    dudy[1] = duydy(x, y, z, c, L, W, US, UD, UT);
    dudy[2] = duzdy(x, y, z, c, L, W, US, UD, UT);
	
	return dudy;
}

// [duxdz,duydz,duzdz]
quakelib::Vec<3> quakelib::Okada::calc_dudz(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument) {
	double	x, y, z;
	Vec<3>	dudz;
    
	if (mu <= 0) throw std::invalid_argument("Mu must be greater than zero.");
    
	x = location[0];
	y = location[1];
	z = location[2];
	
	precalc(dip, lambda, mu);
    
    // strain tensor components
    dudz[0] = duxdz(x, y, z, c, L, W, US, UD, UT);
    dudz[1] = duydz(x, y, z, c, L, W, US, UD, UT);
    dudz[2] = duzdz(x, y, z, c, L, W, US, UD, UT);
	
	return dudz;
}

quakelib::Vec<3> quakelib::Okada::calc_displacement_vector(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument) {
	double x, y, z;
    
	if (mu <= 0) throw std::invalid_argument("Mu must be greater than zero.");
    
	x = location[0];
	y = location[1];
	z = location[2];
	
	precalc(dip, lambda, mu);
    
    return quakelib::Vec<3>(ux(x, y, z, c, L, W, US, UD, UT),
						 uy(x, y, z, c, L, W, US, UD, UT),
						 uz(x, y, z, c, L, W, US, UD, UT));
}

//
// methods that apply to all calculations
void quakelib::Okada::precalc(double dip, double lambda, double mu) {
	op.reset();
	
    alpha = (lambda + mu)/(lambda + 2.0 * mu);
	
	_cos_o_dip = cos_o(dip);
	_sin_o_dip = sin_o(dip);
	_cos_o_2_dip = _cos_o_dip*_cos_o_dip;
	_sin_o_2_dip = _sin_o_dip*_sin_o_dip;
	_one_minus_alpha = 1.0 - alpha;
	_one_minus_alpha_div_two = (1.0 - alpha)/2.0;
	_one_minus_alpha_div_alpha = (1.0 - alpha)/alpha;
	_alpha_div_two = alpha/2.0;

	_nu = 0.5*lambda/(mu + lambda);
	_one_minus_two_nu = 1.0 - 2.0*_nu;

}

double quakelib::Okada::cos_o(double dip) {
    double cosval = cos(dip);
    if (fabs(cosval) < TRIG_TOLERANCE) { return 0.0; }
	else { return cosval; }
}
double quakelib::Okada::sin_o(double dip) {
    double sinval = sin(dip);
    if (fabs(sinval) < TRIG_TOLERANCE) { return 0.0; }
	else { return sinval; }
}

double quakelib::Okada::p(double y, double z, double c) {
	OP_ADD(1); OP_MULT(2);
    return y * _cos_o_dip + d(c, z) * _sin_o_dip;
}
double quakelib::Okada::q(double y, double z, double c) {
	OP_SUB(1); OP_MULT(2);
    return y * _sin_o_dip - d(c, z) * _cos_o_dip;
}
double quakelib::Okada::d(double c, double z) {
	OP_SUB(1);
    return c - z;
}
double quakelib::Okada::R(double xi, double eta, double _q) {
	OP_ADD(2); OP_MULT(3); OP_SQRT(1);
    return sqrt(xi*xi + eta*eta + _q*_q);
}
double quakelib::Okada::ytil(double _q, double eta) {
	OP_ADD(1); OP_MULT(2);
    return eta * _cos_o_dip + _q * _sin_o_dip;
}
double quakelib::Okada::dtil(double _q, double eta) {
	OP_SUB(1); OP_MULT(2);
    return eta * _sin_o_dip - _q * _cos_o_dip;
}
double quakelib::Okada::ctil(double _q, double eta, double z) {
	OP_ADD(1);
    return dtil(_q, eta) + z;
}
double quakelib::Okada::h(double _q, double z) const {
	OP_SUB(1); OP_MULT(1);
    return _q * _cos_o_dip - z;
}

bool quakelib::Okada::on_element_corner(double x, double y, double z, double c, double L, double W) const {
	OP_AND(2); OP_CMP(3);
    if (x == 0.0 && y == 0.0 && z == -c ) {
        return true;
    }
	OP_AND(2); OP_CMP(3);
    if (x == L && y == 0 && z == -c) {
        return true;
    }
	OP_AND(2); OP_CMP(3); OP_MULT(2); OP_ADD(1);
    if (x == L && y == W * _cos_o_dip && z == -c + W * _sin_o_dip) {
        return true;
    }
	OP_AND(2); OP_CMP(3); OP_MULT(2); OP_ADD(1);
    if (x == 0.0 && y == W * _cos_o_dip && z == -c + W * _sin_o_dip) {
        return true;
    }
    return false;
}
bool quakelib::Okada::cos_zero(void) {
    return (fabs(_cos_o_dip) < TRIG_TOLERANCE);
}

double quakelib::Okada::X11(double _R, double xi, double eta, double _q) const {
	OP_CMP(1); OP_ADD(1);
	double _R_plus_xi = _R + xi;
    if (!singularity3(_R_plus_xi)) {
		OP_DIV(1); OP_MULT(1);
        return 1.0 / (_R * _R_plus_xi);
    } else {
        return 0.0;
    }
}
double quakelib::Okada::X32(double _R, double xi) const {
	OP_CMP(1); OP_ADD(1);
	double _R_plus_xi = _R + xi;
    if (!singularity3(_R_plus_xi)) {
		OP_ADD(1); OP_DIV(1); OP_MULT(5);
        return (2.0 * _R + xi ) / ((_R*_R*_R) * (_R_plus_xi*_R_plus_xi));
    } else {
        return 0.0;
    }
}
double quakelib::Okada::X53(double _R, double xi) const {
	OP_CMP(1); OP_ADD(1);
	double _R_plus_xi = _R + xi;
    if (!singularity3(_R_plus_xi)) {
		OP_ADD(2); OP_DIV(1); OP_MULT(11);
		double _R2 = _R*_R;
		double _R5 = _R*_R2*_R2;
        return (8.0 * _R2 + 9.0 * _R * xi + 3.0 * xi*xi ) / (_R5 * (_R_plus_xi*_R_plus_xi*_R_plus_xi));
    } else {
        return 0.0;
    }
}
double quakelib::Okada::Y11(double _R, double eta) const {
	OP_ADD(1); OP_CMP(1);
	double _Rpeta = _R+eta;
    if (!singularity4(_Rpeta)) {
		OP_DIV(1); OP_MULT(1);
        return 1.0 / (_R * (_Rpeta));
    } else {
        return 0.0;
    }
}
double quakelib::Okada::Y32(double _R, double eta) const {
	OP_ADD(1); OP_CMP(1);
	double _Rpeta = _R+eta;
    if (!singularity4(_Rpeta)) {
		OP_ADD(1); OP_DIV(1); OP_MULT(5);
        return (2.0 * _R + eta ) / ((_R*_R*_R) * (_Rpeta*_Rpeta));
    } else {
        return 0.0;
    }
}
double quakelib::Okada::Y53(double _R, double eta) const {
	OP_ADD(1); OP_CMP(1);
	double _Rpeta = _R+eta;
    if (!singularity4(_Rpeta)) {
		OP_ADD(2); OP_DIV(1); OP_MULT(11);
		double _R2 = _R*_R;
		double _R5 = _R*_R2*_R2;
        return (8.0 * _R2 + 9.0 * _R * eta + 3.0 * eta*eta ) / (_R5 * (_Rpeta*_Rpeta*_Rpeta));
    } else {
        return 0.0;
    }
}
double quakelib::Okada::Z32(double _R, double _q, double eta, double z) const {
	OP_SUB(1); OP_DIV(1); OP_MULT(3);
    return _sin_o_dip / (_R*_R*_R) - h(_q,z) * Y32(_R,eta);
}
double quakelib::Okada::Z53(double _R, double _q, double eta, double z) const {
	OP_SUB(1); OP_DIV(1); OP_MULT(6);
    return (3.0 * _sin_o_dip) / (_R*_R*_R*_R*_R) - h(_q,z) * Y53(_R,eta);
}
double quakelib::Okada::Y0(double _R, double xi, double eta) const {
	OP_SUB(1); OP_MULT(2);
    return Y11(_R, eta) - xi*xi * Y32(_R, eta);
}
double quakelib::Okada::Z0(double _R, double xi, double eta, double _q, double z) const {
	OP_SUB(1); OP_MULT(2);
    return Z32(_R,_q,eta,z) - xi*xi * Z53(_R, _q, eta, z);
}

//
//
// strain tensor components
//
//
double quakelib::Okada::get_exx(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    return duxdx(x, y, z, c, L, W, US, UD, UT);
}
double quakelib::Okada::get_eyy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    return duydy(x, y, z, c, L, W, US, UD, UT);
}
double quakelib::Okada::get_ezz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    return duzdz(x, y, z, c, L, W, US, UD, UT);
}
double quakelib::Okada::get_exy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
	OP_ADD(1); OP_MULT(1);
    return 0.5 * ( duxdy(x,y,z,c,L,W,US,UD,UT) + duydx(x,y,z,c,L,W,US,UD,UT) );
}
double quakelib::Okada::get_exz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
	OP_ADD(1); OP_MULT(1);
    return 0.5 * ( duxdz(x,y,z,c,L,W,US,UD,UT) + duzdx(x,y,z,c,L,W,US,UD,UT) );
}
double quakelib::Okada::get_eyz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
	OP_ADD(1); OP_MULT(1);
    return 0.5 * ( duydz(x,y,z,c,L,W,US,UD,UT) + duzdy(x,y,z,c,L,W,US,UD,UT) );
}


//
//
// displacements
//
//
double quakelib::Okada::ux(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
        if (US != 0.0) {
            strike_result  = ( US/(2.0*M_PI) ) * (u1A(x,y,z,c,L,W,M_STRIKE) - u1Ah(x,y,z,c,L,W,M_STRIKE) + u1B(x,y,z,c,L,W,M_STRIKE) + z * u1C(x,y,z,c,L,W,M_STRIKE));
        }
        if (UD != 0.0) {
            dip_result     = ( UD/(2.0*M_PI) ) * (u1A(x,y,z,c,L,W,M_DIP) - u1Ah(x,y,z,c,L,W,M_DIP) + u1B(x,y,z,c,L,W,M_DIP) + z * u1C(x,y,z,c,L,W,M_DIP));
        }
        if (UT != 0.0) {
            tensile_result = ( UT/(2.0*M_PI) ) * (u1A(x,y,z,c,L,W,M_THRUST) - u1Ah(x,y,z,c,L,W,M_THRUST) + u1B(x,y,z,c,L,W,M_THRUST) + z * u1C(x,y,z,c,L,W,M_THRUST));
        }
    }
    
    return strike_result + dip_result + tensile_result;
}
double quakelib::Okada::uy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
        if (US != 0.0) {
            strike_result  = ( US/(2.0*M_PI) ) * ((u2A(x,y,z,c,L,W,M_STRIKE) - u2Ah(x,y,z,c,L,W,M_STRIKE) + u2B(x,y,z,c,L,W,M_STRIKE) + z * u2C(x,y,z,c,L,W,M_STRIKE)) * _cos_o_dip - (u3A(x,y,z,c,L,W,M_STRIKE) - u3Ah(x,y,z,c,L,W,M_STRIKE) + u3B(x,y,z,c,L,W,M_STRIKE) + z * u3C(x,y,z,c,L,W,M_STRIKE)) * _sin_o_dip);
        }
        if (UD != 0.0) {
            dip_result     = ( UD/(2.0*M_PI) ) * ((u2A(x,y,z,c,L,W,M_DIP) - u2Ah(x,y,z,c,L,W,M_DIP) + u2B(x,y,z,c,L,W,M_DIP) + z * u2C(x,y,z,c,L,W,M_DIP)) * _cos_o_dip - (u3A(x,y,z,c,L,W,M_DIP) - u3Ah(x,y,z,c,L,W,M_DIP) + u3B(x,y,z,c,L,W,M_DIP) + z * u3C(x,y,z,c,L,W,M_DIP)) * _sin_o_dip);
        }
        if (UT != 0.0) {
            tensile_result = ( UT/(2.0*M_PI) ) * ((u2A(x,y,z,c,L,W,M_THRUST) - u2Ah(x,y,z,c,L,W,M_THRUST) + u2B(x,y,z,c,L,W,M_THRUST) + z * u2C(x,y,z,c,L,W,M_THRUST)) * _cos_o_dip - (u3A(x,y,z,c,L,W,M_THRUST) - u3Ah(x,y,z,c,L,W,M_THRUST) + u3B(x,y,z,c,L,W,M_THRUST) + z * u3C(x,y,z,c,L,W,M_THRUST)) * _sin_o_dip);
        }
    }
    
    return strike_result + dip_result + tensile_result;    
}
double quakelib::Okada::uz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
        if (US != 0.0) {
            strike_result  = ( US/(2.0*M_PI) ) * ((u2A(x,y,z,c,L,W,M_STRIKE) - u2Ah(x,y,z,c,L,W,M_STRIKE) + u2B(x,y,z,c,L,W,M_STRIKE) - z * u2C(x,y,z,c,L,W,M_STRIKE)) * _sin_o_dip + (u3A(x,y,z,c,L,W,M_STRIKE) - u3Ah(x,y,z,c,L,W,M_STRIKE) + u3B(x,y,z,c,L,W,M_STRIKE) - z * u3C(x,y,z,c,L,W,M_STRIKE)) * _cos_o_dip);
        }
        if (UD != 0.0) {
            dip_result     = ( UD/(2.0*M_PI) ) * ((u2A(x,y,z,c,L,W,M_DIP) - u2Ah(x,y,z,c,L,W,M_DIP) + u2B(x,y,z,c,L,W,M_DIP) - z * u2C(x,y,z,c,L,W,M_DIP)) * _sin_o_dip + (u3A(x,y,z,c,L,W,M_DIP) - u3Ah(x,y,z,c,L,W,M_DIP) + u3B(x,y,z,c,L,W,M_DIP) - z * u3C(x,y,z,c,L,W,M_DIP)) * _cos_o_dip);
        }
        if (UT != 0.0) {
            tensile_result = ( UT/(2.0*M_PI) ) * ((u2A(x,y,z,c,L,W,M_THRUST) - u2Ah(x,y,z,c,L,W,M_THRUST) + u2B(x,y,z,c,L,W,M_THRUST) - z * u2C(x,y,z,c,L,W,M_THRUST)) * _sin_o_dip + (u3A(x,y,z,c,L,W,M_THRUST) - u3Ah(x,y,z,c,L,W,M_THRUST) + u3B(x,y,z,c,L,W,M_THRUST) - z * u3C(x,y,z,c,L,W,M_THRUST)) * _cos_o_dip);
        }
    }
    
    return strike_result + dip_result + tensile_result;    
}
//
// displacement components
// A
double quakelib::Okada::u1A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_p = p(y,z,c);
	_q = q(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return f1AS(x-L,_p-W,_q) - f1AS(x,_p-W,_q) - f1AS(x-L,_p,_q) + f1AS(x,_p,_q);
            break;
        case M_DIP:
            return f1AD(x-L,_p-W,_q) - f1AD(x,_p-W,_q) - f1AD(x-L,_p,_q) + f1AD(x,_p,_q);
            break;
        case M_THRUST:
            return f1AT(x-L,_p-W,_q) - f1AT(x,_p-W,_q) - f1AT(x-L,_p,_q) + f1AT(x,_p,_q);
            break;
        default:
            std::cout << "u1A motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::u2A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return f2AS(x-L,_p-W,_q) - f2AS(x,_p-W,_q) - f2AS(x-L,_p,_q) + f2AS(x,_p,_q);
            break;
        case M_DIP:
            return f2AD(x-L,_p-W,_q) - f2AD(x,_p-W,_q) - f2AD(x-L,_p,_q) + f2AD(x,_p,_q);
            break;
        case M_THRUST:
            return f2AT(x-L,_p-W,_q) - f2AT(x,_p-W,_q) - f2AT(x-L,_p,_q) + f2AT(x,_p,_q);
            break;
        default:
            std::cout << "u2A motion not set!!";
            return 0;
            break;
    }
    
}
double quakelib::Okada::u3A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return f3AS(x-L,_p-W,_q) - f3AS(x,_p-W,_q) - f3AS(x-L,_p,_q) + f3AS(x,_p,_q);
            break;
        case M_DIP:
            return f3AD(x-L,_p-W,_q) - f3AD(x,_p-W,_q) - f3AD(x-L,_p,_q) + f3AD(x,_p,_q);
            break;
        case M_THRUST:
            return f3AT(x-L,_p-W,_q) - f3AT(x,_p-W,_q) - f3AT(x-L,_p,_q) + f3AT(x,_p,_q);
            break;
        default:
            std::cout << "u3A motion not set!!";
            return 0;
            break;
    }
}
// Ah
double quakelib::Okada::u1Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return f1AS(x,_p,_q) - f1AS(x,_p-W,_q) - f1AS(x-L,_p,_q) + f1AS(x-L,_p-W,_q);
            break;
        case M_DIP:
            return f1AD(x,_p,_q) - f1AD(x,_p-W,_q) - f1AD(x-L,_p,_q) + f1AD(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return f1AT(x,_p,_q) - f1AT(x,_p-W,_q) - f1AT(x-L,_p,_q) + f1AT(x-L,_p-W,_q);
            break;
        default:
            std::cout << "u1Ah motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::u2Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return f2AS(x,_p,_q) - f2AS(x,_p-W,_q) - f2AS(x-L,_p,_q) + f2AS(x-L,_p-W,_q);
            break;
        case M_DIP:
            return f2AD(x,_p,_q) - f2AD(x,_p-W,_q) - f2AD(x-L,_p,_q) + f2AD(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return f2AT(x,_p,_q) - f2AT(x,_p-W,_q) - f2AT(x-L,_p,_q) + f2AT(x-L,_p-W,_q);
            break;
        default:
            std::cout << "u2Ah motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::u3Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return f3AS(x,_p,_q) - f3AS(x,_p-W,_q) - f3AS(x-L,_p,_q) + f3AS(x-L,_p-W,_q);
            break;
        case M_DIP:
            return f3AD(x,_p,_q) - f3AD(x,_p-W,_q) - f3AD(x-L,_p,_q) + f3AD(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return f3AT(x,_p,_q) - f3AT(x,_p-W,_q) - f3AT(x-L,_p,_q) + f3AT(x-L,_p-W,_q);
            break;
        default:
            std::cout << "u3Ah motion not set!!";
            return 0;
            break;
    }
}
// B
double quakelib::Okada::u1B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return f1BS(x,_p,_q) - f1BS(x,_p-W,_q) - f1BS(x-L,_p,_q) + f1BS(x-L,_p-W,_q);
            break;
        case M_DIP:
            return f1BD(x,_p,_q) - f1BD(x,_p-W,_q) - f1BD(x-L,_p,_q) + f1BD(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return f1BT(x,_p,_q) - f1BT(x,_p-W,_q) - f1BT(x-L,_p,_q) + f1BT(x-L,_p-W,_q);
            break;
        default:
            std::cout << "u1B motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::u2B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return f2BS(x,_p,_q) - f2BS(x,_p-W,_q) - f2BS(x-L,_p,_q) + f2BS(x-L,_p-W,_q);
            break;
        case M_DIP:
            return f2BD(x,_p,_q) - f2BD(x,_p-W,_q) - f2BD(x-L,_p,_q) + f2BD(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return f2BT(x,_p,_q) - f2BT(x,_p-W,_q) - f2BT(x-L,_p,_q) + f2BT(x-L,_p-W,_q);
            break;
        default:
            std::cout << "u2B motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::u3B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return f3BS(x,_p,_q) - f3BS(x,_p-W,_q) - f3BS(x-L,_p,_q) + f3BS(x-L,_p-W,_q);
            break;
        case M_DIP:
            return f3BD(x,_p,_q) - f3BD(x,_p-W,_q) - f3BD(x-L,_p,_q) + f3BD(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return f3BT(x,_p,_q) - f3BT(x,_p-W,_q) - f3BT(x-L,_p,_q) + f3BT(x-L,_p-W,_q);
            break;
        default:
            std::cout << "u3B motion not set!!";
            return 0;
            break;
    }
}
// C
double quakelib::Okada::u1C (double x, double y, double z, double c, double L, double W, MotionType motion) {
    switch (motion) {
        case M_STRIKE:
            return f1CS(x,p(y,z,c),y,z,c) - f1CS(x,p(y,z,c)-W,y,z,c) - f1CS(x-L,p(y,z,c),y,z,c) + f1CS(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_DIP:
            return f1CD(x,p(y,z,c),y,z,c) - f1CD(x,p(y,z,c)-W,y,z,c) - f1CD(x-L,p(y,z,c),y,z,c) + f1CD(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_THRUST:
            return f1CT(x,p(y,z,c),y,z,c) - f1CT(x,p(y,z,c)-W,y,z,c) - f1CT(x-L,p(y,z,c),y,z,c) + f1CT(x-L,p(y,z,c)-W,y,z,c);
            break;
        default:
            std::cout << "u1C motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::u2C (double x, double y, double z, double c, double L, double W, MotionType motion) {
    switch (motion) {
        case M_STRIKE:
            return f2CS(x,p(y,z,c),y,z,c) - f2CS(x,p(y,z,c)-W,y,z,c) - f2CS(x-L,p(y,z,c),y,z,c) + f2CS(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_DIP:
            return f2CD(x,p(y,z,c),y,z,c) - f2CD(x,p(y,z,c)-W,y,z,c) - f2CD(x-L,p(y,z,c),y,z,c) + f2CD(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_THRUST:
            return f2CT(x,p(y,z,c),y,z,c) - f2CT(x,p(y,z,c)-W,y,z,c) - f2CT(x-L,p(y,z,c),y,z,c) + f2CT(x-L,p(y,z,c)-W,y,z,c);
            break;
        default:
            std::cout << "u2C motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::u3C (double x, double y, double z, double c, double L, double W, MotionType motion) {
    switch (motion) {
        case M_STRIKE:
            return f3CS(x,p(y,z,c),y,z,c) - f3CS(x,p(y,z,c)-W,y,z,c) - f3CS(x-L,p(y,z,c),y,z,c) + f3CS(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_DIP:
            return f3CD(x,p(y,z,c),y,z,c) - f3CD(x,p(y,z,c)-W,y,z,c) - f3CD(x-L,p(y,z,c),y,z,c) + f3CD(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_THRUST:
            return f3CT(x,p(y,z,c),y,z,c) - f3CT(x,p(y,z,c)-W,y,z,c) - f3CT(x-L,p(y,z,c),y,z,c) + f3CT(x-L,p(y,z,c)-W,y,z,c);
            break;
        default:
            std::cout << "u3C motion not set!!";
            return 0;
            break;
    }
}
//
// strike fs
// A
double quakelib::Okada::f1AS(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return Theta(_q,_R,xi,eta)/2.0 + _alpha_div_two * xi * _q * Y11(_R,eta);
}
double quakelib::Okada::f2AS(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return _alpha_div_two * (_q/_R);
}
double quakelib::Okada::f3AS(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
	double _Rpeta = _R+eta;
    if (!singularity4(_Rpeta)) {
        return _one_minus_alpha_div_two * log(_R + eta) - _alpha_div_two * _q*_q * Y11(_R,eta);
    }
    else {
        return -1.0 * _one_minus_alpha_div_two * log(_R - eta) - _alpha_div_two * _q*_q * Y11(_R,eta);
    }
}
// B
double quakelib::Okada::f1BS(double xi, double eta, double _q) {
	OP_SUB(2); OP_MULT(5);
	double _R = R(xi,eta,_q);
    return -1.0 * xi * _q * Y11(_R,eta) - Theta(_q,_R,xi,eta) - _one_minus_alpha_div_alpha * I1(_R,xi,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::f2BS(double xi, double eta, double _q) {
	OP_ADD(2); OP_MULT(3); OP_DIV(2);
	double _R = R(xi,eta,_q);
    return -1.0 * (_q/_R) + _one_minus_alpha_div_alpha * (ytil(_q,eta)/(_R + dtil(_q,eta))) * _sin_o_dip;
}
double quakelib::Okada::f3BS(double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(4);
	double _R = R(xi,eta,_q);
    return _q*_q * Y11(_R,eta) - _one_minus_alpha_div_alpha * I2(_R,xi,eta,_q) * _sin_o_dip;
}
// C
double quakelib::Okada::f1CS(double xi, double eta, double y, double z, double c) {
	OP_SUB(1); OP_MULT(6);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha * xi * Y11(_R,eta) * _cos_o_dip - alpha * xi * _q * Z32(_R,_q,eta,z);
}
double quakelib::Okada::f2CS(double xi, double eta, double y, double z, double c) {
	OP_ADD(1); OP_SUB(1); OP_MULT(8); OP_DIV(2);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return _one_minus_alpha * (_cos_o_dip/_R + 2.0 * _q * Y11(_R,eta) * _sin_o_dip) - alpha * ((ctil(_q,eta,z) * _q)/_R3);
}
double quakelib::Okada::f3CS(double xi, double eta, double y, double z, double c) {
	OP_ADD(1); OP_SUB(2); OP_MULT(10); OP_DIV(1);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return _one_minus_alpha * _q * Y11(_R,eta) * _cos_o_dip - alpha * ((ctil(_q,eta,z) * eta)/_R3 - z * Y11(_R,eta) + (xi*xi) *  Z32(_R,_q,eta,z));
}
//
// dip fs
// A
double quakelib::Okada::f1AD(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return _alpha_div_two * (_q/_R);
}
double quakelib::Okada::f2AD(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return Theta(_q,_R,xi,eta)/2.0 + _alpha_div_two * eta * _q * X11(_R,xi,eta,_q);
}
double quakelib::Okada::f3AD(double xi, double eta, double _q) {
	OP_CMP(1); OP_ADD(1);
	double _R = R(xi,eta,_q);
	double _R_plus_xi = _R + xi;
    if (!singularity3(_R_plus_xi)) {
		OP_SUB(1); OP_MULT(4); OP_LOG(1);
        return _one_minus_alpha_div_two * log(_R_plus_xi) - _alpha_div_two * (_q*_q) * X11(_R,xi,eta,_q);
    }
    else {
		OP_SUB(2); OP_MULT(5); OP_LOG(1);
        return -1.0 * _one_minus_alpha_div_two * log(_R - xi) - _alpha_div_two * (_q*_q) * X11(_R,xi,eta,_q);
    }
}
// B
double quakelib::Okada::f1BD(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return -1.0 * (_q/_R) + _one_minus_alpha_div_alpha * I3(_R,eta,_q) * _sin_o_dip * _cos_o_dip;
}
double quakelib::Okada::f2BD(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return -1.0 * eta * _q * X11(_R,xi,eta,_q) - Theta(_q,_R,xi,eta) - _one_minus_alpha_div_alpha * (xi/(_R + dtil(_q,eta))) * _sin_o_dip * _cos_o_dip;
}
double quakelib::Okada::f3BD(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return (_q*_q) * X11(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * I4(_R,xi,eta,_q) * _sin_o_dip * _cos_o_dip;
}
// C
double quakelib::Okada::f1CD(double xi, double eta, double y, double z, double c) {
	OP_SUB(2); OP_MULT(7); OP_DIV(2);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return _one_minus_alpha * (_cos_o_dip/_R) - _q * Y11(_R,eta) * _sin_o_dip - alpha * ((ctil(_q,eta,z) * _q)/_R3);
}
double quakelib::Okada::f2CD(double xi, double eta, double y, double z, double c) {
	OP_SUB(1); OP_MULT(6);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha * ytil(_q,eta) * X11(_R,xi,eta,_q) - alpha * ctil(_q,eta,z) * eta * _q * X32(_R,xi);
}
double quakelib::Okada::f3CD(double xi, double eta, double y, double z, double c) {
	OP_SUB(3); OP_MULT(8);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return -1.0 * dtil(_q,eta) * X11(_R,xi,eta,_q) - xi * Y11(_R,eta) * _sin_o_dip - alpha * ctil(_q,eta,z) * (X11(_R,xi,eta,_q) - (_q*_q) * X32(_R,xi));
}
//
// tensile fs
// A
double quakelib::Okada::f1AT(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
	double _Rpeta = _R+eta;
    if (!singularity4(_Rpeta)) {
        return -1.0 * _one_minus_alpha_div_two * log(_R + eta) - _alpha_div_two * (_q*_q) * Y11(_R,eta);
    }
    else {
        return _one_minus_alpha_div_two * log(_R - eta) - _alpha_div_two * (_q*_q) * Y11(_R,eta);
    }
}
double quakelib::Okada::f2AT(double xi, double eta, double _q) {
	OP_CMP(1); OP_ADD(1);
	double _R = R(xi,eta,_q);
	double _R_plus_xi = _R + xi;
    if (!singularity3(_R_plus_xi)) {
		OP_SUB(1); OP_MULT(5); OP_LOG(1);
        return -1.0 * _one_minus_alpha_div_two * log(_R_plus_xi) - _alpha_div_two * (_q*_q) * X11(_R,xi,eta,_q);
    } else {
		OP_SUB(2); OP_MULT(4); OP_LOG(1);
        return _one_minus_alpha_div_two * log(_R - xi) - _alpha_div_two * (_q*_q) * X11(_R,xi,eta,_q);
    }
}
double quakelib::Okada::f3AT(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return Theta(_q,_R,xi,eta)/2.0 - _alpha_div_two * _q * (eta * X11(_R,xi,eta,_q) + xi * Y11(_R,eta));
}
// B
double quakelib::Okada::f1BT(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return (_q*_q) * Y11(_R,eta) - _one_minus_alpha_div_alpha * I3(_R,eta,_q) * _sin_o_2_dip;
}
double quakelib::Okada::f2BT(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return (_q*_q) * X11(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * (xi/(_R + dtil(_q,eta))) * _sin_o_2_dip;
}
double quakelib::Okada::f3BT(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return _q * (eta * X11(_R,xi,eta,_q) + xi * Y11(_R,eta)) - Theta(_q,_R,xi,eta) - _one_minus_alpha_div_alpha * I4(_R,xi,eta,_q) * _sin_o_2_dip;
}
// C
double quakelib::Okada::f1CT(double xi, double eta, double y, double z, double c) {
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha * (_sin_o_dip/_R + _q * Y11(_R,eta) * _cos_o_dip) - alpha * (z * Y11(_R,eta) - (_q*_q) * Z32(_R,_q,eta,z));
}
double quakelib::Okada::f2CT(double xi, double eta, double y, double z, double c) {
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha * 2.0 * xi * Y11(_R,eta) * _sin_o_dip + dtil(_q,eta) * X11(_R,xi,eta,_q) - alpha * ctil(_q,eta,z) * (X11(_R,xi,eta,_q) - (_q*_q) * X32(_R,xi));
}
double quakelib::Okada::f3CT(double xi, double eta, double y, double z, double c) {
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha * (ytil(_q,eta) * X11(_R,xi,eta,_q) + xi * Y11(_R,eta) * _cos_o_dip) + alpha * _q * (ctil(_q,eta,z) * eta * X32(_R,xi) + xi * Z32(_R,_q,eta,z));
}
//
// displacement globals
double quakelib::Okada::Theta(double _q, double _R, double xi, double eta) {
    if (!singularity1(_q)) {
        return atan((xi * eta)/(_q * _R));
    } else {
        return 0;
    }
}
double quakelib::Okada::X(double xi, double _q) {
	OP_ADD(1); OP_MULT(2); OP_SQRT(1);
    return sqrt(xi*xi + _q*_q);
}
double quakelib::Okada::I1(double _R, double xi, double eta, double _q) {
    return -1.0 * (xi/(_R + dtil(_q,eta))) * _cos_o_dip - I4(_R,xi,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::I2(double _R, double xi, double eta, double _q) {
    return log(_R + dtil(_q,eta)) + I3(_R,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::I3(double _R, double eta, double _q) {
	double _dtil = dtil(_q,eta);
	double _R_plus_dtil = _R+_dtil;
	double _R_plus_eta = _R+eta;
    if (!cos_zero()) {
        if (!singularity4(_R_plus_eta)) {
            return (1.0/_cos_o_dip) * (ytil(_q,eta)/(_R_plus_dtil)) - (1.0/_cos_o_2_dip) * (log(_R_plus_eta) - _sin_o_dip * log(_R_plus_dtil));
        }
        else {
            return (1.0/_cos_o_dip) * (ytil(_q,eta)/(_R_plus_dtil)) - (1.0/_cos_o_2_dip) * (-1.0 * log(_R - eta) - _sin_o_dip * log(_R_plus_dtil));
        }
    } 
    else {
		double _R_plus_dtil2 = _R_plus_dtil*_R_plus_dtil;
        if (!singularity4(_R_plus_eta)) {
            return 0.5 * (eta/(_R_plus_dtil) + (ytil(_q,eta) * _q)/_R_plus_dtil2 - log(_R_plus_eta));
        }
        else {
            return 0.5 * (eta/(_R_plus_dtil) + (ytil(_q,eta) * _q)/_R_plus_dtil2 + log(_R - eta));
        }
    }
}
double quakelib::Okada::I4(double _R, double xi, double eta, double _q) {
    if (!singularity2(xi)) {
		double _X = X(xi,_q);
		double _dtil = dtil(_q,eta);
		double _R_plus_dtil = _R+_dtil;
        if (!cos_zero()) {
            return (_sin_o_dip/_cos_o_dip) * (xi/(_R_plus_dtil)) + (2.0/_cos_o_2_dip) * atan((eta * (_X + _q * _cos_o_dip) + _X * (_R + _X) * _sin_o_dip)/(xi * (_R + _X) * _cos_o_dip));
        } else {
			double _R_plus_dtil2 = _R_plus_dtil*_R_plus_dtil;
            return 0.5 * ((xi * ytil(_q,eta))/_R_plus_dtil2);
        }
    } else {
        return 0;
    }
    
}


//
//
// displacement derivatives
//
//
//
// dx
//
double quakelib::Okada::duxdx(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
	OP_ADD(2); OP_CMP(1); OP_AND(1);
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
		OP_CMP(3);
        if (US != 0.0) {
			OP_ADD(2); OP_SUB(1); OP_MULT(2); OP_DIV(1);
            strike_result  = ( US/(2.0*M_PI) ) * (j1A(x,y,z,c,L,W,M_STRIKE) - j1Ah(x,y,z,c,L,W,M_STRIKE) + j1B(x,y,z,c,L,W,M_STRIKE) + z * j1C(x,y,z,c,L,W,M_STRIKE));
        }
        if (UD != 0.0) {
			OP_ADD(2); OP_SUB(1); OP_MULT(2); OP_DIV(1);
            dip_result     = ( UD/(2.0*M_PI) ) * (j1A(x,y,z,c,L,W,M_DIP) - j1Ah(x,y,z,c,L,W,M_DIP) + j1B(x,y,z,c,L,W,M_DIP) + z * j1C(x,y,z,c,L,W,M_DIP));
        }
        if (UT != 0.0) {
			OP_ADD(2); OP_SUB(1); OP_MULT(2); OP_DIV(1);
            tensile_result = ( UT/(2.0*M_PI) ) * (j1A(x,y,z,c,L,W,M_THRUST) - j1Ah(x,y,z,c,L,W,M_THRUST) + j1B(x,y,z,c,L,W,M_THRUST) + z * j1C(x,y,z,c,L,W,M_THRUST));
        }
    }
    
    return strike_result + dip_result + tensile_result;
}
double quakelib::Okada::duydx(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
	OP_ADD(2); OP_CMP(1); OP_AND(1);
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
		OP_CMP(3);
        if (US != 0.0) {
			OP_ADD(4); OP_SUB(3); OP_MULT(5); OP_DIV(1);
            strike_result  = ( US/(2.0*M_PI) ) * ((j2A(x,y,z,c,L,W,M_STRIKE) - j2Ah(x,y,z,c,L,W,M_STRIKE) + j2B(x,y,z,c,L,W,M_STRIKE) + z * j2C(x,y,z,c,L,W,M_STRIKE)) * _cos_o_dip - (j3A(x,y,z,c,L,W,M_STRIKE) - j3Ah(x,y,z,c,L,W,M_STRIKE) + j3B(x,y,z,c,L,W,M_STRIKE) + z * j3C(x,y,z,c,L,W,M_STRIKE)) * _sin_o_dip);
        }
        if (UD != 0.0) {
			OP_ADD(4); OP_SUB(3); OP_MULT(5); OP_DIV(1);
            dip_result     = ( UD/(2.0*M_PI) ) * ((j2A(x,y,z,c,L,W,M_DIP) - j2Ah(x,y,z,c,L,W,M_DIP) + j2B(x,y,z,c,L,W,M_DIP) + z * j2C(x,y,z,c,L,W,M_DIP)) * _cos_o_dip - (j3A(x,y,z,c,L,W,M_DIP) - j3Ah(x,y,z,c,L,W,M_DIP) + j3B(x,y,z,c,L,W,M_DIP) + z * j3C(x,y,z,c,L,W,M_DIP)) * _sin_o_dip);
        }
        if (UT != 0.0) {
			OP_ADD(4); OP_SUB(3); OP_MULT(5); OP_DIV(1);
            tensile_result = ( UT/(2.0*M_PI) ) * ((j2A(x,y,z,c,L,W,M_THRUST) - j2Ah(x,y,z,c,L,W,M_THRUST) + j2B(x,y,z,c,L,W,M_THRUST) + z * j2C(x,y,z,c,L,W,M_THRUST)) * _cos_o_dip - (j3A(x,y,z,c,L,W,M_THRUST) - j3Ah(x,y,z,c,L,W,M_THRUST) + j3B(x,y,z,c,L,W,M_THRUST) + z * j3C(x,y,z,c,L,W,M_THRUST)) * _sin_o_dip);
        }
    }
    
    return strike_result + dip_result + tensile_result;    
}
double quakelib::Okada::duzdx(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
	OP_ADD(2); OP_CMP(1); OP_AND(1);
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
		OP_CMP(3);
        if (US != 0.0) {
			OP_ADD(3); OP_SUB(4); OP_MULT(5); OP_DIV(1);
            strike_result  = ( US/(2.0*M_PI) ) * ((j2A(x,y,z,c,L,W,M_STRIKE) - j2Ah(x,y,z,c,L,W,M_STRIKE) + j2B(x,y,z,c,L,W,M_STRIKE) - z * j2C(x,y,z,c,L,W,M_STRIKE)) * _sin_o_dip + (j3A(x,y,z,c,L,W,M_STRIKE) - j3Ah(x,y,z,c,L,W,M_STRIKE) + j3B(x,y,z,c,L,W,M_STRIKE) - z * j3C(x,y,z,c,L,W,M_STRIKE)) * _cos_o_dip);
        }
        if (UD != 0.0) {
			OP_ADD(3); OP_SUB(4); OP_MULT(5); OP_DIV(1);
            dip_result     = ( UD/(2.0*M_PI) ) * ((j2A(x,y,z,c,L,W,M_DIP) - j2Ah(x,y,z,c,L,W,M_DIP) + j2B(x,y,z,c,L,W,M_DIP) - z * j2C(x,y,z,c,L,W,M_DIP)) * _sin_o_dip + (j3A(x,y,z,c,L,W,M_DIP) - j3Ah(x,y,z,c,L,W,M_DIP) + j3B(x,y,z,c,L,W,M_DIP) - z * j3C(x,y,z,c,L,W,M_DIP)) * _cos_o_dip);
        }
        if (UT != 0.0) {
			OP_ADD(3); OP_SUB(4); OP_MULT(5); OP_DIV(1);
            tensile_result = ( UT/(2.0*M_PI) ) * ((j2A(x,y,z,c,L,W,M_THRUST) - j2Ah(x,y,z,c,L,W,M_THRUST) + j2B(x,y,z,c,L,W,M_THRUST) - z * j2C(x,y,z,c,L,W,M_THRUST)) * _sin_o_dip + (j3A(x,y,z,c,L,W,M_THRUST) - j3Ah(x,y,z,c,L,W,M_THRUST) + j3B(x,y,z,c,L,W,M_THRUST) - z * j3C(x,y,z,c,L,W,M_THRUST)) * _cos_o_dip);
        }
    }
    
    return strike_result + dip_result + tensile_result;    
}
//
// dx components
// A
double quakelib::Okada::j1A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df1ASdx(x-L,_p-W,_q) - df1ASdx(x,_p-W,_q) - df1ASdx(x-L,_p,_q) + df1ASdx(x,_p,_q);
            break;
        case M_DIP:
            return df1ADdx(x-L,_p-W,_q) - df1ADdx(x,_p-W,_q) - df1ADdx(x-L,_p,_q) + df1ADdx(x,_p,_q);
            break;
        case M_THRUST:
            return df1ATdx(x-L,_p-W,_q) - df1ATdx(x,_p-W,_q) - df1ATdx(x-L,_p,_q) + df1ATdx(x,_p,_q);
            break;
        default:
            std::cout << "j1A motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::j2A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df2ASdx(x-L,_p-W,_q) - df2ASdx(x,_p-W,_q) - df2ASdx(x-L,_p,_q) + df2ASdx(x,_p,_q);
            break;
        case M_DIP:
            return df2ADdx(x-L,_p-W,_q) - df2ADdx(x,_p-W,_q) - df2ADdx(x-L,_p,_q) + df2ADdx(x,_p,_q);
            break;
        case M_THRUST:
            return df2ATdx(x-L,_p-W,_q) - df2ATdx(x,_p-W,_q) - df2ATdx(x-L,_p,_q) + df2ATdx(x,_p,_q);
            break;
        default:
            std::cout << "j2A motion not set!!";
            return 0;
            break;
    }
    
}
double quakelib::Okada::j3A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df3ASdx(x-L,_p-W,_q) - df3ASdx(x,_p-W,_q) - df3ASdx(x-L,_p,_q) + df3ASdx(x,_p,_q);
            break;
        case M_DIP:
            return df3ADdx(x-L,_p-W,_q) - df3ADdx(x,_p-W,_q) - df3ADdx(x-L,_p,_q) + df3ADdx(x,_p,_q);
            break;
        case M_THRUST:
            return df3ATdx(x-L,_p-W,_q) - df3ATdx(x,_p-W,_q) - df3ATdx(x-L,_p,_q) + df3ATdx(x,_p,_q);
            break;
        default:
            std::cout << "j3A motion not set!!";
            return 0;
            break;
    }
}
// Ah
double quakelib::Okada::j1Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return df1ASdx(x,_p,_q) - df1ASdx(x,_p-W,_q) - df1ASdx(x-L,_p,_q) + df1ASdx(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df1ADdx(x,_p,_q) - df1ADdx(x,_p-W,_q) - df1ADdx(x-L,_p,_q) + df1ADdx(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df1ATdx(x,_p,_q) - df1ATdx(x,_p-W,_q) - df1ATdx(x-L,_p,_q) + df1ATdx(x-L,_p-W,_q);
            break;
        default:
            std::cout << "j1Ah motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::j2Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return df2ASdx(x,_p,_q) - df2ASdx(x,_p-W,_q) - df2ASdx(x-L,_p,_q) + df2ASdx(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df2ADdx(x,_p,_q) - df2ADdx(x,_p-W,_q) - df2ADdx(x-L,_p,_q) + df2ADdx(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df2ATdx(x,_p,_q) - df2ATdx(x,_p-W,_q) - df2ATdx(x-L,_p,_q) + df2ATdx(x-L,_p-W,_q);
            break;
        default:
            std::cout << "j2Ah motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::j3Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return df3ASdx(x,_p,_q) - df3ASdx(x,_p-W,_q) - df3ASdx(x-L,_p,_q) + df3ASdx(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df3ADdx(x,_p,_q) - df3ADdx(x,_p-W,_q) - df3ADdx(x-L,_p,_q) + df3ADdx(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df3ATdx(x,_p,_q) - df3ATdx(x,_p-W,_q) - df3ATdx(x-L,_p,_q) + df3ATdx(x-L,_p-W,_q);
            break;
        default:
            std::cout << "j3Ah motion not set!!";
            return 0;
            break;
    }
}
// B
double quakelib::Okada::j1B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df1BSdx(x,_p,_q) - df1BSdx(x,_p-W,_q) - df1BSdx(x-L,_p,_q) + df1BSdx(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df1BDdx(x,_p,_q) - df1BDdx(x,_p-W,_q) - df1BDdx(x-L,_p,_q) + df1BDdx(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df1BTdx(x,_p,_q) - df1BTdx(x,_p-W,_q) - df1BTdx(x-L,_p,_q) + df1BTdx(x-L,_p-W,_q);
            break;
        default:
            std::cout << "j1B motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::j2B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df2BSdx(x,_p,_q) - df2BSdx(x,_p-W,_q) - df2BSdx(x-L,_p,_q) + df2BSdx(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df2BDdx(x,_p,_q) - df2BDdx(x,_p-W,_q) - df2BDdx(x-L,_p,_q) + df2BDdx(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df2BTdx(x,_p,_q) - df2BTdx(x,_p-W,_q) - df2BTdx(x-L,_p,_q) + df2BTdx(x-L,_p-W,_q);
            break;
        default:
            std::cout << "j2B motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::j3B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df3BSdx(x,_p,_q) - df3BSdx(x,_p-W,_q) - df3BSdx(x-L,_p,_q) + df3BSdx(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df3BDdx(x,_p,_q) - df3BDdx(x,_p-W,_q) - df3BDdx(x-L,_p,_q) + df3BDdx(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df3BTdx(x,_p,_q) - df3BTdx(x,_p-W,_q) - df3BTdx(x-L,_p,_q) + df3BTdx(x-L,_p-W,_q);
            break;
        default:
            std::cout << "j3B motion not set!!";
            return 0;
            break;
    }
}
// C
double quakelib::Okada::j1C (double x, double y, double z, double c, double L, double W, MotionType motion) {
    switch (motion) {
        case M_STRIKE:
            return df1CSdx(x,p(y,z,c),y,z,c) - df1CSdx(x,p(y,z,c)-W,y,z,c) - df1CSdx(x-L,p(y,z,c),y,z,c) + df1CSdx(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_DIP:
            return df1CDdx(x,p(y,z,c),y,z,c) - df1CDdx(x,p(y,z,c)-W,y,z,c) - df1CDdx(x-L,p(y,z,c),y,z,c) + df1CDdx(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_THRUST:
            return df1CTdx(x,p(y,z,c),y,z,c) - df1CTdx(x,p(y,z,c)-W,y,z,c) - df1CTdx(x-L,p(y,z,c),y,z,c) + df1CTdx(x-L,p(y,z,c)-W,y,z,c);
            break;
        default:
            std::cout << "j1C motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::j2C (double x, double y, double z, double c, double L, double W, MotionType motion) {
    switch (motion) {
        case M_STRIKE:
            return df2CSdx(x,p(y,z,c),y,z,c) - df2CSdx(x,p(y,z,c)-W,y,z,c) - df2CSdx(x-L,p(y,z,c),y,z,c) + df2CSdx(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_DIP:
            return df2CDdx(x,p(y,z,c),y,z,c) - df2CDdx(x,p(y,z,c)-W,y,z,c) - df2CDdx(x-L,p(y,z,c),y,z,c) + df2CDdx(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_THRUST:
            return df2CTdx(x,p(y,z,c),y,z,c) - df2CTdx(x,p(y,z,c)-W,y,z,c) - df2CTdx(x-L,p(y,z,c),y,z,c) + df2CTdx(x-L,p(y,z,c)-W,y,z,c);
            break;
        default:
            std::cout << "j2C motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::j3C (double x, double y, double z, double c, double L, double W, MotionType motion) {
    switch (motion) {
        case M_STRIKE:
            return df3CSdx(x,p(y,z,c),y,z,c) - df3CSdx(x,p(y,z,c)-W,y,z,c) - df3CSdx(x-L,p(y,z,c),y,z,c) + df3CSdx(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_DIP:
            return df3CDdx(x,p(y,z,c),y,z,c) - df3CDdx(x,p(y,z,c)-W,y,z,c) - df3CDdx(x-L,p(y,z,c),y,z,c) + df3CDdx(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_THRUST:
            return df3CTdx(x,p(y,z,c),y,z,c) - df3CTdx(x,p(y,z,c)-W,y,z,c) - df3CTdx(x-L,p(y,z,c),y,z,c) + df3CTdx(x-L,p(y,z,c)-W,y,z,c);
            break;
        default:
            std::cout << "j3C motion not set!!";
            return 0;
            break;
    }
}
//
// dx strike dfs
// A
double quakelib::Okada::df1ASdx(double xi, double eta, double _q) {
	OP_SUB(2); OP_MULT(7);
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha_div_two * _q * Y11(_R,eta) - _alpha_div_two * xi*xi * _q * Y32(_R,eta);
}
double quakelib::Okada::df2ASdx(double xi, double eta, double _q) {
	OP_MULT(5); OP_DIV(1);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return -1.0 * _alpha_div_two * (xi * _q)/_R3;
}
double quakelib::Okada::df3ASdx(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(6);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * xi * Y11(_R,eta) + _alpha_div_two * xi * (_q*_q) * Y32(_R,eta);
}
// B
double quakelib::Okada::df1BSdx(double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(5);
	double _R = R(xi,eta,_q);
    return xi*xi * _q * Y32(_R,eta) - _one_minus_alpha_div_alpha * J1(_R,xi,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::df2BSdx(double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(5); OP_DIV(1);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return (xi * _q)/_R3 - _one_minus_alpha_div_alpha * J2(_R,xi,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::df3BSdx(double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(6);
	double _R = R(xi,eta,_q);
    return -1.0 * xi * (_q*_q) * Y32(_R,eta) - _one_minus_alpha_div_alpha * J3(_R,xi,eta,_q) * _sin_o_dip;
}
// C
double quakelib::Okada::df1CSdx(double xi, double eta, double y, double z, double c) {
	OP_SUB(1); OP_MULT(4);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha * Y0(_R,xi,eta) * _cos_o_dip - alpha * _q * Z0(_R,xi,eta,_q,z);
}
double quakelib::Okada::df2CSdx(double xi, double eta, double y, double z, double c) {
	OP_ADD(2); OP_SUB(1); OP_MULT(14); OP_DIV(2);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return -1.0 * _one_minus_alpha * xi * (_cos_o_dip/_R3 + 2.0 * _q * Y32(_R,eta) * _sin_o_dip) + alpha * (3.0 * ctil(_q,eta,z) * xi * _q)/_R5;
}
double quakelib::Okada::df3CSdx(double xi, double eta, double y, double z, double c) {
	OP_ADD(1); OP_SUB(3); OP_MULT(14); OP_DIV(1);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R5 = _R*_R*_R*_R*_R;
    return -1.0 * _one_minus_alpha * xi * _q * Y32(_R,eta) * _cos_o_dip + alpha * xi * ((3.0 * ctil(_q,eta,z) * eta)/_R5 - z * Y32(_R,eta) - Z32(_R,_q,eta,z) - Z0(_R,xi,eta,_q,z));
}
//
// dx dip dfs
// A
double quakelib::Okada::df1ADdx(double xi, double eta, double _q) {
	OP_MULT(5); OP_DIV(1);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return -1.0 * _alpha_div_two * (xi * _q)/_R3;
}
double quakelib::Okada::df2ADdx(double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(6); OP_DIV(2);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return -1.0 * (_q/2.0) * Y11(_R,eta) - _alpha_div_two * (eta * _q)/_R3;
}
double quakelib::Okada::df3ADdx(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(5); OP_DIV(2);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return _one_minus_alpha_div_two * 1.0/_R + _alpha_div_two * (_q*_q)/_R3;
}
// B
double quakelib::Okada::df1BDdx(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(6); OP_DIV(1);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return (xi * _q)/_R3 + _one_minus_alpha_div_alpha * J4(_R,xi,eta,_q) * _sin_o_dip * _cos_o_dip;
}
double quakelib::Okada::df2BDdx(double xi, double eta, double _q) {
	OP_ADD(2); OP_MULT(7); OP_DIV(1);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return (eta * _q)/_R3 + _q * Y11(_R,eta) + _one_minus_alpha_div_alpha * J5(_R,eta,_q) * _sin_o_dip * _cos_o_dip;
}
double quakelib::Okada::df3BDdx(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(7); OP_DIV(1);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return -1.0 * ((_q*_q)/_R3) + _one_minus_alpha_div_alpha * J6(_R,xi,eta,_q) * _sin_o_dip * _cos_o_dip;
}
// C
double quakelib::Okada::df1CDdx(double xi, double eta, double y, double z, double c) {
	OP_ADD(2); OP_MULT(14); OP_DIV(2);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return -1.0 * _one_minus_alpha * (xi/_R3) * _cos_o_dip + xi * _q * Y32(_R,eta) * _sin_o_dip + alpha * ((3.0 * ctil(_q,eta,z) * xi * _q)/_R5);
}
double quakelib::Okada::df2CDdx(double xi, double eta, double y, double z, double c) {
	OP_ADD(1); OP_MULT(10); OP_DIV(2);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return -1.0 * _one_minus_alpha * (ytil(_q,eta)/_R3) + alpha * ((3.0 * ctil(_q,eta,z) * eta * _q)/_R5);
}
double quakelib::Okada::df3CDdx(double xi, double eta, double y, double z, double c) {
	OP_ADD(1); OP_MULT(7); OP_SUB(2); OP_DIV(3);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R2 = _R*_R;
	double _R3 = _R2*_R;
    return dtil(_q,eta)/_R3 - Y0(_R,xi,eta) * _sin_o_dip + alpha * (ctil(_q,eta,z)/_R3) * (1.0 - (3.0 * (_q*_q))/_R2);
}
//
// dx tensile dfs
// A
double quakelib::Okada::df1ATdx(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha_div_two * xi * Y11(_R,eta) + _alpha_div_two * xi * (_q*_q) * Y32(_R,eta);
}
double quakelib::Okada::df2ATdx(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return -1.0 * _one_minus_alpha_div_two * 1.0/_R + _alpha_div_two * (_q*_q)/_R3;
}
double quakelib::Okada::df3ATdx(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha_div_two * _q * Y11(_R,eta) - _alpha_div_two * (_q*_q*_q) * Y32(_R,eta);
}
// B
double quakelib::Okada::df1BTdx(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return -xi * (_q*_q) * Y32(_R,eta) - _one_minus_alpha_div_alpha * J4(_R,xi,eta,_q) * _sin_o_2_dip;
}
double quakelib::Okada::df2BTdx(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
    return -1.0 * ((_q*_q)/_R3) - _one_minus_alpha_div_alpha * J5(_R,eta,_q) * _sin_o_2_dip;
}
double quakelib::Okada::df3BTdx(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return (_q*_q*_q) * Y32(_R,eta) - _one_minus_alpha_div_alpha * J6(_R,xi,eta,_q) * _sin_o_2_dip;
}
// C
double quakelib::Okada::df1CTdx(double xi, double eta, double y, double z, double c) {
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return _one_minus_alpha * (xi/_R3) * _sin_o_dip + xi * _q * Y32(_R,eta) * _cos_o_dip + alpha * xi *((3.0 * ctil(_q,eta,z) * eta)/_R5 - 2.0 * Z32(_R,_q,eta,z) - Z0(_R,xi,eta,_q,z));
}
double quakelib::Okada::df2CTdx(double xi, double eta, double y, double z, double c) {
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R2 = _R*_R;
	double _R3 = _R2*_R;
    return _one_minus_alpha * 2.0 * Y0(_R,xi,eta) * _sin_o_dip - dtil(_q,eta)/_R3 + alpha * (ctil(_q,eta,z)/_R3) * (1.0 - (3.0 * (_q*_q))/_R2);
}
double quakelib::Okada::df3CTdx(double xi, double eta, double y, double z, double c) {
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return -1.0 * _one_minus_alpha * (ytil(_q,eta)/_R3 - Y0(_R,xi,eta) * _cos_o_dip) - alpha * ((3.0 * ctil(_q,eta,z) * eta * _q)/_R5 - _q * Z0(_R,xi,eta,_q,z));
}
//
// dx globals
double quakelib::Okada::J1(double _R, double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(2);
    return J5(_R,eta,_q) * _cos_o_dip - J6(_R,xi,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::J2(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(2); OP_DIV(1);
	double _dtil = dtil(_q,eta);
    return ((xi * ytil(_q,eta))/(_R + _dtil)) * D11(_R,_dtil);
}
double quakelib::Okada::J3(double _R, double xi, double eta, double _q) {
    if (!cos_zero()) {
		OP_SUB(1); OP_MULT(2); OP_DIV(1); OP_CMP(1);
        return (1.0/_cos_o_dip) * (K1(_R,xi,eta,_q) - J2(_R,xi,eta,_q) * _sin_o_dip);
    } else {
		OP_SUB(1); OP_ADD(1); OP_MULT(4); OP_DIV(1); OP_CMP(1);
		double _dtil = dtil(_q,eta);
        return -1.0 * (xi/pow(_R+_dtil,2.0)) * ((_q*_q) * D11(_R,_dtil) - 0.5);
    }
}

double quakelib::Okada::J4(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(4); OP_SUB(1);
    return -1.0 * xi * Y11(_R,eta) - J2(_R,xi,eta,_q) * _cos_o_dip + J3(_R,xi,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::J5(double _R, double eta, double _q) {
	OP_ADD(2); OP_MULT(3); OP_DIV(1);
	double _dtil = dtil(_q,eta);
	double _ytil = ytil(_q,eta);
    return -1.0 * (_dtil + (_ytil*_ytil)/(_R + _dtil)) * D11(_R,_dtil);
}
double quakelib::Okada::J6(double _R, double xi, double eta, double _q) {
    if (!cos_zero()) {
		OP_SUB(1); OP_MULT(2); OP_DIV(1); OP_CMP(1);
        return (1.0/_cos_o_dip) * (K3(_R,xi,eta,_q) - J5(_R,eta,_q) * _sin_o_dip);
    }
    else {
		OP_ADD(1); OP_SUB(1); OP_MULT(5); OP_DIV(1); OP_CMP(1);
		double _dtil = dtil(_q,eta);
		double _Rpdtil = _R+_dtil;
        return -1.0 * (ytil(_q,eta)/(_Rpdtil*_Rpdtil)) * (xi*xi * D11(_R,_dtil) - 0.5);
    } 
}
double quakelib::Okada::K1(double _R, double xi, double eta, double _q) {
	double _dtil = dtil(_q,eta);
    if (!cos_zero()) {
		OP_SUB(1); OP_MULT(2); OP_DIV(1); OP_CMP(1);
        return (xi/_cos_o_dip) * (D11(_R,_dtil) - Y11(_R,eta) * _sin_o_dip);
    }
    else {
		OP_ADD(1); OP_MULT(2); OP_DIV(1); OP_CMP(1);
        return ((xi * _q)/(_R+_dtil)) * D11(_R,_dtil);
    } 
}
double quakelib::Okada::K2(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(1); OP_DIV(1);
    return 1.0/_R + K3(_R,xi,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::K3(double _R, double xi, double eta, double _q) {
	double _dtil = dtil(_q,eta);
    if (!cos_zero()) {
		OP_SUB(1); OP_MULT(3); OP_DIV(1); OP_CMP(1);
        return (1.0/_cos_o_dip) * (_q * Y11(_R,eta) - ytil(_q,eta) * D11(_R,_dtil));
    }
    else {
		OP_ADD(1); OP_SUB(1); OP_MULT(3); OP_DIV(1); OP_CMP(1);
        return (_sin_o_dip/(_R + _dtil)) * (xi*xi * D11(_R,_dtil) - 1);
    }  
}
double quakelib::Okada::K4(double _R, double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(3);
    return xi * Y11(_R,eta) * _cos_o_dip - K1(_R,xi,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::D11(double _R, double _dtil) {
	OP_ADD(1); OP_MULT(1); OP_DIV(1);
    return 1.0/(_R * (_R + _dtil));
}

//
// dy
//
double quakelib::Okada::duxdy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
        if (US != 0.0) {
            strike_result  = ( US/(2.0*M_PI) ) * (k1A(x,y,z,c,L,W,M_STRIKE) - k1Ah(x,y,z,c,L,W,M_STRIKE) + k1B(x,y,z,c,L,W,M_STRIKE) + z * k1C(x,y,z,c,L,W,M_STRIKE));
        }
        if (UD != 0.0) {
            dip_result     = ( UD/(2.0*M_PI) ) * (k1A(x,y,z,c,L,W,M_DIP) - k1Ah(x,y,z,c,L,W,M_DIP) + k1B(x,y,z,c,L,W,M_DIP) + z * k1C(x,y,z,c,L,W,M_DIP));
        }
        if (UT != 0.0) {
            tensile_result = ( UT/(2.0*M_PI) ) * (k1A(x,y,z,c,L,W,M_THRUST) - k1Ah(x,y,z,c,L,W,M_THRUST) + k1B(x,y,z,c,L,W,M_THRUST) + z * k1C(x,y,z,c,L,W,M_THRUST));
        }
    }
    
    return strike_result + dip_result + tensile_result;
}
double quakelib::Okada::duydy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
        if (US != 0.0) {
            strike_result  = ( US/(2.0*M_PI) ) * ((k2A(x,y,z,c,L,W,M_STRIKE) - k2Ah(x,y,z,c,L,W,M_STRIKE) + k2B(x,y,z,c,L,W,M_STRIKE) + z * k2C(x,y,z,c,L,W,M_STRIKE)) * _cos_o_dip - (k3A(x,y,z,c,L,W,M_STRIKE) - k3Ah(x,y,z,c,L,W,M_STRIKE) + k3B(x,y,z,c,L,W,M_STRIKE) + z * k3C(x,y,z,c,L,W,M_STRIKE)) * _sin_o_dip);
        }
        if (UD != 0.0) {
            dip_result     = ( UD/(2.0*M_PI) ) * ((k2A(x,y,z,c,L,W,M_DIP) - k2Ah(x,y,z,c,L,W,M_DIP) + k2B(x,y,z,c,L,W,M_DIP) + z * k2C(x,y,z,c,L,W,M_DIP)) * _cos_o_dip - (k3A(x,y,z,c,L,W,M_DIP) - k3Ah(x,y,z,c,L,W,M_DIP) + k3B(x,y,z,c,L,W,M_DIP) + z * k3C(x,y,z,c,L,W,M_DIP)) * _sin_o_dip);
        }
        if (UT != 0.0) {
            tensile_result = ( UT/(2.0*M_PI) ) * ((k2A(x,y,z,c,L,W,M_THRUST) - k2Ah(x,y,z,c,L,W,M_THRUST) + k2B(x,y,z,c,L,W,M_THRUST) + z * k2C(x,y,z,c,L,W,M_THRUST)) * _cos_o_dip - (k3A(x,y,z,c,L,W,M_THRUST) - k3Ah(x,y,z,c,L,W,M_THRUST) + k3B(x,y,z,c,L,W,M_THRUST) + z * k3C(x,y,z,c,L,W,M_THRUST)) * _sin_o_dip);
        }
    }
    
    return strike_result + dip_result + tensile_result; 
}
double quakelib::Okada::duzdy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
        if (US != 0.0) {
            strike_result  = ( US/(2.0*M_PI) ) * ((k2A(x,y,z,c,L,W,M_STRIKE) - k2Ah(x,y,z,c,L,W,M_STRIKE) + k2B(x,y,z,c,L,W,M_STRIKE) - z * k2C(x,y,z,c,L,W,M_STRIKE)) * _sin_o_dip + (k3A(x,y,z,c,L,W,M_STRIKE) - k3Ah(x,y,z,c,L,W,M_STRIKE) + k3B(x,y,z,c,L,W,M_STRIKE) - z * k3C(x,y,z,c,L,W,M_STRIKE)) * _cos_o_dip);
        }
        if (UD != 0.0) {
            dip_result     = ( UD/(2.0*M_PI) ) * ((k2A(x,y,z,c,L,W,M_DIP) - k2Ah(x,y,z,c,L,W,M_DIP) + k2B(x,y,z,c,L,W,M_DIP) - z * k2C(x,y,z,c,L,W,M_DIP)) * _sin_o_dip + (k3A(x,y,z,c,L,W,M_DIP) - k3Ah(x,y,z,c,L,W,M_DIP) + k3B(x,y,z,c,L,W,M_DIP) - z * k3C(x,y,z,c,L,W,M_DIP)) * _cos_o_dip);
        }
        if (UT != 0.0) {
            tensile_result = ( UT/(2.0*M_PI) ) * ((k2A(x,y,z,c,L,W,M_THRUST) - k2Ah(x,y,z,c,L,W,M_THRUST) + k2B(x,y,z,c,L,W,M_THRUST) - z * k2C(x,y,z,c,L,W,M_THRUST)) * _sin_o_dip + (k3A(x,y,z,c,L,W,M_THRUST) - k3Ah(x,y,z,c,L,W,M_THRUST) + k3B(x,y,z,c,L,W,M_THRUST) - z * k3C(x,y,z,c,L,W,M_THRUST)) * _cos_o_dip);
        }
    }
    
    return strike_result + dip_result + tensile_result;
}
//
// dy components
// A
double quakelib::Okada::k1A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df1ASdy(x-L,_p-W,_q) - df1ASdy(x,_p-W,_q) - df1ASdy(x-L,_p,_q) + df1ASdy(x,_p,_q);
            break;
        case M_DIP:
            return df1ADdy(x-L,_p-W,_q) - df1ADdy(x,_p-W,_q) - df1ADdy(x-L,_p,_q) + df1ADdy(x,_p,_q);
            break;
        case M_THRUST:
            return df1ATdy(x-L,_p-W,_q) - df1ATdy(x,_p-W,_q) - df1ATdy(x-L,_p,_q) + df1ATdy(x,_p,_q);
            break;
        default:
            std::cout << "k1A motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::k2A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df2ASdy(x-L,_p-W,_q) - df2ASdy(x,_p-W,_q) - df2ASdy(x-L,_p,_q) + df2ASdy(x,_p,_q);
            break;
        case M_DIP:
            return df2ADdy(x-L,_p-W,_q) - df2ADdy(x,_p-W,_q) - df2ADdy(x-L,_p,_q) + df2ADdy(x,_p,_q);
            break;
        case M_THRUST:
            return df2ATdy(x-L,_p-W,_q) - df2ATdy(x,_p-W,_q) - df2ATdy(x-L,_p,_q) + df2ATdy(x,_p,_q);
            break;
        default:
            std::cout << "k2A motion not set!!";
            return 0;
            break;
    }
    
}
double quakelib::Okada::k3A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df3ASdy(x-L,_p-W,_q) - df3ASdy(x,_p-W,_q) - df3ASdy(x-L,_p,_q) + df3ASdy(x,_p,_q);
            break;
        case M_DIP:
            return df3ADdy(x-L,_p-W,_q) - df3ADdy(x,_p-W,_q) - df3ADdy(x-L,_p,_q) + df3ADdy(x,_p,_q);
            break;
        case M_THRUST:
            return df3ATdy(x-L,_p-W,_q) - df3ATdy(x,_p-W,_q) - df3ATdy(x-L,_p,_q) + df3ATdy(x,_p,_q);
            break;
        default:
            std::cout << "k3A motion not set!!";
            return 0;
            break;
    }
}
// Ah
double quakelib::Okada::k1Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return df1ASdy(x,_p,_q) - df1ASdy(x,_p-W,_q) - df1ASdy(x-L,_p,_q) + df1ASdy(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df1ADdy(x,_p,_q) - df1ADdy(x,_p-W,_q) - df1ADdy(x-L,_p,_q) + df1ADdy(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df1ATdy(x,_p,_q) - df1ATdy(x,_p-W,_q) - df1ATdy(x-L,_p,_q) + df1ATdy(x-L,_p-W,_q);
            break;
        default:
            std::cout << "k1Ah motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::k2Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return df2ASdy(x,_p,_q) - df2ASdy(x,_p-W,_q) - df2ASdy(x-L,_p,_q) + df2ASdy(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df2ADdy(x,_p,_q) - df2ADdy(x,_p-W,_q) - df2ADdy(x-L,_p,_q) + df2ADdy(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df2ATdy(x,_p,_q) - df2ATdy(x,_p-W,_q) - df2ATdy(x-L,_p,_q) + df2ATdy(x-L,_p-W,_q);
            break;
        default:
            std::cout << "k2Ah motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::k3Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return df3ASdy(x,_p,_q) - df3ASdy(x,_p-W,_q) - df3ASdy(x-L,_p,_q) + df3ASdy(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df3ADdy(x,_p,_q) - df3ADdy(x,_p-W,_q) - df3ADdy(x-L,_p,_q) + df3ADdy(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df3ATdy(x,_p,_q) - df3ATdy(x,_p-W,_q) - df3ATdy(x-L,_p,_q) + df3ATdy(x-L,_p-W,_q);
            break;
        default:
            std::cout << "k3Ah motion not set!!";
            return 0;
            break;
    }
}
// B
double quakelib::Okada::k1B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df1BSdy(x,_p,_q) - df1BSdy(x,_p-W,_q) - df1BSdy(x-L,_p,_q) + df1BSdy(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df1BDdy(x,_p,_q) - df1BDdy(x,_p-W,_q) - df1BDdy(x-L,_p,_q) + df1BDdy(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df1BTdy(x,_p,_q) - df1BTdy(x,_p-W,_q) - df1BTdy(x-L,_p,_q) + df1BTdy(x-L,_p-W,_q);
            break;
        default:
            std::cout << "k1B motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::k2B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df2BSdy(x,_p,_q) - df2BSdy(x,_p-W,_q) - df2BSdy(x-L,_p,_q) + df2BSdy(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df2BDdy(x,_p,_q) - df2BDdy(x,_p-W,_q) - df2BDdy(x-L,_p,_q) + df2BDdy(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df2BTdy(x,_p,_q) - df2BTdy(x,_p-W,_q) - df2BTdy(x-L,_p,_q) + df2BTdy(x-L,_p-W,_q);
            break;
        default:
            std::cout << "k2B motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::k3B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df3BSdy(x,_p,_q) - df3BSdy(x,_p-W,_q) - df3BSdy(x-L,_p,_q) + df3BSdy(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df3BDdy(x,_p,_q) - df3BDdy(x,_p-W,_q) - df3BDdy(x-L,_p,_q) + df3BDdy(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df3BTdy(x,_p,_q) - df3BTdy(x,_p-W,_q) - df3BTdy(x-L,_p,_q) + df3BTdy(x-L,_p-W,_q);
            break;
        default:
            std::cout << "k3B motion not set!!";
            return 0;
            break;
    }
}
// C
double quakelib::Okada::k1C (double x, double y, double z, double c, double L, double W, MotionType motion) {
    switch (motion) {
        case M_STRIKE:
            return df1CSdy(x,p(y,z,c),y,z,c) - df1CSdy(x,p(y,z,c)-W,y,z,c) - df1CSdy(x-L,p(y,z,c),y,z,c) + df1CSdy(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_DIP:
            return df1CDdy(x,p(y,z,c),y,z,c) - df1CDdy(x,p(y,z,c)-W,y,z,c) - df1CDdy(x-L,p(y,z,c),y,z,c) + df1CDdy(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_THRUST:
            return df1CTdy(x,p(y,z,c),y,z,c) - df1CTdy(x,p(y,z,c)-W,y,z,c) - df1CTdy(x-L,p(y,z,c),y,z,c) + df1CTdy(x-L,p(y,z,c)-W,y,z,c);
            break;
        default:
            std::cout << "k1C motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::k2C (double x, double y, double z, double c, double L, double W, MotionType motion) {
    switch (motion) {
        case M_STRIKE:
            return df2CSdy(x,p(y,z,c),y,z,c) - df2CSdy(x,p(y,z,c)-W,y,z,c) - df2CSdy(x-L,p(y,z,c),y,z,c) + df2CSdy(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_DIP:
            return df2CDdy(x,p(y,z,c),y,z,c) - df2CDdy(x,p(y,z,c)-W,y,z,c) - df2CDdy(x-L,p(y,z,c),y,z,c) + df2CDdy(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_THRUST:
            return df2CTdy(x,p(y,z,c),y,z,c) - df2CTdy(x,p(y,z,c)-W,y,z,c) - df2CTdy(x-L,p(y,z,c),y,z,c) + df2CTdy(x-L,p(y,z,c)-W,y,z,c);
            break;
        default:
            std::cout << "k2C motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::k3C (double x, double y, double z, double c, double L, double W, MotionType motion) {
    switch (motion) {
        case M_STRIKE:
            return df3CSdy(x,p(y,z,c),y,z,c) - df3CSdy(x,p(y,z,c)-W,y,z,c) - df3CSdy(x-L,p(y,z,c),y,z,c) + df3CSdy(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_DIP:
            return df3CDdy(x,p(y,z,c),y,z,c) - df3CDdy(x,p(y,z,c)-W,y,z,c) - df3CDdy(x-L,p(y,z,c),y,z,c) + df3CDdy(x-L,p(y,z,c)-W,y,z,c);
            break;
        case M_THRUST:
            return df3CTdy(x,p(y,z,c),y,z,c) - df3CTdy(x,p(y,z,c)-W,y,z,c) - df3CTdy(x-L,p(y,z,c),y,z,c) + df3CTdy(x-L,p(y,z,c)-W,y,z,c);
            break;
        default:
            std::cout << "k3C motion not set!!";
            return 0;
            break;
    }
}
//
// dy strike dfs
// A
double quakelib::Okada::df1ASdy(double xi, double eta, double _q) {
	OP_ADD(2); OP_MULT(6); OP_DIV(1);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * xi * Y11(_R,eta) * _sin_o_dip + (dtil(_q,eta)/2.0) * X11(_R,xi,eta,_q) + _alpha_div_two * xi * F(_R,xi,eta,_q);
}
double quakelib::Okada::df2ASdy(double xi, double eta, double _q) {
	OP_MULT(1);
	double _R = R(xi,eta,_q);
    return _alpha_div_two * E(_R,xi,eta,_q);
}
double quakelib::Okada::df3ASdy(double xi, double eta, double _q) {
	OP_ADD(1); OP_SUB(1); OP_MULT(5); OP_DIV(1);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * (_cos_o_dip/_R + _q * Y11(_R,eta) * _sin_o_dip) - _alpha_div_two * _q * F(_R,xi,eta,_q);
}
// B
double quakelib::Okada::df1BSdy(double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(6);
	double _R = R(xi,eta,_q);
    return -1.0 * xi * F(_R,xi,eta,_q) - dtil(_q,eta) * X11(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * (xi * Y11(_R,eta) + J4(_R,xi,eta,_q)) * _sin_o_dip;
}
double quakelib::Okada::df2BSdy(double xi, double eta, double _q) {
	OP_ADD(2); OP_SUB(1); OP_MULT(2); OP_DIV(1);
	double _R = R(xi,eta,_q);
    return -E(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * (1.0/_R + J5(_R,eta,_q)) * _sin_o_dip;
}
double quakelib::Okada::df3BSdy(double xi, double eta, double _q) {
	OP_SUB(2); OP_MULT(4);
	double _R = R(xi,eta,_q);
    return _q * F(_R,xi,eta,_q) - _one_minus_alpha_div_alpha * (_q * Y11(_R,eta) - J6(_R,xi,eta,_q)) * _sin_o_dip;
}
// C
double quakelib::Okada::df1CSdy(double xi, double eta, double y, double z, double c) {
	OP_SUB(1); OP_MULT(6);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha * xi * P(_R,xi,eta,_q) * _cos_o_dip - alpha * xi * Q(_R,xi,eta,y,z,c);
}
double quakelib::Okada::df2CSdy(double xi, double eta, double y, double z, double c) {
	OP_ADD(1); OP_SUB(5); OP_MULT(14); OP_DIV(5);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return 2.0 * _one_minus_alpha * (dtil(_q,eta)/_R3 - Y0(_R,xi,eta) * _sin_o_dip) * _sin_o_dip - (ytil(_q,eta)/_R3) * _cos_o_dip - alpha * (((ctil(_q,eta,z) + dtil(_q,eta))/_R3) * _sin_o_dip - eta/_R3 - (3.0 * ctil(_q,eta,z) * ytil(_q,eta) * _q)/_R5);
}
double quakelib::Okada::df3CSdy(double xi, double eta, double y, double z, double c) {
	OP_ADD(5); OP_SUB(2); OP_MULT(16); OP_DIV(4);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return -1.0 * _one_minus_alpha * _q/_R3 + (ytil(_q,eta)/_R3 - Y0(_R,xi,eta) * _cos_o_dip) * _sin_o_dip + alpha * (((ctil(_q,eta,z) + dtil(_q,eta))/_R3) * _cos_o_dip + (3.0 * ctil(_q,eta,z) * dtil(_q,eta) * _q)/_R5 - (Y0(_R,xi,eta) * _cos_o_dip + _q * Z0(_R,xi,eta,_q,z)) * _sin_o_dip);
}
//
// dy dip dfs
// A
double quakelib::Okada::df1ADdy(double xi, double eta, double _q) {
	OP_MULT(1);
	double _R = R(xi,eta,_q);
    return _alpha_div_two * E(_R,xi,eta,_q);
}
double quakelib::Okada::df2ADdy(double xi, double eta, double _q) {
	OP_ADD(2); OP_MULT(6); OP_DIV(1);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * dtil(_q,eta) * X11(_R,xi,eta,_q) + (xi/2.0) * Y11(_R,eta) * _sin_o_dip + _alpha_div_two * eta * G(_R,xi,eta,_q);
}
double quakelib::Okada::df3ADdy(double xi, double eta, double _q) {
	OP_MULT(4); OP_SUB(1);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * ytil(_q,eta) * X11(_R,xi,eta,_q) - _alpha_div_two * _q * G(_R,xi,eta,_q);
}
// B
double quakelib::Okada::df1BDdy(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(4);
	double _R = R(xi,eta,_q);
    return -1.0 * E(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * J1(_R,xi,eta,_q) * _sin_o_dip * _cos_o_dip;
}
double quakelib::Okada::df2BDdy(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(7); OP_SUB(1);
	double _R = R(xi,eta,_q);
    return -1.0 * eta * G(_R,xi,eta,_q) - xi * Y11(_R,eta) * _sin_o_dip + _one_minus_alpha_div_alpha * J2(_R,xi,eta,_q) * _sin_o_dip * _cos_o_dip;
}
double quakelib::Okada::df3BDdy(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(4);
	double _R = R(xi,eta,_q);
    return _q * G(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * J3(_R,xi,eta,_q) * _sin_o_dip * _cos_o_dip;
}
// C
double quakelib::Okada::df1CDdy(double xi, double eta, double y, double z, double c) {
	OP_ADD(2); OP_MULT(13); OP_SUB(2); OP_DIV(3);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return -1.0 * _one_minus_alpha * (eta/_R3) + Y0(_R,xi,eta) * (_sin_o_dip*_sin_o_dip) - alpha * (((ctil(_q,eta,z) + dtil(_q,eta))/_R3) * _sin_o_dip - (3.0 * ctil(_q,eta,z) * ytil(_q,eta) * _q)/_R5);
}
double quakelib::Okada::df2CDdy(double xi, double eta, double y, double z, double c) {
	OP_ADD(1); OP_MULT(11); OP_SUB(3);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _ytil = ytil(_q,eta);
    return _one_minus_alpha * (X11(_R,xi,eta,_q) - (_ytil*_ytil) * X32(_R,xi)) - alpha * ctil(_q,eta,z) * ((dtil(_q,eta) + 2.0 * _q * _cos_o_dip) * X32(_R,xi) - _ytil * eta * _q * X53(_R,xi));
}
double quakelib::Okada::df3CDdy(double xi, double eta, double y, double z, double c) {
	OP_ADD(3); OP_MULT(12); OP_SUB(1);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return xi * P(_R,xi,eta,_q) * _sin_o_dip + ytil(_q,eta) * dtil(_q,eta) * X32(_R,xi) + alpha * ctil(_q,eta,z) * ((ytil(_q,eta) + 2.0 * _q * _sin_o_dip) * X32(_R,xi) - ytil(_q,eta) * (_q*_q) * X53(_R,xi));
}
//
// dy tensile dfs
// A
double quakelib::Okada::df1ATdy(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha_div_two * (_cos_o_dip/_R + _q * Y11(_R,eta) * _sin_o_dip) - _alpha_div_two * _q * F(_R,xi,eta,_q);
}
double quakelib::Okada::df2ATdy(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha_div_two * ytil(_q,eta) * X11(_R,xi,eta,_q) - _alpha_div_two * _q * G(_R,xi,eta,_q);
}
double quakelib::Okada::df3ATdy(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * (dtil(_q,eta) * X11(_R,xi,eta,_q) + xi * Y11(_R,eta) * _sin_o_dip) + _alpha_div_two * _q * H(_R,xi,eta,_q);
}
// B
double quakelib::Okada::df1BTdy(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return _q * F(_R,xi,eta,_q) - _one_minus_alpha_div_alpha * J1(_R,xi,eta,_q) * (_sin_o_dip*_sin_o_dip);
}
double quakelib::Okada::df2BTdy(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return _q * G(_R,xi,eta,_q) - _one_minus_alpha_div_alpha * J2(_R,xi,eta,_q) * (_sin_o_dip*_sin_o_dip);
}
double quakelib::Okada::df3BTdy(double xi, double eta, double _q) {
	double _R = R(xi,eta,_q);
    return -1.0 * _q * H(_R,xi,eta,_q) - _one_minus_alpha_div_alpha * J3(_R,xi,eta,_q) * (_sin_o_dip*_sin_o_dip);
}
// C
double quakelib::Okada::df1CTdy(double xi, double eta, double y, double z, double c) {
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return _one_minus_alpha * (_q/_R3 + Y0(_R,xi,eta) * _sin_o_dip * _cos_o_dip) + alpha * ((z/_R3) * _cos_o_dip + (3.0 * ctil(_q,eta,z) * dtil(_q,eta) * _q)/_R5 - _q * Z0(_R,xi,eta,_q,z) * _sin_o_dip);
}
double quakelib::Okada::df2CTdy(double xi, double eta, double y, double z, double c) {
	OP_ADD(2); OP_MULT(15); OP_SUB(2);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha * 2.0 * xi * P(_R,xi,eta,_q) * _sin_o_dip - ytil(_q,eta) * dtil(_q,eta) * X32(_R,xi) + alpha * ctil(_q,eta,z) * ((ytil(_q,eta) + 2.0 * _q * _sin_o_dip) * X32(_R,xi) - ytil(_q,eta) * (_q*_q) * X53(_R,xi));
}
double quakelib::Okada::df3CTdy(double xi, double eta, double y, double z, double c) {
	OP_ADD(4); OP_MULT(16); OP_SUB(2);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _ytil = ytil(_q,eta);
    return -1.0 * _one_minus_alpha * (xi * P(_R,xi,eta,_q) * _cos_o_dip - X11(_R,xi,eta,_q) + (_ytil*_ytil) * X32(_R,xi)) + alpha * ctil(_q,eta,z) * ((dtil(_q,eta) + 2.0 * _q * _cos_o_dip) * X32(_R,xi) - _ytil * eta * _q * X53(_R,xi)) + alpha * xi * Q(_R,xi,eta,y,z,c);
}
//
// dy globals
double quakelib::Okada::E(double _R, double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(3); OP_DIV(2);
	double _R3 = _R*_R*_R;
    return _sin_o_dip/_R - (ytil(_q,eta) * _q)/_R3;
}
double quakelib::Okada::F(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(5); OP_DIV(1);
    return dtil(_q,eta)/(_R*_R*_R) + xi*xi * Y32(_R,eta) * _sin_o_dip;
}
double quakelib::Okada::G(double _R, double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(4);
    return 2.0 * X11(_R,xi,eta,_q) * _sin_o_dip - ytil(_q,eta) * _q * X32(_R,xi);
}
double quakelib::Okada::H(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(5);
    return dtil(_q,eta) * _q * X32(_R,xi) + xi * _q * Y32(_R,eta) * _sin_o_dip;
}
double quakelib::Okada::P(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(4);
	double _R3 = _R*_R*_R;
    return _cos_o_dip/_R3 + _q * Y32(_R,eta) * _sin_o_dip;
}
double quakelib::Okada::Q(double _R, double xi, double eta, double y, double z, double c) {
	OP_ADD(2); OP_MULT(8); OP_SUB(2); OP_DIV(1);
	double _q = q(y,z,c);
	double _R5 = _R*_R*_R*_R*_R;
    return (3.0 * ctil(_q,eta,z) * dtil(_q,eta))/_R5 - (z * Y32(_R,eta) + Z32(_R,_q,eta,z) + Z0(_R,xi,eta,_q,z)) * _sin_o_dip;
}

//
// dz
//
double quakelib::Okada::duxdz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
        if (US != 0.0) {
            strike_result  = ( US/(2.0*M_PI) ) * (l1A(x,y,z,c,L,W,M_STRIKE) + l1Ah(x,y,z,c,L,W,M_STRIKE) + l1B(x,y,z,c,L,W,M_STRIKE) + u1C(x,y,z,c,L,W,M_STRIKE) + z * l1C(x,y,z,c,L,W,M_STRIKE));
        }
        if (UD != 0.0) {
            dip_result     = ( UD/(2.0*M_PI) ) * (l1A(x,y,z,c,L,W,M_DIP) + l1Ah(x,y,z,c,L,W,M_DIP) + l1B(x,y,z,c,L,W,M_DIP) + u1C(x,y,z,c,L,W,M_DIP) + z * l1C(x,y,z,c,L,W,M_DIP));
        }
        if (UT != 0.0) {
            tensile_result = ( UT/(2.0*M_PI) ) * (l1A(x,y,z,c,L,W,M_THRUST) + l1Ah(x,y,z,c,L,W,M_THRUST) + l1B(x,y,z,c,L,W,M_THRUST) + u1C(x,y,z,c,L,W,M_THRUST) + z * l1C(x,y,z,c,L,W,M_THRUST));
        }
    }
    
    return strike_result + dip_result + tensile_result;
}
double quakelib::Okada::duydz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
        if (US != 0.0) {
            strike_result  = ( US/(2.0*M_PI) ) * ((l2A(x,y,z,c,L,W,M_STRIKE) + l2Ah(x,y,z,c,L,W,M_STRIKE) + l2B(x,y,z,c,L,W,M_STRIKE) + u2C(x,y,z,c,L,W,M_STRIKE) + z * l2C(x,y,z,c,L,W,M_STRIKE)) * _cos_o_dip - (l3A(x,y,z,c,L,W,M_STRIKE) + l3Ah(x,y,z,c,L,W,M_STRIKE) + l3B(x,y,z,c,L,W,M_STRIKE) + u3C(x,y,z,c,L,W,M_STRIKE) + z * l3C(x,y,z,c,L,W,M_STRIKE)) * _sin_o_dip);
        }
        if (UD != 0.0) {
            dip_result     = ( UD/(2.0*M_PI) ) * ((l2A(x,y,z,c,L,W,M_DIP) + l2Ah(x,y,z,c,L,W,M_DIP) + l2B(x,y,z,c,L,W,M_DIP) + u2C(x,y,z,c,L,W,M_DIP) + z * l2C(x,y,z,c,L,W,M_DIP)) * _cos_o_dip - (l3A(x,y,z,c,L,W,M_DIP) + l3Ah(x,y,z,c,L,W,M_DIP) + l3B(x,y,z,c,L,W,M_DIP) + u3C(x,y,z,c,L,W,M_DIP) + z * l3C(x,y,z,c,L,W,M_DIP)) * _sin_o_dip);
        }
        if (UT != 0.0) {
            tensile_result = ( UT/(2.0*M_PI) ) * ((l2A(x,y,z,c,L,W,M_THRUST) + l2Ah(x,y,z,c,L,W,M_THRUST) + l2B(x,y,z,c,L,W,M_THRUST) + u2C(x,y,z,c,L,W,M_THRUST) + z * l2C(x,y,z,c,L,W,M_THRUST)) * _cos_o_dip - (l3A(x,y,z,c,L,W,M_THRUST) + l3Ah(x,y,z,c,L,W,M_THRUST) + l3B(x,y,z,c,L,W,M_THRUST) + u3C(x,y,z,c,L,W,M_THRUST) + z * l3C(x,y,z,c,L,W,M_THRUST)) * _sin_o_dip);
        }
    }
    
    return strike_result + dip_result + tensile_result;
}
double quakelib::Okada::duzdz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT) {
    double strike_result = 0.0, dip_result = 0.0, tensile_result = 0.0;
    
    if (z <= 0.0 && !on_element_corner(x, y, z, c, L, W)) {
        if (US != 0.0) {
            strike_result  = ( US/(2.0*M_PI) ) * ((l2A(x,y,z,c,L,W,M_STRIKE) + l2Ah(x,y,z,c,L,W,M_STRIKE) + l2B(x,y,z,c,L,W,M_STRIKE) - u2C(x,y,z,c,L,W,M_STRIKE) - z * l2C(x,y,z,c,L,W,M_STRIKE)) * _sin_o_dip + (l3A(x,y,z,c,L,W,M_STRIKE) + l3Ah(x,y,z,c,L,W,M_STRIKE) + l3B(x,y,z,c,L,W,M_STRIKE) - u3C(x,y,z,c,L,W,M_STRIKE) - z * l3C(x,y,z,c,L,W,M_STRIKE)) * _cos_o_dip);
        }
        if (UD != 0.0) {
            dip_result     = ( UD/(2.0*M_PI) ) * ((l2A(x,y,z,c,L,W,M_DIP) + l2Ah(x,y,z,c,L,W,M_DIP) + l2B(x,y,z,c,L,W,M_DIP) - u2C(x,y,z,c,L,W,M_DIP) - z * l2C(x,y,z,c,L,W,M_DIP)) * _sin_o_dip + (l3A(x,y,z,c,L,W,M_DIP) + l3Ah(x,y,z,c,L,W,M_DIP) + l3B(x,y,z,c,L,W,M_DIP) - u3C(x,y,z,c,L,W,M_DIP) - z * l3C(x,y,z,c,L,W,M_DIP)) * _cos_o_dip);
        }
        if (UT != 0.0) {
            tensile_result = ( UT/(2.0*M_PI) ) * ((l2A(x,y,z,c,L,W,M_THRUST) + l2Ah(x,y,z,c,L,W,M_THRUST) + l2B(x,y,z,c,L,W,M_THRUST) - u2C(x,y,z,c,L,W,M_THRUST) - z * l2C(x,y,z,c,L,W,M_THRUST)) * _sin_o_dip + (l3A(x,y,z,c,L,W,M_THRUST) + l3Ah(x,y,z,c,L,W,M_THRUST) + l3B(x,y,z,c,L,W,M_THRUST) - u3C(x,y,z,c,L,W,M_THRUST) - z * l3C(x,y,z,c,L,W,M_THRUST)) * _cos_o_dip);
        }
    }
    
    return strike_result + dip_result + tensile_result;
}
//
// dz components
// A
double quakelib::Okada::l1A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df1ASdz(x-L,_p-W,_q) - df1ASdz(x,_p-W,_q) - df1ASdz(x-L,_p,_q) + df1ASdz(x,_p,_q);
            break;
        case M_DIP:
            return df1ADdz(x-L,_p-W,_q) - df1ADdz(x,_p-W,_q) - df1ADdz(x-L,_p,_q) + df1ADdz(x,_p,_q);
            break;
        case M_THRUST:
            return df1ATdz(x-L,_p-W,_q) - df1ATdz(x,_p-W,_q) - df1ATdz(x-L,_p,_q) + df1ATdz(x,_p,_q);
            break;
        default:
            std::cout << "l1A motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::l2A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df2ASdz(x-L,_p-W,_q) - df2ASdz(x,_p-W,_q) - df2ASdz(x-L,_p,_q) + df2ASdz(x,_p,_q);
            break;
        case M_DIP:
            return df2ADdz(x-L,_p-W,_q) - df2ADdz(x,_p-W,_q) - df2ADdz(x-L,_p,_q) + df2ADdz(x,_p,_q);
            break;
        case M_THRUST:
            return df2ATdz(x-L,_p-W,_q) - df2ATdz(x,_p-W,_q) - df2ATdz(x-L,_p,_q) + df2ATdz(x,_p,_q);
            break;
        default:
            std::cout << "l2A motion not set!!";
            return 0;
            break;
    }
    
}
double quakelib::Okada::l3A (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df3ASdz(x-L,_p-W,_q) - df3ASdz(x,_p-W,_q) - df3ASdz(x-L,_p,_q) + df3ASdz(x,_p,_q);
            break;
        case M_DIP:
            return df3ADdz(x-L,_p-W,_q) - df3ADdz(x,_p-W,_q) - df3ADdz(x-L,_p,_q) + df3ADdz(x,_p,_q);
            break;
        case M_THRUST:
            return df3ATdz(x-L,_p-W,_q) - df3ATdz(x,_p-W,_q) - df3ATdz(x-L,_p,_q) + df3ATdz(x,_p,_q);
            break;
        default:
            std::cout << "l3A motion not set!!";
            return 0;
            break;
    }
}
// Ah
double quakelib::Okada::l1Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return df1ASdz(x,_p,_q) - df1ASdz(x,_p-W,_q) - df1ASdz(x-L,_p,_q) + df1ASdz(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df1ADdz(x,_p,_q) - df1ADdz(x,_p-W,_q) - df1ADdz(x-L,_p,_q) + df1ADdz(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df1ATdz(x,_p,_q) - df1ATdz(x,_p-W,_q) - df1ATdz(x-L,_p,_q) + df1ATdz(x-L,_p-W,_q);
            break;
        default:
            std::cout << "l1Ah motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::l2Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return df2ASdz(x,_p,_q) - df2ASdz(x,_p-W,_q) - df2ASdz(x-L,_p,_q) + df2ASdz(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df2ADdz(x,_p,_q) - df2ADdz(x,_p-W,_q) - df2ADdz(x-L,_p,_q) + df2ADdz(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df2ATdz(x,_p,_q) - df2ATdz(x,_p-W,_q) - df2ATdz(x-L,_p,_q) + df2ATdz(x-L,_p-W,_q);
            break;
        default:
            std::cout << "l2Ah motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::l3Ah(double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,-z,c);
	_p = p(y,-z,c);
    switch (motion) {
        case M_STRIKE:
            return df3ASdz(x,_p,_q) - df3ASdz(x,_p-W,_q) - df3ASdz(x-L,_p,_q) + df3ASdz(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df3ADdz(x,_p,_q) - df3ADdz(x,_p-W,_q) - df3ADdz(x-L,_p,_q) + df3ADdz(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df3ATdz(x,_p,_q) - df3ATdz(x,_p-W,_q) - df3ATdz(x-L,_p,_q) + df3ATdz(x-L,_p-W,_q);
            break;
        default:
            std::cout << "l3Ah motion not set!!";
            return 0;
            break;
    }
}
// B
double quakelib::Okada::l1B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df1BSdz(x,_p,_q) - df1BSdz(x,_p-W,_q) - df1BSdz(x-L,_p,_q) + df1BSdz(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df1BDdz(x,_p,_q) - df1BDdz(x,_p-W,_q) - df1BDdz(x-L,_p,_q) + df1BDdz(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df1BTdz(x,_p,_q) - df1BTdz(x,_p-W,_q) - df1BTdz(x-L,_p,_q) + df1BTdz(x-L,_p-W,_q);
            break;
        default:
            std::cout << "l1B motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::l2B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df2BSdz(x,_p,_q) - df2BSdz(x,_p-W,_q) - df2BSdz(x-L,_p,_q) + df2BSdz(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df2BDdz(x,_p,_q) - df2BDdz(x,_p-W,_q) - df2BDdz(x-L,_p,_q) + df2BDdz(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df2BTdz(x,_p,_q) - df2BTdz(x,_p-W,_q) - df2BTdz(x-L,_p,_q) + df2BTdz(x-L,_p-W,_q);
            break;
        default:
            std::cout << "l2B motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::l3B (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p, _q;
	_q = q(y,z,c);
	_p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df3BSdz(x,_p,_q) - df3BSdz(x,_p-W,_q) - df3BSdz(x-L,_p,_q) + df3BSdz(x-L,_p-W,_q);
            break;
        case M_DIP:
            return df3BDdz(x,_p,_q) - df3BDdz(x,_p-W,_q) - df3BDdz(x-L,_p,_q) + df3BDdz(x-L,_p-W,_q);
            break;
        case M_THRUST:
            return df3BTdz(x,_p,_q) - df3BTdz(x,_p-W,_q) - df3BTdz(x-L,_p,_q) + df3BTdz(x-L,_p-W,_q);
            break;
        default:
            std::cout << "l3B motion not set!!";
            return 0;
            break;
    }
}
// C
double quakelib::Okada::l1C (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df1CSdz(x,_p,y,z,c) - df1CSdz(x,_p-W,y,z,c) - df1CSdz(x-L,_p,y,z,c) + df1CSdz(x-L,_p-W,y,z,c);
            break;
        case M_DIP:
            return df1CDdz(x,_p,y,z,c) - df1CDdz(x,_p-W,y,z,c) - df1CDdz(x-L,_p,y,z,c) + df1CDdz(x-L,_p-W,y,z,c);
            break;
        case M_THRUST:
            return df1CTdz(x,_p,y,z,c) - df1CTdz(x,_p-W,y,z,c) - df1CTdz(x-L,_p,y,z,c) + df1CTdz(x-L,_p-W,y,z,c);
            break;
        default:
            std::cout << "l1C motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::l2C (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df2CSdz(x,_p,y,z,c) - df2CSdz(x,_p-W,y,z,c) - df2CSdz(x-L,_p,y,z,c) + df2CSdz(x-L,_p-W,y,z,c);
            break;
        case M_DIP:
            return df2CDdz(x,_p,y,z,c) - df2CDdz(x,_p-W,y,z,c) - df2CDdz(x-L,_p,y,z,c) + df2CDdz(x-L,_p-W,y,z,c);
            break;
        case M_THRUST:
            return df2CTdz(x,_p,y,z,c) - df2CTdz(x,_p-W,y,z,c) - df2CTdz(x-L,_p,y,z,c) + df2CTdz(x-L,_p-W,y,z,c);
            break;
        default:
            std::cout << "k2C motion not set!!";
            return 0;
            break;
    }
}
double quakelib::Okada::l3C (double x, double y, double z, double c, double L, double W, MotionType motion) {
	double _p = p(y,z,c);
    switch (motion) {
        case M_STRIKE:
            return df3CSdz(x,_p,y,z,c) - df3CSdz(x,_p-W,y,z,c) - df3CSdz(x-L,_p,y,z,c) + df3CSdz(x-L,_p-W,y,z,c);
            break;
        case M_DIP:
            return df3CDdz(x,_p,y,z,c) - df3CDdz(x,_p-W,y,z,c) - df3CDdz(x-L,_p,y,z,c) + df3CDdz(x-L,_p-W,y,z,c);
            break;
        case M_THRUST:
            return df3CTdz(x,_p,y,z,c) - df3CTdz(x,_p-W,y,z,c) - df3CTdz(x-L,_p,y,z,c) + df3CTdz(x-L,_p-W,y,z,c);
            break;
        default:
            std::cout << "l3C motion not set!!";
            return 0;
            break;
    }
}
//
// dz strike dfs
// A
double quakelib::Okada::df1ASdz(double xi, double eta, double _q) {
	OP_ADD(2); OP_MULT(6); OP_DIV(1);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * xi * Y11(_R,eta) * _cos_o_dip + (ytil(_q,eta)/2.0) * X11(_R,xi,eta,_q) + _alpha_div_two * xi * Fp(_R,xi,eta,_q);
}
double quakelib::Okada::df2ASdz(double xi, double eta, double _q) {
	OP_MULT(1);
	double _R = R(xi,eta,_q);
    return _alpha_div_two * Ep(_R,xi,eta,_q);
}
double quakelib::Okada::df3ASdz(double xi, double eta, double _q) {
	OP_SUB(2); OP_MULT(6); OP_DIV(1);
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha_div_two * (_sin_o_dip/_R - _q * Y11(_R,eta) * _cos_o_dip) - _alpha_div_two * _q * Fp(_R,xi,eta,_q);
}
// B
double quakelib::Okada::df1BSdz(double xi, double eta, double _q) {
	OP_ADD(1); OP_SUB(1); OP_MULT(5);
	double _R = R(xi,eta,_q);
    return -1.0 * xi * Fp(_R,xi,eta,_q) - ytil(_q,eta) * X11(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * K1(_R,xi,eta,_q) * _sin_o_dip;
}
double quakelib::Okada::df2BSdz(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(4);
	double _R = R(xi,eta,_q);
	double _dtil = dtil(_q,eta);
    return -1.0 * Ep(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * ytil(_q,eta) * D11(_R,_dtil) * _sin_o_dip;
}
double quakelib::Okada::df3BSdz(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(3);
	double _R = R(xi,eta,_q);
    return _q * Fp(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * K2(_R,xi,eta,_q) * _sin_o_dip;
}
// C
double quakelib::Okada::df1CSdz(double xi, double eta, double y, double z, double c) {
	OP_MULT(5); OP_SUB(1);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha * xi * Pp(_R,xi,eta,_q) * _cos_o_dip - alpha * xi * Qp(_R,xi,eta,_q,z);
}
double quakelib::Okada::df2CSdz(double xi, double eta, double y, double z, double c) {
	OP_ADD(3); OP_SUB(2); OP_MULT(14); OP_DIV(4);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return 2.0 * _one_minus_alpha * (ytil(_q,eta)/_R3 - Y0(_R,xi,eta) * _cos_o_dip) * _sin_o_dip + (dtil(_q,eta)/_R3) * _cos_o_dip - alpha * (((ctil(_q,eta,z) + dtil(_q,eta))/_R3) * _cos_o_dip + (3.0 * ctil(_q,eta,z) * dtil(_q,eta) * _q)/_R5);
}
double quakelib::Okada::df3CSdz(double xi, double eta, double y, double z, double c) {
	OP_ADD(2); OP_SUB(4); OP_MULT(15); OP_DIV(3);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return (ytil(_q,eta)/_R3 - Y0(_R,xi,eta) * _cos_o_dip) * _cos_o_dip - alpha * (((ctil(_q,eta,z) + dtil(_q,eta))/_R3) * _sin_o_dip - (3.0 * ctil(_q,eta,z) * ytil(_q,eta) * _q)/_R5 - Y0(_R,xi,eta) * (_sin_o_dip*_sin_o_dip) + _q * Z0(_R,xi,eta,_q,z) * _cos_o_dip);
}
//
// dz dip dfs
// A
double quakelib::Okada::df1ADdz(double xi, double eta, double _q) {
	OP_MULT(1);
	double _R = R(xi,eta,_q);
    return _alpha_div_two * Ep(_R,xi,eta,_q);
}
double quakelib::Okada::df2ADdz(double xi, double eta, double _q) {
	OP_ADD(2); OP_MULT(6); OP_DIV(1);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * ytil(_q,eta) * X11(_R,xi,eta,_q) + (xi/2.0) * Y11(_R,eta) * _cos_o_dip + _alpha_div_two * eta * Gp(_R,xi,eta,_q);
}
double quakelib::Okada::df3ADdz(double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(5);
	double _R = R(xi,eta,_q);
    return -1.0 * _one_minus_alpha_div_two * dtil(_q,eta) * X11(_R,xi,eta,_q) - _alpha_div_two * _q * Gp(_R,xi,eta,_q);
}
// B
double quakelib::Okada::df1BDdz(double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(4);
	double _R = R(xi,eta,_q);
    return -1.0 * Ep(_R,xi,eta,_q) - _one_minus_alpha_div_alpha * K3(_R,xi,eta,_q) * _sin_o_dip * _cos_o_dip;
}
double quakelib::Okada::df2BDdz(double xi, double eta, double _q) {
	OP_SUB(2); OP_MULT(8);
	double _R = R(xi,eta,_q);
	double _dtil = dtil(_q,eta);
    return -1.0 * eta * Gp(_R,xi,eta,_q) - xi * Y11(_R,eta) * _cos_o_dip - _one_minus_alpha_div_alpha * xi * D11(_R,_dtil) * _sin_o_dip * _cos_o_dip;
}
double quakelib::Okada::df3BDdz(double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(4);
	double _R = R(xi,eta,_q);
    return _q * Gp(_R,xi,eta,_q) - _one_minus_alpha_div_alpha * K4(_R,xi,eta,_q) * _sin_o_dip * _cos_o_dip;
}
// C
double quakelib::Okada::df1CDdz(double xi, double eta, double y, double z, double c) {
	OP_ADD(3); OP_SUB(1); OP_MULT(12); OP_DIV(3);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return -1.0 * (_q/_R3) + Y0(_R,xi,eta) * _sin_o_dip * _cos_o_dip - alpha * (((ctil(_q,eta,z) + dtil(_q,eta))/_R3) * _cos_o_dip + (3.0 * ctil(_q,eta,z) * dtil(_q,eta) * _q)/_R5);
}
double quakelib::Okada::df2CDdz(double xi, double eta, double y, double z, double c) {
	OP_ADD(1); OP_SUB(2); OP_MULT(11);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha * ytil(_q,eta) * dtil(_q,eta) * X32(_R,xi) - alpha * ctil(_q,eta,z) * ((ytil(_q,eta) - 2.0 * _q * _sin_o_dip) * X32(_R,xi) + dtil(_q,eta) * eta * _q * X53(_R,xi));
}
double quakelib::Okada::df3CDdz(double xi, double eta, double y, double z, double c) {
	OP_ADD(1); OP_SUB(4); OP_MULT(13);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _dtil = dtil(_q,eta);
	double _dtil2 = _dtil*_dtil;
    return -1.0 * xi * Pp(_R,xi,eta,_q) * _sin_o_dip + X11(_R,xi,eta,_q) - _dtil2 * X32(_R,xi) - alpha * ctil(_q,eta,z) * ((_dtil - 2.0 * _q * _cos_o_dip) * X32(_R,xi) - _dtil * (_q*_q) * X53(_R,xi));
}
//
// dz tensile dfs
// A
double quakelib::Okada::df1ATdz(double xi, double eta, double _q) {
	OP_MULT(5); OP_SUB(2); OP_DIV(1);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * (_sin_o_dip/_R - _q * Y11(_R,eta) * _cos_o_dip) - _alpha_div_two * _q * Fp(_R,xi,eta,_q);
}
double quakelib::Okada::df2ATdz(double xi, double eta, double _q) {
	OP_MULT(4); OP_SUB(1);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * dtil(_q,eta) * X11(_R,xi,eta,_q) - _alpha_div_two * _q * Gp(_R,xi,eta,_q);
}
double quakelib::Okada::df3ATdz(double xi, double eta, double _q) {
	OP_ADD(2); OP_MULT(6);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha_div_two * (ytil(_q,eta) * X11(_R,xi,eta,_q) + xi * Y11(_R,eta) * _cos_o_dip) + _alpha_div_two * _q * Hp(_R,xi,eta,_q);
}
// B
double quakelib::Okada::df1BTdz(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(4);
	double _R = R(xi,eta,_q);
    return _q * Fp(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * K3(_R,xi,eta,_q) * (_sin_o_dip*_sin_o_dip);
}
double quakelib::Okada::df2BTdz(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(5);
	double _R = R(xi,eta,_q);
	double _dtil = dtil(_q,eta);
    return _q * Gp(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * xi * D11(_R,_dtil) * (_sin_o_dip*_sin_o_dip);
}
double quakelib::Okada::df3BTdz(double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(5);
	double _R = R(xi,eta,_q);
    return -1.0 * _q * Hp(_R,xi,eta,_q) + _one_minus_alpha_div_alpha * K4(_R,xi,eta,_q) * (_sin_o_dip*_sin_o_dip);
}
// C
double quakelib::Okada::df1CTdz(double xi, double eta, double y, double z, double c) {
	OP_ADD(2); OP_MULT(16); OP_SUB(3); OP_DIV(3);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _R3 = _R*_R*_R;
	double _R5 = _R3*_R*_R;
    return -1.0 * (eta/_R3) + Y0(_R,xi,eta) * (_cos_o_dip*_cos_o_dip) - alpha * ((z/_R3) * _sin_o_dip - (3.0 * ctil(_q,eta,z) * ytil(_q,eta) * _q)/_R5 - Y0(_R,xi,eta) * (_sin_o_dip*_sin_o_dip) + _q * Z0(_R,xi,eta,_q,z) * _cos_o_dip);
}
double quakelib::Okada::df2CTdz(double xi, double eta, double y, double z, double c) {
	OP_ADD(4); OP_SUB(1); OP_MULT(15);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
	double _dtil = dtil(_q,eta);
	double _dtil2 = _dtil*_dtil;
    return _one_minus_alpha * 2.0 * xi * Pp(_R,xi,eta,_q) * _sin_o_dip - X11(_R,xi,eta,_q) + _dtil2 * X32(_R,xi) - alpha * ctil(_q,eta,z) * ((_dtil - 2.0 * _q * _cos_o_dip) * X32(_R,xi) - _dtil * (_q*_q) * X53(_R,xi));
}
double quakelib::Okada::df3CTdz(double xi, double eta, double y, double z, double c) {
	OP_ADD(4); OP_SUB(1); OP_MULT(15);
	double _q = q(y,z,c);
	double _R = R(xi,eta,_q);
    return _one_minus_alpha * (xi * Pp(_R,xi,eta,_q) * _cos_o_dip + ytil(_q,eta) * dtil(_q,eta) * X32(_R,xi)) + alpha * ctil(_q,eta,z) * ((ytil(_q,eta) - 2.0 * _q * _sin_o_dip) * X32(_R,xi) + dtil(_q,eta) * eta * _q * X53(_R,xi)) + alpha * xi * Qp(_R,xi,eta,_q,z);
}
//
// dz globals
double quakelib::Okada::Ep(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(3); OP_DIV(2);
	double _R3 = _R*_R*_R;
    return _cos_o_dip/_R + (dtil(_q,eta) * _q)/_R3;
}
double quakelib::Okada::Fp(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(6); OP_DIV(1);
	double _R3 = _R*_R*_R;
    return ytil(_q,eta)/_R3 + (xi*xi) * Y32(_R,eta) * _cos_o_dip;
}
double quakelib::Okada::Gp(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(4);
    return 2.0 * X11(_R,xi,eta,_q) * _cos_o_dip + dtil(_q,eta) * _q * X32(_R,xi);
}
double quakelib::Okada::Hp(double _R, double xi, double eta, double _q) {
	OP_ADD(1); OP_MULT(5);
    return ytil(_q,eta) * _q * X32(_R,xi) + xi * _q * Y32(_R,eta) * _cos_o_dip;
}
double quakelib::Okada::Pp(double _R, double xi, double eta, double _q) {
	OP_SUB(1); OP_MULT(4); OP_DIV(1);
	double _R3 = _R*_R*_R;
    return _sin_o_dip/_R3 - _q * Y32(_R,eta) * _cos_o_dip;
}
double quakelib::Okada::Qp(double _R, double xi, double eta, double _q, double z) {
	OP_ADD(3); OP_MULT(9); OP_SUB(1); OP_DIV(1);
	double _R5 = _R*_R*_R*_R*_R;
    return (3.0 * ctil(_q,eta,z) * ytil(_q,eta))/_R5 + _q * Y32(_R,eta) - (z * Y32(_R,eta) + Z32(_R,_q,eta,z) + Z0(_R,xi,eta,_q,z)) * _cos_o_dip;
}

//------------------- Not 100% verified   ---Kasey---
// gravity change on the free surface (z=0)
double quakelib::Okada::dg(double x, double y, double c, double dip, double L, double W, double US, double UD, double UT, double lambda, double mu) {
    OP_CMP(3); OP_MULT(3); OP_ADD(3); OP_SUB(1);
    // Evaluating delta_g at free surface so z=0

    if (mu <= 0) throw std::invalid_argument("Mu must be greater than zero.");

	precalc(dip, lambda, mu);

    // Everything is in M-K-S units
    double G   = 0.000000000066738; //Big G gravitation constant
    double RHO = 2670.0;   //mean crustal density (rough estimate)
    double B   = 0.00000309; //free-air gravity gradient (taken from Okubo '92)

	double _p = p(y,0.0,c);
	double _q = q(y,0.0,c);
    double dgS= 0.0; //contribution from Strike
    double dgD= 0.0; //Dip
    double dgT= 0.0; //Tensile
    double dgC= 0.0; //Cavitation, currently set such that no matter fills the cavity.
                     //If adding this feature, multiply dgC instead by DENSITY_DIFF*BIG_G
                     //in the return call. DENSITY_DIFF is mean crustal density minus the
                     //density of the cavity filling material.

    if (US != 0.0) {
        OP_MULT(1);
        dgS = US*dSg(x,_p,_q,L,W);
    }
    if (UD != 0.0) {
        OP_MULT(1);
        dgD = UD*dDg(x,_p,_q,L,W);
    }
    if (UT != 0.0) {
        OP_MULT(2);
        dgT = UT*dTg(x,_p,_q,L,W);
        dgC = UT*dCg(x,_p,_q,L,W);
    }

    //Free-air effect of surface elevation change (not sure what this means, taken from Okubo '92)
    // here uz used to get the height change on surface of halfspace at (x,y)
    Vec<3> location;
    Vec<3> displace;

    location[0] = x;
    location[1] = y;
    location[2] = 0.0;

    displace = calc_displacement_vector(location,c,dip,L,W,US,UD,UT,lambda,mu);

    double dgFree = B*displace[2];

    return RHO*G*(dgS + dgD + dgT + dgC) - dgFree;
}
//
// double bar evaluation (chinnery)
//
double quakelib::Okada::dSg(double x, double _p, double _q, double L, double W) {
    OP_ADD(1); OP_SUB(6);
    return Sg(x,_p,_q) - Sg(x,_p-W,_q) - Sg(x-L,_p,_q) + Sg(x-L,_p-W,_q);
}
double quakelib::Okada::dDg(double x, double _p, double _q, double L, double W) {
    OP_ADD(1); OP_SUB(6);
    return Dg(x,_p,_q) - Dg(x,_p-W,_q) - Dg(x-L,_p,_q) + Dg(x-L,_p-W,_q);
}
double quakelib::Okada::dTg(double x, double _p, double _q, double L, double W) {
    OP_ADD(1); OP_SUB(6);
    return Tg(x,_p,_q) - Tg(x,_p-W,_q) - Tg(x-L,_p,_q) + Tg(x-L,_p-W,_q);
}
double quakelib::Okada::dCg(double x, double _p, double _q, double L, double W) {
    OP_ADD(1); OP_SUB(6);
    return Cg(x,_p,_q) - Cg(x,_p-W,_q) - Cg(x-L,_p,_q) + Cg(x-L,_p-W,_q);
}
//
// gravity globals
//
double quakelib::Okada::Sg(double xi, double eta, double _q){
	double _R = R(xi,eta,_q);
	double _Rpeta = _R+eta;
    if (!singularity4(_Rpeta)){
        OP_ADD(1); OP_MULT(4); OP_DIV(2); OP_SUB(1);
        double _RRpeta = _R*_Rpeta;
        double _q2 = _q*_q;
        return _q2*_cos_o_dip/_RRpeta - _q*_sin_o_dip/_R;
    } else {
        OP_ADD(1); OP_MULT(2); OP_DIV(1);
        return -1.0*_q*_sin_o_dip/_R;
    }
}
double quakelib::Okada::Dg(double xi, double eta, double _q){
	double _R = R(xi,eta,_q);
	double _Rpxi = _R+xi;
    if (!singularity3(_Rpxi)) {
        OP_ADD(1); OP_MULT(4); OP_DIV(1); OP_SUB(1);
        double _RRpxi = _R*_Rpxi;
        double _dtil = dtil(_q,eta);
        return 2.0*I2g(_R,xi,eta,_q)*_sin_o_dip - _q*_dtil/_RRpxi ;
    } else {
        OP_ADD(1); OP_MULT(2);
        return 2.0*I2g(_R,xi,eta,_q)*_sin_o_dip;
    }
}
double quakelib::Okada::Tg(double xi, double eta, double _q){
    OP_ADD(2); OP_MULT(2);
	double _R = R(xi,eta,_q);
    double _ytil = ytil(_q,eta);
	double _Rpeta = _R+eta;
	double _Rpxi = _R+xi;
    double ret_value = 2.0*I2g(_R,xi,eta,_q)*_cos_o_dip;
    if (!singularity3(_Rpxi)) {
        OP_MULT(2); OP_DIV(1); OP_ADD(1);
        double _RRpxi = _R*_Rpxi;
        ret_value += _q*_ytil/_RRpxi;
    }
    if (!singularity4(_Rpeta)) {
        OP_MULT(3); OP_ADD(1); OP_DIV(1);
        double _RRpeta = _R*_Rpeta;
        ret_value += _q*xi*_cos_o_dip/_RRpeta;
    }
    return ret_value;
}
//Cg is the contribution from matter infall into newly created cavity
//Currently this is set such that no matter fills the opening
double quakelib::Okada::Cg(double xi, double eta, double _q){
    OP_ADD(1); OP_MULT(2);
	double _R = R(xi,eta,_q);
    double _Rpxi = _R+xi;
    double ret_value = 2.0*I2g(_R,xi,eta,_q)*_cos_o_dip;
    if (!singularity3(_Rpxi)) {
        OP_MULT(1); OP_SUB(1); OP_LOG(1);
        ret_value -= _sin_o_dip*log(_Rpxi);
    } else {
        OP_ADD(1); OP_MULT(1); OP_LOG(1); OP_SUB(1);
        ret_value += _sin_o_dip*log(_R-xi);
    }
    return ret_value;
}
double quakelib::Okada::I2g(double _R, double xi, double eta, double _q){
    if (!singularity1(_q)) {
        OP_ADD(2); OP_DIV(1);
        return atan((_R+xi+eta)/_q);
    } else {
        return 0.0;
    }
}

//Below is testing if Okubo's implementation of delta_h matches the current Okada displacement z
double quakelib::Okada::dg2(double x, double y, double c, double dip, double L, double W, double US, double UD, double UT, double lambda, double mu){

	if (mu <= 0) throw std::invalid_argument("Mu must be greater than zero.");

	precalc(dip, lambda, mu);

    // Everything is in M-K-S units
    double G   = 0.000000000066738; //Big G gravitation constant
    double RHO = 2670.0;   //mean crustal density (rough estimate)
    double B   = 0.00000309; //free-air gravity gradient (taken from Okubo '92)

	double _p = p(y,0.0,c);
	double _q = q(y,0.0,c);
    double dgS= 0.0; //contribution from Strike
    double dgD= 0.0; //Dip
    double dgT= 0.0; //Tensile
    double dgC= 0.0; //Cavitation, currently set such that no matter fills the cavity.
                     //If adding this feature, multiply dgC instead by DENSITY_DIFF*BIG_G
                     //in the return call. DENSITY_DIFF is mean crustal density minus the
                     //density of the cavity filling material.

    if (US != 0.0) {
        OP_MULT(1);
        dgS = US*dSg(x,_p,_q,L,W);
    }
    if (UD != 0.0) {
        OP_MULT(1);
        dgD = UD*dDg(x,_p,_q,L,W);
    }
    if (UT != 0.0) {
        OP_MULT(2);
        dgT = UT*dTg(x,_p,_q,L,W);
        dgC = UT*dCg(x,_p,_q,L,W);
    }

    //Free-air effect of surface elevation change (not sure what this means, taken from Okubo '92)
    // here uz used to get the height change on surface of halfspace at (x,y)

    double dgFree = B*dH(x,y,c,dip,L,W,US,UD,UT,lambda,mu);

    return RHO*G*(dgS + dgD + dgT + dgC) - dgFree;
}
////-------------------------
double quakelib::Okada::dH(double x, double y, double c, double dip, double L, double W, double US, double UD, double UT, double lambda, double mu) {
	double dHs = 0.0;
	double dHd = 0.0;
	double dHt = 0.0;
	double _p  = p(y,0.0,c);
	double _q  = q(y,0.0,c);

	precalc(dip, lambda, mu);

	if (US != 0.0) {
		dHs = US*dSh(x,_p,_q,L,W)/(2.0*M_PI);
	}
	if (UD != 0.0) {
		dHd = UD*dDh(x,_p,_q,L,W)/(2.0*M_PI);
	}
	if (UT != 0.0) {
		dHt = UT*dTh(x,_p,_q,L,W)/(2.0*M_PI);
	}

	return dHs + dHd + dHt;
}
////-------------------------chinnery
double quakelib::Okada::dSh(double x, double _p, double _q, double L, double W){
    OP_ADD(1); OP_SUB(6);
    return Sh(x,_p,_q) - Sh(x,_p-W,_q) - Sh(x-L,_p,_q) + Sh(x-L,_p-W,_q);
}
double quakelib::Okada::dDh(double x, double _p, double _q, double L, double W){
    OP_ADD(1); OP_SUB(6);
    return Dh(x,_p,_q) - Dh(x,_p-W,_q) - Dh(x-L,_p,_q) + Dh(x-L,_p-W,_q);
}
double quakelib::Okada::dTh(double x, double _p, double _q, double L, double W){
    OP_ADD(1); OP_SUB(6);
    return Th(x,_p,_q) - Th(x,_p-W,_q) - Th(x-L,_p,_q) + Th(x-L,_p-W,_q);
}
////-------------------------
double quakelib::Okada::Sh(double xi, double eta, double _q){
	double _R 			= R(xi,eta,_q);
	double _Rpeta 		= _R+eta;
	double _dtil  		= dtil(_q,eta);
	double non_singular = -1.0*I4h(_R,xi,eta,_q)*_sin_o_dip;
	if (!singularity4(_Rpeta)){
		double _RRpeta  = _R*_Rpeta;
		return non_singular - _dtil*_q/_RRpeta - _q*_sin_o_dip/_Rpeta;
	} else {
		return non_singular;
	}
}
double quakelib::Okada::Dh(double xi, double eta, double _q){
	double _R 			= R(xi,eta,_q);
	double _qR			= _q*_R;
    double _Rpxi 		= _R+xi;
	double _dtil  		= dtil(_q,eta);
	double non_singular = I5h(_R,xi,eta,_q)*_sin_o_dip*_cos_o_dip -_sin_o_dip*atan(xi*eta/_qR);
	if (!singularity3(_Rpxi)){
		double _RRpxi   = _R*_Rpxi;
		return non_singular - _dtil*_q/_RRpxi;
	} else{
		return non_singular;
	}
}
double quakelib::Okada::Th(double xi, double eta, double _q){
	double _R 			= R(xi,eta,_q);
    double _ytil		= ytil(_q,eta);
    double _Rpxi 		= _R+xi;
	double _Rpeta 		= _R+eta;
	double _qR			= _q*_R;
	double non_singular = -1.0*_sin_o_2_dip*I5h(_R,xi,eta,_q) - _cos_o_dip*atan(xi*eta/_qR);
	if (!singularity4(_Rpeta)){
		double _RRpeta  = _R*_Rpeta;
		non_singular 	+= _cos_o_dip*xi*_q/_RRpeta;
	}
	if (!singularity3(_Rpxi)){
		double _RRpxi   = _R*_Rpxi;
		non_singular    += _ytil*_q/_RRpxi;
	}
	return non_singular;
}
double quakelib::Okada::I4h(double _R, double xi, double eta, double _q){
	double _Rpeta = _R+eta;
	double _dtil  = dtil(_q,eta);
	if (!cos_zero()){
	    if (!singularity4(_Rpeta)){
	    	return _one_minus_two_nu*(log(_R+_dtil)-_sin_o_dip*log(_Rpeta))/_cos_o_dip;
	    } else {
	    	return _one_minus_two_nu*(log(_R+_dtil)+ _sin_o_dip*log(_R-eta))/_cos_o_dip;
	    }
	} else {
		return -1.0*_one_minus_two_nu*_q/(_R+_dtil);
	}
}
double quakelib::Okada::I5h(double _R, double xi, double eta, double _q){
	double _dtil  = dtil(_q,eta);
	if (!cos_zero()){
    	return 2.0*_one_minus_two_nu*I1h(_R,xi,eta,_q)/_cos_o_dip;
	} else {
		return -1.0*_one_minus_two_nu*xi*_sin_o_dip/(_R+_dtil);
	}
}
double quakelib::Okada::I1h(double _R, double xi, double eta, double _q){
	double _Rpeta = _R+eta;
	double xi_cos = xi*_cos_o_dip;
	if (!singularity2(xi)){
		return atan((_Rpeta*(1.0+_sin_o_dip)-_q*_cos_o_dip)/xi_cos);
	} else {
		return 0.0;
	}
}




















