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

#include "QuakeLibUtil.h"

#ifndef _QUAKELIB_OKADA_H_
#define _QUAKELIB_OKADA_H_

//#define COUNT_FLOPS
#define TOLERANCE 0.0001
#define TRIG_TOLERANCE 0.0001

//KWS added rough constants to be refined later
#define EARTH_CRUST_DENSITY     2800.0              // in kg/m^3
#define BIG_G                   0.0000000000667     // in m^3/kg/s^2
#define GRAV_BETA               0.00000309          // in m/s^2 free air gravity gradient

namespace quakelib {
	enum MotionType {
		M_UNDEFINED,
		M_STRIKE,
		M_DIP,
		M_THRUST,
	};
	
	//! Tracks operation counts in complex calculations.
	class OpCount {
	private:
		int _add, _sub, _mult, _div, _sqrt, _cmp, _and, _abs, _log;
	public:
		OpCount(void) : _add(0), _sub(0), _mult(0), _div(0), _sqrt(0), _cmp(0), _and(0), _abs(0), _log(0) {};
		//! Reset operation counts to 0.
		void reset(void) { _add = _sub = _mult = _div = _sqrt = _cmp = _and = _abs = _log = 0; };
		//! Increment the count of add operations.
		void add(int count) { _add += count; };
		//! Increment the count of subtraction operations.
		void sub(int count) { _sub += count; };
		//! Increment the count of multiplication operations.
		void mult(int count) { _mult += count; };
		//! Increment the count of division operations.
		void div(int count) { _div += count; };
		//! Increment the count of square root operations.
		void sqrt(int count) { _sqrt += count; };
		//! Increment the count of comparison operations.
		void cmp(int count) { _cmp += count; };
		//! Increment the count of logical AND operations.
		void land(int count) { _and += count; };
		//! Increment the count of abs() operations.
		void abs(int count) { _abs += count; };
		//! Increment the count of logarithm operations.
		void log(int count) { _log += count; };
		//! Null operation.
		void null(void) const {};
		
		friend std::ostream& operator<<(std::ostream& os, const OpCount& opc);
	};
	
	//! Output the operation counts for this object.
	std::ostream& operator<<(std::ostream& os, const OpCount& opc);
	
	extern OpCount op;
	
#ifdef COUNT_FLOPS
	// Define operation counter functions
#define OP_ADD(x) (op.add(x))
#define OP_SUB(x) (op.sub(x))
#define OP_MULT(x) (op.mult(x))
#define OP_DIV(x) (op.div(x))
#define OP_SQRT(x) (op.sqrt(x))
#define OP_CMP(x) (op.cmp(x))
#define OP_AND(x) (op.land(x))
#define OP_ABS(x) (op.abs(x))
#define OP_LOG(x) (op.log(x))
#else
	// If COUNT_FLOPS is off, these functions are empty
#define OP_ADD(x) (op.null())
#define OP_SUB(x) (op.null())
#define OP_MULT(x) (op.null())
#define OP_DIV(x) (op.null())
#define OP_SQRT(x) (op.null())
#define OP_CMP(x) (op.null())
#define OP_AND(x) (op.null())
#define OP_ABS(x) (op.null())
#define OP_LOG(x) (op.null())
#endif
	
	//! Calculates Okada stress and displacement functions for fault.
	class Okada {
	public:
		// [sxx,syy,szz,sxy,sxz,syz]
		//! Calculate stress tensor at a given location from a rectanglular fault
		//! moving with parameters as specified in Okada's paper.
		Tensor<3,3> calc_stress_tensor(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument);
		// [ux,uy,uz]
		//! Calculate displacement vector at a given location from a rectanglular fault
		//! with parameters as specified in Okada's paper.
		Vec<3> calc_displacement_vector(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument);
		//! Calculate dus at a given location from a rectanglular fault
		//! moving with parameters as specified in Okada's paper.
		// [duxdx,duydx,duzdx]
		Vec<3> calc_dudx(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument);
		// [duxdy,duydy,duzdy]
		Vec<3> calc_dudy(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument);
		// [duxdz,duydz,duzdz]
		Vec<3> calc_dudz(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument);
		
	public:
		//
		// Precalculated values to improve performance
		double _cos_o_dip, _sin_o_dip;
		double _cos_o_2_dip, _sin_o_2_dip;
		double _one_minus_alpha, _one_minus_alpha_div_two, _one_minus_alpha_div_alpha, _alpha_div_two;
		
		void precalc(double dip, double lambda, double mu);
		
		// methods and variables that apply to all calculations
		double cos_o(double dip);
		double sin_o(double dip);
		double alpha;
		double p(double y, double z, double c);
		double q(double y, double z, double c);
		double d(double c, double z);
		double R(double xi, double eta, double _q);
		double ytil(double _q, double eta);
		double dtil(double _q, double eta);
		double ctil(double _q, double eta, double z);
		double h(double _q, double z) const;
		
		bool on_element_corner(double x, double y, double z, double c, double L, double W) const;
		bool singularity1(double _q) const {
			OP_CMP(1); OP_ABS(1);
			return (fabs(_q) < TOLERANCE);
		};
		bool singularity2(double xi) const {
			OP_CMP(1); OP_ABS(1);
			return (fabs(xi) < TOLERANCE);
		};
		bool singularity3(double _R_plus_xi) const {
			OP_CMP(1); OP_ABS(1);
			return (fabs(_R_plus_xi) < TOLERANCE);
		};
		bool singularity4(double _R_plus_eta) const {
			OP_CMP(1); OP_ABS(1);
			return (fabs(_R_plus_eta) < TOLERANCE);
		};
		bool cos_zero(void);
		
		double X11(double _R, double xi, double eta, double _q) const;
		double X32(double _R, double xi) const;
		double X53(double _R, double xi) const;
		double Y11(double _R, double eta) const;
		double Y32(double _R, double eta) const;
		double Y53(double _R, double eta) const;
		double Z32(double _R, double _q, double eta, double z) const;
		double Z53(double _R, double _q, double eta, double z) const;
		double Y0(double _R, double xi, double eta) const;
		double Z0(double _R, double xi, double eta, double _q, double z) const;
		
		//
		// strain tensor components
		double get_exx(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double get_eyy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double get_ezz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double get_exy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double get_exz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double get_eyz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		
		//
		//
		// displacements
		//
		//
		double ux(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double uy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double uz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		//
		// displacement components
		// A
		double u1A (double x, double y, double z, double c, double L, double W, MotionType motion);
		double u2A (double x, double y, double z, double c, double L, double W, MotionType motion);
		double u3A (double x, double y, double z, double c, double L, double W, MotionType motion);
		// Ah
		double u1Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		double u2Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		double u3Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		// B
		double u1B (double x, double y, double z, double c, double L, double W, MotionType motion);
		double u2B (double x, double y, double z, double c, double L, double W, MotionType motion);
		double u3B (double x, double y, double z, double c, double L, double W, MotionType motion);
		// C
		double u1C (double x, double y, double z, double c, double L, double W, MotionType motion);
		double u2C (double x, double y, double z, double c, double L, double W, MotionType motion);
		double u3C (double x, double y, double z, double c, double L, double W, MotionType motion);
		//
		// strike fs
		// A
		double f1AS(double xi, double eta, double _q);
		double f2AS(double xi, double eta, double _q);
		double f3AS(double xi, double eta, double _q);
		// B
		double f1BS(double xi, double eta, double _q);
		double f2BS(double xi, double eta, double _q);
		double f3BS(double xi, double eta, double _q);
		// C
		double f1CS(double xi, double eta, double y, double z, double c);
		double f2CS(double xi, double eta, double y, double z, double c);
		double f3CS(double xi, double eta, double y, double z, double c);
		//
		// dip fs
		// A
		double f1AD(double xi, double eta, double _q);
		double f2AD(double xi, double eta, double _q);
		double f3AD(double xi, double eta, double _q);
		// B
		double f1BD(double xi, double eta, double _q);
		double f2BD(double xi, double eta, double _q);
		double f3BD(double xi, double eta, double _q);
		// C
		double f1CD(double xi, double eta, double y, double z, double c);
		double f2CD(double xi, double eta, double y, double z, double c);
		double f3CD(double xi, double eta, double y, double z, double c);
		//
		// tensile fs
		// A
		double f1AT(double xi, double eta, double _q);
		double f2AT(double xi, double eta, double _q);
		double f3AT(double xi, double eta, double _q);
		// B
		double f1BT(double xi, double eta, double _q);
		double f2BT(double xi, double eta, double _q);
		double f3BT(double xi, double eta, double _q);
		// C
		double f1CT(double xi, double eta, double y, double z, double c);
		double f2CT(double xi, double eta, double y, double z, double c);
		double f3CT(double xi, double eta, double y, double z, double c);
		//
		// displacement globals
		double Theta(double _q, double _R, double xi, double eta);
		double X(double xi, double _q);
		double I1(double _R, double xi, double eta, double _q);
		double I2(double _R, double xi, double eta, double _q);
		double I3(double _R, double eta, double _q);
		double I4(double _R, double xi, double eta, double _q);
		
		
		
		//
		//
		// displacement derivatives
		//
		//
		//
		// dx
		//    
		double duxdx(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double duydx(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double duzdx(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		//
		// dx components
		// A
		double j1A (double x, double y, double z, double c, double L, double W, MotionType motion);
		double j2A (double x, double y, double z, double c, double L, double W, MotionType motion);
		double j3A (double x, double y, double z, double c, double L, double W, MotionType motion);
		// Ah
		double j1Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		double j2Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		double j3Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		// B
		double j1B (double x, double y, double z, double c, double L, double W, MotionType motion);
		double j2B (double x, double y, double z, double c, double L, double W, MotionType motion);
		double j3B (double x, double y, double z, double c, double L, double W, MotionType motion);
		// C
		double j1C (double x, double y, double z, double c, double L, double W, MotionType motion);
		double j2C (double x, double y, double z, double c, double L, double W, MotionType motion);
		double j3C (double x, double y, double z, double c, double L, double W, MotionType motion);
		//
		// dx strike dfs
		// A
		double df1ASdx(double xi, double eta, double _q);
		double df2ASdx(double xi, double eta, double _q);
		double df3ASdx(double xi, double eta, double _q);
		// B
		double df1BSdx(double xi, double eta, double _q);
		double df2BSdx(double xi, double eta, double _q);
		double df3BSdx(double xi, double eta, double _q);
		// C
		double df1CSdx(double xi, double eta, double y, double z, double c);
		double df2CSdx(double xi, double eta, double y, double z, double c);
		double df3CSdx(double xi, double eta, double y, double z, double c);
		//
		// dx dip dfs
		// A
		double df1ADdx(double xi, double eta, double _q);
		double df2ADdx(double xi, double eta, double _q);
		double df3ADdx(double xi, double eta, double _q);
		// B
		double df1BDdx(double xi, double eta, double _q);
		double df2BDdx(double xi, double eta, double _q);
		double df3BDdx(double xi, double eta, double _q);
		// C
		double df1CDdx(double xi, double eta, double y, double z, double c);
		double df2CDdx(double xi, double eta, double y, double z, double c);
		double df3CDdx(double xi, double eta, double y, double z, double c);
		//
		// dx tensile dfs
		// A
		double df1ATdx(double xi, double eta, double _q);
		double df2ATdx(double xi, double eta, double _q);
		double df3ATdx(double xi, double eta, double _q);
		// B
		double df1BTdx(double xi, double eta, double _q);
		double df2BTdx(double xi, double eta, double _q);
		double df3BTdx(double xi, double eta, double _q);
		// C
		double df1CTdx(double xi, double eta, double y, double z, double c);
		double df2CTdx(double xi, double eta, double y, double z, double c);
		double df3CTdx(double xi, double eta, double y, double z, double c);
		//
		// dx globals
		double J1(double _R, double xi, double eta, double _q);
		double J2(double _R, double xi, double eta, double _q);
		double J3(double _R, double xi, double eta, double _q);
		double J4(double _R, double xi, double eta, double _q);
		double J5(double _R, double eta, double _q);
		double J6(double _R, double xi, double eta, double _q);
		double K1(double _R, double xi, double eta, double _q);
		double K2(double _R, double xi, double eta, double _q);
		double K3(double _R, double xi, double eta, double _q);
		double K4(double _R, double xi, double eta, double _q);
		double D11(double _R, double _dtil);
		
		//
		// dy
		//
		double duxdy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double duydy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double duzdy(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		//
		// dy components
		// A
		double k1A (double x, double y, double z, double c, double L, double W, MotionType motion);
		double k2A (double x, double y, double z, double c, double L, double W, MotionType motion);
		double k3A (double x, double y, double z, double c, double L, double W, MotionType motion);
		// Ah
		double k1Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		double k2Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		double k3Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		// B
		double k1B (double x, double y, double z, double c, double L, double W, MotionType motion);
		double k2B (double x, double y, double z, double c, double L, double W, MotionType motion);
		double k3B (double x, double y, double z, double c, double L, double W, MotionType motion);
		// C
		double k1C (double x, double y, double z, double c, double L, double W, MotionType motion);
		double k2C (double x, double y, double z, double c, double L, double W, MotionType motion);
		double k3C (double x, double y, double z, double c, double L, double W, MotionType motion);
		//
		// dy strike dfs
		// A
		double df1ASdy(double xi, double eta, double _q);
		double df2ASdy(double xi, double eta, double _q);
		double df3ASdy(double xi, double eta, double _q);
		// B
		double df1BSdy(double xi, double eta, double _q);
		double df2BSdy(double xi, double eta, double _q);
		double df3BSdy(double xi, double eta, double _q);
		// C
		double df1CSdy(double xi, double eta, double y, double z, double c);
		double df2CSdy(double xi, double eta, double y, double z, double c);
		double df3CSdy(double xi, double eta, double y, double z, double c);
		//
		// dy dip dfs
		// A
		double df1ADdy(double xi, double eta, double _q);
		double df2ADdy(double xi, double eta, double _q);
		double df3ADdy(double xi, double eta, double _q);
		// B
		double df1BDdy(double xi, double eta, double _q);
		double df2BDdy(double xi, double eta, double _q);
		double df3BDdy(double xi, double eta, double _q);
		// C
		double df1CDdy(double xi, double eta, double y, double z, double c);
		double df2CDdy(double xi, double eta, double y, double z, double c);
		double df3CDdy(double xi, double eta, double y, double z, double c);
		//
		// dy tensile dfs
		// A
		double df1ATdy(double xi, double eta, double _q);
		double df2ATdy(double xi, double eta, double _q);
		double df3ATdy(double xi, double eta, double _q);
		// B
		double df1BTdy(double xi, double eta, double _q);
		double df2BTdy(double xi, double eta, double _q);
		double df3BTdy(double xi, double eta, double _q);
		// C
		double df1CTdy(double xi, double eta, double y, double z, double c);
		double df2CTdy(double xi, double eta, double y, double z, double c);
		double df3CTdy(double xi, double eta, double y, double z, double c);
		//
		// dy globals
		double E(double _R, double xi, double eta, double _q);
		double F(double _R, double xi, double eta, double _q);
		double G(double _R, double xi, double eta, double _q);
		double H(double _R, double xi, double eta, double _q);
		double P(double _R, double xi, double eta, double _q);
		double Q(double _R, double xi, double eta, double y, double z, double c);
		
		//
		// dz
		//
		double duxdz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double duydz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		double duzdz(double x, double y, double z, double c, double L, double W, double US, double UD, double UT);
		//
		// dz components
		// A
		double l1A (double x, double y, double z, double c, double L, double W, MotionType motion);
		double l2A (double x, double y, double z, double c, double L, double W, MotionType motion);
		double l3A (double x, double y, double z, double c, double L, double W, MotionType motion);
		// Ah
		double l1Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		double l2Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		double l3Ah(double x, double y, double z, double c, double L, double W, MotionType motion);
		// B
		double l1B (double x, double y, double z, double c, double L, double W, MotionType motion);
		double l2B (double x, double y, double z, double c, double L, double W, MotionType motion);
		double l3B (double x, double y, double z, double c, double L, double W, MotionType motion);
		// C
		double l1C (double x, double y, double z, double c, double L, double W, MotionType motion);
		double l2C (double x, double y, double z, double c, double L, double W, MotionType motion);
		double l3C (double x, double y, double z, double c, double L, double W, MotionType motion);
		//
		// dz strike dfs
		// A
		double df1ASdz(double xi, double eta, double _q);
		double df2ASdz(double xi, double eta, double _q);
		double df3ASdz(double xi, double eta, double _q);
		// B
		double df1BSdz(double xi, double eta, double _q);
		double df2BSdz(double xi, double eta, double _q);
		double df3BSdz(double xi, double eta, double _q);
		// C
		double df1CSdz(double xi, double eta, double y, double z, double c);
		double df2CSdz(double xi, double eta, double y, double z, double c);
		double df3CSdz(double xi, double eta, double y, double z, double c);
		//
		// dz dip dfs
		// A
		double df1ADdz(double xi, double eta, double _q);
		double df2ADdz(double xi, double eta, double _q);
		double df3ADdz(double xi, double eta, double _q);
		// B
		double df1BDdz(double xi, double eta, double _q);
		double df2BDdz(double xi, double eta, double _q);
		double df3BDdz(double xi, double eta, double _q);
		// C
		double df1CDdz(double xi, double eta, double y, double z, double c);
		double df2CDdz(double xi, double eta, double y, double z, double c);
		double df3CDdz(double xi, double eta, double y, double z, double c);
		//
		// dz tensile dfs
		// A
		double df1ATdz(double xi, double eta, double _q);
		double df2ATdz(double xi, double eta, double _q);
		double df3ATdz(double xi, double eta, double _q);
		// B
		double df1BTdz(double xi, double eta, double _q);
		double df2BTdz(double xi, double eta, double _q);
		double df3BTdz(double xi, double eta, double _q);
		// C
		double df1CTdz(double xi, double eta, double y, double z, double c);
		double df2CTdz(double xi, double eta, double y, double z, double c);
		double df3CTdz(double xi, double eta, double y, double z, double c);
		//
		// dy globals
		double Ep(double _R, double xi, double eta, double _q);
		double Fp(double _R, double xi, double eta, double _q);
		double Gp(double _R, double xi, double eta, double _q);
		double Hp(double _R, double xi, double eta, double _q);
		double Pp(double _R, double xi, double eta, double _q);
        double Qp(double _R, double xi, double eta, double _q, double z);
        // ===================================================================
        //Added by KWS (untested), below is for change in gravity functions
        //
        // dg
        //
        double dg(double x, double y, double c, double L, double W, double US, double UD, double UT);
        //
        // dg components
        //
        double dSg(double x, double _p, double _q, double L, double W);
        double dDg(double x, double _p, double _q, double L, double W);
        double dTg(double x, double _p, double _q, double L, double W);
        double dCg(double x, double _p, double _q, double L, double W);
        //
        // dg globals
        double Sg(double xi, double eta, double _q);
        double Dg(double xi, double eta, double _q);
        double Tg(double xi, double eta, double _q);
        double Cg(double xi, double eta, double _q);
        double I2g(double _R, double xi, double eta, double _q);
	};
}

#endif
