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

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <stdlib.h>
#include <climits>

//#define HAVE_GEOGRAPHIC_LIB

#define assertThrow(COND, ERR_MSG) assert(COND);

#ifdef HAVE_GEOGRAPHIC_LIB
#include </opt/local/include/GeographicLib/Geodesic.hpp>
#include </opt/local/include/GeographicLib/Constants.hpp>
#endif

#ifndef _QUAKELIB_UTIL_H_
#define _QUAKELIB_UTIL_H_

#define EARTH_MEAN_RADIUS		6371000.0			// in m
#define EARTH_EQUATORIAL_RADIUS	6378137.0			// in m
#define EARTH_POLAR_RADIUS		6356752.314245		// in m
#define WGS_84_FLATTENING		1/298.257223563

#define CROSS_TOLERANCE 0.000001 

namespace quakelib {
	// Provides general setup information regarding the library.
	std::string SetupInfo(void);
	
	/*!
	 Represents a vector with associated operations.
	 */
	template <unsigned int dim>
	class Vec {
	private:
		double _x[dim];
		
		void rotation_matrix(double rot_matrix[dim*dim], const double &theta) throw(std::domain_error) {
			if (dim != 3) {
				throw std::domain_error("Vec::rotation_matrix");
			} else {
				//[r11,r12,r13,r21,r22,r23,r31,r32,r33]
				Vec<dim>		rot_axis;
				double			w[dim];
				unsigned int	i;
				
				rot_axis = this->unit_vector();
				for (i=0;i<dim;++i) w[i] = rot_axis[i];
				
				rot_matrix[0] = cos(theta) + pow(w[0],2.0) * (1.0 - cos(theta));
				rot_matrix[1] = w[0] * w[1] * (1.0 - cos(theta)) + w[2] * sin(theta);
				rot_matrix[2] = w[0] * w[2] * (1.0 - cos(theta)) - w[1] * sin(theta);
				rot_matrix[3] = w[1] * w[0] * (1.0 - cos(theta)) - w[2] * sin(theta);
				rot_matrix[4] = cos(theta) + pow(w[1],2.0) * (1.0 - cos(theta));
				rot_matrix[5] = w[1] * w[2] * (1.0 - cos(theta)) + w[0] * sin(theta);
				rot_matrix[6] = w[2] * w[0] * (1.0 - cos(theta)) + w[1] * sin(theta);
				rot_matrix[7] = w[2] * w[1] * (1.0 - cos(theta)) - w[0] * sin(theta);
				rot_matrix[8] = cos(theta) + pow(w[2],2.0) * (1.0 - cos(theta));
			}
		}
		
	public:
		//! Initialize a vector with all 0 elements.
		Vec(void) { for(unsigned int i=0;i<dim;++i) _x[i] = 0; };
		//! Initialize a vector of size dim with element values specified by the input array.
		Vec(const double vals[dim]) { for(unsigned int i=0;i<dim;++i) _x[i] = vals[i]; };
		//! Initialize a vector of dim>=1 as (x,0,...0). If dim < 1 an exception is thrown.
		Vec(const double &x) throw(std::out_of_range) { if (dim < 1) throw std::out_of_range("Vec"); _x[0] = x; };
		//! Initialize a vector of dim>=2 as (x,y,0,...0). If dim < 2 an exception is thrown.
		Vec(const double &x, const double &y) throw(std::out_of_range) { if (dim < 2) throw std::out_of_range("Vec"); _x[0] = x; _x[1] = y; };
		//! Initialize a vector of dim>=3 as (x,y,z,0...0). If dim < 3 an exception is thrown.
		Vec(const double &x, const double &y, const double &z) throw(std::out_of_range) { if (dim < 3) throw std::out_of_range("Vec"); _x[0] = x; _x[1] = y; _x[2] = z; };
		
		static Vec<dim> nan_vec(void) {
			Vec<dim> new_vec;
			for (unsigned int i=0;i<dim;++i) new_vec[i] = nan("");
			return new_vec;
		};
		
		//! Returns the cross product between this vector and the specified vector.
		//! TODO: Keep this?  If the result is below a specified tolerance, returns null vector
		Vec<dim> cross(const Vec<dim> &vec) const throw(std::domain_error) {
			if (dim != 3) throw std::domain_error("Vec::cross");
			double vals[dim];
			vals[0] = _x[1] * vec._x[2] - _x[2] * vec._x[1];
			vals[1] = _x[2] * vec._x[0] - _x[0] * vec._x[2];
			vals[2] = _x[0] * vec._x[1] - _x[1] * vec._x[0];
			Vec<dim> res(vals);
			if (res.mag() >= CROSS_TOLERANCE) return res;
			else return Vec<dim>();
		};
		
		//! Calculate the angle (in radians) between this vector and another.
		double vector_angle(const Vec<dim> &vec) const {
			double cos_r = this->dot_product(vec)/(this->mag() * vec.mag());
			if (cos_r > 1.0) cos_r = 1.0; if (cos_r < -1.0) cos_r = -1.0;
			return acos(cos_r);
		};
		
		//! Calculate the distance between this point and another.
		double dist(const Vec<dim> &a) const { return (*this-a).mag(); };
		
		//! Calculate the magnitude of this vector.
		double mag(void) const { double v=0; for(unsigned int i=0;i<dim;++i) v += pow(_x[i], 2); return sqrt(v); };
		
		//! Calculate the dot product of this vector with another.
		double dot_product(const Vec<dim> &a) const { double v=0; for (unsigned int i=0;i<dim;++i) v += _x[i]*a._x[i]; return v; };
		
		//! Get the unit vector along this vector.
		Vec<dim> unit_vector(void) const { double m = mag(), vals[dim]; for (unsigned int i=0;i<dim;++i) vals[i] = _x[i]/m; return Vec<dim>(vals); };
		
		//! Add a vector to this vector.
		Vec<dim> &operator+= (const Vec<dim> &a) { for (unsigned int i=0;i<dim;++i) _x[i] += a._x[i]; return *this; };
		//! Add a vector to this vector.
		const Vec<dim> operator+ (const Vec<dim> &a) const { Vec<dim> res(*this); res += a; return res; };
		
		//! Subtract a vector from this vector.
		Vec<dim> &operator-= (const Vec<dim> &a) { for (unsigned int i=0;i<dim;++i) _x[i] -= a._x[i]; return *this; };
		//! Subtract a vector from this vector.
		const Vec<dim> operator- (const Vec<dim> &a) const { Vec<dim> res(*this); res -= a; return res; };
		//! Unary negative vector operation.
		const Vec<dim> operator- (void) const { Vec<dim> res(*this); res *= -1.0; return res; };
		
		//! Multiply all elements of this vector by a constant.
		void operator*=(const double &a) { for (unsigned int i=0;i<dim;++i) _x[i] *= a; };
		//! Multiply all elements of this vector by a constant.
		const Vec<dim> operator* (const double &a) const { Vec<dim> res(*this); res *= a; return res; };
		
		//! Divide all elements of this vector by a constant.
		void operator/=(const double &a) { for (unsigned int i=0;i<dim;++i) _x[i] /= a; };
		//! Divide all elements of this vector by a constant.
		const Vec<dim> operator/ (const double &a) const { Vec<dim> res(*this); res /= a; return res; };
		
		//! Test equality of vectors.
		bool operator==(const Vec<dim> &a) { for (unsigned int i=0;i<dim;++i) if (_x[i] != a._x[i]) return false; return true; };
		//! Test inequality of vectors.
		bool operator!=(const Vec<dim> &a) { return !(*this==a); };

		//! Subscript operator to read vector elements.
		double operator[](const unsigned int &d) const throw(std::out_of_range) { if (dim <= d) throw std::out_of_range("Vec[]"); return _x[d]; };
		//! Subscript operator to read/modify vector elements.
		double &operator[](const unsigned int &d) throw(std::out_of_range) { if (dim <= d) throw std::out_of_range("Vec[]"); return _x[d]; };
		
		//! Vector rotation related functions
		Vec<dim> rotate_around_axis(Vec<dim> axis, const double &theta) throw(std::domain_error) {
			if (dim != 3) throw std::domain_error("Vec::rotate_around_axis");
			unsigned int	i, n;
			double			r[dim*dim], outv[dim];
			axis.rotation_matrix(r, theta);
			for (i=0;i<dim;++i) { outv[i] = 0; for (n=0;n<dim;++n) outv[i] += r[i*dim+n]*_x[n]; }
			return Vec<dim>(outv);
		}
		
		// Return number of bytes used by this object
		unsigned long mem_bytes(void) const { return sizeof(double)*dim; };
		
		//! Ordering operator for vectors, only for use with STL.
		bool operator< (const Vec<dim> &a) const { for (unsigned int i=0;i<dim;++i) if (_x[i] < a._x[i]) return true; return false; };
	};
	
	std::ostream& operator<<(std::ostream& os, const Vec<2>& pt);
	std::ostream& operator<<(std::ostream& os, const Vec<3>& pt);
	
	template <unsigned int ncols>
	class TensorRow {
	private:
		double			_entries[ncols];
		
	public:
		//! Multiply the elements of this row by a vector.
		double operator* (const Vec<ncols> &vec) const { double res; unsigned int i; for (i=0,res=0;i<ncols;++i) res += vec[i]*_entries[i]; return res; };
		//! Subscript operator to read vector elements.
		const double operator[](const unsigned int &c) const throw(std::out_of_range) { if (ncols <= c) throw std::out_of_range("TensorRow[]"); return _entries[c]; };
		//! Subscript operator to read/modify vector elements.
		double &operator[](const unsigned int &c) throw(std::out_of_range) { if (ncols <= c) throw std::out_of_range("TensorRow[]"); return _entries[c]; };
	};
	
	template <unsigned int ncols, unsigned int nrows>
	class Tensor {
	private:
		TensorRow<ncols>		_rows[nrows];
		
	public:
		//! Multiply the elements of this tensor by a vector.
		Vec<nrows> operator* (const Vec<ncols> &vec) const { Vec<nrows> res; unsigned int i; for (i=0;i<nrows;++i) res[i] = _rows[i]*vec; return res; };
		//! Subscript operator to read vector elements.
		const TensorRow<ncols> operator[](const unsigned int &r) const throw(std::out_of_range) { if (nrows <= r) throw std::out_of_range("Tensor[]"); return _rows[r]; };
		//! Subscript operator to read/modify vector elements.
		TensorRow<ncols> &operator[](const unsigned int &r) throw(std::out_of_range) { if (nrows <= r) throw std::out_of_range("Tensor[]"); return _rows[r]; };
	};
	
	std::ostream& operator<<(std::ostream& os, const Tensor<3,3>& tensor);
	std::ostream& operator<<(std::ostream& os, const TensorRow<3>& tensor);
	
	/*!
	 A latitude/longitude and depth under or on the Earth used in calculations.
	 Latitude must be in [-90, 90] degrees and longitude must be within [-180, 180] degrees.
	 Altitude indicates meters above ground (negative is underground).
	 */
	class LatLonDepth {
	private:
		double		_lat, _lon, _altitude;
		
	public:
		//! Default constructor - sets latitude, longitude and altitude to be 0.
		LatLonDepth(void) : _lat(0), _lon(0), _altitude(0) {};
		//! Constructor with specified latitude and longitude. Altitude defaults to 0 unless specified.
		LatLonDepth(const double &lat, const double &lon, const double &altitude=0) throw(std::invalid_argument) : _lat(lat), _lon(lon), _altitude(altitude) {
			if (fabs(lat)>90) throw std::invalid_argument("LatLonDepth::lat must be in [-90,90].");
			if (fabs(lon)>180) throw std::invalid_argument("LatLonDepth::lon must be in [-180,180].");
		}
		
		//! Get the latitude in degrees of this point.
		double lat(void) const { return _lat; };
		//! Get the longitude in degrees of this point.
		double lon(void) const { return _lon; };
		//! Get the altitude in meters of this point.
		double altitude(void) const { return _altitude; };
		
		//! Set the latitude in degrees of this point.
		void set_lat(const double &lat) throw(std::invalid_argument) { if (fabs(lat)>90) throw std::invalid_argument("LatLonDepth::lat must be in [-90,90]."); _lat = lat; };
		//! Set the longitude in degrees of this point.
		void set_lon(const double &lon) throw(std::invalid_argument) { if (fabs(lon)>180) throw std::invalid_argument("LatLonDepth::lon must be in [-180,180]."); _lon = lon; };
		//! Set the altitude in meters of this point.
		void set_altitude(const double &altitude) throw(std::invalid_argument) { _altitude = altitude; };
		
		// Return number of bytes used by this object
		unsigned long mem_bytes(void) const { return sizeof(double)*3; };
		
		//! Tests point equality (identical latitude, longitude and altitude).
		bool operator==(LatLonDepth &pt) const { return (_lat==pt._lat && _lon==pt._lon && _altitude==pt._altitude); };
		//! Tests point inequality.
		bool operator!=(LatLonDepth &pt) const { return (!(*this == pt)); };
	};
	
	std::ostream& operator<<(std::ostream& os, const LatLonDepth& pt);
	
	/*!
	 \brief A class to perform conversions between units.
	 
	 This class isn't needed because the conversions are particularly complicated,
	 but rather to remind the user what units are being converted to what.
	 The class also performs longitude/latitude conversion to Cartesian coordinates.
	 */
	class Conversion {
	private:
		//! Base coordinate for conversion.
		LatLonDepth			_base;
		
		void dist_vincenty(double &distance, double &start_azimuth, double &end_azimuth, const LatLonDepth &p1, const LatLonDepth &p2) const;
		
	public:
		//! Create a Conversion class with a base LatLonDepth of (0,0,0).
		Conversion(void) : _base() {};
		//! Create a Conversion class with the specified base LatLonDepth point.
		Conversion(const LatLonDepth &base) : _base(base) {};
		
		//! Set the base LatLonDepth used to convert to Cartesian coordinates.
		void set_base_lat_lon_depth(const LatLonDepth &new_base) { _base = new_base; };
		//! Get the base LatLonDepth.
		LatLonDepth get_base_lat_lon_depth(void) const { return _base; };
		
		//! Convert the specified point to a Cartesian coordinate using the base as (0,0,0).
		Vec<3> convert2xyz(const LatLonDepth &in_pt) const;
		
		//! Convert the specified Cartesian coordinate to latitude/longitude using the base as (0,0,0).
		LatLonDepth convert2LatLon(const Vec<3> &in_pt) const;
		
		//! Convert degrees to radians.
		double deg2rad(const double &degrees) const { return degrees*M_PI/180; };
		
		//! Convert radians to degrees.
		double rad2deg(const double &radians) const { return 180*radians/M_PI; };
		
		//! Convert years to seconds.
		double year2sec(const double &years) const { return years*365.25*24*60*60; };
		
		//! Convert seconds to years.
		double sec2year(const double &seconds) const { return seconds*(1.0/60.0)*(1.0/60.0)*(1.0/24.0)*(1.0/365.25); };
		
		//! Convert meters per second into centimeters per year.
		double m_per_sec2cm_per_yr(const double &mps) const { return year2sec(mps)*100; };
		
		//! Convert centimeters per year into meters per second.
		double cm_per_yr2m_per_sec(const double &cpy) const { return sec2year(cpy)/100; };
		
		//! Convert meters into kilometers.
		double m2km(const double &meters) const { return meters*1e-3; };
		
		//! Convert kilometers into meters.
		double km2m(const double &km) const { return km*1e3; };
		
		//! Convert centimeters into meters.
		double cm2m(const double &cm) const { return cm*1e-2; };
		
		//! Convert square kilometers into square meters.
		double sqkm2sqm(const double &sqkm) const { return sqkm*1e6; };
		
		//! Convert square meters into square kilometers.
		double sqm2sqkm(const double &sqm) const { return sqm*1e-6; };
		
		//! Convert pascals to bars.
		double pascal2bar(const double &pascals) const { return pascals*1e-5; };
		
		//! Convert bars to pascals.
		double bar2pascal(const double &bars) const { return bars*1e5; };
	};

	
	
	
	
	/*
	 Top class representing a dense matrix for computation.
	 */
	template <class CELL_TYPE>
	class DenseMatrix {
	protected:
		unsigned int	_ncols, _nrows;
	public:
		DenseMatrix(const unsigned int &ncols, const unsigned int &nrows) : _ncols(ncols), _nrows(nrows) {};
		virtual ~DenseMatrix(void) {};
		virtual void allocateRow(const unsigned int &row) = 0;
		virtual CELL_TYPE val(const unsigned int &row, const unsigned int &col) const = 0;
		virtual void setVal(const unsigned int &row, const unsigned int &col, const CELL_TYPE &new_val) = 0;
		virtual bool transpose(void) const = 0;
		virtual bool compressed(void) const = 0;
		virtual bool compressRow(const unsigned int &row, const float &ratio) = 0;
		virtual bool decompressRow(const unsigned int &row) = 0;
		virtual CELL_TYPE *getRow(CELL_TYPE *buf, const unsigned int &row) const = 0;
		virtual CELL_TYPE *getCol(CELL_TYPE *buf, const unsigned int &col) const = 0;
		virtual unsigned long mem_bytes(void) const = 0;
	};
	
	/*
	 Subclass representing a non-compressed dense matrix.
	 */
	template <class CELL_TYPE>
	class DenseStd : public DenseMatrix<CELL_TYPE> {
	protected:
		CELL_TYPE		*_data;
	public:
		DenseStd(const unsigned int &ncols, const unsigned int &nrows);
		virtual ~DenseStd(void);
		void allocateRow(const unsigned int &row) {};
		bool compressed(void) const { return false; };
		bool compressRow(const unsigned int &row, const float &ratio) { return false; };
		bool decompressRow(const unsigned int &row) { return false; };
		unsigned long mem_bytes(void) const;
	};
	
	/*
	 Sub-subclass representing a normally oriented non-compressed dense matrix.
	 */
	template <class CELL_TYPE>
	class DenseStdStraight : public DenseStd<CELL_TYPE> {
	public:
		DenseStdStraight(const unsigned int &ncols, const unsigned int &nrows) : DenseStd<CELL_TYPE>(ncols, nrows) {};
		virtual ~DenseStdStraight(void) {};
		bool transpose(void) const { return false; };
		CELL_TYPE val(const unsigned int &row, const unsigned int &col) const;
		void setVal(const unsigned int &row, const unsigned int &col, const CELL_TYPE &new_val);
		CELL_TYPE *getRow(CELL_TYPE *buf, const unsigned int &row) const;
		CELL_TYPE *getCol(CELL_TYPE *buf, const unsigned int &col) const;
	};
	
	/*
	 Sub-subclass representing a transposed non-compressed dense matrix.
	 */
	template <class CELL_TYPE>
	class DenseStdTranspose : public DenseStd<CELL_TYPE> {
	public:
		DenseStdTranspose(const unsigned int &ncols, const unsigned int &nrows) : DenseStd<CELL_TYPE>(nrows, ncols) {};
		virtual ~DenseStdTranspose(void) {};
		bool transpose(void) const { return true; };
		CELL_TYPE val(const unsigned int &row, const unsigned int &col) const;
		void setVal(const unsigned int &row, const unsigned int &col, const CELL_TYPE &new_val);
		CELL_TYPE *getRow(CELL_TYPE *buf, const unsigned int &row) const;
		CELL_TYPE *getCol(CELL_TYPE *buf, const unsigned int &col) const;
	};
	
	/*
	 Structure of compressed row run.
	 */
	template <class CELL_TYPE>
	struct RowRun {
		int			_length;
		CELL_TYPE	_val;
	};
	
	template <class CELL_TYPE>
	class CompressedRow {
	private:
		RowRun<CELL_TYPE>	*runs;
		CELL_TYPE			*raw_data;
		bool				compressed;
		unsigned int		data_dim;		// Dimension of the data, either number of runs (when compressed == true) or number of values (when compressed == false)
	public:
		CompressedRow(const unsigned int &ncols);
		void init(const unsigned int &dim);
		unsigned int getRowLen(void);
		bool compressRow(const float &ratio);
		bool decompressRow(void);
		CELL_TYPE val(const unsigned int &col) const;
		void setVal(const unsigned int &col, const CELL_TYPE &new_val);
		void copyRowContents(CELL_TYPE *dest) const;
		unsigned long mem_bytes(void) const;
	};
	
	template <class CELL_TYPE>
	class CompressedRowMatrix : public DenseMatrix<CELL_TYPE> {
	protected:
		CompressedRow<CELL_TYPE>	**_rows;
	public:
		CompressedRowMatrix(const unsigned int &ncols, const unsigned int &nrows);
		virtual ~CompressedRowMatrix(void) {};
		void allocateRow(const unsigned int &row);
		bool compressed(void) const { return true; };
		bool compressRow(const unsigned int &row, const float &ratio);
		bool decompressRow(const unsigned int &row);
		unsigned long mem_bytes(void) const;
	};
	
	template <class CELL_TYPE>
	class CompressedRowMatrixStraight : public CompressedRowMatrix<CELL_TYPE> {
	public:
		CompressedRowMatrixStraight(const unsigned int &ncols, const unsigned int &nrows) : CompressedRowMatrix<CELL_TYPE>(ncols, nrows) {};
		bool transpose(void) const { return false; };
		CELL_TYPE val(const unsigned int &row, const unsigned int &col) const;
		void setVal(const unsigned int &row, const unsigned int &col, const CELL_TYPE &new_val);
		CELL_TYPE *getRow(CELL_TYPE *buf, const unsigned int &row) const;
		CELL_TYPE *getCol(CELL_TYPE *buf, const unsigned int &col) const;
	};
	
	template <class CELL_TYPE>
	class CompressedRowMatrixTranspose : public CompressedRowMatrix<CELL_TYPE> {
	public:
		CompressedRowMatrixTranspose(const unsigned int &ncols, const unsigned int &nrows) : CompressedRowMatrix<CELL_TYPE>(ncols, nrows) {};
		bool transpose(void) const { return true; };
		CELL_TYPE val(const unsigned int &row, const unsigned int &col) const;
		void setVal(const unsigned int &row, const unsigned int &col, const CELL_TYPE &new_val);
		CELL_TYPE *getRow(CELL_TYPE *buf, const unsigned int &row) const;
		CELL_TYPE *getCol(CELL_TYPE *buf, const unsigned int &col) const;
	};
	
	template <class CELL_TYPE>
	class HierarchicalMatrix {
	protected:
		unsigned int	_rows, _cols;
		
	public:
		HierarchicalMatrix(void) : _rows(0), _cols(0) {};
		virtual ~HierarchicalMatrix(void) {};
		virtual void set(unsigned int row, unsigned int col, CELL_TYPE new_val) = 0;
		virtual CELL_TYPE get(unsigned int row, unsigned int col) const = 0;
	};
	
	template <class CELL_TYPE>
	class HFullMatrix : public HierarchicalMatrix<CELL_TYPE> {
	private:
		CELL_TYPE		*_vals;
		
	public:
		HFullMatrix(void) : HierarchicalMatrix<CELL_TYPE>(), _vals(NULL) {};
		
		HFullMatrix(unsigned int nrows, unsigned int ncols) {
			this->_rows = nrows;
			this->_cols = ncols;
			_vals = (CELL_TYPE*)malloc(this->_rows*this->_cols*sizeof(CELL_TYPE));
		}
		
		virtual ~HFullMatrix(void) {
			if (_vals) free(_vals);
		}
		
		virtual void set(unsigned int row, unsigned int col, CELL_TYPE new_val) {
			assertThrow(row < this->_rows, "Out of bounds");
			assertThrow(col < this->_cols, "Out of bounds");
			_vals[row*this->_cols+col] = new_val;
		}
		
		virtual CELL_TYPE get(unsigned int row, unsigned int col) const {
			assertThrow(row < this->_rows, "Out of bounds");
			assertThrow(col < this->_cols, "Out of bounds");
			return _vals[row*this->_cols+col];
		}
	};
	
	
	template <class CELL_TYPE>
	class HSuperMatrix : public HierarchicalMatrix<CELL_TYPE> {
	private:
		unsigned int					_mat_rows, _mat_cols;
		HierarchicalMatrix<CELL_TYPE>	**_matrices;
		
	public:
		virtual void set(unsigned int row, unsigned int col, CELL_TYPE new_val) {
			unsigned int sub_matrix_row, sub_matrix_col, new_row, new_col;
			
			assertThrow(row < this->_rows, "Out of bounds");
			assertThrow(col < this->_cols, "Out of bounds");
			
			sub_matrix_row = row/(this->_rows/_mat_rows);
			sub_matrix_col = col/(this->_cols/_mat_cols);
			
			assertThrow(sub_matrix_row < this->_mat_rows, "Out of bounds");
			assertThrow(sub_matrix_col < this->_mat_cols, "Out of bounds");
			
			new_row = row-sub_matrix_row*(this->_rows/_mat_rows);
			new_col = col-sub_matrix_col*(this->_cols/_mat_cols);
			
			_matrices[sub_matrix_row*_mat_rows+sub_matrix_col]->set(new_row, new_col, new_val);
		}
		
		virtual CELL_TYPE get(unsigned int row, unsigned int col) const {
			unsigned int sub_matrix_row, sub_matrix_col, new_row, new_col;
			
			assertThrow(row < this->_rows, "Out of bounds");
			assertThrow(col < this->_cols, "Out of bounds");
			
			sub_matrix_row = row/_mat_rows;
			sub_matrix_col = col/_mat_cols;
			
			assertThrow(sub_matrix_row < this->_mat_rows, "Out of bounds");
			assertThrow(sub_matrix_col < this->_mat_cols, "Out of bounds");
			
			new_row = row-sub_matrix_row*(this->_rows/_mat_rows);
			new_col = col-sub_matrix_col*(this->_cols/_mat_cols);
			
			return _matrices[sub_matrix_row*_mat_rows+sub_matrix_col]->get(new_row, new_col);
		}
	};
	
	
	// Rectangular shaped bound of Cartesian space
	// The bound consists of all points in [_min_bound[i], _max_bound[i]) for i in [0,dim)
	template <unsigned int dim>
	class RectBound {
	private:
		bool			_valid;
		Vec<dim>		_min_bound, _max_bound;
		
	public:
		RectBound(void) : _min_bound(), _max_bound(), _valid(false) {};
		RectBound(const Vec<dim> &new_bound1, const Vec<dim> &new_bound2) : _valid(true) {
			for (unsigned int i=0;i<dim;++i) {
				_min_bound[i] = fmin(new_bound1[i], new_bound2[i]);
				_max_bound[i] = fmax(new_bound1[i], new_bound2[i]);
			}
		};
		
		//! Test equality of bounds.
		bool operator==(const RectBound<dim> &a) { if (!_valid) return false; else return (_min_bound == a._min_bound && _max_bound == a._max_bound); };
		//! Test inequality of vectors.
		bool operator!=(const RectBound<dim> &a) { return !(*this==a); };
		
		// Whether this bound is valid or not
		bool valid(void) const { return _valid; };
		
		// Get the maximum length among all dimensions
		double max_length(void) const;
		
		// Get the center point of the bound
		Vec<dim> center() const;
		
		// Get the minimum and maximum of the bound
		Vec<dim> min_bound() const { return _min_bound; };
		Vec<dim> max_bound() const { return _max_bound; };
		
		// Return a bit flag indicating which subdivision the given point is in.
		// This flag will have bit N set if point dimension N is greater than the midpoint of the bound in dimension N.
		// This function does not check if the point is actually within the bound.
		unsigned int get_child_subdivision(const Vec<dim> &pt) const;
		
		// Returns a new RectBound specified by the subdivision bit flag.
		// For a bound of [a, b), the left part of the subbound is [a, (a+b)/2) and the right is [(a+b)/2, b)
		RectBound<dim> get_child_bound(const unsigned int &subdivision) const;
		
		// Returns true if test_pt is contained in this bound, false otherwise
		bool in_bound(const Vec<dim> &test_pt);
		
		// Extend the bound to include the specified point
		// If the bound is currently empty, the new bound will consist only of this point
		void extend_bound(const Vec<dim> &new_pt);
		
		// Return number of bytes used by this object
		unsigned long mem_bytes(void) const { return sizeof(bool)+_min_bound.mem_bytes()+_max_bound.mem_bytes(); };
	};
	
	std::ostream& operator<<(std::ostream& os, const RectBound<2>& pt);
	std::ostream& operator<<(std::ostream& os, const RectBound<3>& pt);
	
	enum TraverseCommand {
		CONTINUE_TRAVERSE,
		BACK_TO_PARENT,
		ALL_FINISH
	};
    
	template <unsigned int dim>
    class Octree;
    
	template <unsigned int dim>
    class OctreeTraverser {
    public:
        virtual TraverseCommand process_node(const quakelib::Octree<dim> &cur_node) = 0;
    };
    
	template <unsigned int dim>
	class Octree {
	private:
		// Constant describing the number of children for a branch (2^dim)
		static const unsigned int	_num_children = (0x01 << dim);
		
		// Pointer to array of children
		Octree<dim>		**_children;
		
		// Rectangular bound containing this node
		RectBound<dim>	_bound;
		
		// Whether this node is a leaf (true) or a branch (false)
		bool			_is_leaf;
		
		// If the node is a leaf, record the corresponding point and ID
		unsigned int	_pt_id;
		Vec<dim>		_pt;
		
		std::string vtk_string_internal(std::vector<int> &depth_rec, const int &cur_depth);
		
	public:
		// Initialize this octree branch to contain 2^N children for an N-dimensional tree.
		// All children start as NULL indicating this node is a leaf.
		Octree(const RectBound<dim> &bound);
		
		// Deallocate the octree and all children
		~Octree(void);
		
		// Adds the ID at the specified point in the tree, creating new nodes and moving other IDs as necessary
		// Returns true if the point was successfully added, false if the point was outside the tree bounds
		bool add_point(const Vec<dim> &pt, const unsigned int &pt_id);
		
		// Whether this node is a leaf or not
		bool is_leaf(void) const { return _is_leaf; };
		
		// Return the ID associated with this node
		unsigned int id(void) const { return _pt_id; };
		
		// Return the rectangular bound associated with this node
		RectBound<dim> bound(void) const { return _bound; };
		
		// Find out which leaf in the octree this point belongs in
		Octree *get_leaf_containing_point(const Vec<dim> &pt);
		
		// Get the maximum depth from this point in the octree
		unsigned int max_depth(void) const;
		
		// Get the number of active children only at this node
		unsigned int num_active_children(void) const;
		
		// Get the number of leaves for all children below this node
		unsigned int num_leaves(void) const;
		
		// Get the total number of descendents below this node
		unsigned int num_descendents(void) const;
		
		// Returns a string of the VTK ASCII representation for this tree
		std::string vtk_string(void);
		
		TraverseCommand traverse(OctreeTraverser<dim> *traverser) const;
        
		// Calculates the number of bytes in memory used by this node and all descendents
		unsigned long mem_bytes(void) const;
	};
	
    // Class to traverse an octree and get a set of all IDs within
	template <unsigned int dim>
    class GetOctreeIDs : public OctreeTraverser<dim> {
    public:
        std::set<unsigned int>	id_set;
        
        virtual TraverseCommand process_node(const quakelib::Octree<dim> &cur_node) {
            if (cur_node.is_leaf()) id_set.insert(cur_node.id());
            return CONTINUE_TRAVERSE;
        }
    };
    
	template <unsigned int dim>
    class BHAnalyzer : public OctreeTraverser<dim> {
    private:
        RectBound<dim>          _target_bound;
        double                  _theta;
        unsigned int            _num_checks;
        
    public:
        std::set<std::pair<unsigned int, unsigned int> >	run_bounds;
        
        BHAnalyzer(RectBound<dim> target_bound, double theta) : _target_bound(target_bound), _theta(theta), _num_checks(0) {};
        unsigned int num_checks(void) const { return _num_checks; };
        
        virtual TraverseCommand process_node(const quakelib::Octree<dim> &cur_node) {
            GetOctreeIDs<dim>           octree_ids;
            std::set<unsigned int>::const_iterator	it;
            unsigned int                first_block, prev_block;
            
            if (cur_node.bound().max_length()/cur_node.bound().center().dist(_target_bound.center()) < _theta || cur_node.num_active_children() == 0) {
                cur_node.traverse(&octree_ids);
                for (it=octree_ids.id_set.begin();it!=octree_ids.id_set.end();) {
                    first_block = prev_block = *it;
                    it++;
                    while (it != octree_ids.id_set.end() && prev_block+1 == *it) {
                        prev_block = *it;
                        it++;
                        _num_checks++;
                    }
                    run_bounds.insert(std::make_pair(first_block, prev_block));
                }
                return BACK_TO_PARENT;
            } else {
                return CONTINUE_TRAVERSE;
            }
        }
    };
}

#endif
