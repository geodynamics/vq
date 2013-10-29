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

#include "QuakeLibUtil.h"

std::string quakelib::SetupInfo(void) {
	std::stringstream		ss;
	
	ss << "QuakeLib v1.0" << std::endl;
	ss << "Using geographic lib: ";
#ifdef HAVE_GEOGRAPHIC_LIB
	ss << "TRUE" << std::endl;
#else
	ss << "FALSE" << std::endl;
#endif
	
	return ss.str();
}
std::ostream& quakelib::operator<<(std::ostream& os, const LatLonDepth& pt) {
	os << "[" << pt.lat() << "," << pt.lon();
	if (pt.altitude() != 0) os << "," << pt.altitude();
	os << "]";
	return os;
}

std::ostream& quakelib::operator<<(std::ostream& os, const Vec<2>& v) {
	for (int i=0;i<2;++i) os << (i == 0 ? "(" : ",") << v[i] << (i == 1 ? ")" : "");
	return os;
}

std::ostream& quakelib::operator<<(std::ostream& os, const Vec<3>& v) {
	for (int i=0;i<3;++i) os << (i == 0 ? "(" : ",") << v[i] << (i == 2 ? ")" : "");
	return os;
}

std::ostream& quakelib::operator<<(std::ostream& os, const Tensor<3,3>& tensor) {
	for (int i=0;i<3;++i) os << (i == 0 ? "[" : ",") << tensor[i] << (i == 2 ? "]" : "");
	return os;
}

std::ostream& quakelib::operator<<(std::ostream& os, const TensorRow<3>& tr) {
	for (int i=0;i<3;++i) os << (i == 0 ? "[" : ",") << tr[i] << (i == 2 ? "]" : "");
	return os;
}

std::ostream& quakelib::operator<<(std::ostream& os, const RectBound<2>& rb) {
	if (!rb.valid()) os << "[invalid bound]";
	os << "[" << rb.min_bound() << "; " << rb.max_bound() << "]";
	return os;
}

std::ostream& quakelib::operator<<(std::ostream& os, const RectBound<3>& rb) {
	if (!rb.valid()) os << "[invalid bound]";
	os << "[" << rb.min_bound() << "; " << rb.max_bound() << "]";
	return os;
}

/*!
 Converts an (x, y) value (given in km) to lattitude/longitude (in degrees decimal) relative to a set origin.
 This just inverts the conversion done in convert2xy
 */
#ifdef HAVE_GEOGRAPHIC_LIB
quakelib::LatLonDepth quakelib::Conversion::convert2LatLon(const Vec<3> &in_pt) const {
	double	a = GeographicLib::Constants::WGS84_a(),	// major radius
			f = GeographicLib::Constants::WGS84_f();	// flattening
	const	GeographicLib::Geodesic geod(a, f);
	double	s12, azi1, a12;
	double	new_lat, new_lon, new_altitude;
	
	s12 = sqrt( in_pt[0] * in_pt[0] + in_pt[1] * in_pt[1] );
	azi1 = rad2deg(atan2(in_pt[0],in_pt[1]));
	a12 = geod.Direct(_base.lat(),_base.lon(), azi1, s12, new_lat, new_lon);
	new_altitude = in_pt[2];
	
	return LatLonDepth(new_lat, new_lon, new_altitude);
}
#else
/*!
 Converts XYZ coordinate to latitude/longitude off _base coordinate using 
 Vincenty inverse formula for ellipsoids.
 
 Vincenty Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2011
 http://www.movable-type.co.uk/scripts/latlong.html
 
 from: Vincenty inverse formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the 
 Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975    
 http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf                                             
 */
quakelib::LatLonDepth quakelib::Conversion::convert2LatLon(const Vec<3> &in_pt) const {
	double a = EARTH_EQUATORIAL_RADIUS,  f = WGS_84_FLATTENING, b = (1-f)*a;
	double s = sqrt( in_pt[0] * in_pt[0] + in_pt[1] * in_pt[1] );
	double alpha1 = atan2(in_pt[0], in_pt[1]);
	double sinAlpha1 = sin(alpha1);
	double cosAlpha1 = cos(alpha1);
	
	double tanU1 = (1-f) * tan(deg2rad(_base.lat()));
	double cosU1 = 1.0 / sqrt((1 + tanU1*tanU1)), sinU1 = tanU1*cosU1;
	double sigma1 = atan2(tanU1, cosAlpha1);
	double sinAlpha = cosU1 * sinAlpha1;
	double cosSqAlpha = 1 - sinAlpha*sinAlpha;
	double uSq = cosSqAlpha * (a*a - b*b) / (b*b);
	double A = 1 + uSq/16384.0*(4096.0+uSq*(-768.0+uSq*(320.0-175.0*uSq)));
	double B = uSq/1024.0 * (256.0+uSq*(-128.0+uSq*(74.0-47.0*uSq)));
	
	double sigma = s / (b*A), sigmaP = 2*M_PI;
	double cos2SigmaM, sinSigma, cosSigma, deltaSigma;
	while (fabs(sigma-sigmaP) > 1e-14) {
		cos2SigmaM = cos(2*sigma1 + sigma);
		sinSigma = sin(sigma);
		cosSigma = cos(sigma);
		deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
		sigmaP = sigma;
		sigma = s / (b*A) + deltaSigma;
	}
	
	double tmp = sinU1*sinSigma - cosU1*cosSigma*cosAlpha1;
	double lat2 = atan2(sinU1*cosSigma + cosU1*sinSigma*cosAlpha1,
						(1-f)*sqrt(sinAlpha*sinAlpha + tmp*tmp));
	double lambda = atan2(sinSigma*sinAlpha1, cosU1*cosSigma - sinU1*sinSigma*cosAlpha1);
	double C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
	double L = lambda - (1-C) * f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
	double lon2 = fmod((deg2rad(_base.lon())+L+3*M_PI),(2*M_PI)) - M_PI;  // normalise to -180...+180
	
	return LatLonDepth(rad2deg(lat2), rad2deg(lon2), in_pt[2]);
}
#endif
/*!
 Converts a latitude/longitude point to (x, y) in meters from a set origin.
 */
#ifdef HAVE_GEOGRAPHIC_LIB
quakelib::Vec<3> quakelib::Conversion::convert2xyz(const LatLonDepth &in_pt) const {
	double	a = GeographicLib::Constants::WGS84_a(),  // major radius
	f = GeographicLib::Constants::WGS84_f();  // flattening
	const	GeographicLib::Geodesic geod(a, f);
	double	new_vals[3];
	double	s12, azi1, azi2, a12;
	
	// s12 is the distance from (latitude0, longitude0) to (in_pt.lat,in_pt.lon)
	// azi1 is the direction (angle in degrees clockwise from north) the point (in_pt.lat,in_pt.lon) is in
	a12 = geod.Inverse(_base.lat(),_base.lon(),in_pt.lat(),in_pt.lon(),s12,azi1,azi2);
	
	new_vals[0] = s12 * sin(deg2rad(azi1));
	new_vals[1] = s12 * cos(deg2rad(azi1));
	new_vals[2] = in_pt.altitude();
	
	return Vec<3>(new_vals);
}
#else
quakelib::Vec<3> quakelib::Conversion::convert2xyz(const LatLonDepth &in_pt) const {
	double		new_vals[3];
	double		dist, start_azimuth, end_azimuth;
	
	dist_vincenty(dist, start_azimuth, end_azimuth, LatLonDepth(_base.lat(), _base.lon()), LatLonDepth(in_pt.lat(), in_pt.lon()));

	new_vals[0] = dist*sin(start_azimuth);
	new_vals[1] = dist*cos(start_azimuth);
	new_vals[2] = in_pt.altitude();
	
	return Vec<3>(new_vals);
}
#endif

quakelib::VectorList quakelib::Conversion::convertArray2xyz(const FloatList &lats, const FloatList &lons) const {
	quakelib::VectorList conversions;
	
	for(FloatList::size_type lat_id = 0; lat_id != lats.size(); lat_id++) {
		for(FloatList::size_type lon_id = 0; lon_id != lons.size(); lon_id++) {
			conversions.push_back(convert2xyz(LatLonDepth(lats[lat_id],lons[lon_id],0.0)));
		}
	}
	
	return conversions;
};

/*!
 Calculates geodetic distance between two points specified by latitude/longitude using 
 Vincenty inverse formula for ellipsoids.
 
 Vincenty Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2011
 http://www.movable-type.co.uk/scripts/latlong.html
                                                                                                
 from: Vincenty inverse formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the 
       Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975    
       http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf                                             
 */
void quakelib::Conversion::dist_vincenty(double &distance, double &start_azimuth, double &end_azimuth, const LatLonDepth &p1, const LatLonDepth &p2) const {
	double		p1lon = deg2rad(p1.lon()), p1lat = deg2rad(p1.lat());
	double		p2lon = deg2rad(p2.lon()), p2lat = deg2rad(p2.lat());
	double		a = EARTH_EQUATORIAL_RADIUS;
	double		f = WGS_84_FLATTENING;
	double		b = (1-f)*a;
	double		L = p2lon-p1lon;
	double		U1 = atan((1-f)*tan(p1lat));
	double		U2 = atan((1-f)*tan(p2lat));
	double		sinU1 = sin(U1), cosU1 = cos(U1);
	double		sinU2 = sin(U2), cosU2 = cos(U2);
	double		lambda = L, lambdaP;
	int			max_iter = 100;
	double		sinLambda, cosLambda, sinSigma, cosSigma;
	double		sigma, sinAlpha, cosSqAlpha, cos2SigmaM;
	double		C, uSq, A, B, deltaSigma;
	
	do {
		sinLambda = sin(lambda);
		cosLambda = cos(lambda);
		sinSigma = sqrt(pow((cosU2*sinLambda),2) + pow((cosU1*sinU2-sinU1*cosU2*cosLambda),2));
		if (sinSigma == 0) {		// coincident points
			distance = start_azimuth = end_azimuth = 0;
			return;
		}
		cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
		sigma = atan2(sinSigma, cosSigma);
		sinAlpha = cosU1*cosU2*sinLambda/sinSigma;
		cosSqAlpha = 1.0 - sinAlpha*sinAlpha;
		cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha;
		if (isnan(cos2SigmaM)) cos2SigmaM = 0;	// equatorial line: cosSqAlpha = 0
		C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
		lambdaP = lambda;
		lambda = L + (1-C) * f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
	} while (fabs(lambda-lambdaP) > 1e-14 && --max_iter>0);
	
	if (max_iter == 0) {
		distance = start_azimuth = end_azimuth = nan("");
		return;
	};
	
	uSq = cosSqAlpha * (a*a - b*b)/(b*b);
	A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
	B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
	deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)- B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
	distance = b*A*(sigma-deltaSigma);
	start_azimuth = atan2(cosU2*sinLambda, cosU1*sinU2-sinU1*cosU2*cosLambda);
	end_azimuth = atan2(cosU1*sinLambda, -sinU1*cosU2+cosU1*sinU2*cosLambda);
}

template <class CELL_TYPE>
quakelib::DenseStd<CELL_TYPE>::DenseStd(const unsigned int &ncols, const unsigned int &nrows) : DenseMatrix<CELL_TYPE>(ncols, nrows) {
	_data = (CELL_TYPE*)valloc(sizeof(CELL_TYPE)*ncols*nrows);
	//assertThrow(_data, "Not enough memory to allocate matrix.");
	for (unsigned int i=0;i<ncols*nrows;++i) _data[i] = 0;
}

template <class CELL_TYPE>
quakelib::DenseStd<CELL_TYPE>::~DenseStd(void) {
	if (_data) free(_data);
	_data = NULL;
}

template <class CELL_TYPE>
quakelib::CompressedRowMatrix<CELL_TYPE>::CompressedRowMatrix(const unsigned int &ncols, const unsigned int &nrows) : DenseMatrix<CELL_TYPE>(ncols, nrows) {
	_rows = (CompressedRow<CELL_TYPE>**)valloc(nrows*sizeof(CompressedRow<CELL_TYPE>*));
	//assertThrow(_rows, "Not enough memory to allocate matrix.");
	// Wait to allocate rows until they are actually needed
	for (unsigned int i=0;i<nrows;++i) _rows[i] = NULL;
}

template <class CELL_TYPE>
void quakelib::CompressedRowMatrix<CELL_TYPE>::allocateRow(const unsigned int &row) {
	_rows[row] = new CompressedRow<CELL_TYPE>(this->_ncols);
	//assertThrow(_rows[row], "Not enough memory to allocate matrix.");
}

template <class CELL_TYPE>
bool quakelib::CompressedRowMatrix<CELL_TYPE>::compressRow(const unsigned int &row, const float &ratio) {
	//assertThrow(row >= 0 && row < this->_nrows, "Outside row bounds.");
	return _rows[row]->compressRow(ratio);
}

template <class CELL_TYPE>
bool quakelib::CompressedRowMatrix<CELL_TYPE>::decompressRow(const unsigned int &row) {
	//assertThrow(row >= 0 && row < this->_nrows, "Outside row bounds.");
	return _rows[row]->decompressRow();
}

template <class CELL_TYPE>
quakelib::CompressedRow<CELL_TYPE>::CompressedRow(const unsigned int &ncols) {
	data_dim = ncols;
	compressed = false;
	raw_data = NULL;
	runs = NULL;
	init(ncols);
}

template <class CELL_TYPE>
CELL_TYPE quakelib::DenseStdStraight<CELL_TYPE>::val(const unsigned int &row, const unsigned int &col) const { return this->_data[row*this->_ncols+col]; };

template <class CELL_TYPE>
unsigned long quakelib::DenseStd<CELL_TYPE>::mem_bytes(void) const {
	return sizeof(CELL_TYPE)*this->_ncols*this->_nrows;
}

template <class CELL_TYPE>
unsigned long quakelib::CompressedRow<CELL_TYPE>::mem_bytes(void) const {
	if (compressed) return sizeof(RowRun<CELL_TYPE>*)+sizeof(RowRun<CELL_TYPE>)*data_dim+sizeof(CELL_TYPE*)+sizeof(bool)+sizeof(unsigned int);
	else return sizeof(RowRun<CELL_TYPE>*)+sizeof(CELL_TYPE*)+sizeof(CELL_TYPE)*data_dim+sizeof(bool)+sizeof(unsigned int);
}

template <class CELL_TYPE>
unsigned long quakelib::CompressedRowMatrix<CELL_TYPE>::mem_bytes(void) const {
	unsigned long sum=0;
	for (int i=0;i<this->_nrows;++i) sum += _rows[i]->mem_bytes();
	return sum;
}

template <class CELL_TYPE>
void quakelib::CompressedRow<CELL_TYPE>::init(const unsigned int &dim) {
	raw_data = (CELL_TYPE*)valloc(sizeof(CELL_TYPE)*dim);
	for (unsigned int i=0;i<dim;++i) raw_data[i] = 0;
}

template <class CELL_TYPE>
unsigned int quakelib::CompressedRow<CELL_TYPE>::getRowLen(void) {
	unsigned int	i, uncompressed_len;
	if (compressed) {
		uncompressed_len = 0;
		for (i=0;i<data_dim;++i) uncompressed_len += runs[i]._length;
	} else {
		uncompressed_len = data_dim;
	}
	return uncompressed_len;
}

template <class CELL_TYPE>
CELL_TYPE quakelib::DenseStdTranspose<CELL_TYPE>::val(const unsigned int &row, const unsigned int &col) const { return this->_data[col*this->_nrows+row]; };

template <class CELL_TYPE>
void quakelib::DenseStdTranspose<CELL_TYPE>::setVal(const unsigned int &row, const unsigned int &col, const CELL_TYPE &new_val) { /*assert(col<width); assert(row<height);*/ this->_data[col*this->_nrows+row] = new_val; };

template <class CELL_TYPE>
void quakelib::DenseStdStraight<CELL_TYPE>::setVal(const unsigned int &row, const unsigned int &col, const CELL_TYPE &new_val) { this->_data[row*this->_ncols+col] = new_val; };

template <class CELL_TYPE>
CELL_TYPE quakelib::CompressedRowMatrixStraight<CELL_TYPE>::val(const unsigned int &row, const unsigned int &col) const { return this->_rows[row]->val(col); };

template <class CELL_TYPE>
void quakelib::CompressedRowMatrixStraight<CELL_TYPE>::setVal(const unsigned int &row, const unsigned int &col, const CELL_TYPE &new_val) { this->_rows[row]->setVal(col, new_val); };

template <class CELL_TYPE>
CELL_TYPE quakelib::CompressedRowMatrixTranspose<CELL_TYPE>::val(const unsigned int &row, const unsigned int &col) const { return this->_rows[col]->val(row); };

template <class CELL_TYPE>
void quakelib::CompressedRowMatrixTranspose<CELL_TYPE>::setVal(const unsigned int &row, const unsigned int &col, const CELL_TYPE &new_val) { this->_rows[col]->setVal(row, new_val); };

template <class CELL_TYPE>
CELL_TYPE *quakelib::DenseStdStraight<CELL_TYPE>::getRow(CELL_TYPE *buf, const unsigned int &row) const {
	return &this->_data[row*this->_ncols];
}

template <class CELL_TYPE>
CELL_TYPE *quakelib::DenseStdStraight<CELL_TYPE>::getCol(CELL_TYPE *buf, const unsigned int &col) const {
	//assertThrow(false, "Not implemented.");
}

template <class CELL_TYPE>
CELL_TYPE *quakelib::DenseStdTranspose<CELL_TYPE>::getRow(CELL_TYPE *buf, const unsigned int &row) const {
	//assertThrow(false, "Not implemented.");
}

template <class CELL_TYPE>
CELL_TYPE *quakelib::DenseStdTranspose<CELL_TYPE>::getCol(CELL_TYPE *buf, const unsigned int &col) const {
	return &this->_data[col*this->_nrows];
}

template <class CELL_TYPE>
CELL_TYPE *quakelib::CompressedRowMatrixStraight<CELL_TYPE>::getRow(CELL_TYPE *buf, const unsigned int &row) const {
	this->_rows[row]->copyRowContents(buf);
	return buf;
}

template <class CELL_TYPE>
CELL_TYPE *quakelib::CompressedRowMatrixStraight<CELL_TYPE>::getCol(CELL_TYPE *buf, const unsigned int &col) const {
	//assertThrow(false, "Not implemented.");
}

template <class CELL_TYPE>
CELL_TYPE *quakelib::CompressedRowMatrixTranspose<CELL_TYPE>::getRow(CELL_TYPE *buf, const unsigned int &row) const {
	//assertThrow(false, "Not implemented.");
}

template <class CELL_TYPE>
CELL_TYPE *quakelib::CompressedRowMatrixTranspose<CELL_TYPE>::getCol(CELL_TYPE *buf, const unsigned int &col) const {
	this->_rows[col]->copyRowContents(buf);
	return buf;
}

// This shouldn't be used for compressed rows in performance critical sections
template <class CELL_TYPE>
CELL_TYPE quakelib::CompressedRow<CELL_TYPE>::val(const unsigned int &col) const {
	unsigned int		i, cur_col;
	if (compressed) {
		cur_col = 0;
		for (i=0;i<data_dim;++i) {
			if (cur_col >= col) break;
			cur_col += runs[i]._length;
		}
		return runs[i]._val;
	} else {
		return raw_data[col];
	}
}

template <class CELL_TYPE>
void quakelib::CompressedRow<CELL_TYPE>::setVal(const unsigned int &col, const CELL_TYPE &new_val) {
	if (compressed) {
		//assertThrow(false, "Shouldn't access compressed row.");
	} else {
		raw_data[col] = new_val;
	}
}

// Decompress from runs to raw data
template <class CELL_TYPE>
bool quakelib::CompressedRow<CELL_TYPE>::decompressRow(void) {
	unsigned int new_data_dim;
	// If not already compressed we can return
	if (!compressed) return false;
	
	// Create the new raw_data array and copy the compressed data into it
	new_data_dim = getRowLen();
	init(new_data_dim);
	copyRowContents(raw_data);
	
	// Delete the old compressed data
	delete runs;
	runs = NULL;
	compressed = false;
	data_dim = new_data_dim;
	
	return true;
}

/*
 Try compressing the raw data of this row.  If the compression ratio (new size/old size)
 is greater than the specified ratio, the row remains uncompressed.
 */
template <class CELL_TYPE>
bool quakelib::CompressedRow<CELL_TYPE>::compressRow(const float &ratio) {
	double			last_val, compress_ratio;
	unsigned int	i, num_runs, cur_run, run_start;
	
	// If already compressed we can return
	if (compressed) return false;
	
	num_runs = 1;
	last_val = raw_data[0];
	// Count the number of runs in the data
	for (i=1;i<data_dim;++i) {
		if (raw_data[i] != last_val) {
			num_runs++;
			last_val = raw_data[i];
		}
	}
	
	// Compression ratio is new_size/old_size
	compress_ratio = (num_runs*sizeof(RowRun<CELL_TYPE>))/(data_dim*sizeof(CELL_TYPE));
	
	// If we achieve an acceptable compression ratio for this row, compress it
	if (compress_ratio <= ratio) {
		runs = (RowRun<CELL_TYPE>*)valloc(sizeof(RowRun<CELL_TYPE>)*num_runs);
		cur_run = 0;
		last_val = raw_data[0];
		run_start = 0;
		for (i=1;i<data_dim;++i) {
			if (raw_data[i] != last_val || i == data_dim-1) {
				runs[cur_run]._val = last_val;
				runs[cur_run]._length = i-run_start;
				cur_run++;
				last_val = raw_data[i];
				run_start = i;
			}
		}
		data_dim = num_runs;
		compressed = true;
		delete raw_data;
		raw_data = NULL;
	}
	return compressed;
}

template <class CELL_TYPE>
void quakelib::CompressedRow<CELL_TYPE>::copyRowContents(CELL_TYPE *dest) const {
	unsigned int		i, n, p;
	
	if (compressed) {
		for (p=0,i=0;i<data_dim;++i) {
			for (n=0;n<runs[i]._length;++n,++p) {
				dest[p] = runs[i]._val;
			}
		}
	} else {
		for (i=0;i<data_dim;++i) dest[i] = raw_data[i];
	}
}

template <unsigned int dim>
double quakelib::RectBound<dim>::max_length() const {
	double			max_len;
	unsigned int	i;
	
	if (!_valid) return nan("");
	
	max_len = _max_bound[0] - _min_bound[0];
	for (i=1;i<dim;++i) max_len = fmax(max_len, _max_bound[i] - _min_bound[i]);
	return max_len;
}

template <unsigned int dim>
quakelib::Vec<dim> quakelib::RectBound<dim>::center(void) const {
	if (!_valid) return Vec<dim>::nan_vec();
	
	return (_min_bound+_max_bound)/2.0;
}

template <unsigned int dim>
unsigned int quakelib::RectBound<dim>::get_child_subdivision(const Vec<dim> &pt) const {
	unsigned int		i, bitflag, subdivision = 0;
	
	if (!_valid) return 0;
	
	for (i=0,bitflag=1;i<dim;++i,bitflag*=2) if (pt[i] > (_min_bound[i]+_max_bound[i])/2.0) subdivision |= bitflag;
	
	return subdivision;
}

template <unsigned int dim>
quakelib::RectBound<dim> quakelib::RectBound<dim>::get_child_bound(const unsigned int &subdivision) const {
	double			new_min_bound[dim], new_max_bound[dim];
	unsigned int	i;
	
	if (!_valid) return RectBound();
	
	for (i=0;i<dim;++i) {
		new_min_bound[i] = ((subdivision >> i) & 0x01 ? (_min_bound[i]+_max_bound[i])/2.0 : _min_bound[i]);
		new_max_bound[i] = ((subdivision >> i) & 0x01 ? _max_bound[i] : (_min_bound[i]+_max_bound[i])/2.0);
	}
	
	return RectBound<dim>(new_min_bound, new_max_bound);
}

template <unsigned int dim>
bool quakelib::RectBound<dim>::in_bound(const Vec<dim> &test_pt) {
	unsigned int	i;
	
	if (!_valid) return false;
	
	for (i=0;i<dim;++i) {
		if (test_pt[i] < _min_bound[i] || test_pt[i] >= _max_bound[i]) return false;
	}
	
	return true;
}

template <unsigned int dim>
void quakelib::RectBound<dim>::extend_bound(const Vec<dim> &new_pt) {
	unsigned int	i;
	
	for (i=0;i<dim;++i) {
		if (!_valid) {
			_min_bound[i] = new_pt[i];
			_max_bound[i] = nextafter(new_pt[i], 2*new_pt[i]);
		} else {
			if (new_pt[i] < _min_bound[i]) _min_bound[i] = new_pt[i];
			else if (new_pt[i] >= _max_bound[i]) _max_bound[i] = nextafter(new_pt[i], 2*new_pt[i]);
		}
	}
	
	_valid = true;
}

template <unsigned int dim>
quakelib::Octree<dim>::Octree(const RectBound<dim> &bound) : _pt_id(UINT_MAX), _pt(), _is_leaf(true), _bound(bound) {
	unsigned int	i;
	
	_children = new Octree*[_num_children];
	for (i=0;i<_num_children;++i) _children[i] = NULL;
}

template <unsigned int dim>
quakelib::Octree<dim>::~Octree(void) {
	unsigned int	i;
	
	for (i=0;i<_num_children;++i) if (_children[i]) delete _children[i];
	delete _children;
}

template <unsigned int dim>
quakelib::Octree<dim> *quakelib::Octree<dim>::get_leaf_containing_point(const Vec<dim> &pt) {
	unsigned int				subdiv;
	
	subdiv = _bound.get_child_subdivision(pt);
	if (_children[subdiv] != NULL) return _children[subdiv]->get_leaf_containing_point(pt);
	else return this;
}

template <unsigned int dim>
bool quakelib::Octree<dim>::add_point(const Vec<dim> &pt, const unsigned int &pt_id) {
	unsigned int		subdiv, saved_id;
	Vec<dim>			saved_pt;
	bool				added;
	
	// Ensure the point is in the bound
	if (!_bound.in_bound(pt)) return false;
	
	// Get the subdivision this point will belong to
	subdiv = _bound.get_child_subdivision(pt);
	
	// If this node doesn't have any ID assigned
	if (_pt_id == UINT_MAX) {
		// And it is a leaf, then we can assign the ID
		if (_is_leaf) {
			_pt = pt;
			_pt_id = pt_id;
			return true;
		} else {
			// Otherwise, it is a branch and we need to go down
			subdiv = _bound.get_child_subdivision(pt);
			if (_children[subdiv] == NULL) _children[subdiv] = new Octree<dim>(_bound.get_child_bound(subdiv));
			return _children[subdiv]->add_point(pt, pt_id);
		}
	} else {
		// Otherwise, there is already an ID here so it must be a leaf
		// TODO: assert this is a leaf?
		
		// Mark this node as a branch now
		_is_leaf = false;
		
		// Save the current ID and pt for reinsertion and clear the current values
		saved_pt = _pt;
		saved_id = _pt_id;
		_pt = Vec<dim>();
		_pt_id = UINT_MAX;
		
		// Continue inserting the given pt and id
		subdiv = _bound.get_child_subdivision(pt);
		if (_children[subdiv] == NULL) _children[subdiv] = new Octree<dim>(_bound.get_child_bound(subdiv));
		added = _children[subdiv]->add_point(pt, pt_id);
		
		// And reinsert the saved pt and id
		subdiv = _bound.get_child_subdivision(saved_pt);
		if (_children[subdiv] == NULL) _children[subdiv] = new Octree<dim>(_bound.get_child_bound(subdiv));
		return added && _children[subdiv]->add_point(saved_pt, saved_id);
	}
}

template <unsigned int dim>
quakelib::TraverseCommand quakelib::Octree<dim>::traverse(OctreeTraverser<dim> *traverser) const {
	TraverseCommand		res;
	
    res = traverser->process_node(*this);
	if (res == ALL_FINISH) return ALL_FINISH;
	else if (res == BACK_TO_PARENT) return BACK_TO_PARENT;
	
	for (int i=0;i<_num_children;++i) {
		// If child doesn't exist, skip it
		if (!_children[i]) continue;
		
		// Otherwise, try traversing the child
		res = _children[i]->traverse(traverser);
		if (res == ALL_FINISH) return ALL_FINISH;
		else if (res == BACK_TO_PARENT) return CONTINUE_TRAVERSE;
	}
	
	return CONTINUE_TRAVERSE;
}

template <unsigned int dim>
unsigned int quakelib::Octree<dim>::max_depth(void) const {
	unsigned int	i, cur_max = 0, child_depth;
	
	for (i=0;i<_num_children;++i) {
		if (_children[i]) {
			child_depth = 1+_children[i]->max_depth();
			if (child_depth > cur_max) cur_max = child_depth;
		}
	}
	
	return cur_max;
}

template <unsigned int dim>
unsigned int quakelib::Octree<dim>::num_active_children(void) const {
	unsigned int	i, child_count = 0;
	
	for (i=0;i<_num_children;++i) if (_children[i]) child_count++;
	
	return child_count;
}

template <unsigned int dim>
unsigned int quakelib::Octree<dim>::num_leaves(void) const {
	unsigned int	i, total=0;
	
	if (_is_leaf) return 1;
	
	for (i=0;i<_num_children;++i) if (_children[i]) total += _children[i]->num_leaves();
	
	return total;
}

template <unsigned int dim>
unsigned int quakelib::Octree<dim>::num_descendents(void) const {
	unsigned int	i, total=0;
	
	for (i=0;i<_num_children;++i) if (_children[i]) total += _children[i]->num_descendents();
	
	return total+1;
}

template <unsigned int dim>
std::string quakelib::Octree<dim>::vtk_string(void) {
	std::vector<int>	rec;
	
	return vtk_string_internal(rec, 0);
}

template <unsigned int dim>
std::string quakelib::Octree<dim>::vtk_string_internal(std::vector<int> &depth_rec, const int &cur_depth) {
	unsigned int		i;
	std::stringstream	ss;
	Vec<dim>			bounds[2];
	
	if (cur_depth == 0) {
		ss << "# vtk DataFile Version 3.0" << std::endl;
		ss << "vtk output" << std::endl;
		ss << "ASCII" << std::endl;
		ss << "DATASET POLYDATA" << std::endl;
		ss << "POINTS " << _num_children*num_descendents() << " float" << std::endl;
	}
	
	bounds[0] = _bound.min_bound();
	bounds[1] = _bound.max_bound();
	for (i=0;i<dim;++i) bounds[1][i] = nextafter(bounds[1][i], -bounds[1][i]);
	
	for (i=0;i<_num_children;++i) ss << bounds[(i>>0)%2][0] << " " << bounds[(i>>1)%2][1] << " " << bounds[(i>>2)%2][2] << std::endl;
	depth_rec.push_back(cur_depth);
	
	for (i=0;i<_num_children;++i) if (_children[i]) ss << _children[i]->vtk_string_internal(depth_rec, cur_depth+1);
	
	if (cur_depth==0) {
		ss << "POLYGONS " << 2*dim*num_descendents() << " " << 2*dim*5*num_descendents() << std::endl;
		for (i=0;i<num_descendents();++i) {
			ss << "4 " << i*_num_children+0 << " " << i*_num_children+1 << " " << i*_num_children+3 << " " << i*_num_children+2 << std::endl;	// top
			ss << "4 " << i*_num_children+1 << " " << i*_num_children+5 << " " << i*_num_children+7 << " " << i*_num_children+3 << std::endl;	// left
			ss << "4 " << i*_num_children+5 << " " << i*_num_children+4 << " " << i*_num_children+6 << " " << i*_num_children+7 << std::endl;	// bottom
			ss << "4 " << i*_num_children+0 << " " << i*_num_children+2 << " " << i*_num_children+6 << " " << i*_num_children+4 << std::endl;	// right
			ss << "4 " << i*_num_children+3 << " " << i*_num_children+2 << " " << i*_num_children+6 << " " << i*_num_children+7 << std::endl;	// front
			ss << "4 " << i*_num_children+0 << " " << i*_num_children+4 << " " << i*_num_children+5 << " " << i*_num_children+1 << std::endl;	// back
		}
		ss << "CELL_DATA " << 2*dim*num_descendents() << std::endl;
		ss << "SCALARS element float" << std::endl;
		ss << "LOOKUP_TABLE default" << std::endl;
		for (i=0;i<2*dim*num_descendents();++i) {
			ss << depth_rec[i/(2*dim)] << std::endl;
		}
	}
	
	return ss.str();
}

template <unsigned int dim>
unsigned long quakelib::Octree<dim>::mem_bytes(void) const {
	unsigned long		i, total = 0;
	
	total += _num_children*sizeof(Octree<dim>*);
	for (i=0;i<_num_children;++i) if (_children[i]) total += _children[i]->mem_bytes();
	total += _bound.mem_bytes();
	total += sizeof(bool);
	total += sizeof(unsigned int);
	total += _pt.mem_bytes();
	
	return total;
}

template class quakelib::DenseStd<float>;
template class quakelib::DenseStdStraight<float>;
template class quakelib::DenseStdTranspose<float>;
template class quakelib::CompressedRowMatrix<float>;
template class quakelib::CompressedRowMatrixStraight<float>;
template class quakelib::CompressedRowMatrixTranspose<float>;

template class quakelib::DenseStd<double>;
template class quakelib::DenseStdStraight<double>;
template class quakelib::DenseStdTranspose<double>;
template class quakelib::CompressedRowMatrix<double>;
template class quakelib::CompressedRowMatrixStraight<double>;
template class quakelib::CompressedRowMatrixTranspose<double>;

template class quakelib::HFullMatrix<float>;
template class quakelib::HSuperMatrix<float>;

template class quakelib::HFullMatrix<double>;
template class quakelib::HSuperMatrix<double>;

template class quakelib::RectBound<2>;
template class quakelib::RectBound<3>;

template class quakelib::Octree<2>;
template class quakelib::Octree<3>;
