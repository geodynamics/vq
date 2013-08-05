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

void quakelib::EQSimGeometry::write(std::ostream &out_stream) const {
	quakelib::EQSimGeomSectionMap::const_iterator		sit;
	
	// TODO: write the metadata
	out_stream	<< "200"
				<< " " << num_sections()
				<< " " << num_vertices()
				<< " " << num_triangles()
				<< " " << num_rectangles()
				<< " " << lat_lo()
				<< " " << lat_hi()
				<< " " << lon_lo()
				<< " " << lon_hi()
				<< " " << depth_lo()
				<< " " << depth_hi()
				<< std::endl;
	
	// Write the sections
	for (sit=sections.begin();sit!=sections.end();++sit) {
		sit->second.write(out_stream);
	}
}

void quakelib::EQSimGeometrySection::write(std::ostream &out_stream) const {
	quakelib::EQSimGeomVertexMap::const_iterator			vit;
	quakelib::EQSimGeomTriangleMap::const_iterator		tit;
	quakelib::EQSimGeomRectangleMap::const_iterator		rit;
	
	out_stream	<< "201"
				<< " " << _sid
				<< " " << _name
				<< " " << num_vertices()
				<< " " << num_triangles()
				<< " " << num_rectangles()
				<< " " << lat_lo()
				<< " " << lat_hi()
				<< " " << lon_lo()
				<< " " << lon_hi()
				<< " " << depth_lo()
				<< " " << depth_hi()
				<< " " << das_lo()
				<< " " << das_hi()
				<< " " << fid()
				<< std::endl;
	
	// Write the vertices
	for (vit=vertices.begin();vit!=vertices.end();++vit) {
		vit->second.write(out_stream);
	}

	// Write the triangles
	for (tit=triangles.begin();tit!=triangles.end();++tit) {
		tit->second.write(out_stream);
	}
	
	// Write the rectangles
	for (rit=rectangles.begin();rit!=rectangles.end();++rit) {
		rit->second.write(out_stream);
	}
}

void quakelib::EQSimGeometryVertex::write(std::ostream &out_stream) const {
	out_stream	<< "202"
	<< " " << _index
	<< " " << _loc.lat()
	<< " " << _loc.lon()
	<< " " << _loc.altitude()
	<< " " << _das
	<< " " << _trace_flag
	<< std::endl;
}

void quakelib::EQSimGeometryTriangle::write(std::ostream &out_stream) const {
	out_stream	<< "203"
	<< " " << _index
	<< " " << _vertex[0]
	<< " " << _vertex[1]
	<< " " << _vertex[2]
	<< " " << _rake
	<< " " << _slip_rate
	<< " " << _aseis_factor
	<< " " << _strike
	<< " " << _dip
	<< std::endl;
}

void quakelib::EQSimGeometryRectangle::write(std::ostream &out_stream) const {
	out_stream	<< "204"
	<< " " << _index
	<< " " << _vertex[0]
	<< " " << _vertex[1]
	<< " " << _vertex[2]
	<< " " << _vertex[3]
	<< " " << _rake
	<< " " << _slip_rate
	<< " " << _aseis_factor
	<< " " << _strike
	<< " " << _dip
	<< " " << _perfect_flag
	<< std::endl;
}

void quakelib::EQSimGeometrySection::fill_vertex_finder(const Conversion &conv, VertexFinder &finder) {
	quakelib::EQSimGeomVertexMap::const_iterator		it;
	Vec<3>								coord;
	
	// Fill the vertex finder with the vertex locations
	for (it=vertices.begin();it!=vertices.end();++it) {
		coord = conv.convert2xyz(it->second.loc());
		finder.add_point(it->first, coord);
	}
}

void quakelib::VertexFinder::add_point(const UIndex &index, const Vec<3> &loc) {
	x_vertex.insert(std::make_pair(loc[0], index));
	y_vertex.insert(std::make_pair(loc[1], index));
	z_vertex.insert(std::make_pair(loc[2], index));
}

// Returns a list of vertex indices for vertices within the specified Euclidean distance of the specified location.
void quakelib::VertexFinder::find_vertices_in_distance(std::set<UIndex> &vertex_indices, const Vec<3> &loc, const double &distance) {
	std::multimap<double, UIndex>::const_iterator	it, x_min, x_max, y_min, y_max, z_min, z_max;
	std::map<UIndex, double>						x_pos, y_pos;
	Vec<3>												test_pos;
	
	// Get the iterators to the specified min/max vertices
	// This helps narrow down the range before we do a direct distance calculation
	x_min = x_vertex.lower_bound(loc[0]-distance);
	x_max = x_vertex.upper_bound(loc[0]+distance);
	y_min = y_vertex.lower_bound(loc[1]-distance);
	y_max = y_vertex.upper_bound(loc[1]+distance);
	z_min = z_vertex.lower_bound(loc[2]-distance);
	z_max = z_vertex.upper_bound(loc[2]+distance);
	
	// Record the vertices within the specified ranges
	for (it=x_min;it!=x_max;++it) x_pos.insert(std::make_pair(it->second, it->first));
	for (it=y_min;it!=y_max;++it) y_pos.insert(std::make_pair(it->second, it->first));
	
	// Find which vertices are in all three ranges and are within the required distance
	for (it=z_min;it!=z_max;++it) {
		if (x_pos.count(it->second) > 0 && y_pos.count(it->second) > 0) {
			test_pos[0] = x_pos[it->second];
			test_pos[1] = y_pos[it->second];
			test_pos[2] = it->first;
			if (loc.dist(test_pos) <= distance) vertex_indices.insert(it->second);
		}
	}
}

quakelib::EQSimGeometryVertex &quakelib::EQSimGeometrySection::new_vertex(void) throw(std::out_of_range) {
	UIndex											new_index;
	EQSimGeometryVertex								new_vertex;
	std::pair<EQSimGeomVertexMap::iterator, bool>	res;
	
	if (vertices.size() == 0) new_index = 1;
	else new_index = vertices.rbegin()->first + 1;
	
	new_vertex.set_index(new_index);
	res = vertices.insert(std::make_pair(new_index, new_vertex));
	if (!res.second) throw std::out_of_range("new_vertex");
	return res.first->second;
}

quakelib::EQSimGeometryTriangle &quakelib::EQSimGeometrySection::new_triangle(void) throw(std::out_of_range) {
	UIndex											new_index;
	EQSimGeometryTriangle							new_triangle;
	std::pair<EQSimGeomTriangleMap::iterator, bool>	res;
	
	if (triangles.size() == 0) new_index = 1;
	else new_index = triangles.rbegin()->first + 1;
	
	new_triangle.set_index(new_index);
	res = triangles.insert(std::make_pair(new_index, new_triangle));
	if (!res.second) throw std::out_of_range("new_triangle");
	return res.first->second;
}

quakelib::EQSimGeometryRectangle &quakelib::EQSimGeometrySection::new_rectangle(void) throw(std::out_of_range) {
	UIndex										new_index;
	EQSimGeometryRectangle						new_rectangle;
	std::pair<EQSimGeomRectangleMap::iterator, bool>	res;
	
	if (rectangles.size() == 0) new_index = 1;
	else new_index = rectangles.rbegin()->first + 1;
	
	new_rectangle.set_index(new_index);
	res = rectangles.insert(std::make_pair(new_index, new_rectangle));
	if (!res.second) throw std::out_of_range("new_rectangle");
	return res.first->second;
}

quakelib::EQSimGeometrySection &quakelib::EQSimGeometry::new_section(void) throw(std::out_of_range) {
	UIndex									new_sid;
	EQSimGeometrySection					new_section;
	std::pair<EQSimGeomSectionMap::iterator, bool>	res;
	
	if (sections.size() == 0) new_sid = 1;
	else new_sid = sections.rbegin()->first + 1;
	
	new_section.set_sid(new_sid);
	res = sections.insert(std::make_pair(new_sid, new_section));
	if (!res.second) throw std::out_of_range("new_section");
	return res.first->second;
}

bool quakelib::EQSimGeometrySection::is_perfect_rect(const Vec<3> v[4], const double &tolerance) {
	Vec<3>		v01, v02, v03, v12, v23, v30;
	double		vol, rel_diff1, rel_diff2, rel_vol;
	
	// Check if the points are coplanar by getting the volume
	// of the parallelepiped formed by vectors between them
	v01 = v[1] - v[0]; v02 = v[2] - v[0]; v03 = v[3] - v[0];
	vol = v01.dot_product(v02.cross(v03));
	// Check that the volume is relatively small compared to the edge size
	rel_vol = fabs(vol / (v01.mag() * v02.mag() * v03.mag()));
	if (rel_vol < tolerance) {
		// Ensure that parallel sides are of the same length (within a specified relative tolerance)
		v12 = v[2] - v[1]; v23 = v[3] - v[2]; v30 = v[0] - v[3];
		rel_diff1 = fabs((v01.mag()-v23.mag())/v01.mag());
		rel_diff2 = fabs((v12.mag()-v30.mag())/v12.mag());
		if (rel_diff1 < tolerance && rel_diff2 < tolerance) return true;
	}
	return false;
}

quakelib::EQSimGeometryRectangle &quakelib::EQSimGeometrySection::rect_from_lat_lon_depth(const std::vector<LatLonDepth> &locs,
																	  const Vec<2> &slip,
																	  const double &tolerance) throw (std::invalid_argument) {
	EQSimGeometryVertex			*verts[4];
	EQSimGeometryRectangle		*new_rect;
	unsigned int				i;
	Vec<3>						v[4];
	Conversion					conv(locs[0]);		// base Cartesian conversion on the first point
	double						high_dep, low_dep;
	double						deps[4], *min_dep1, *min_dep2, *max_dep, strike, dip, rake;
	Vec<3>						base_vec, strike_vec;
	int							ind1, ind2, ind3;
	Vec<3>						dip_vec, horiz_vec;
	Vec<2>						horiz_rake_vec, slip_vec;
	
	if (locs.size() != 4) throw std::invalid_argument("Must specify 4 points for rectangle.");
	
	// Ensure there are only two unique depths
	for (i=0;i<4;++i) deps[i] = locs[i].altitude();
	high_dep = fmax(deps[0], fmax(deps[1], fmax(deps[2], deps[3])));
	low_dep = fmin(deps[0], fmin(deps[1], fmin(deps[2], deps[3])));
	for (i=0;i<4;++i) if (deps[i] != high_dep && deps[i] != low_dep) throw std::invalid_argument("Points must have at most two unique depths.");
	
	// Make new vertices with the specified lat/lon/depth points
	for (i=0;i<4;++i) {
		verts[i] = &(this->new_vertex());
		verts[i]->set_loc(locs[i]);
	}
	
	// Make a new rectangle with the created vertices
	new_rect = &(this->new_rectangle());
	for (i=0;i<4;++i) new_rect->set_vertex(i, verts[i]->index());
	
	// Convert the lat/lon to Cartesian coordinates
	for (i=0;i<4;++i) v[i] = conv.convert2xyz(locs[i]);
	// Determine if this is a perfect rectangle
	new_rect->set_perfect_flag(is_perfect_rect(v, tolerance));
	
	// Calculate the dip, strike and rake of the fault
	// Start by finding the two shallowest points which will define the strike angle
	// Get the first shallowest point, then get the second after removing the first
	min_dep1 = std::max_element(deps, deps+4);
	*min_dep1 = -INFINITY;
	min_dep2 = std::max_element(deps, deps+4);
	max_dep = std::min_element(deps, deps+4);
	ind1 = min_dep1 - deps;
	ind2 = min_dep2 - deps;
	ind3 = max_dep - deps;
	
	// Calculate strike
	// TODO: project vertices to surface(?)
	// TODO: handle > 180 degree cases
	double vals[3] = {0, 1, 0};
	base_vec = Vec<3>(vals);
	strike_vec = v[ind2] - v[ind1];
	strike = acos(base_vec.dot_product(strike_vec)/(strike_vec.mag()));
	new_rect->set_strike(conv.rad2deg(strike));
	
	// Calculate dip
	dip_vec = v[ind2] - v[(ind2+1)%4];
	horiz_vec = dip_vec;
	horiz_vec[2] = 0;		// set Z value
	dip = acos(horiz_vec.dot_product(dip_vec)/(dip_vec.mag()*horiz_vec.mag()));
	new_rect->set_dip(conv.rad2deg(dip));
	
	// Calculate rake
	double rake_vals[2] = {1, 0};
	horiz_rake_vec = Vec<2>(rake_vals);
	if (slip.mag() == 0) slip_vec = Vec<2>(rake_vals);
	else slip_vec = slip;
	rake = acos(horiz_rake_vec.dot_product(slip_vec)/(slip_vec.mag()));
	new_rect->set_rake(conv.rad2deg(rake));
	
	// Calculate the slip rate
	new_rect->set_slip_rate(slip.mag());
	
	return *new_rect;
}

// Calculates the distance along strike and trace values for the different vertices in this section
// TODO: implement this
void quakelib::EQSimGeometrySection::calculate_trace_das(void) {
	double							min_depth;
	EQSimGeomRectangleMap::iterator	it;
	
	// Trace points will only be on the minimum section depth
	min_depth = depth_hi();
	
	// Travel along the rectangles and record which points are closest to each other
	for (it=rectangles.begin();it!=rectangles.end();++it) {
		
		if (it->second.vertex(0) != min_depth) {
			//it->second.set_trace_flag(NOT_ON_TRACE);
		} else {
			
		}
	}
}

quakelib::UIndex quakelib::EQSimGeometrySection::add_element(const ElementRect &elem, const Conversion &conv) {
	EQSimGeometryRectangle	&rect = new_rectangle();
	Vec<3>					v[4];
	unsigned int			i;
	
	// Convert the element vertex values to the EqSim file format
	for (i=0;i<4;++i) {
		EQSimGeometryVertex		&vert = new_vertex();
		v[i] = elem.vert(i);
		vert.set_loc(conv.convert2LatLon(v[i]));
		vert.set_das(elem.das(i));
		vert.set_trace_flag(elem.trace_flag(i));
		rect.set_vertex(i, vert.index());
	}
	
	// Set rectangle values
	rect.set_dip(conv.rad2deg(elem.dip()));
	rect.set_rake(conv.rad2deg(elem.rake()));
	rect.set_slip_rate(elem.slip_rate());
	rect.set_perfect_flag(is_perfect_rect(v, 1.0e-2));
	rect.set_aseismic(elem.aseismic());
	
	return rect.index();
}

// Finds vertices which are duplicates of each other within a specified distance.
// The distance parameter indicates how close two vertices must be (in XYZ coordinates) to be considered a duplicate.
void quakelib::EQSimGeometrySection::remove_duplicate_vertices(const double &distance) {
	EQSimGeomVertexMap::const_iterator		it;
	std::set<UIndex>						matching_vertices;
	std::set<UIndex>::const_iterator		vit;
	quakelib::VertexFinder					finder;
	Conversion								conv((vertices.begin())->second.loc());			// do x, y, z conversions based on the first point
	Vec<3>									vert_loc;
	quakelib::EQSimGeomRemapping			dup_map;
	
	// Create an index of vertex locations
	fill_vertex_finder(conv, finder);
	
	// Go through all the vertices
	for (it=vertices.begin();it!=vertices.end();++it) {
		// Check that we haven't already merged this vertex
		if (dup_map.vert_remap.count(it->first) == 0) {
			matching_vertices.clear();
			vert_loc = conv.convert2xyz(it->second.loc());
			finder.find_vertices_in_distance(matching_vertices, vert_loc, distance);
			// Remove the original vertex
			matching_vertices.erase(it->first);
			for (vit=matching_vertices.begin();vit!=matching_vertices.end();++vit) {
				dup_map.vert_remap.insert(std::make_pair(*vit, it->first));
			}
		}
	}
	
	// Apply the remapping
	apply_remap(dup_map);
	
	// Remove vertices no longer in use
	remove_unused_vertices();
}

void quakelib::EQSimGeometrySection::remove_unused_vertices(void) {
	std::set<UIndex>						used_verts;
	EQSimGeomTriangleMap::const_iterator	tit;
	EQSimGeomRectangleMap::const_iterator	rit;
	EQSimGeomVertexMap::const_iterator		vit;
	UIndex									ind;
	
	// Determine which vertices are used by triangles
	for (tit=triangles.begin();tit!=triangles.end();++rit) {
		used_verts.insert(tit->second.vertex(0));
		used_verts.insert(tit->second.vertex(1));
		used_verts.insert(tit->second.vertex(2));
	}
	
	// Determine which vertices are used by rectangles
	for (rit=rectangles.begin();rit!=rectangles.end();++rit) {
		used_verts.insert(rit->second.vertex(0));
		used_verts.insert(rit->second.vertex(1));
		used_verts.insert(rit->second.vertex(2));
		used_verts.insert(rit->second.vertex(3));
	}
	
	// Go through the vertices - if one is not being used then remove it
	for (vit=vertices.begin();vit!=vertices.end();) {
		ind = vit->first;
		vit++;
		if (used_verts.count(ind) == 0) {
			vertices.erase(ind);
		}
	}
}

// Creates a remapping of vertex indices so that the indices are contiguous.
// Input is the index number to start remapping with, and the function returns the index number that mapping ended on
quakelib::UIndex quakelib::EQSimGeometrySection::remap_continguous(const UIndex &start_ind) {
	EQSimGeomVertexMap::const_iterator		vit;
	UIndex									cur_ind = start_ind;
	EQSimGeomRemapping						remap;
	
	for (vit=vertices.begin();vit!=vertices.end();++vit) {
		remap.vert_remap.insert(std::make_pair(vit->first, cur_ind));
		cur_ind++;
	}
	
	apply_remap(remap);
	
	return cur_ind;
}

void quakelib::EQSimGeometryTriangle::apply_remap(const EQSimGeomRemapping &remap) {
	EQSimGeomIndexMap::const_iterator		it;
	
	// Remap the triangle index
	it = remap.tri_remap.find(_index);
	if (it != remap.tri_remap.end()) _index = it->second;
	
	// Remap the vertices
	for (int i=0;i<3;++i) {
		if (remap.vert_remap.count(_vertex[i]) > 0) {
			_vertex[i] = remap.vert_remap.find(_vertex[i])->second;
		}
	}
}

void quakelib::EQSimGeometryRectangle::apply_remap(const EQSimGeomRemapping &remap) {
	EQSimGeomIndexMap::const_iterator		it;
	
	// Remap the rectangle index
	it = remap.rect_remap.find(_index);
	if (it != remap.rect_remap.end()) _index = it->second;
	
	// Remap the vertices
	for (int i=0;i<4;++i) {
		if (remap.vert_remap.count(_vertex[i]) > 0) {
			std::cerr << "Changed vert " << i << " of rectangle " << _index << " from " << _vertex[i] << " to " << remap.vert_remap.find(_vertex[i])->second << "." << std::endl;
			_vertex[i] = remap.vert_remap.find(_vertex[i])->second;
		}
	}
}

void quakelib::EQSimGeometrySection::apply_remap(const EQSimGeomRemapping &remap) {
	EQSimGeomTriangleMap::iterator		tit;
	EQSimGeomRectangleMap::iterator		rit;
	EQSimGeomVertexMap::iterator		vit;
	EQSimGeomIndexMap::const_iterator	it;
	EQSimGeomVertexMap					new_vert_map;
	
	// Renumber vertex indices
	for (it=remap.vert_remap.begin();it!=remap.vert_remap.end();++it) {
		new_vert_map[it->second] = vertices[it->first];
		new_vert_map[it->second].set_index(it->second);
	}
	vertices.swap(new_vert_map);
	
	// Remap triangle vertices
	for (tit=triangles.begin();tit!=triangles.end();++tit) {
		tit->second.apply_remap(remap);
	}

	// Remap rectangle vertices
	for (rit=rectangles.begin();rit!=rectangles.end();++rit) {
		rit->second.apply_remap(remap);
	}
}

void quakelib::EQSimGeometrySection::erase_vertex(const UIndex &ind) {
	// TODO: WRITE
}

size_t quakelib::EQSimGeometry::num_vertices(void) const {
	EQSimGeomSectionMap::const_iterator it;
	unsigned int		total = 0;
	
	for (it=sections.begin();it!=sections.end();++it) total += it->second.num_vertices();
	
	return total;
}

size_t quakelib::EQSimGeometry::num_triangles(void) const {
	EQSimGeomSectionMap::const_iterator it;
	unsigned int		total = 0;
	
	for (it=sections.begin();it!=sections.end();++it) total += it->second.num_triangles();
	
	return total;
}

size_t quakelib::EQSimGeometry::num_rectangles(void) const {
	EQSimGeomSectionMap::const_iterator it;
	unsigned int		total = 0;
	
	for (it=sections.begin();it!=sections.end();++it) total += it->second.num_rectangles();
	
	return total;
}

double quakelib::EQSimGeometry::lat_hi(void) const {
	EQSimGeomSectionMap::const_iterator	it;
	double		cur_high = -INFINITY;
	for (it=sections.begin();it!=sections.end();++it) cur_high = fmax(cur_high, it->second.lat_hi());
	return cur_high;
}

double quakelib::EQSimGeometrySection::lat_hi(void) const {
	EQSimGeomVertexMap::const_iterator	it;
	double		cur_high = -INFINITY;
	for (it=vertices.begin();it!=vertices.end();++it) cur_high = fmax(cur_high, it->second.loc().lat());
	return cur_high;
}

double quakelib::EQSimGeometry::lat_lo(void) const {
	EQSimGeomSectionMap::const_iterator	it;
	double		cur_low = INFINITY;
	for (it=sections.begin();it!=sections.end();++it) cur_low = fmin(cur_low, it->second.lat_lo());
	return cur_low;
}

double quakelib::EQSimGeometrySection::lat_lo(void) const {
	EQSimGeomVertexMap::const_iterator	it;
	double		cur_low = INFINITY;
	for (it=vertices.begin();it!=vertices.end();++it) cur_low = fmin(cur_low, it->second.loc().lat());
	return cur_low;
}

double quakelib::EQSimGeometry::lon_hi(void) const {
	EQSimGeomSectionMap::const_iterator	it;
	double		cur_high = -INFINITY;
	for (it=sections.begin();it!=sections.end();++it) cur_high = fmax(cur_high, it->second.lon_hi());
	return cur_high;
}

double quakelib::EQSimGeometrySection::lon_hi(void) const {
	EQSimGeomVertexMap::const_iterator	it;
	double		cur_high = -INFINITY;
	for (it=vertices.begin();it!=vertices.end();++it) cur_high = fmax(cur_high, it->second.loc().lon());
	return cur_high;
}

double quakelib::EQSimGeometry::lon_lo(void) const {
	EQSimGeomSectionMap::const_iterator	it;
	double		cur_low = INFINITY;
	for (it=sections.begin();it!=sections.end();++it) cur_low = fmin(cur_low, it->second.lon_lo());
	return cur_low;
}

double quakelib::EQSimGeometrySection::lon_lo(void) const {
	EQSimGeomVertexMap::const_iterator	it;
	double		cur_low = INFINITY;
	for (it=vertices.begin();it!=vertices.end();++it) cur_low = fmin(cur_low, it->second.loc().lon());
	return cur_low;
}

double quakelib::EQSimGeometry::depth_hi(void) const {
	EQSimGeomSectionMap::const_iterator	it;
	double		cur_high = -INFINITY;
	for (it=sections.begin();it!=sections.end();++it) cur_high = fmax(cur_high, it->second.depth_hi());
	return cur_high;
}

double quakelib::EQSimGeometrySection::depth_hi(void) const {
	EQSimGeomVertexMap::const_iterator	it;
	double		cur_high = -INFINITY;
	for (it=vertices.begin();it!=vertices.end();++it) cur_high = fmax(cur_high, it->second.loc().altitude());
	return cur_high;
}

double quakelib::EQSimGeometry::depth_lo(void) const {
	EQSimGeomSectionMap::const_iterator	it;
	double		cur_low = INFINITY;
	for (it=sections.begin();it!=sections.end();++it) cur_low = fmin(cur_low, it->second.depth_lo());
	return cur_low;
}

double quakelib::EQSimGeometrySection::depth_lo(void) const {
	EQSimGeomVertexMap::const_iterator	it;
	double		cur_low = INFINITY;
	for (it=vertices.begin();it!=vertices.end();++it) cur_low = fmin(cur_low, it->second.loc().altitude());
	return cur_low;
}

double quakelib::EQSimGeometry::das_hi(void) const {
	EQSimGeomSectionMap::const_iterator	it;
	double		cur_high = -INFINITY;
	for (it=sections.begin();it!=sections.end();++it) cur_high = fmax(cur_high, it->second.das_hi());
	return cur_high;
}

double quakelib::EQSimGeometrySection::das_hi(void) const {
	EQSimGeomVertexMap::const_iterator	it;
	double		cur_high = -INFINITY;
	for (it=vertices.begin();it!=vertices.end();++it) cur_high = fmax(cur_high, it->second.das());
	return cur_high;
}

double quakelib::EQSimGeometry::das_lo(void) const {
	EQSimGeomSectionMap::const_iterator	it;
	double		cur_low = INFINITY;
	for (it=sections.begin();it!=sections.end();++it) cur_low = fmin(cur_low, it->second.das_lo());
	return cur_low;
}

double quakelib::EQSimGeometrySection::das_lo(void) const {
	EQSimGeomVertexMap::const_iterator	it;
	double		cur_low = INFINITY;
	for (it=vertices.begin();it!=vertices.end();++it) cur_low = fmin(cur_low, it->second.das());
	return cur_low;
}

bool quakelib::EQSimGeometryReader::parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream) {
	switch (rec_num) {
		case 200: parse_fault_summary_record(line_num, line_stream); return true;
		case 201: parse_section_record(line_num, line_stream); return true;
		case 202: parse_vertex_record(line_num, line_stream); return true;
		case 203: parse_triangle_record(line_num, line_stream); return true;
		case 204: parse_rectangle_record(line_num, line_stream); return true;
		default: return false;
	}
}

void quakelib::EQSimGeometryReader::parse_fault_summary_record(const int &line_num, std::istringstream &line_stream) {
	_line_num = line_num;
	line_stream >> _n_section	// Field 1: Total number of fault sections in the file
				>> _n_vertex	// Field 2: Total number of vertices in the file
				>> _n_triangle	// Field 3: Total number of triangles in the file
				>> _n_rectangle	// Field 4: Total number of rectangles in the file
				>> _lat_lo		// Field 5: Lowest value of latitude in the file (decimal degrees, positive north)
				>> _lat_hi		// Field 6: Highest value of latitude in the file (decimal degrees, positive north)
				>> _lon_lo		// Field 7: Lowest value of longitude in the file (decimal degrees, positive east)
				>> _lon_hi		// Field 8: Highest value of longitude in the file (decimal degrees, positive east)
				>> _depth_lo	// Field 9: Lowest value of depth in the file (meters, negative underground)
				>> _depth_hi;	// Field 10: Highest value of depth in the file (meters, negative underground)
}

void quakelib::EQSimParsedGeometrySection::parse(const int &line_num, std::istringstream &line_stream) {
	_line_num = line_num;
	line_stream >> _sid			// Field 1: Section identification number (positive integer, may not be consecutive)
				>> _name		// Field 2: Section name
				>> _n_vertex	// Field 3: Total number of vertices in the section
				>> _n_triangle	// Field 4: Total number of triangles in the section
				>> _n_rectangle	// Field 5: Total number of rectangles in the section
				>> _lat_lo		// Field 6: Lowest value of latitude in the section (decimal degrees, positive north)
				>> _lat_hi		// Field 7: Highest value of latitude in the section (decimal degrees, positive north)
				>> _lon_lo		// Field 8: Lowest value of longitude in the section (decimal degrees, positive east)
				>> _lon_hi		// Field 9: Highest value of longitude in the section (decimal degrees, positive east)
				>> _depth_lo	// Field 10: Lowest value of depth in the section (meters, negative underground)
				>> _depth_hi	// Field 11: Highest value of depth in the section (meters, negative underground)
				>> _das_lo		// Field 12: Lowest value of distance-along-strike in the section (meters)
				>> _das_hi		// Field 13: Highest value of distance-along-strike in the section (meters)
				>> _fid;		// Field 14: Fault identification number (positive integer)
}

void quakelib::EQSimGeometryReader::parse_section_record(const int &line_num, std::istringstream &line_stream) {
	EQSimParsedGeometrySection	pgs;
	EQSimGeometrySection		gs;
	
	pgs.parse(line_num, line_stream);
	gs.set_sid(pgs._sid);
	gs.set_name(pgs._name);
	gs.set_fid(pgs._fid);
	
	cur_section = gs.sid();
	
	_parsed_sections.insert(std::make_pair(cur_section, pgs));
	sections.insert(std::make_pair(cur_section, gs));
}

void quakelib::EQSimGeometryVertex::parse(const int &line_num, std::istringstream &line_stream) {
	double		new_lat, new_lon, new_depth;
	int			trace_flag;
	
	_line_num = line_num;
	line_stream >> _index		// Field 1: Vertex index number
				>> new_lat		// Field 2: Latitude (decimal degrees, positive north)
				>> new_lon		// Field 3: Longitude (decimal degrees, positive east)
				>> new_depth	// Field 4: Depth (meters, negative underground)
				>> _das			// Field 5: Distance-along-strike (meters)
				>> trace_flag;	// Field 6: Trace flag (0 = not on trace, 1 = on trace but not initial or final, 2 = initial point on trace, 3 = final point on trace)
	
	switch (trace_flag) {
		case 0: this->set_trace_flag(NOT_ON_TRACE); break;
		case 1: this->set_trace_flag(MIDDLE_TRACE); break;
		case 2: this->set_trace_flag(BEGINNING_TRACE); break;
		case 3: this->set_trace_flag(END_TRACE); break;
		default: this->set_trace_flag(UNDEFINED_TRACE_STATUS); break;
	}
	this->set_loc(LatLonDepth(new_lat, new_lon, new_depth));
}

void quakelib::EQSimGeometryReader::parse_vertex_record(const int &line_num, std::istringstream &line_stream) {
	quakelib::EQSimGeometryVertex		gv;
	
	gv.parse(line_num, line_stream);
	
	sections[cur_section].vertices.insert(std::make_pair(gv.index(), gv));
}

void quakelib::EQSimGeometryTriangle::parse(const int &line_num, std::istringstream &line_stream) {
	_line_num = line_num;
	line_stream >> _index			// Field 1: Vertex index number
				>> _vertex[0]		// Field 2: Vertex index number for corner #1 (counting counterclockwise as viewed from positive side of element)
				>> _vertex[1]		// Field 3: Vertex index number for corner #2 (counting counterclockwise as viewed from positive side of element)
				>> _vertex[2]		// Field 4: Vertex index number for corner #3 (counting counterclockwise as viewed from positive side of element)
				>> _rake			// Field 5: Rake angle (decimal degrees, 0.0 = left lateral, 90.0 = positive side moves up)
				>> _slip_rate		// Field 6: Element slip rate (meters/second)
				>> _aseis_factor;	// Field 7: Element aseismicity factor
}

void quakelib::EQSimGeometryReader::parse_triangle_record(const int &line_num, std::istringstream &line_stream) {
	quakelib::EQSimGeometryTriangle		gt;
	
	gt.parse(line_num, line_stream);
	
	sections[cur_section].triangles.insert(std::make_pair(gt.index(), gt));
}

void quakelib::EQSimGeometryRectangle::parse(const int &line_num, std::istringstream &line_stream) {
	_line_num = line_num;
	
	line_stream >> _index		// Field 1: Rectangle index number
				>> _vertex[0]	// Field 2: Vertex index number for corner #1 (counting counterclockwise as viewed from positive side of element)
				>> _vertex[1]	// Field 3: Vertex index number for corner #2 (counting counterclockwise as viewed from positive side of element)
				>> _vertex[2]	// Field 4: Vertex index number for corner #3 (counting counterclockwise as viewed from positive side of element)
				>> _vertex[3]	// Field 5: Vertex index number for corner #4 (counting counterclockwise as viewed from positive side of element)
				>> _rake		// Field 6: Rake angle (decimal degrees, 0.0 = left lateral, 90.0 = positive side moves up)
				>> _slip_rate	// Field 7: Element slip rate (meters/second)
				>> _aseis_factor// Field 8: Element aseismicity factor
				>> _strike		// Field 9: Strike angle (decimal degrees)
				>> _dip			// Field 10: Dip angle (decimal degrees)
				>> _perfect_flag;// Field 11: Perfect flag (0 = not perfect rectangle, 1 = perfect rectangle)
}

void quakelib::EQSimGeometryReader::parse_rectangle_record(const int &line_num, std::istringstream &line_stream) {
	quakelib::EQSimGeometryRectangle		gr;
	
	gr.parse(line_num, line_stream);
	
	sections[cur_section].rectangles.insert(std::make_pair(gr.index(), gr));
}

void quakelib::EQSimGeometryVertex::validate(EQSimErrors &errors) const {
	// Check that the trace flag is from 0 to 3
	if (_trace_flag < 0 || _trace_flag > 3) {
		std::stringstream		error_msg;
		error_msg << "Trace flag of vertex (" << _index
		<< ") must be 0, 1, 2 or 3 (currently " << _trace_flag << ").";
		errors.report(error_msg.str(), _line_num);
	}
}

void quakelib::EQSimGeometryTriangle::validate(EQSimErrors &errors) const {
	// TODO: assert vertices are in correct order
	
	// Check that the rectangle aseismicity factor is in acceptable bounds
	if (_aseis_factor < 0 || _aseis_factor > 1) {
		std::stringstream		error_msg;
		error_msg << "Aseismicity factor " << _aseis_factor << " of triangle (" << _index
		<< ") is outside acceptable bounds.";
		errors.report(error_msg.str(), _line_num);
	}
}

void quakelib::EQSimGeometryRectangle::validate(EQSimErrors &errors) const {
	// TODO: assert vertices are in correct order
	
	// Check that the rectangle aseismicity factor is in acceptable bounds
	if (_aseis_factor < 0 || _aseis_factor > 1) {
		std::stringstream		error_msg;
		error_msg << "Aseismicity factor " << _aseis_factor << " of rectangle (" << _index
		<< ") is outside acceptable bounds.";
		errors.report(error_msg.str(), _line_num);
	}
	
	// Check that the rectangle perfect flag is either 0 or 1
	if (_perfect_flag != 0 && _perfect_flag != 1) {
		std::stringstream		error_msg;
		error_msg << "Perfect flag of rectangle (" << _index
		<< ") must be either 0 or 1 (currently " << _perfect_flag << ").";
		errors.report(error_msg.str(), _line_num);
	}
}

void quakelib::EQSimGeometrySection::validate(EQSimErrors &errors) const {
	EQSimGeomTriangleMap::const_iterator		tit;
	EQSimGeomRectangleMap::const_iterator		rit;
	EQSimGeomVertexMap::const_iterator			vit;
	double									das1, das2, das3, das4, dep1, dep2, dep3, dep4;
	double									min_das, max_das, min_dep, max_dep;
	
	for (vit=vertices.begin();vit!=vertices.end();++vit) {
		// Perform internal vertex correctness check
		vit->second.validate(errors);
		
		// Check that vertex index matches the map index
		if (vit->second.index() != vit->first) {
			std::stringstream		error_msg;
			error_msg << "Index mismatch for vertex " << vit->first << " (internal index is " << vit->second.index() << ").";
			errors.report(error_msg.str(), vit->second.line_num());
		}
	}
	
	for (tit=triangles.begin();tit!=triangles.end();++tit) {
		// Perform internal triangle correctness check
		tit->second.validate(errors);
		
		// Check that triangle index matches the map index
		if (tit->second.index() != tit->first) {
			std::stringstream		error_msg;
			error_msg << "Index mismatch for triangle " << tit->first << " (internal index is " << tit->second.index() << ").";
			errors.report(error_msg.str(), tit->second.line_num());
		}
		
		// Check that the vertices used by triangles are within their own sections
		// Check vertex 1
		if (vertices.count(tit->second.vertex(0)) == 0) {
			std::stringstream		error_msg;
			error_msg << "Vertex 1 of triangle " << tit->first << " not found in section " << _sid << ".";
			errors.report(error_msg.str(), tit->second.line_num());
		}
		// Check vertex 2
		if (vertices.count(tit->second.vertex(1)) == 0) {
			std::stringstream		error_msg;
			error_msg << "Vertex 2 of triangle " << tit->first << " not found in section " << _sid << ".";
			errors.report(error_msg.str(), tit->second.line_num());
		}
		// Check vertex 3
		if (vertices.count(tit->second.vertex(2)) == 0) {
			std::stringstream		error_msg;
			error_msg << "Vertex 3 of triangle " << tit->first << " not found in section " << _sid << ".";
			errors.report(error_msg.str(), tit->second.line_num());
		}
	}
	
	for (rit=rectangles.begin();rit!=rectangles.end();++rit) {
		// Perform internal rectangle correctness check
		rit->second.validate(errors);
		
		// Check that rectangle index matches the map index
		if (rit->second.index() != rit->first) {
			std::stringstream		error_msg;
			error_msg << "Index mismatch for rectangle " << rit->first << " (internal index is " << rit->second.index() << ").";
			errors.report(error_msg.str(), rit->second.line_num());
		}
		
		// Check that the vertices used by rectangles are within their own sections
		// Check vertex 1
		if (vertices.count(rit->second.vertex(0)) == 0) {
			std::stringstream		error_msg;
			error_msg << "Vertex 1 of rectangle " << rit->first << " not found in section " << _sid << ".";
			errors.report(error_msg.str(), rit->second.line_num());
		}
		// Check vertex 2
		if (vertices.count(rit->second.vertex(1)) == 0) {
			std::stringstream		error_msg;
			error_msg << "Vertex 2 of rectangle " << rit->first << " not found in section " << _sid << ".";
			errors.report(error_msg.str(), rit->second.line_num());
		}
		// Check vertex 3
		if (vertices.count(rit->second.vertex(2)) == 0) {
			std::stringstream		error_msg;
			error_msg << "Vertex 3 of rectangle " << rit->first << " not found in section " << _sid << ".";
			errors.report(error_msg.str(), rit->second.line_num());
		}
		// Check vertex 4
		if (vertices.count(rit->second.vertex(3)) == 0) {
			std::stringstream		error_msg;
			error_msg << "Vertex 4 of rectangle " << rit->first << " not found in section " << _sid << ".";
			errors.report(error_msg.str(), rit->second.line_num());
		}
		
		// Check that fault is rectangular in DAS/depth coordinates
		das1 = vertices.find(rit->second.vertex(0))->second.das();
		das2 = vertices.find(rit->second.vertex(1))->second.das();
		das3 = vertices.find(rit->second.vertex(2))->second.das();
		das4 = vertices.find(rit->second.vertex(3))->second.das();
		dep1 = vertices.find(rit->second.vertex(0))->second.loc().altitude();
		dep2 = vertices.find(rit->second.vertex(1))->second.loc().altitude();
		dep3 = vertices.find(rit->second.vertex(2))->second.loc().altitude();
		dep4 = vertices.find(rit->second.vertex(3))->second.loc().altitude();
		
		min_das = fmin(das1, fmin(das2, fmin(das3, das4)));
		max_das = fmax(das1, fmax(das2, fmax(das3, das4)));
		min_dep = fmin(dep1, fmin(dep2, fmin(dep3, dep4)));
		max_dep = fmax(dep1, fmax(dep2, fmax(dep3, dep4)));
		
		if ((das1 != min_das && das1 != max_das) || (das2 != min_das && das2 != max_das) ||
			(das3 != min_das && das3 != max_das) || (das4 != min_das && das4 != max_das)) {
			std::stringstream		error_msg;
			error_msg << "Rectangle " << rit->first << " is not rectangular in DAS.";
			errors.report(error_msg.str(), rit->second.line_num());
		}
		if ((dep1 != min_dep && dep1 != max_dep) || (dep2 != min_dep && dep2 != max_dep) ||
			(dep3 != min_dep && dep3 != max_dep) || (dep4 != min_dep && dep4 != max_dep)) {
			std::stringstream		error_msg;
			error_msg << "Rectangle " << rit->first << " is not rectangular in depth.";
			errors.report(error_msg.str(), rit->second.line_num());
		}
	}
}

// TODO: add check for consecutively increasing indices
void quakelib::EQSimGeometry::validate(EQSimErrors &errors) const {
	EQSimGeomSectionMap::const_iterator				sit;
	EQSimParsedGeomSectionMap::const_iterator		pit;
	EQSimGeomVertexMap::const_iterator				vit;
	
	// Assert correctness of each parsed section
	for (pit=_parsed_sections.begin();pit!=_parsed_sections.end();++pit) {
		sit = sections.find(pit->second._sid);
		
		// Check that the section has the correct number of vertices
		if (sit->second.num_vertices() != pit->second._n_vertex) {
			std::stringstream		error_msg;
			error_msg << "Specified number of vertices (" << pit->second._n_vertex << ") for section "
			<< pit->second._sid << " does not match actual count in file (" << sit->second.num_vertices() << ").";
			errors.report(error_msg.str(), pit->second._line_num);
		}
		
		// Check that the section has the correct number of triangles
		if (sit->second.num_triangles() != pit->second._n_triangle) {
			std::stringstream		error_msg;
			error_msg << "Specified number of triangles (" << pit->second._n_triangle << ") for section "
			<< pit->second._sid << " does not match actual count in file (" << sit->second.num_triangles() << ").";
			errors.report(error_msg.str(), pit->second._line_num);
		}
		
		// Check that the section has the correct number of rectangles
		if (sit->second.num_rectangles() != pit->second._n_rectangle) {
			std::stringstream		error_msg;
			error_msg << "Specified number of rectangles (" << pit->second._n_rectangle << ") for section "
			<< pit->second._sid << " does not match actual count in file (" << sit->second.num_rectangles() << ").";
			errors.report(error_msg.str(), pit->second._line_num);
		}
		
		// Check that the bounds of each vertex is within the section bounds
		for (vit=sit->second.vertices.begin();vit!=sit->second.vertices.end();++vit) {
			// Check that the vertex latitude is within the section bounds
			if (vit->second.loc().lat() < pit->second._lat_lo || vit->second.loc().lat() > pit->second._lat_hi) {
				std::stringstream		error_msg;
				error_msg << "Vertex " << vit->first << " latitude (" << vit->second.loc().lat() << ") is outside section bounds.";
				errors.report(error_msg.str(), vit->second.line_num());
			}
			
			// Check that the vertex longitude is within the section bounds
			if (vit->second.loc().lon() < pit->second._lon_lo || vit->second.loc().lon() > pit->second._lon_hi) {
				std::stringstream		error_msg;
				error_msg << "Vertex " << vit->first << " longitude (" << vit->second.loc().lon() << ") is outside section bounds.";
				errors.report(error_msg.str(), vit->second.line_num());
			}
			
			// Check that the vertex depth is within the section bounds
			if (vit->second.loc().altitude() < pit->second._depth_lo || vit->second.loc().altitude() > pit->second._depth_hi) {
				std::stringstream		error_msg;
				error_msg << "Vertex " << vit->first << " depth (" << vit->second.loc().altitude() << ") is outside section bounds.";
				errors.report(error_msg.str(), vit->second.line_num());
			}
			
			// Check that the vertex DAS is within the section bounds
			if (vit->second.das() < pit->second._das_lo || vit->second.das() > pit->second._das_hi) {
				std::stringstream		error_msg;
				error_msg << "Vertex " << vit->first << " DAS (" << vit->second.das() << ") is outside section bounds.";
				errors.report(error_msg.str(), vit->second.line_num());
			}
		}
	}
	
	// Assert correctness of each section
	for (sit=sections.begin();sit!=sections.end();++sit) {
		sit->second.validate(errors);
	}
}

void quakelib::EQSimGeometryReader::validate(EQSimErrors &errors) const {
	quakelib::EQSimMetadataReader::validate(errors);
	
	// Check that the number of sections was the same as specified
	if (sections.size() != _n_section) {
		std::stringstream		error_msg;
		error_msg << "Specified number of fault sections (" << _n_section
		<< ") does not match actual count in file (" << sections.size() << ").";
		errors.report(error_msg.str(), _line_num);
	}
	
	// Check that the number of vertices was the same as specified
	if (num_vertices() != _n_vertex) {
		std::stringstream		error_msg;
		error_msg << "Specified number of vertices (" << _n_vertex
		<< ") does not match actual count in file (" << num_vertices() << ").";
		errors.report(error_msg.str(), _line_num);
	}
	
	// Check that the number of triangles was the same as specified
	if (num_triangles() != _n_triangle) {
		std::stringstream		error_msg;
		error_msg << "Specified number of triangles (" << _n_triangle
		<< ") does not match actual count in file (" << num_triangles() << ").";
		errors.report(error_msg.str(), _line_num);
	}
	
	// Check that the number of triangles was the same as specified
	if (num_rectangles() != _n_rectangle) {
		std::stringstream		error_msg;
		error_msg << "Specified number of rectangles (" << _n_rectangle
		<< ") does not match actual count in file (" << num_rectangles() << ").";
		errors.report(error_msg.str(), _line_num);
	}
	
	// Check that the latitude, longitude and depths are within the specified limits
	if (lat_hi() > _lat_hi || lat_lo() < _lat_lo) {
		std::stringstream		error_msg;
		error_msg << "Vertex latitude limits in file (" << lat_lo() << "," << lat_hi()
		<< ") are outside specified latitude limits (" << _lat_lo << "," << _lat_hi << ").";
		errors.report(error_msg.str(), _line_num);
	}
	if (lon_hi() > _lon_hi || lon_lo() < _lon_lo) {
		std::stringstream		error_msg;
		error_msg << "Vertex longitude limits in file (" << lon_lo() << "," << lon_hi()
		<< ") are outside specified longitude limits (" << _lon_lo << "," << _lon_hi << ").";
		errors.report(error_msg.str(), _line_num);
	}
	if (depth_hi() > _depth_hi || depth_lo() < _depth_lo) {
		std::stringstream		error_msg;
		error_msg << "Vertex depth limits in file (" << depth_lo() << "," << depth_hi()
		<< ") are outside specified depth limits (" << _depth_lo << "," << _depth_hi << ").";
		errors.report(error_msg.str(), _line_num);
	}
	
	// TODO: assert section parsed bounds are within model parse bounds
	
	EQSimGeometry::validate(errors);
}

quakelib::EQSimGeometryWriter::EQSimGeometryWriter(void) : EQSimGeometry(), EQSimMetadataWriter() {
	internal::RecordDesc	geom_summary_rec, geom_section_rec, geom_vertex_rec, geom_triangle_rec, geom_rectangle_rec;
	
	set_spec_level(2);
	meta_add_record(META_SIGNATURE, "EQSim_Input_Geometry_2");
	
	geom_summary_rec = internal::RecordDesc(0, "summary", 10, "Record 200: Fault system summary");
	geom_summary_rec.add_field(1, internal::FieldDesc("n_section", 1, 1, "Total number of fault sections in the file"));
	geom_summary_rec.add_field(2, internal::FieldDesc("n_vertex", 1, 2, "Total number of vertices in the file"));
	geom_summary_rec.add_field(3, internal::FieldDesc("n_triangle", 1, 3, "Total number of triangles in the file"));
	geom_summary_rec.add_field(4, internal::FieldDesc("n_rectangle", 1, 4, "Total number of rectangles in the file"));
	geom_summary_rec.add_field(5, internal::FieldDesc("lat_lo", 2, 5, "Lowest value of latitude in the file (decimal degrees, positive north)"));
	geom_summary_rec.add_field(6, internal::FieldDesc("lat_hi", 2, 6, "Highest value of latitude in the file (decimal degrees, positive north)"));
	geom_summary_rec.add_field(7, internal::FieldDesc("lon_lo", 2, 7, "Lowest value of longitude in the file (decimal degrees, positive east)"));
	geom_summary_rec.add_field(8, internal::FieldDesc("lon_hi", 2, 8, "Highest value of longitude in the file (decimal degrees, positive east)"));
	geom_summary_rec.add_field(9, internal::FieldDesc("depth_lo", 2, 9, "Lowest value of depth in the file (meters, negative underground)"));
	geom_summary_rec.add_field(10, internal::FieldDesc("depth_hi", 2, 10, "Highest value of depth in the file (meters, negative underground)"));
	
	geom_section_rec = internal::RecordDesc(0, "section", 14, "Record 201: Fault section information");
	geom_section_rec.add_field(1, internal::FieldDesc("sid", 1, 1, "Section identification number (positive integer, may not be consecutive)"));
	geom_section_rec.add_field(2, internal::FieldDesc("name", 3, 2, "Section name"));
	geom_section_rec.add_field(3, internal::FieldDesc("n_vertex", 1, 3, "Total number of vertices in the section"));
	geom_section_rec.add_field(4, internal::FieldDesc("n_triangle", 1, 4, "Total number of triangles in the section"));
	geom_section_rec.add_field(5, internal::FieldDesc("n_rectangle", 1, 5, "Total number of rectangles in the section"));
	geom_section_rec.add_field(6, internal::FieldDesc("lat_lo", 2, 6, "Lowest value of latitude in the section (decimal degrees, positive north)"));
	geom_section_rec.add_field(7, internal::FieldDesc("lat_hi", 2, 7, "Highest value of latitude in the section (decimal degrees, positive north)"));
	geom_section_rec.add_field(8, internal::FieldDesc("lon_lo", 2, 8, "Lowest value of longitude in the section (decimal degrees, positive east)"));
	geom_section_rec.add_field(9, internal::FieldDesc("lon_hi", 2, 9, "Highest value of longitude in the section (decimal degrees, positive east)"));
	geom_section_rec.add_field(10, internal::FieldDesc("depth_lo", 2, 10, "Lowest value of depth in the section (meters, negative underground)"));
	geom_section_rec.add_field(11, internal::FieldDesc("depth_hi", 2, 11, "Highest value of depth in the section (meters, negative underground)"));
	geom_section_rec.add_field(12, internal::FieldDesc("das_lo", 2, 12, "Lowest value of distance-along-strike in the section (meters)"));
	geom_section_rec.add_field(13, internal::FieldDesc("das_hi", 2, 13, "Highest value of distance-along-strike in the section (meters)"));
	geom_section_rec.add_field(14, internal::FieldDesc("fault_id", 1, 14, "Fault identification number (positive integer)"));
	
	geom_vertex_rec = internal::RecordDesc(0, "vertex", 6, "Record 202: Vertex");
	geom_vertex_rec.add_field(1, internal::FieldDesc("index", 1, 1, "Vertex index number (consecutive integers, starting with 1)"));
	geom_vertex_rec.add_field(2, internal::FieldDesc("lat", 2, 2, "Latitude (decimal degrees, positive north)"));
	geom_vertex_rec.add_field(3, internal::FieldDesc("lon", 2, 3, "Longitude (decimal degrees, positive east)"));
	geom_vertex_rec.add_field(4, internal::FieldDesc("depth", 2, 4, "Depth (meters, negative underground)"));
	geom_vertex_rec.add_field(5, internal::FieldDesc("das", 2, 5, "Distance-along-strike (meters)"));
	geom_vertex_rec.add_field(6, internal::FieldDesc("trace_flag", 1, 6, "Trace flag (0 = not on trace, 1 = on trace but not initial or final, 2 = initial point on trace, 3 = final point on trace)"));
	
	geom_triangle_rec = internal::RecordDesc(0, "triangle", 9, "Record 203: Triangle");
	geom_triangle_rec.add_field(1, internal::FieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
	geom_triangle_rec.add_field(2, internal::FieldDesc("vertex_1", 1, 2, "Vertex index number for corner #1 (counting counterclockwise as viewed from positive side of element)"));
	geom_triangle_rec.add_field(3, internal::FieldDesc("vertex_2", 1, 3, "Vertex index number for corner #2 (counting counterclockwise as viewed from positive side of element)"));
	geom_triangle_rec.add_field(4, internal::FieldDesc("vertex_3", 1, 4, "Vertex index number for corner #3 (counting counterclockwise as viewed from positive side of element)"));
	geom_triangle_rec.add_field(5, internal::FieldDesc("rake", 2, 5, "Rake angle (decimal degrees)"));
	geom_triangle_rec.add_field(6, internal::FieldDesc("slip_rate", 2, 6, "Element slip rate (meters/second)"));
	geom_triangle_rec.add_field(7, internal::FieldDesc("aseis_factor", 2, 7, "Element aseismicity factor"));
	geom_triangle_rec.add_field(8, internal::FieldDesc("strike", 2, 8, "Strike angle (decimal degrees)"));
	geom_triangle_rec.add_field(9, internal::FieldDesc("dip", 2, 9, "Dip angle (decimal degrees)"));
	
	geom_rectangle_rec = internal::RecordDesc(0, "rectangle", 11, "Record 204: Rectangle");
	geom_rectangle_rec.add_field(1, internal::FieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
	geom_rectangle_rec.add_field(2, internal::FieldDesc("vertex_1", 1, 2, "Vertex index number for corner #1 (counting counterclockwise as viewed from positive side of element)"));
	geom_rectangle_rec.add_field(3, internal::FieldDesc("vertex_2", 1, 3, "Vertex index number for corner #2 (counting counterclockwise as viewed from positive side of element)"));
	geom_rectangle_rec.add_field(4, internal::FieldDesc("vertex_3", 1, 4, "Vertex index number for corner #3 (counting counterclockwise as viewed from positive side of element)"));
	geom_rectangle_rec.add_field(5, internal::FieldDesc("vertex_4", 1, 5, "Vertex index number for corner #4 (counting counterclockwise as viewed from positive side of element)"));
	geom_rectangle_rec.add_field(6, internal::FieldDesc("rake", 2, 6, "Rake angle (decimal degrees)"));
	geom_rectangle_rec.add_field(7, internal::FieldDesc("slip_rate", 2, 7, "Element slip rate (meters/second)"));
	geom_rectangle_rec.add_field(8, internal::FieldDesc("aseis_factor", 2, 8, "Element aseismicity factor"));
	geom_rectangle_rec.add_field(9, internal::FieldDesc("strike", 2, 9, "Strike angle (decimal degrees)"));
	geom_rectangle_rec.add_field(10, internal::FieldDesc("dip", 2, 10, "Dip angle (decimal degrees)"));
	geom_rectangle_rec.add_field(11, internal::FieldDesc("perfect_flag", 1, 11, "Perfect flag (0 = not perfect rectangle, 1 = perfect rectangle)"));
	
	add_record_desc_record(200, geom_summary_rec);
	add_record_desc_record(201, geom_section_rec);
	add_record_desc_record(202, geom_vertex_rec);
	add_record_desc_record(203, geom_triangle_rec);
	add_record_desc_record(204, geom_rectangle_rec);
}
