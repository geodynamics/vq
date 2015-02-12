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

#include "QuakeLibIO.h"
#include "QuakeLibEQSim.h"

quakelib::ModelSection &quakelib::ModelWorld::section(const UIndex &ind) throw(std::domain_error) {
    std::map<UIndex, ModelSection>::iterator it = _sections.find(ind);

    if (it == _sections.end()) throw std::domain_error("quakelib::ModelWorld::section");
    else return it->second;
}

quakelib::ModelElement &quakelib::ModelWorld::element(const UIndex &ind) throw(std::domain_error) {
    std::map<UIndex, ModelElement>::iterator it = _elements.find(ind);

    if (it == _elements.end()) throw std::domain_error("quakelib::ModelWorld::element");
    else return it->second;
}

quakelib::ModelVertex &quakelib::ModelWorld::vertex(const UIndex &ind) throw(std::domain_error) {
    std::map<UIndex, ModelVertex>::iterator it = _vertices.find(ind);

    if (it == _vertices.end()) throw std::domain_error("quakelib::ModelWorld::vertex");
    else return it->second;
}

std::string quakelib::ModelIO::next_line(std::istream &in_stream) {
    std::string line = "";
    size_t      pos;

    do {
        std::getline(in_stream, line);
        _comment = "";
        // Cut off any initial whitespace
        pos = line.find_first_not_of(" \t");

        if (pos != std::string::npos) line = line.substr(pos, std::string::npos);

        // Comment consists of hash mark until the end of the line
        pos = line.find("#");

        if (pos != std::string::npos) _comment = line.substr(pos, std::string::npos);

        // Extract the non-comment part of the line
        line = line.substr(0, line.find("#"));

        // If the line is empty, we keep going
        if (line.length() > 0) break;
    } while (in_stream && !in_stream.eof());

    return line;
}

void quakelib::ModelIO::next_line(std::ostream &out_stream) const {
    if (!_comment.empty()) out_stream << " # " << _comment;

    out_stream << "\n";
}

void quakelib::ModelSection::get_field_descs(std::vector<FieldDesc> &descs) {
    FieldDesc       field_desc;
#ifdef HDF5_FOUND
    hid_t           section_name_datatype;

    // Create the datatype for the section name strings
    section_name_datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(section_name_datatype, (size_t)NAME_MAX_LEN);
#endif

    field_desc.name = "id";
    field_desc.details = "Unique ID of the section.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SectionData, _id);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(UIndex);
#endif
    descs.push_back(field_desc);

    field_desc.name = "fault_id";
    field_desc.details = "ID of the parent fault.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SectionData, _fault_id);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(UIndex);
#endif
    descs.push_back(field_desc);

    field_desc.name = "name";
    field_desc.details = "Name of the fault.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SectionData, _name);
    field_desc.type = section_name_datatype;
    field_desc.size = sizeof(char)*NAME_MAX_LEN;
#endif
    descs.push_back(field_desc);

    // TODO: handle this release properly
    //H5Tclose(section_name_datatype);
}

void quakelib::ModelSection::read_data(const SectionData &in_data) {
    memcpy(&_data, &in_data, sizeof(SectionData));
}

void quakelib::ModelSection::write_data(SectionData &out_data) const {
    memcpy(&out_data, &_data, sizeof(SectionData));
}

void quakelib::ModelSection::read_ascii(std::istream &in_stream) {
    std::stringstream   ss(next_line(in_stream));

    ss >> _data._id;
    ss >> _data._fault_id;
    ss >> _data._name;
}

void quakelib::ModelSection::write_ascii(std::ostream &out_stream) const {
    out_stream << _data._id << " ";
    out_stream << _data._fault_id << " ";
    out_stream << _data._name;

    next_line(out_stream);
}

void quakelib::ModelElement::get_field_descs(std::vector<FieldDesc> &descs) {
    FieldDesc       field_desc;

    field_desc.name = "id";
    field_desc.details = "Unique ID of the element.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _id);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(UIndex);
#endif
    descs.push_back(field_desc);

    field_desc.name = "section_id";
    field_desc.details = "ID of the section associated with the element.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _section_id);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(UIndex);
#endif
    descs.push_back(field_desc);

    field_desc.name = "vertex_0";
    field_desc.details = "ID of vertex 0.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _vertices[0]);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(UIndex);
#endif
    descs.push_back(field_desc);

    field_desc.name = "vertex_1";
    field_desc.details = "ID of vertex 1.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _vertices[1]);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(UIndex);
#endif
    descs.push_back(field_desc);

    field_desc.name = "vertex_2";
    field_desc.details = "ID of vertex 2.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _vertices[2]);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(UIndex);
#endif
    descs.push_back(field_desc);

    field_desc.name = "is_quad";
    field_desc.details = "Whether the vertices constitute 3 points of a triangle (zero) or 3 points of a parallelogram (non-zero).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _is_quad);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "slip_rate";
    field_desc.details = "Long term slip rate of element in meters per second.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _slip_rate);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "aseismic";
    field_desc.details = "Fraction of slip on element that is aseismic.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _aseismic);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "rake";
    field_desc.details = "Rake angle of element in radians.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _rake);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "lame_mu";
    field_desc.details = "Lame's parameter describing the shear modulus of the material for this element (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _lame_mu);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "lame_lambda";
    field_desc.details = "Lame's lambda parameter of the material for this element, in Pascals.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _lame_lambda);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "max_slip";
    field_desc.details = "Maximum slip distance for this element, in meters.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(ElementData, _max_slip);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);
}

void quakelib::ModelElement::read_data(const ElementData &in_data) {
    memcpy(&_data, &in_data, sizeof(ElementData));
}

void quakelib::ModelElement::write_data(ElementData &out_data) const {
    memcpy(&out_data, &_data, sizeof(ElementData));
}

void quakelib::ModelElement::read_ascii(std::istream &in_stream) {
    unsigned int        i;
    std::stringstream   ss(next_line(in_stream));

    ss >> _data._id;
    ss >> _data._section_id;

    for (i=0; i<3; ++i) ss >> _data._vertices[i];

    ss >> _data._is_quad;
    ss >> _data._slip_rate;
    ss >> _data._aseismic;
    ss >> _data._rake;
    ss >> _data._lame_mu;
    ss >> _data._lame_lambda;
    ss >> _data._max_slip;
}

void quakelib::ModelElement::write_ascii(std::ostream &out_stream) const {
    unsigned int        i;

    out_stream << _data._id << " ";
    out_stream << _data._section_id << " ";

    for (i=0; i<3; ++i) out_stream << _data._vertices[i] << " ";

    out_stream << _data._is_quad << " ";
    out_stream << _data._slip_rate << " ";
    out_stream << _data._aseismic << " ";
    out_stream << _data._rake << " ";
    out_stream << _data._lame_mu << " ";
    out_stream << _data._lame_lambda << " ";
    out_stream << _data._max_slip;

    next_line(out_stream);
}

void quakelib::FaultTracePoint::read_ascii(std::istream &in_stream) {
    float      lat, lon, alt;
    std::stringstream   ss(next_line(in_stream));

    ss >> lat;
    ss >> lon;
    ss >> alt;
    _pos = LatLonDepth(lat, lon, alt);
    ss >> _depth_along_dip;
    ss >> _slip_rate;
    ss >> _aseismic;
    ss >> _rake;
    ss >> _dip;
    ss >> _lame_mu;
    ss >> _lame_lambda;
}

void quakelib::FaultTracePoint::write_ascii(std::ostream &out_stream) const {
    out_stream << _pos.lat() << " ";
    out_stream << _pos.lon() << " ";
    out_stream << _pos.altitude() << " ";
    out_stream << _depth_along_dip << " ";
    out_stream << _slip_rate << " ";
    out_stream << _aseismic << " ";
    out_stream << _rake << " ";
    out_stream << _dip << " ";
    out_stream << _lame_mu << " ";
    out_stream << _lame_lambda;

    next_line(out_stream);
}

/*
 Provides a spline to help mesh faults along trace points. This spline is referenced by t (in [0,1])
 representing how long along the trace a point is. The spline is a linear interpolation between successive points.
 */
class TraceSpline {
    private:
        //! The points comprising the spline
        std::vector<quakelib::Vec<3> >  _pts;
        //! Total length of the distances between successive points
        double                          _spline_len;
        //! Individual lengths between successive points
        //! _point_dists[i] is the distance between _pts[i] and _pts[i+1]
        std::vector<double>             _point_dists;

    public:
        TraceSpline(void) : _spline_len(0) {};

        // Add another point to the spline
        void add_point(const quakelib::Vec<3> &new_pt) {
            double add_dist = (_pts.size() > 0 ? _pts.back().dist(new_pt) : 0);
            _spline_len += add_dist;

            if (_pts.size() > 0) _point_dists.push_back(add_dist);

            _pts.push_back(new_pt);
        }

        // Return the element index and inner t corresponding to parameter t
        void get_element(const double t, unsigned int &index, double &inner_t) {
            double spline_dist = t * _spline_len;

            // Ensure t is in [0,1]
            assertThrow(t >= 0 && t <= 1, std::out_of_range("TraceSpline::interpolate"));

            // Go through each point
            for (unsigned int i=0; i<_pts.size()-1; ++i) {
                // If we're between the points for this t, interpolate and return
                if (spline_dist <= _point_dists[i]) {
                    index = i;

                    if (_point_dists[i] != 0) inner_t = spline_dist / _point_dists[i];
                    else inner_t = 0;

                    return;
                }

                spline_dist -= _point_dists[i];
            }

            // If we reach the end, we return the final point and inner_t = 0
            index = _pts.size()-1;
            inner_t = 0;
        }

        // Return the point on this spline at t
        quakelib::Vec<3> interpolate(const double t) {
            unsigned int    ind;
            double          inner_t;

            assertThrow(t >= 0 && t <= 1, std::out_of_range("TraceSpline::interpolate"));
            get_element(t, ind, inner_t);

            if (ind == _pts.size() - 1) return _pts.back();
            else return _pts.at(ind) * (1-inner_t) + _pts.at(ind+1) * inner_t;
        }

        // Given a starting point t (in [0,1]) and an element size,
        // returns the next t which is linear elem_size away (not distance along the spline)
        double advance_element(const double t, const double elem_size) {
            unsigned int    ind;
            double          inner_t;

            // Find the starting point on the spline at t
            quakelib::Vec<3> start_pt = interpolate(t);

            // Get the element that the starting point is associated iwth
            get_element(t, ind, inner_t);

            // Keep going until we find a trace point which would create
            // an element greater than our target size
            double cur_dist = t * _spline_len;

            while (ind+1 < _pts.size() && start_pt.dist(_pts.at(ind+1)) < elem_size) {
                ind++;

                if (ind < _point_dists.size()) cur_dist += (1-inner_t) * _point_dists.at(ind);
                else cur_dist += elem_size;

                inner_t = 0;
            }

            // If we're past the end of the trace, return our best guess
            // for t based on the size of the last segment
            // This is needed to adjust the element size during meshing
            if (ind+1 == _pts.size()) return cur_dist/_spline_len;

            // Now we know the points between which our element must exist (ind and ind+1)
            double      x_0 = _pts.at(ind)[0], x_1 = _pts.at(ind+1)[0], x_s = start_pt[0];
            double      y_0 = _pts.at(ind)[1], y_1 = _pts.at(ind+1)[1], y_s = start_pt[1], l = elem_size;

            // Calculate the inner t between these two points
            double next_t = fabs((sqrt(pow(-2*x_0*x_0+2*x_0*x_1+2*x_0*x_s-2*x_1*x_s-2*y_0*y_0+2*y_0*y_1+2*y_0*y_s-2*y_1*y_s, 2)-
                                       4*(x_0*x_0-2*x_0*x_1+x_1*x_1+y_0*y_0-2*y_0*y_1+y_1*y_1)*
                                       (x_0*x_0-2*x_0*x_s+x_s*x_s+y_0*y_0-2*y_0*y_s+y_s*y_s-l*l))+
                                  2*x_0*x_0-2*x_0*x_1-2*x_0*x_s+2*x_1*x_s+2*y_0*y_0-2*y_0*y_1-2*y_0*y_s+2*y_1*y_s)/
                                 (2*(x_0*x_0-2*x_0*x_1+x_1*x_1+y_0*y_0-2*y_0*y_1+y_1*y_1)));

            // Given this point, recalculate t and return
            cur_dist += (next_t-inner_t) * _point_dists[ind];
            return cur_dist/_spline_len;
        };
};

void quakelib::ModelWorld::create_section(std::vector<unsigned int> &unused_trace_segments, const std::vector<FaultTracePoint> &trace, const LatLonDepth &base_coord, const UIndex &fault_id, const float &element_size, const std::string &section_name, const std::string &taper_method, const bool resize_trace_elements) {
    Vec<3>              cur_trace_point, next_trace_point, element_end, element_step_vec, vert_step;
    std::vector<UIndex> elem_ids;
    std::set<unsigned int> unused_trace_pts;
    double              elem_depth, elem_slip_rate, elem_aseismic;
    double              elem_rake, elem_dip;
    double              elem_lame_mu, elem_lame_lambda;
    unsigned int        num_vert_elems, ve, elem_count;
    double              taper_t;
    Conversion          conv(base_coord);
    unsigned int        i, num_trace_pts;
    TraceSpline         spline;

    if (element_size <= 0) return;

    num_trace_pts = trace.size();

    if (num_trace_pts == 0 || num_trace_pts == 1) return;

    ModelSection &section = new_section();
    section.set_name(section_name);
    section.set_fault_id(fault_id);

    // Create a spline with the trace points
    for (i=0; i<num_trace_pts; ++i) {
        Vec<3> pt = conv.convert2xyz(trace.at(i).pos());
        spline.add_point(pt);
        unused_trace_pts.insert(i);
    }

    // Initially we don't know the appropriate element size to exactly mesh
    // the trace points. This is solved by starting with the user suggested
    // element size and slowly shrinking it down to half the original element size until the mesh fits the
    // trace exactly.
    double cur_elem_size_guess = element_size;
    double step_size = element_size/(2*1000);
    double best_step = element_size, best_t = 0;
    unsigned int best_elem_count = 0;

    while (cur_elem_size_guess > element_size/2.0) {
        double cur_t, next_t, sum_t=0;
        elem_count = 0;
        cur_t = 0;

        while (cur_t < 1) {
            next_t = spline.advance_element(cur_t, cur_elem_size_guess);

            if (next_t < 1) {
                sum_t += next_t-cur_t;
                elem_count++;
                cur_t = next_t;
            } else {
                break;
            }
        }

        // If we used a fixed element size, one time through is enough
        if (!resize_trace_elements) {
            best_step = cur_elem_size_guess;
            best_elem_count = elem_count;
            break;
        }

        // Record which element size got us closest to the end of the trace
        if (cur_t > best_t) {
            best_t = cur_t;
            best_step = cur_elem_size_guess;
            best_elem_count = elem_count;
        }

        cur_elem_size_guess -= step_size;
    }

    double  horiz_elem_size = best_step;
    double  vert_elem_size;
    unsigned int cur_elem_ind, next_elem_ind;
    double cur_t = 0, next_t;

    double dist_along_strike = 0;
    double taper_flow = 0;
    double taper_full = 0;
    double fault_area = 0;
    cur_trace_point = spline.interpolate(cur_t);
    cur_elem_ind = 0;
    unused_trace_pts.erase(0);

    for (i=0; i<best_elem_count; ++i) {
        // Get the next t value along the trace by advancing the element size
        next_t = spline.advance_element(cur_t, horiz_elem_size);

        // If we go past the end of the spline, align ourselves with the end
        if (next_t > 1) next_t = 1;

        // And get the actual point corresponding to this t
        next_trace_point = spline.interpolate(next_t);

        double          inner_t;
        // And the element and inner_t corresponding to it
        spline.get_element(cur_t, next_elem_ind, inner_t);

        // Mark the element corresponding to this point as having been used
        unused_trace_pts.erase(next_elem_ind);

        element_step_vec = next_trace_point-cur_trace_point;

        // Use the t value between the trace points for interpolation
        elem_depth = inner_t *trace.at(next_elem_ind).depth_along_dip()+(1.0-inner_t)*trace.at(cur_elem_ind).depth_along_dip();

        // Change vertical element size to exactly match the required depth
        num_vert_elems = round(elem_depth/element_size);
        vert_elem_size = elem_depth / num_vert_elems;

        // TODO: change this to an assertion throw
        if (num_vert_elems == 0) std::cerr << "WARNING: Depth is smaller than element size in trace segment "
                                               << next_elem_ind << ". Element size may be too big." << std::endl;

        elem_slip_rate = conv.cm_per_yr2m_per_sec(inner_t *trace.at(next_elem_ind).slip_rate()+(1.0-inner_t)*trace.at(cur_elem_ind).slip_rate());
        elem_aseismic = inner_t *trace.at(next_elem_ind).aseismic()+(1.0-inner_t)*trace.at(cur_elem_ind).aseismic();
        elem_dip = conv.deg2rad(inner_t *trace.at(next_elem_ind).dip()+(1.0-inner_t)*trace.at(cur_elem_ind).dip());
        elem_rake = conv.deg2rad(inner_t *trace.at(next_elem_ind).rake()+(1.0-inner_t)*trace.at(cur_elem_ind).rake());
        elem_lame_mu = inner_t *trace.at(next_elem_ind).lame_mu()+(1.0-inner_t)*trace.at(cur_elem_ind).lame_mu();
        elem_lame_lambda = inner_t *trace.at(next_elem_ind).lame_lambda()+(1.0-inner_t)*trace.at(cur_elem_ind).lame_lambda();

        // Set up the vertical step to go along the dip
        vert_step = element_step_vec.rotate_around_axis(Vec<3>(0,0,-1), M_PI/2);
        vert_step = vert_step.rotate_around_axis(element_step_vec, elem_dip);

        // Create each of the elements along dip
        for (ve=0; ve<num_vert_elems; ++ve) {
            // Calculate values for element based on interpolation between trace points
            taper_t = 1;

            double mid_t = (cur_t+next_t)/2.0;

            if (taper_method == "taper_full" || taper_method == "taper_renorm") {
                double x = mid_t;
                double z = (float(ve)+0.5)/num_vert_elems;
                taper_t *= 4*(x-x*x)*sqrt(1-z);
            } else if (taper_method == "taper") {
                double inside_dist = (0.5 - fabs(0.5-mid_t));

                if (inside_dist <= elem_depth) {
                    double x = inside_dist/elem_depth;
                    double z = (float(ve)+0.5)/num_vert_elems;
                    taper_t *= sqrt(x)*sqrt(1-z);
                }
            }

            taper_flow += taper_t *elem_slip_rate*(horiz_elem_size*vert_elem_size);
            taper_full += elem_slip_rate*(horiz_elem_size*vert_elem_size);

            // Create the new vertices
            ModelVertex &v0 = new_vertex();
            ModelVertex &v1 = new_vertex();
            ModelVertex &v2 = new_vertex();

            // Set xyz for the vertices
            v0.set_xyz(cur_trace_point+vert_step*ve, base_coord);
            v1.set_xyz(cur_trace_point+vert_step*(ve+1), base_coord);
            v2.set_xyz(cur_trace_point+vert_step*ve+element_step_vec, base_coord);

            // Set distance along strike for each vertex
            v0.set_das(dist_along_strike);
            v1.set_das(dist_along_strike);
            v2.set_das(dist_along_strike+horiz_elem_size);

            // Create element and set up with the created vertices and values
            ModelElement &elem = new_element();
            elem_ids.push_back(elem.id());
            elem.set_section_id(section.id());
            elem.set_vertex(0, v0.id());
            elem.set_vertex(1, v1.id());
            elem.set_vertex(2, v2.id());
            elem.set_is_quad(true);
            elem.set_slip_rate(elem_slip_rate*taper_t);
            elem.set_aseismic(elem_aseismic);
            elem.set_rake(elem_rake);
            elem.set_lame_mu(elem_lame_mu);
            elem.set_lame_lambda(elem_lame_lambda);

            fault_area += horiz_elem_size*vert_elem_size;
        }

        cur_t = next_t;
        dist_along_strike += horiz_elem_size;
        cur_trace_point = next_trace_point;
    }

    if (taper_method == "taper_renorm") {
        double renorm_factor = taper_full/taper_flow, cur_slip_rate;

        for (unsigned int i=0; i<elem_ids.size(); ++i) {
            cur_slip_rate = element(elem_ids[i]).slip_rate();
            element(elem_ids[i]).set_slip_rate(renorm_factor*cur_slip_rate);
        }
    }

    // Go through the created elements and assign maximum slip based on the total fault area
    // From Table 2A in Wells Coppersmith 1994
    double moment_magnitude = 4.07+0.98*log10(conv.sqm2sqkm(fault_area));

    for (unsigned int i=0; i<elem_ids.size(); ++i) {
        double max_slip = pow(10, (3.0/2.0)*(moment_magnitude+10.7))/(1e7*element(elem_ids[i]).lame_mu()*fault_area);
        element(elem_ids[i]).set_max_slip(max_slip);
    }
}

quakelib::ModelSection &quakelib::ModelWorld::new_section(void) {
    UIndex  max_ind = next_section_index();
    _sections.insert(std::make_pair(max_ind, ModelSection()));
    _sections.find(max_ind)->second.set_id(max_ind);
    return _sections.find(max_ind)->second;
}

quakelib::ModelElement &quakelib::ModelWorld::new_element(void) {
    UIndex  max_ind = next_element_index();
    _elements.insert(std::make_pair(max_ind, ModelElement()));
    _elements.find(max_ind)->second.set_id(max_ind);
    return _elements.find(max_ind)->second;
}

quakelib::ModelVertex &quakelib::ModelWorld::new_vertex(void) {
    UIndex  max_ind = next_vertex_index();
    _vertices.insert(std::make_pair(max_ind, ModelVertex()));
    _vertices.find(max_ind)->second.set_id(max_ind);
    return _vertices.find(max_ind)->second;
}

void quakelib::ModelVertex::get_field_descs(std::vector<FieldDesc> &descs) {
    FieldDesc       field_desc;

    field_desc.name = "id";
    field_desc.details = "Unique ID of the vertex.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(VertexData, _id);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(UIndex);
#endif
    descs.push_back(field_desc);

    field_desc.name = "latitude";
    field_desc.details = "Latitude of the vertex.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(VertexData, _lat);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "longitude";
    field_desc.details = "Longitude of the vertex.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(VertexData, _lon);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "altitude";
    field_desc.details = "Altitude of the vertex in meters (negative is below ground).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(VertexData, _alt);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "das";
    field_desc.details = "Vertex distance along fault strike in meters.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(VertexData, _das);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);
}

void quakelib::ModelVertex::read_data(const VertexData &in_data) {
    memcpy(&_data, &in_data, sizeof(VertexData));
}

void quakelib::ModelVertex::write_data(VertexData &out_data) const {
    memcpy(&out_data, &_data, sizeof(VertexData));
}

void quakelib::ModelVertex::read_ascii(std::istream &in_stream) {
    std::stringstream   ss(next_line(in_stream));

    ss >> _data._id;
    ss >> _data._lat;
    ss >> _data._lon;
    ss >> _data._alt;
    ss >> _data._das;
}

void quakelib::ModelVertex::write_ascii(std::ostream &out_stream) const {
    out_stream << _data._id << " ";
    out_stream << _data._lat << " ";
    out_stream << _data._lon << " ";
    out_stream << _data._alt << " ";
    out_stream << _data._das;

    next_line(out_stream);
}

void quakelib::ModelWorld::clear(void) {
    _sections.clear();
    _elements.clear();
    _vertices.clear();
}

int quakelib::ModelWorld::read_file_ascii(const std::string &file_name) {
    std::ifstream       in_file;
    unsigned int        i, num_sections, num_elements, num_vertices;
    LatLonDepth         min_latlon, max_latlon;

    // Clear the world first to avoid incorrectly mixing indices
    clear();

    in_file.open(file_name.c_str());

    if (!in_file.is_open()) return -1;

    // Read the first line describing the number of sections, etc
    std::stringstream desc_line(next_line(in_file));
    desc_line >> num_sections;
    desc_line >> num_elements;
    desc_line >> num_vertices;

    // Read sections
    for (i=0; i<num_sections; ++i) {
        ModelSection     new_section;
        new_section.read_ascii(in_file);
        _sections.insert(std::make_pair(new_section.id(), new_section));
    }

    // Read elements
    for (i=0; i<num_elements; ++i) {
        ModelElement     new_elem;
        new_elem.read_ascii(in_file);
        _elements.insert(std::make_pair(new_elem.id(), new_elem));
    }

    // Read vertices
    for (i=0; i<num_vertices; ++i) {
        ModelVertex     new_vert;
        new_vert.read_ascii(in_file);
        _vertices.insert(std::make_pair(new_vert.id(), new_vert));
    }

    in_file.close();

    // Reset the internal Cartesian coordinate system
    get_bounds(min_latlon, max_latlon);
    min_latlon.set_altitude(0);
    reset_base_coord(min_latlon);

    return 0;
}

int quakelib::ModelWorld::read_file_trace_latlon(std::vector<unsigned int> &unused_trace_segments, const std::string &file_name, const float &elem_size, const std::string &taper_method, const bool resize_trace_elements) {
    std::ifstream                   in_file;
    std::vector<FaultTracePoint>    trace_pts;
    std::string                     cur_section_name;
    unsigned int                    i, num_trace_pts;
    UIndex                          fault_id;
    double                          min_lat, min_lon;

    // Clear the world first to avoid incorrectly mixing indices
    clear();

    in_file.open(file_name.c_str());

    if (!in_file.is_open()) return -1;

    while (in_file) {
        std::stringstream           ss(next_line(in_file));
        ModelWorld                  new_world;

        cur_section_name = "";
        num_trace_pts = 0;

        ss >> fault_id;
        ss >> num_trace_pts;
        ss >> cur_section_name;

        if (cur_section_name.empty() || num_trace_pts == 0) break;

        min_lat = min_lon = DBL_MAX;

        for (i=0; i<num_trace_pts; ++i) {
            FaultTracePoint new_trace_pt(in_file);
            trace_pts.push_back(new_trace_pt);
            min_lat = fmin(min_lat, new_trace_pt.pos().lat());
            min_lon = fmin(min_lon, new_trace_pt.pos().lon());
        }

        new_world.create_section(unused_trace_segments, trace_pts, LatLonDepth(min_lat, min_lon), fault_id, elem_size, cur_section_name, taper_method, resize_trace_elements);
        this->insert(new_world);
    }

    in_file.close();

    return 0;
}

void quakelib::ModelWorld::reset_base_coord(const LatLonDepth &new_base) {
    std::map<UIndex, ModelVertex>::iterator         it;

    for (it=_vertices.begin(); it!=_vertices.end(); ++it) {
        it->second.set_lld(it->second.lld(), new_base);
    }
}

// TODO: Currently only supports sections where top element is at the same depth, change to support more complex faults
// Also assumes elements will be in order along the trace
// TODO: add element comments to output
int quakelib::ModelWorld::write_file_trace_latlon(const std::string &file_name, const float &depth_along_dip) {
    std::ofstream       out_file;
    eiterator           eit, last_element;
    siterator           sit;
    UIndex              sid;
    unsigned int        i;
    double              max_alt;
    bool                element_on_trace;
    Conversion          c;

    out_file.open(file_name.c_str());

    // Write traces by section
    for (sit=begin_section(); sit!=end_section(); ++sit) {
        std::vector<FaultTracePoint>    trace_pts;

        trace_pts.clear();
        sid = sit->id();

        // Start by going through all elements
        max_alt = -DBL_MAX;

        for (eit=begin_element(sid); eit!=end_element(sid); ++eit) {
            for (i=0; i<3; ++i) max_alt = fmax(max_alt, vertex(eit->vertex(i)).xyz()[2]);
        }

        // Go through the elements again and grab those which have vertices at the correct depth
        // TODO: interpolate between elements
        for (eit=begin_element(sid); eit!=end_element(sid); ++eit) {
            element_on_trace = (vertex(eit->vertex(0)).xyz()[2] == max_alt);

            // If the element is on the trace, print it out
            if (element_on_trace) {
                Vec<3>      a, b;
                double      dip_angle;
                a = vertex(eit->vertex(1)).xyz() - vertex(eit->vertex(0)).xyz();
                b = vertex(eit->vertex(2)).xyz() - vertex(eit->vertex(0)).xyz();
                dip_angle = a.cross(b).unit_vector().vector_angle(Vec<3>(0,0,1));

                FaultTracePoint trace_pt(vertex(eit->vertex(0)).lld(),
                                         depth_along_dip,
                                         c.m_per_sec2cm_per_yr(eit->slip_rate()),
                                         eit->aseismic(),
                                         c.rad2deg(eit->rake()),
                                         c.rad2deg(dip_angle),
                                         eit->lame_mu(),
                                         eit->lame_lambda());
                trace_pts.push_back(trace_pt);
                last_element = eit;
            }
        }

        Vec<3>      a, b;
        double      dip_angle;
        a = vertex(last_element->vertex(1)).xyz() - vertex(last_element->vertex(0)).xyz();
        b = vertex(last_element->vertex(2)).xyz() - vertex(last_element->vertex(0)).xyz();
        dip_angle = a.cross(b).unit_vector().vector_angle(Vec<3>(0,0,1));
        FaultTracePoint trace_pt(vertex(last_element->vertex(2)).lld(),
                                 depth_along_dip,
                                 c.m_per_sec2cm_per_yr(last_element->slip_rate()),
                                 last_element->aseismic(),
                                 c.rad2deg(last_element->rake()),
                                 c.rad2deg(dip_angle),
                                 last_element->lame_mu(),
                                 last_element->lame_lambda());
        trace_pts.push_back(trace_pt);

        // Write the fault header
        out_file << "# fault_id: ID number of the parent fault of this section\n";
        out_file << "# num_points: Number of trace points comprising this section\n";
        out_file << "# section_name: Name of the section\n";

        // Write out the recorded trace for this fault
        out_file << sit->fault_id() << " " << trace_pts.size() << " " << sit->name() << "\n";

        // Write out the trace point header
        out_file << "# latitude: Latitude of trace point\n";
        out_file << "# longitude: Longitude of trace point\n";
        out_file << "# altitude: Altitude of trace point (meters)\n";
        out_file << "# depth_along_dip: Depth along dip (meters)\n";
        out_file << "# slip_rate: Slip rate at trace point (centimeters/year)\n";
        out_file << "# aseismic: Fraction of slip that is aseismic at point\n";
        out_file << "# rake: Fault rake at trace point (degrees)\n";
        out_file << "# dip: Fault dip at trace point (degrees)\n";
        out_file << "# lame_mu: Lame's mu parameter at trace point (Pascals)\n";
        out_file << "# lame_lambda: Lame's lambda parameter at trace point (Pascals)\n";

        // And each of the trace points
        for (i=0; i<trace_pts.size(); ++i) trace_pts[i].write_ascii(out_file);
    }

    out_file.close();

    return 0;
}

int quakelib::ModelWorld::write_file_ascii(const std::string &file_name) const {
    std::ofstream                                   out_file;
    std::vector<FieldDesc>                          descs;
    std::vector<FieldDesc>::iterator                dit;
    std::map<UIndex, ModelVertex>::const_iterator   vit;
    std::map<UIndex, ModelElement>::const_iterator  eit;
    std::map<UIndex, ModelSection>::const_iterator  fit;

    out_file.open(file_name.c_str());
    out_file << "# Number of sections\n";
    out_file << "# Number of elements\n";
    out_file << "# Number of vertices\n";
    out_file << _sections.size() << " " << _elements.size() << " " << _vertices.size();
    next_line(out_file);

    // Write section header
    descs.clear();
    ModelSection::get_field_descs(descs);

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_file << "# " << dit->name << ": " << dit->details << "\n";
    }

    out_file << "# ";

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_file << dit->name << " ";
    }

    out_file << "\n";

    // Write sections
    for (fit=_sections.begin(); fit!=_sections.end(); ++fit) {
        fit->second.write_ascii(out_file);
    }

    // Write element header
    descs.clear();
    ModelElement::get_field_descs(descs);

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_file << "# " << dit->name << ": " << dit->details << "\n";
    }

    out_file << "# ";

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_file << dit->name << " ";
    }

    out_file << "\n";

    // Write elements
    for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
        eit->second.write_ascii(out_file);
    }

    // Write vertex header
    descs.clear();
    ModelVertex::get_field_descs(descs);

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_file << "# " << dit->name << ": " << dit->details << "\n";
    }

    out_file << "# ";

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_file << dit->name << " ";
    }

    out_file << "\n";

    // Write vertices
    for (vit=_vertices.begin(); vit!=_vertices.end(); ++vit) {
        vit->second.write_ascii(out_file);
    }

    out_file.close();

    return 0;
}

int quakelib::ModelWorld::read_file_hdf5(const std::string &file_name) {
#ifdef HDF5_FOUND
    hid_t       plist_id, data_file;
    herr_t      res;
    LatLonDepth min_latlon, max_latlon;

    // Clear the world first to avoid incorrectly mixing indices
    clear();

    if (!H5Fis_hdf5(file_name.c_str())) return -1;

    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    if (plist_id < 0) exit(-1);

    data_file = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, plist_id);

    if (data_file < 0) exit(-1);

    read_section_hdf5(data_file);
    read_element_hdf5(data_file);
    read_vertex_hdf5(data_file);

    // Release HDF5 handles
    res = H5Pclose(plist_id);

    if (res < 0) exit(-1);

    res = H5Fclose(data_file);

    if (res < 0) exit(-1);

    // Reset the internal Cartesian coordinate system
    get_bounds(min_latlon, max_latlon);
    min_latlon.set_altitude(0);
    reset_base_coord(min_latlon);
#else
    // TODO: Error out
#endif
    return 0;
}

#ifdef HDF5_FOUND
void quakelib::ModelWorld::read_section_hdf5(const hid_t &data_file) {
    std::vector<FieldDesc>                          descs;
    std::map<UIndex, ModelSection>::const_iterator  fit;
    hsize_t                     num_fields, num_sections;
    unsigned int                i;
    SectionData                 *section_data;
    size_t                      *field_offsets;
    size_t                      *field_sizes;
    herr_t                      res;

    descs.clear();
    ModelSection::get_field_descs(descs);
    num_fields = descs.size();
    field_offsets = new size_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_offsets[i] = descs[i].offset;
        field_sizes[i] = descs[i].size;
    }

    res = H5TBget_table_info(data_file, ModelSection::hdf5_table_name().c_str(), &num_fields, &num_sections);

    if (res < 0) exit(-1);

    // TODO: check that num_fields matches the descs

    section_data = new SectionData[num_sections];
    res = H5TBread_records(data_file, ModelSection::hdf5_table_name().c_str(), 0, num_sections, sizeof(SectionData), field_offsets, field_sizes, section_data);

    if (res < 0) exit(-1);

    // Read section data into the World
    for (i=0; i<num_sections; ++i) {
        ModelSection  new_section;
        new_section.read_data(section_data[i]);
        _sections.insert(std::make_pair(new_section.id(), new_section));
    }

    // Free memory for HDF5 related data
    delete section_data;
    delete field_offsets;
    delete field_sizes;
}

void quakelib::ModelWorld::read_element_hdf5(const hid_t &data_file) {
    std::vector<FieldDesc>                          descs;
    std::map<UIndex, ModelElement>::const_iterator  fit;
    hsize_t                     num_fields, num_elements;
    unsigned int                i;
    ElementData                 *element_data;
    size_t                      *field_offsets;
    size_t                      *field_sizes;
    herr_t                      res;

    descs.clear();
    ModelElement::get_field_descs(descs);
    num_fields = descs.size();
    field_offsets = new size_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_offsets[i] = descs[i].offset;
        field_sizes[i] = descs[i].size;
    }

    res = H5TBget_table_info(data_file, ModelElement::hdf5_table_name().c_str(), &num_fields, &num_elements);

    if (res < 0) exit(-1);

    // TODO: check that num_fields matches the descs

    element_data = new ElementData[num_elements];
    res = H5TBread_records(data_file, ModelElement::hdf5_table_name().c_str(), 0, num_elements, sizeof(ElementData), field_offsets, field_sizes, element_data);

    if (res < 0) exit(-1);

    // Read element data into the World
    for (i=0; i<num_elements; ++i) {
        ModelElement  new_element;
        new_element.read_data(element_data[i]);
        _elements.insert(std::make_pair(new_element.id(), new_element));
    }

    // Free memory for HDF5 related data
    delete element_data;
    delete field_offsets;
    delete field_sizes;
}

void quakelib::ModelWorld::read_vertex_hdf5(const hid_t &data_file) {
    std::vector<FieldDesc>                          descs;
    std::map<UIndex, ModelVertex>::const_iterator  fit;
    hsize_t                     num_fields, num_vertices;
    unsigned int                i;
    VertexData                  *vertex_data;
    size_t                      *field_offsets;
    size_t                      *field_sizes;
    herr_t                      res;

    descs.clear();
    ModelVertex::get_field_descs(descs);
    num_fields = descs.size();
    field_offsets = new size_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_offsets[i] = descs[i].offset;
        field_sizes[i] = descs[i].size;
    }

    res = H5TBget_table_info(data_file, ModelVertex::hdf5_table_name().c_str(), &num_fields, &num_vertices);

    if (res < 0) exit(-1);

    // TODO: check that num_fields matches the descs

    vertex_data = new VertexData[num_vertices];
    res = H5TBread_records(data_file, ModelVertex::hdf5_table_name().c_str(), 0, num_vertices, sizeof(VertexData), field_offsets, field_sizes, vertex_data);

    if (res < 0) exit(-1);

    // Read vertex data into the World
    for (i=0; i<num_vertices; ++i) {
        ModelVertex  new_vertex;
        new_vertex.read_data(vertex_data[i]);
        _vertices.insert(std::make_pair(new_vertex.id(), new_vertex));
    }

    // Free memory for HDF5 related data
    delete vertex_data;
    delete field_offsets;
    delete field_sizes;
}

void quakelib::ModelWorld::write_section_hdf5(const hid_t &data_file) const {
    std::vector<FieldDesc>                          descs;
    std::map<UIndex, ModelSection>::const_iterator  fit;
    size_t                      num_fields, num_sections;
    unsigned int                i;
    SectionData                 blank_section, *section_data;
    char                        **field_names, **field_details;
    size_t                      *field_offsets;
    hid_t                       *field_types;
    size_t                      *field_sizes;
    herr_t                      res;

    // Set up the section table definition
    descs.clear();
    ModelSection::get_field_descs(descs);
    num_fields = descs.size();
    num_sections = _sections.size();
    field_names = new char *[num_fields];
    field_details = new char *[num_fields];
    field_offsets = new size_t[num_fields];
    field_types = new hid_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_names[i] = new char[descs[i].name.length()+1];
        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
        field_names[i][descs[i].name.length()] = '\0';
        field_details[i] = new char[descs[i].details.length()+1];
        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
        field_details[i][descs[i].details.length()] = '\0';
        field_offsets[i] = descs[i].offset;
        field_types[i] = descs[i].type;
        field_sizes[i] = descs[i].size;
    }

    // TODO: factor this out?
    blank_section = ModelSection().data();

    // Fill in the data for the sections
    section_data = new SectionData[num_sections];

    for (i=0,fit=_sections.begin(); fit!=_sections.end(); ++i,++fit) {
        fit->second.write_data(section_data[i]);
    }

    // Create the section table
    res = H5TBmake_table("Fault Sections",
                         data_file,
                         ModelSection::hdf5_table_name().c_str(),
                         num_fields,
                         num_sections,
                         sizeof(SectionData),
                         (const char **)field_names,
                         field_offsets,
                         field_types,
                         num_sections,
                         &blank_section,
                         0,
                         section_data);

    if (res < 0) exit(-1);

    // Add the details of each field as an attribute
    for (i=0; i<num_fields; ++i) {
        std::stringstream   ss;
        ss << "FIELD_" << i << "_DETAILS";
        res = H5LTset_attribute_string(data_file,
                                       ModelSection::hdf5_table_name().c_str(),
                                       ss.str().c_str(),
                                       field_details[i]);

        if (res < 0) exit(-1);
    }

    // Free memory for HDF5 related data
    delete section_data;

    for (i=0; i<num_fields; ++i) delete field_names[i];

    delete field_names;

    for (i=0; i<num_fields; ++i) delete field_details[i];

    delete field_details;
    delete field_offsets;
    delete field_types;
    delete field_sizes;
}

void quakelib::ModelWorld::write_element_hdf5(const hid_t &data_file) const {
    std::vector<FieldDesc>                          descs;
    std::map<UIndex, ModelElement>::const_iterator  eit;
    size_t                      num_fields, num_elements;
    unsigned int                i;
    ElementData                 blank_element, *element_data;
    char                        **field_names, **field_details;
    size_t                      *field_offsets;
    hid_t                       *field_types;
    size_t                      *field_sizes;
    herr_t                      res;

    // Set up the element table definition
    descs.clear();
    ModelElement::get_field_descs(descs);
    num_fields = descs.size();
    num_elements = _elements.size();
    field_names = new char *[num_fields];
    field_details = new char *[num_fields];
    field_offsets = new size_t[num_fields];
    field_types = new hid_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_names[i] = new char[descs[i].name.length()+1];
        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
        field_names[i][descs[i].name.length()] = '\0';
        field_details[i] = new char[descs[i].details.length()+1];
        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
        field_details[i][descs[i].details.length()] = '\0';
        field_offsets[i] = descs[i].offset;
        field_types[i] = descs[i].type;
        field_sizes[i] = descs[i].size;
    }

    // Get a blank element for table filling
    blank_element = ModelElement().data();

    // Fill in the data for the elements
    element_data = new ElementData[num_elements];

    for (i=0,eit=_elements.begin(); eit!=_elements.end(); ++i,++eit) {
        eit->second.write_data(element_data[i]);
    }

    // Create the elements table
    res = H5TBmake_table("Elements",
                         data_file,
                         ModelElement::hdf5_table_name().c_str(),
                         num_fields,
                         num_elements,
                         sizeof(ElementData),
                         (const char **)field_names,
                         field_offsets,
                         field_types,
                         num_elements,
                         &blank_element,
                         0,
                         element_data);

    if (res < 0) exit(-1);

    // Add the details of each field as an attribute
    for (i=0; i<num_fields; ++i) {
        std::stringstream   ss;
        ss << "FIELD_" << i << "_DETAILS";
        res = H5LTset_attribute_string(data_file,
                                       ModelElement::hdf5_table_name().c_str(),
                                       ss.str().c_str(),
                                       field_details[i]);

        if (res < 0) exit(-1);
    }

    // Free memory for HDF5 related data
    delete element_data;

    for (i=0; i<num_fields; ++i) delete field_names[i];

    delete field_names;

    for (i=0; i<num_fields; ++i) delete field_details[i];

    delete field_details;
    delete field_offsets;
    delete field_types;
    delete field_sizes;
}

void quakelib::ModelWorld::write_vertex_hdf5(const hid_t &data_file) const {
    std::vector<FieldDesc>                          descs;
    std::map<UIndex, ModelVertex>::const_iterator   vit;
    size_t                      num_fields, num_vertices;
    unsigned int                i;
    VertexData                  blank_vertex, *vertex_data;
    char                        **field_names, **field_details;
    size_t                      *field_offsets;
    hid_t                       *field_types;
    size_t                      *field_sizes;
    herr_t                      res;

    // Set up the vertex table definition
    descs.clear();
    ModelVertex::get_field_descs(descs);
    num_fields = descs.size();
    num_vertices = _vertices.size();
    field_names = new char *[num_fields];
    field_details = new char *[num_fields];
    field_offsets = new size_t[num_fields];
    field_types = new hid_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_names[i] = new char[descs[i].name.length()+1];
        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
        field_names[i][descs[i].name.length()] = '\0';
        field_details[i] = new char[descs[i].details.length()+1];
        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
        field_details[i][descs[i].details.length()] = '\0';
        field_offsets[i] = descs[i].offset;
        field_types[i] = descs[i].type;
        field_sizes[i] = descs[i].size;
    }

    blank_vertex = ModelVertex().data();

    // Fill in the data for the vertices
    vertex_data = new VertexData[num_vertices];

    for (i=0,vit=_vertices.begin(); vit!=_vertices.end(); ++i,++vit) {
        vit->second.write_data(vertex_data[i]);
    }

    // Create the vertices table
    res = H5TBmake_table("Vertices",
                         data_file,
                         ModelVertex::hdf5_table_name().c_str(),
                         num_fields,
                         num_vertices,
                         sizeof(VertexData),
                         (const char **)field_names,
                         field_offsets,
                         field_types,
                         num_vertices,
                         &blank_vertex,
                         0,
                         vertex_data);

    if (res < 0) exit(-1);

    // Add the details of each field as an attribute
    for (i=0; i<num_fields; ++i) {
        std::stringstream   ss;
        ss << "FIELD_" << i << "_DETAILS";
        res = H5LTset_attribute_string(data_file,
                                       ModelVertex::hdf5_table_name().c_str(),
                                       ss.str().c_str(),
                                       field_details[i]);

        if (res < 0) exit(-1);
    }

    // Free memory for HDF5 related data
    delete vertex_data;

    for (i=0; i<num_fields; ++i) delete field_names[i];

    delete field_names;

    for (i=0; i<num_fields; ++i) delete field_details[i];

    delete field_details;
    delete field_offsets;
    delete field_types;
    delete field_sizes;
}
#endif

int quakelib::ModelWorld::write_file_hdf5(const std::string &file_name) const {
#ifdef HDF5_FOUND
    hid_t       plist_id, data_file;
    herr_t      res;

    // Create access properties
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    if (plist_id < 0) exit(-1);

    // Create the data file, overwriting any old files
    data_file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    if (data_file < 0) exit(-1);

    write_section_hdf5(data_file);
    write_element_hdf5(data_file);
    write_vertex_hdf5(data_file);

    // Release HDF5 handles
    res = H5Pclose(plist_id);

    if (res < 0) exit(-1);

    res = H5Fclose(data_file);

    if (res < 0) exit(-1);

    return 0;
#else
    return 1;
#endif
}

void quakelib::ModelWorld::get_bounds(LatLonDepth &minimum, LatLonDepth &maximum) const {
    std::map<UIndex, ModelVertex>::const_iterator    it;
    double      min_lat, min_lon, min_alt;
    double      max_lat, max_lon, max_alt;

    min_lat = min_lon = min_alt = DBL_MAX;
    max_lat = max_lon = max_alt = -DBL_MAX;

    for (it=_vertices.begin(); it!=_vertices.end(); ++it) {
        min_lat = fmin(min_lat, it->second.lld().lat());
        max_lat = fmax(max_lat, it->second.lld().lat());
        min_lon = fmin(min_lon, it->second.lld().lon());
        max_lon = fmax(max_lon, it->second.lld().lon());
        min_alt = fmin(min_alt, it->second.lld().altitude());
        max_alt = fmax(max_alt, it->second.lld().altitude());
    }

    if (min_lat == DBL_MAX || min_lon == DBL_MAX || min_alt == DBL_MAX) {
        minimum = LatLonDepth();
    } else {
        minimum = LatLonDepth(min_lat, min_lon, min_alt);
    }

    if (max_lat == -DBL_MAX || max_lon == -DBL_MAX || max_alt == -DBL_MAX) {
        maximum = LatLonDepth();
    } else {
        maximum = LatLonDepth(max_lat, max_lon, max_alt);
    }
}

quakelib::SimElement quakelib::ModelWorld::create_sim_element(const UIndex &element_id) const {
    SimElement          new_element;
    std::map<UIndex, ModelElement>::const_iterator  eit;
    std::map<UIndex, ModelVertex>::const_iterator   vit;
    unsigned int        i;

    eit = _elements.find(element_id);

    for (i=0; i<3; ++i) {
        vit = _vertices.find(eit->second.vertex(i));
        new_element.set_vert(i, vit->second.xyz());
    }

    new_element.set_is_quad(eit->second.is_quad());
    new_element.set_rake(eit->second.rake());
    new_element.set_slip_rate(eit->second.slip_rate());
    new_element.set_aseismic(eit->second.aseismic());
    new_element.set_lame_mu(eit->second.lame_mu());
    new_element.set_lame_lambda(eit->second.lame_lambda());
    new_element.set_max_slip(eit->second.max_slip());

    return new_element;
}

int quakelib::ModelWorld::read_files_eqsim(const std::string &geom_file_name, const std::string &cond_file_name, const std::string &fric_file_name) {
    quakelib::ModelWorld            eqsim_world;
    quakelib::EQSimGeometryReader   geometry_data;
    quakelib::EQSimConditionReader  condition_data;
    quakelib::EQSimFrictionReader   friction_data;
    quakelib::EQSimErrors           error_report;
    bool                            read_cond_file;
    quakelib::EQSimGeomSectionMap::const_iterator   sit;
    quakelib::EQSimGeomRectangleMap::const_iterator it;
    quakelib::eiterator eit;
    quakelib::LatLonDepth           base;
    std::map<UIndex, double> fault_areas;

    // Clear the world first to avoid incorrectly mixing indices
    clear();

    geometry_data.parse_file(geom_file_name);
    geometry_data.validate(error_report);
    error_report.write(std::cerr);
    error_report.clear();

    friction_data.parse_file(fric_file_name);
    friction_data.validate(error_report);
    error_report.write(std::cerr);
    error_report.clear();

    read_cond_file = condition_data.parse_file(cond_file_name);

    if (read_cond_file) {
        condition_data.validate(error_report);
        error_report.write(std::cerr);
        error_report.clear();
    }

    // Take the conversion base as the middle of the section map
    base = quakelib::LatLonDepth(geometry_data.lat_lo(), geometry_data.lon_lo());

    // Initiate converter
    Conversion          conv(base);

    // Triangle elements are currently not supported
    if (geometry_data.num_triangles() > 0) {
        std::cerr << "ERROR: Currently cannot handle EQSim triangle elements. These elements will be ignored." << std::endl;
    }

    // Go through the geometry and create sections for each EQSim section
    for (sit=geometry_data.sections.begin(); sit!=geometry_data.sections.end(); ++sit) {
        eqsim_world.insert(sit->second.create_model_section());

        // Assuming aligned rectangular elements
        for (it=sit->second.rectangles.begin(); it!=sit->second.rectangles.end(); ++it) {
            quakelib::ModelElement   new_element;
            unsigned int             i;

            new_element = it->second.create_model_element();
            new_element.set_section_id(sit->second.sid());

            for (i=0; i<4; ++i) eqsim_world.insert(sit->second.vertices.find(it->second.vertex(i))->second.create_model_vertex(base));

            new_element.set_lame_lambda(friction_data.get_lame_lambda());
            new_element.set_lame_mu(friction_data.get_lame_mu());
            new_element.set_max_slip(0);    // Set a temporary maximum slip of 0 (this will be changed below)

            // Insert partially finished element into the eqsim_world
            eqsim_world.insert(new_element);

            // Compute area of the current element, add it to the total for this section
            fault_areas[sit->second.sid()] += eqsim_world.create_sim_element(new_element.id()).area();
        }

        // Go through the created elements and assign maximum slip based on fault section area
        for (eit=eqsim_world.begin_element(); eit!=eqsim_world.end_element(); ++eit) {
            // From Table 2A in Wells Coppersmith 1994
            double moment_magnitude = 4.07+0.98*log10(conv.sqm2sqkm(fault_areas[eit->section_id()]));

            // Need to document where this scaling law comes from
            double max_slip = pow(10, (3.0/2.0)*(moment_magnitude+10.7))/(1e7*(eit->lame_mu())*fault_areas[eit->section_id()]);

            // Set the max slip for the current element
            eit->set_max_slip(max_slip);
        }
    }

    insert(eqsim_world);

    return 0;
}

int quakelib::ModelWorld::write_files_eqsim(const std::string &geom_file_name, const std::string &cond_file_name, const std::string &fric_file_name) {
    EQSimGeometryWriter     geometry_data;
    EQSimFrictionWriter     friction_data;
    eiterator               eit;
    siterator               sit;
    UIndex                  sid;
    UIndex                  vind, eind;
    LatLonDepth             base = min_bound();
    Conversion              c(base);

    vind = eind = 1;
    friction_data.set_lame_lambda_mu(3.2e10, 3.0e10);

    for (sit=begin_section(); sit!=end_section(); ++sit) {
        EQSimGeometrySection &section = geometry_data.new_section();

        // Set section properties
        section.set_name(sit->name());
        section.set_fid(sit->fault_id());

        sid = sit->id();

        for (eit=begin_element(sid); eit!=end_element(sid); ++eit) {
            // Create SimElement to allow dip calculation
            SimElement      sim_elem = create_sim_element(eit->id());

            // Create EQSim rectangle and vertices
            EQSimGeometryRectangle  &elem = section.new_rectangle();
            EQSimGeometryVertex &v0 = section.new_vertex();
            v0.set_index(vind++);
            EQSimGeometryVertex &v1 = section.new_vertex();
            v1.set_index(vind++);
            EQSimGeometryVertex &v2 = section.new_vertex();
            v2.set_index(vind++);
            EQSimGeometryVertex &v3 = section.new_vertex();
            v3.set_index(vind++);

            // Set element properties
            friction_data.set_strengths(eind, 2.261600e+007, 0);
            elem.set_index(eind++);
            elem.set_vertex(0, v0.index());
            elem.set_vertex(1, v1.index());
            elem.set_vertex(2, v2.index());
            elem.set_vertex(3, v3.index());
            elem.set_rake(c.rad2deg(eit->rake()));
            elem.set_slip_rate(eit->slip_rate());
            elem.set_aseismic(eit->aseismic());
            // TODO: set strike
            elem.set_dip(c.rad2deg(sim_elem.dip()));
            elem.set_perfect_flag(1);

            // Set vertex properties
            v0.set_loc(vertex(eit->vertex(0)).lld());
            v0.set_das(vertex(eit->vertex(0)).das());
            v0.set_trace_flag(NOT_ON_TRACE);
            v1.set_loc(vertex(eit->vertex(1)).lld());
            v1.set_das(vertex(eit->vertex(1)).das());
            v1.set_trace_flag(NOT_ON_TRACE);
            v2.set_loc(vertex(eit->vertex(2)).lld());
            v2.set_das(vertex(eit->vertex(2)).das());
            v2.set_trace_flag(NOT_ON_TRACE);
            // Fourth vertex is implicit so we calculate it
            v3.set_loc(c.convert2LatLon(sim_elem.implicit_vert()));
            v3.set_das(vertex(eit->vertex(2)).das());
            v3.set_trace_flag(NOT_ON_TRACE);
        }
    }

    geometry_data.open(geom_file_name);
    geometry_data.write();
    geometry_data.close();

    friction_data.open(fric_file_name);
    friction_data.write();
    friction_data.close();

    return 0;
}

int quakelib::ModelWorld::write_file_kml(const std::string &file_name) {
    std::ofstream                                   out_file;
    std::map<UIndex, ModelSection>::const_iterator  fit;
    std::map<UIndex, ModelElement>::const_iterator  eit;
    LatLonDepth                                     min_bound, max_bound, center;
    Vec<3>                                          min_xyz, max_xyz;
    double                                          dx, dy, range;
    unsigned int                                    i;

    out_file.open(file_name.c_str());

    get_bounds(min_bound, max_bound);
    center = LatLonDepth(max_bound.lat() - (max_bound.lat()-min_bound.lat())/2,
                         max_bound.lon() - (max_bound.lon()-min_bound.lon())/2);
    Conversion c(center);
    min_xyz = c.convert2xyz(min_bound);
    max_xyz = c.convert2xyz(max_bound);
    dx = max_xyz[0]-min_xyz[0];
    dy = max_xyz[1]-min_xyz[1];
    range = fmax(dx, dy) * (1.0/tan(c.deg2rad(30)));

    out_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out_file << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    out_file << "<Document>\n";
    out_file << "<LookAt>\n";
    out_file << "\t<latitude>" << center.lat() << "</latitude>\n";
    out_file << "\t<longitude>" << center.lon() << "</longitude>\n";
    out_file << "\t<altitude>0</altitude>\n";
    out_file << "\t<range>" << range << "</range>\n";
    out_file << "\t<tilt>0</tilt>\n";
    out_file << "\t<heading>0</heading>\n";
    out_file << "\t<altitudeMode>absolute</altitudeMode>\n";
    out_file << "</LookAt>\n";
    out_file << "<Style id=\"sectionLabel\">\n";
    out_file << "\t<IconStyle>\n";
    out_file << "\t\t<Icon>\n";
    out_file << "\t\t\t<href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href>\n";
    out_file << "\t\t</Icon>\n";
    out_file << "\t</IconStyle>\n";
    out_file << "</Style>\n";
    out_file << "<Folder id=\"fault_names\">\n";

    for (fit=_sections.begin(); fit!=_sections.end(); ++fit) {
        out_file << "\t<Placemark id=\"section_" << fit->second.id() << "_label\">\n";
        out_file << "\t\t<name>" << fit->second.id() << " " << fit->second.name() << "</name>\n";
        out_file << "\t\t<styleUrl>#sectionLabel</styleUrl>\n";
        out_file << "\t\t<Point>\n";
        out_file << "\t\t\t<extrude>1</extrude>\n";
        out_file << "\t\t\t<altitudeMode>relativeToGround</altitudeMode>\n";
        // Find the deepest element for this section
        UIndex          best_vertex;
        double          min_altitude = DBL_MAX, cur_alt;
        ModelVertex     vert;

        for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
            if (eit->second.section_id() == fit->second.id()) {
                for (i=0; i<3; ++i) {
                    cur_alt = _vertices.find(eit->second.vertex(i))->second.lld().altitude();

                    if (cur_alt < min_altitude) {
                        min_altitude = cur_alt;
                        best_vertex = eit->second.vertex(i);
                    }
                }
            }
        }

        vert = _vertices.find(best_vertex)->second;
        out_file << "\t\t\t<coordinates>" << vert.lld().lon() << "," << vert.lld().lat() << "," << fabs(vert.lld().altitude()) << "</coordinates>\n";
        out_file << "\t\t</Point>\n";
        out_file << "\t</Placemark>\n";
    }

    out_file << "</Folder>\n";
    out_file << "<Folder id=\"faults\">\n";

    // Go through the sections
    for (fit=_sections.begin(); fit!=_sections.end(); ++fit) {
        // And output the elements for each section
        out_file << "\t<Folder id=\"section_" << fit->first << "\">\n";
        out_file << "\t\t<name>" << fit->first << " " << fit->second.name() << "</name>\n";

        for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
            if (fit->first == eit->second.section_id()) {
                LatLonDepth         lld[4];
                unsigned int        i, npoints;

                for (i=0; i<3; ++i) {
                    std::map<UIndex, ModelVertex>::const_iterator   it;
                    it = _vertices.find(eit->second.vertex(i));
                    lld[i] = it->second.lld();
                }

                // If this is a quad element, calculate the 4th implicit point
                if (eit->second.is_quad()) {
                    Vec<3>              xyz[3];

                    for (i=0; i<3; ++i) xyz[i] = c.convert2xyz(lld[i]);

                    lld[3] = lld[2];
                    lld[2] = c.convert2LatLon(xyz[2]+(xyz[1]-xyz[0]));
                }

                // Output the KML format polygon for this section
                out_file << "\t\t<Placemark>\n";
                out_file << "\t\t<description>\n";
                out_file << "Fault name: " << fit->second.name() << "\n";
                out_file << "Slip rate: " << c.m_per_sec2cm_per_yr(eit->second.slip_rate()) << " cm/year\n";
                out_file << "Rake: " << c.rad2deg(eit->second.rake()) << " degrees\n";
                out_file << "Aseismic: " << eit->second.aseismic() << "\n";
                out_file << "\t\t</description>\n";
                out_file << "\t\t\t<styleUrl>#baseStyle</styleUrl>\n";
                out_file << "\t\t\t<Polygon>\n";
                out_file << "\t\t\t\t<extrude>0</extrude>\n";
                out_file << "\t\t\t\t<altitudeMode>relativeToGround</altitudeMode>\n";
                out_file << "\t\t\t\t<outerBoundaryIs>\n";
                out_file << "\t\t\t\t\t<LinearRing>\n";
                out_file << "\t\t\t\t\t\t<coordinates>\n";
                npoints = (eit->second.is_quad() ? 4 : 3);

                for (i=0; i<npoints+1; ++i) out_file << "\t\t\t\t\t\t\t" << lld[i%npoints].lon() << "," << lld[i%npoints].lat() << "," << fabs(lld[i%npoints].altitude()) << "\n";

                out_file << "\t\t\t\t\t\t</coordinates>\n";
                out_file << "\t\t\t\t\t</LinearRing>\n";
                out_file << "\t\t\t\t</outerBoundaryIs>\n";
                out_file << "\t\t\t</Polygon>\n";
                out_file << "\t\t</Placemark>\n";
            }
        }

        out_file << "\t</Folder>\n";
    }

    out_file << "</Folder>\n";
    out_file << "</Document>\n";
    out_file << "</kml>\n";

    out_file.close();

    return 0;
}

void quakelib::ModelSection::apply_remap(const ModelRemapping &remap) {
    _data._id = remap.get_section_id(_data._id);
}

void quakelib::ModelElement::apply_remap(const ModelRemapping &remap) {
    unsigned int            i;

    _data._id = remap.get_element_id(_data._id);
    _data._section_id = remap.get_section_id(_data._section_id);

    for (i=0; i<3; ++i) _data._vertices[i] = remap.get_vertex_id(_data._vertices[i]);
}

void quakelib::ModelVertex::apply_remap(const ModelRemapping &remap) {
    _data._id = remap.get_vertex_id(_data._id);
}

quakelib::UIndex quakelib::ModelRemapping::get_section_id(const UIndex &orig_id) const {
    IndexRemap::const_iterator it=section_remap.find(orig_id);
    return (it==section_remap.end() ? orig_id : it->second);
};

quakelib::UIndex quakelib::ModelRemapping::get_element_id(const UIndex &orig_id) const {
    IndexRemap::const_iterator it=element_remap.find(orig_id);
    return (it==element_remap.end() ? orig_id : it->second);
};

quakelib::UIndex quakelib::ModelRemapping::get_vertex_id(const UIndex &orig_id) const {
    IndexRemap::const_iterator it=vertex_remap.find(orig_id);
    return (it==vertex_remap.end() ? orig_id : it->second);
};

// Returns true if the remapping will overwrite (and lose) some sections/elements/vertices, false otherwise
bool quakelib::ModelWorld::overwrite(const ModelRemapping &remap) {
    std::map<UIndex, ModelSection>::const_iterator  fit;
    std::map<UIndex, ModelElement>::const_iterator  eit;
    std::map<UIndex, ModelVertex>::const_iterator   vit;
    std::set<UIndex>    remap_inds;
    UIndex              remap_ind;

    remap_inds.clear();

    for (fit=_sections.begin(); fit!=_sections.end(); ++fit) {
        remap_ind = remap.get_section_id(fit->first);

        if (remap_ind == INVALID_INDEX) remap_ind = fit->first;

        if (remap_inds.count(remap_ind) > 0) return true;

        remap_inds.insert(remap_ind);
    }

    remap_inds.clear();

    for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
        remap_ind = remap.get_element_id(eit->first);

        if (remap_ind == INVALID_INDEX) remap_ind = eit->first;

        if (remap_inds.count(remap_ind) > 0) return true;

        remap_inds.insert(remap_ind);
    }

    remap_inds.clear();

    for (vit=_vertices.begin(); vit!=_vertices.end(); ++vit) {
        remap_ind = remap.get_section_id(vit->first);

        if (remap_ind == INVALID_INDEX) remap_ind = vit->first;

        if (remap_inds.count(remap_ind) > 0) return true;

        remap_inds.insert(remap_ind);
    }

    return false;
}

// Apply the specified remapping to the sections, elements and vertices in this world
// Depending on the remapping it is entirely possible for things to be overwritten, it is the responsibility of the user to check this
void quakelib::ModelWorld::apply_remap(const ModelRemapping &remap) {
    std::map<UIndex, ModelSection>::iterator    fit;
    std::map<UIndex, ModelElement>::iterator    eit;
    std::map<UIndex, ModelVertex>::iterator     vit;
    std::vector<UIndex>                         erased_sections, erased_elements, erased_vertices;
    std::vector<UIndex>::iterator               it;
    std::vector<ModelSection>                   section_set;
    std::vector<ModelElement>                   element_set;
    std::vector<ModelVertex>                    vertex_set;
    unsigned int                                i;

    // Pull out any sections that will be remapped
    for (fit=_sections.begin(); fit!=_sections.end(); ++fit) {
        if (remap.get_section_id(fit->first) != INVALID_INDEX) {
            section_set.push_back(fit->second);
            section_set.back().apply_remap(remap);
            erased_sections.push_back(fit->first);
        } else {
            fit->second.apply_remap(remap);
        }
    }

    // Remove the old sections from the map
    for (it=erased_sections.begin(); it!=erased_sections.end(); ++it) _sections.erase(*it);

    erased_sections.clear();

    // Pull out any elements that will be remapped
    for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
        if (remap.get_element_id(eit->first) != INVALID_INDEX) {
            element_set.push_back(eit->second);
            element_set.back().apply_remap(remap);
            erased_elements.push_back(eit->first);
        } else {
            eit->second.apply_remap(remap);
        }
    }

    // Remove the old elements from the map
    for (it=erased_elements.begin(); it!=erased_elements.end(); ++it) _elements.erase(*it);

    erased_elements.clear();

    // Pull out any vertices that will be remapped
    for (vit=_vertices.begin(); vit!=_vertices.end(); ++vit) {
        if (remap.get_vertex_id(vit->first) != INVALID_INDEX) {
            vertex_set.push_back(vit->second);
            vertex_set.back().apply_remap(remap);
            erased_vertices.push_back(vit->first);
        } else {
            vit->second.apply_remap(remap);
        }
    }

    // Remove the old vertices from the map
    for (it=erased_vertices.begin(); it!=erased_vertices.end(); ++it) _vertices.erase(*it);

    erased_vertices.clear();

    // Go through each stored section and put it back in the world
    for (i=0; i<section_set.size(); ++i) _sections.insert(std::make_pair(section_set[i].id(),section_set[i]));

    section_set.clear();

    // Go through each element, change the indices and put it back in the world
    for (i=0; i<element_set.size(); ++i) _elements.insert(std::make_pair(element_set[i].id(),element_set[i]));

    element_set.clear();

    // Go through each vertex, change the indices and put it back in the world
    for (i=0; i<vertex_set.size(); ++i) _vertices.insert(std::make_pair(vertex_set[i].id(),vertex_set[i]));

    vertex_set.clear();
}

quakelib::ModelRemapping quakelib::ModelWorld::remap_indices_contiguous(const UIndex &start_section_index,
                                                                        const UIndex &start_element_index,
                                                                        const UIndex &start_vertex_index) const {
    ModelRemapping                                  remap;
    std::map<UIndex, ModelSection>::const_iterator  fit;
    std::map<UIndex, ModelElement>::const_iterator  eit;
    std::map<UIndex, ModelVertex>::const_iterator   vit;
    UIndex                                          cur_ind;

    // Change section indices to be in contiguous ascending order
    for (cur_ind=start_section_index,fit=_sections.begin(); fit!=_sections.end(); ++cur_ind,++fit) remap.remap_section(fit->first, cur_ind);

    // Change element indices to be in contiguous ascending order
    for (cur_ind=start_element_index,eit=_elements.begin(); eit!=_elements.end(); ++cur_ind,++eit) remap.remap_element(eit->first, cur_ind);

    // Change vertex indices to be in contiguous ascending order
    for (cur_ind=start_vertex_index,vit=_vertices.begin(); vit!=_vertices.end(); ++cur_ind,++vit) remap.remap_vertex(vit->first, cur_ind);

    return remap;
}

quakelib::ModelRemapping quakelib::ModelWorld::remove_duplicate_vertices_remap(void) const {
    ModelRemapping                                  remap;
    std::map<UIndex, ModelVertex>::const_iterator   vit;
    std::map<Vec<3>, UIndex>                        vert_map;
    std::map<Vec<3>, UIndex>::const_iterator        mit;

    for (vit=_vertices.begin(); vit!=_vertices.end(); ++vit) {
        Vec<3>  pos = vit->second.xyz();
        mit = vert_map.find(pos);

        if (mit!= vert_map.end()) {
            remap.remap_vertex(vit->first, mit->second);
        } else {
            vert_map.insert(std::make_pair(pos, vit->first));
        }
    }

    return remap;
}

quakelib::LatLonDepth quakelib::ModelWorld::min_bound(const UIndex &fid) const {
    std::map<UIndex, ModelElement>::const_iterator  eit;
    bool            check_elem_verts;
    unsigned int    i;
    double          min_lat, min_lon, min_alt;

    min_lat = min_lon = min_alt = DBL_MAX;

    for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
        check_elem_verts = (fid == INVALID_INDEX) || (eit->second.section_id() == fid);

        if (check_elem_verts) {
            // TODO: for quads determine 4th implicit point and check it
            for (i=0; i<3; ++i) {
                const ModelVertex &v = _vertices.find(eit->second.vertex(i))->second;
                min_lat = fmin(min_lat, v.lld().lat());
                min_lon = fmin(min_lon, v.lld().lon());
                min_alt = fmin(min_alt, v.lld().altitude());
            }
        }
    }

    return LatLonDepth(min_lat, min_lon, min_alt);
}

quakelib::LatLonDepth quakelib::ModelWorld::max_bound(const UIndex &fid) const {
    std::map<UIndex, ModelElement>::const_iterator  eit;
    bool            check_elem_verts;
    unsigned int    i;
    double          max_lat, max_lon, max_alt;

    max_lat = max_lon = max_alt = DBL_MAX;

    for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
        check_elem_verts = (fid == INVALID_INDEX) || (eit->second.section_id() == fid);

        if (check_elem_verts) {
            // TODO: for quads determine 4th implicit point and check it
            for (i=0; i<3; ++i) {
                const ModelVertex &v = _vertices.find(eit->second.vertex(i))->second;
                max_lat = fmax(max_lat, v.lld().lat());
                max_lon = fmax(max_lon, v.lld().lon());
                max_alt = fmax(max_alt, v.lld().altitude());
            }
        }
    }

    return LatLonDepth(max_lat, max_lon, max_alt);
}

// TODO: rebase vertices
void quakelib::ModelWorld::insert(const quakelib::ModelWorld &other_world) {
    UIndex                      next_section_ind, next_element_ind, next_vertex_ind;
    ModelRemapping              remap;
    std::vector<ModelSection>   section_list;
    std::vector<ModelElement>   element_list;
    std::vector<ModelVertex>    vertex_list;
    unsigned int                i;
    std::map<UIndex, ModelSection>::const_iterator  fit;
    std::map<UIndex, ModelElement>::const_iterator  eit;
    std::map<UIndex, ModelVertex>::const_iterator   vit;
    LatLonDepth                 this_lld_min, this_lld_max;

    // Get the maximum indices in this model
    next_section_ind = next_section_index();
    next_element_ind = next_element_index();
    next_vertex_ind = next_vertex_index();

    // Create a remapping of the sections/elements/vertices in the
    // other model to allow it to be inserted in this model
    remap = other_world.remap_indices_contiguous(next_section_ind, next_element_ind, next_vertex_ind);

    // Make a copy of sections from other world, apply remapping, then insert into this world
    // Don't directly insert, since there's no guarantee that this world and the other aren't the same
    for (fit=other_world._sections.begin(); fit!=other_world._sections.end(); ++fit) {
        ModelSection  other_section = fit->second;
        other_section.apply_remap(remap);
        section_list.push_back(other_section);
    }

    for (i=0; i<section_list.size(); ++i) _sections.insert(std::make_pair(section_list[i].id(), section_list[i]));

    section_list.clear();

    // Insert elements from other world into this one, applying remapping
    for (eit=other_world._elements.begin(); eit!=other_world._elements.end(); ++eit) {
        ModelElement  other_element = eit->second;
        other_element.apply_remap(remap);
        element_list.push_back(other_element);
    }

    for (i=0; i<element_list.size(); ++i) _elements.insert(std::make_pair(element_list[i].id(), element_list[i]));

    element_list.clear();

    // Insert vertices from other world into this one, applying remapping
    for (vit=other_world._vertices.begin(); vit!=other_world._vertices.end(); ++vit) {
        ModelVertex     other_vertex = vit->second;
        other_vertex.apply_remap(remap);
        vertex_list.push_back(other_vertex);
    }

    for (i=0; i<vertex_list.size(); ++i) _vertices.insert(std::make_pair(vertex_list[i].id(), vertex_list[i]));

    // Get the new minimum latlon for the updated model and reset the vertices to match
    this->get_bounds(this_lld_min, this_lld_max);
    reset_base_coord(this_lld_min);
}

// TODO: correctly handle preexisting map entries
void quakelib::ModelWorld::insert(const quakelib::ModelSection &new_section) {
    _sections.insert(std::make_pair(new_section.id(), new_section));
}

// TODO: correctly handle preexisting map entries
void quakelib::ModelWorld::insert(const quakelib::ModelElement &new_element) {
    _elements.insert(std::make_pair(new_element.id(), new_element));
}

// TODO: correctly handle preexisting map entries
void quakelib::ModelWorld::insert(const quakelib::ModelVertex &new_vertex) {
    _vertices.insert(std::make_pair(new_vertex.id(), new_vertex));
}

size_t quakelib::ModelWorld::num_sections(void) const {
    return _sections.size();
}

size_t quakelib::ModelWorld::num_faults(void) const {
    std::set<UIndex>    fault_ids;
    std::map<UIndex, ModelSection>::const_iterator  sit;

    for (sit=_sections.begin(); sit!=_sections.end(); ++sit) {
        fault_ids.insert(sit->second.fault_id());
    }

    return fault_ids.size();
}

size_t quakelib::ModelWorld::num_elements(const quakelib::UIndex &fid) const {
    std::map<UIndex, ModelElement>::const_iterator  eit;
    size_t  num_elements = 0;

    if (fid == INVALID_INDEX) return _elements.size();

    for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
        if (eit->second.section_id() == fid) num_elements++;
    }

    return num_elements;
}

size_t quakelib::ModelWorld::num_vertices(const quakelib::UIndex &fid) const {
    std::set<UIndex>        section_vertices;
    std::map<UIndex, ModelElement>::const_iterator  eit;
    unsigned int            i;

    if (fid == INVALID_INDEX) return _vertices.size();

    for (eit=_elements.begin(); eit!=_elements.end(); ++eit) {
        if (eit->second.section_id() == fid) {
            for (i=0; i<3; ++i) section_vertices.insert(eit->second.vertex(i));
        }
    }

    return section_vertices.size();
}

void quakelib::ModelEvent::read_ascii(std::istream &in_stream) {
    std::stringstream   ss(next_line(in_stream));

    ss >> _data._event_number;
    ss >> _data._event_year;
    ss >> _data._event_trigger;
    ss >> _data._event_magnitude;
    ss >> _data._shear_stress_init;
    ss >> _data._shear_stress_final;
    ss >> _data._normal_stress_init;
    ss >> _data._normal_stress_final;
    ss >> _data._start_sweep_rec;
    ss >> _data._end_sweep_rec;
}

void quakelib::ModelEvent::write_ascii(std::ostream &out_stream) const {
    out_stream << _data._event_number << " ";
    out_stream << _data._event_year << " ";
    out_stream << _data._event_trigger << " ";
    out_stream << _data._event_magnitude << " ";
    out_stream << _data._shear_stress_init << " ";
    out_stream << _data._shear_stress_final << " ";
    out_stream << _data._normal_stress_init << " ";
    out_stream << _data._normal_stress_final << " ";
    out_stream << _data._start_sweep_rec << " ";
    out_stream << _data._end_sweep_rec;

    next_line(out_stream);
}

void quakelib::ModelSweeps::read_ascii(std::istream &in_stream, const unsigned int num_records) {
    for (unsigned int i=0; i<num_records; ++i) {
        std::stringstream   ss(next_line(in_stream));
        SweepData   new_sweep;
        ss >> new_sweep._event_number;
        ss >> new_sweep._sweep_number;
        ss >> new_sweep._element_id;
        ss >> new_sweep._slip;
        ss >> new_sweep._area;
        ss >> new_sweep._mu;
        ss >> new_sweep._shear_init;
        ss >> new_sweep._shear_final;
        ss >> new_sweep._normal_init;
        ss >> new_sweep._normal_final;
        // Record the sweep/element in the mapping
        std::pair<UIndex, UIndex> sweep_elem = std::make_pair(new_sweep._sweep_number, new_sweep._element_id);
        _rel.insert(std::make_pair(sweep_elem, _sweeps.size()));
        // Put the sweep on the list
        _sweeps.push_back(new_sweep);
    }
}

void quakelib::ModelSweeps::write_ascii(std::ostream &out_stream) const {
    std::vector<SweepData>::const_iterator it;

    for (it=_sweeps.begin(); it!=_sweeps.end(); ++it) {
        out_stream << it->_event_number << " ";
        out_stream << it->_sweep_number << " ";
        out_stream << it->_element_id << " ";
        out_stream << it->_slip << " ";
        out_stream << it->_area << " ";
        out_stream << it->_mu << " ";
        out_stream << it->_shear_init << " ";
        out_stream << it->_shear_final << " ";
        out_stream << it->_normal_init << " ";
        out_stream << it->_normal_final;

        next_line(out_stream);
    }
}

void quakelib::ModelSweeps::read_data(const SweepData &in_data) {
    // Record the sweep/element in the mapping
    std::pair<UIndex, UIndex> sweep_elem = std::make_pair(in_data._sweep_number, in_data._element_id);
    _rel.insert(std::make_pair(sweep_elem, _sweeps.size()));
    // Put the sweep on the list
    _sweeps.push_back(in_data);
}

void quakelib::ModelSweeps::write_data(SweepData &out_data) const {
    //memcpy(&out_data, &_data, sizeof(SweepData));
}

void quakelib::ModelEvent::read_data(const EventData &in_data) {
    memcpy(&_data, &in_data, sizeof(EventData));
}

void quakelib::ModelEvent::write_data(EventData &out_data) const {
    memcpy(&out_data, &_data, sizeof(EventData));
}

void quakelib::ModelSweeps::get_field_descs(std::vector<quakelib::FieldDesc> &descs) {
    FieldDesc       field_desc;

    // Sweep table definition
    field_desc.name = "event_number";
    field_desc.details = "Event number corresponding to this sweep.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _event_number);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "sweep_number";
    field_desc.details = "Sweep number within the event.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _sweep_number);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "block_id";
    field_desc.details = "Element ID.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _element_id);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "block_slip";
    field_desc.details = "Slip on element in this sweep (meters).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _slip);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "block_area";
    field_desc.details = "Area of element (square meters).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _area);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "block_mu";
    field_desc.details = "Element Lame mu parameter (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _mu);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "shear_init";
    field_desc.details = "Shear stress of element before sweep (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _shear_init);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "shear_final";
    field_desc.details = "Shear stress of element after sweep (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _shear_final);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "normal_init";
    field_desc.details = "Normal stress of element before sweep (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _normal_init);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "normal_final";
    field_desc.details = "Normal stress of element after sweep (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(SweepData, _normal_final);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);
}

void quakelib::ModelSweeps::write_ascii_header(std::ostream &out_stream) {
    std::vector<FieldDesc>                  descs;
    std::vector<FieldDesc>::const_iterator  dit;

    // Write section header
    ModelSweeps::get_field_descs(descs);

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_stream << "# " << dit->name << ": " << dit->details << "\n";
    }

    out_stream << "# ";

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_stream << dit->name << " ";
    }

    out_stream << "\n";
}

#ifdef HDF5_FOUND
void quakelib::ModelSweeps::setup_sweeps_hdf5(const hid_t &data_file) {
    std::vector<FieldDesc>  descs;
    size_t                  num_fields;
    unsigned int            i;
    SweepData               blank_sweep;
    char                    **field_names, **field_details;
    size_t                  *field_offsets;
    hid_t                   *field_types;
    size_t                  *field_sizes;
    herr_t                  res;

    // Set up the section table definition
    descs.clear();
    ModelSweeps::get_field_descs(descs);
    num_fields = descs.size();
    field_names = new char *[num_fields];
    field_details = new char *[num_fields];
    field_offsets = new size_t[num_fields];
    field_types = new hid_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_names[i] = new char[descs[i].name.length()+1];
        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
        field_names[i][descs[i].name.length()] = '\0';
        field_details[i] = new char[descs[i].details.length()+1];
        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
        field_details[i][descs[i].details.length()] = '\0';
        field_offsets[i] = descs[i].offset;
        field_types[i] = descs[i].type;
        field_sizes[i] = descs[i].size;
    }

    blank_sweep._event_number = blank_sweep._sweep_number = UNDEFINED_EVENT_ID;
    blank_sweep._element_id = UNDEFINED_ELEMENT_ID;
    blank_sweep._slip = blank_sweep._area = blank_sweep._mu = std::numeric_limits<float>::quiet_NaN();
    blank_sweep._shear_init = blank_sweep._shear_final = std::numeric_limits<float>::quiet_NaN();
    blank_sweep._normal_init = blank_sweep._normal_final = std::numeric_limits<float>::quiet_NaN();

    // Create the sweep table
    res = H5TBmake_table("Sweeps Table",
                         data_file,
                         ModelSweeps::hdf5_table_name().c_str(),
                         num_fields,
                         0,
                         sizeof(SweepData),
                         (const char **)field_names,
                         field_offsets,
                         field_types,
                         100,
                         &blank_sweep,
                         0,
                         NULL);

    if (res < 0) exit(-1);

    // Add the details of each field as an attribute
    for (i=0; i<num_fields; ++i) {
        std::stringstream   ss;
        ss << "FIELD_" << i << "_DETAILS";
        res = H5LTset_attribute_string(data_file, ModelSweeps::hdf5_table_name().c_str(), ss.str().c_str(), field_details[i]);

        if (res < 0) exit(-1);
    }

    // Free memory for HDF5 related data
    for (i=0; i<num_fields; ++i) delete field_names[i];

    delete field_names;

    for (i=0; i<num_fields; ++i) delete field_details[i];

    delete field_details;
    delete field_offsets;
    delete field_types;
    delete field_sizes;
}

void quakelib::ModelSweeps::append_sweeps_hdf5(const hid_t &data_file) const {
    std::vector<FieldDesc>                  descs;
    std::vector<SweepData>::const_iterator  it;
    size_t                                  num_fields, num_sweeps;
    unsigned int                i;
    SweepData                   *sweep_data;
    size_t                      *field_offsets;
    size_t                      *field_sizes;
    herr_t                      res;

    // Set up the section table definition
    descs.clear();
    ModelSweeps::get_field_descs(descs);
    num_fields = descs.size();
    num_sweeps = _sweeps.size();
    field_offsets = new size_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_offsets[i] = descs[i].offset;
        field_sizes[i] = descs[i].size;
    }

    // Fill in the data for the sections
    sweep_data = new SweepData[num_sweeps];

    for (i=0,it=_sweeps.begin(); it!=_sweeps.end(); ++i,++it) {
        memcpy(&(sweep_data[i]), &(*it), sizeof(SweepData));
    }

    // Create the section table
    res = H5TBappend_records(data_file,
                             ModelSweeps::hdf5_table_name().c_str(),
                             num_sweeps,
                             sizeof(SweepData),
                             field_offsets,
                             field_sizes,
                             sweep_data);

    if (res < 0) exit(-1);

    // Free memory for HDF5 related data
    delete sweep_data;

    delete field_offsets;
    delete field_sizes;
}
#endif

void quakelib::ModelEvent::get_field_descs(std::vector<FieldDesc> &descs) {
    FieldDesc       field_desc;

    field_desc.name = "event_number";
    field_desc.details = "Event number.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _event_number);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "event_year";
    field_desc.details = "Event year.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _event_year);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "event_trigger";
    field_desc.details = "Event trigger element ID.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _event_trigger);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(UIndex);
#endif
    descs.push_back(field_desc);

    field_desc.name = "event_magnitude";
    field_desc.details = "Event magnitude.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _event_magnitude);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "event_shear_init";
    field_desc.details = "Total initial shear stress of elements involved in event (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _shear_stress_init);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "event_shear_final";
    field_desc.details = "Total final shear stress of elements involved in event (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _shear_stress_final);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "event_normal_init";
    field_desc.details = "Total initial normal stress of elements involved in event (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _normal_stress_init);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "event_normal_final";
    field_desc.details = "Total final normal stress of elements involved in event (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _normal_stress_final);
    field_desc.type = H5T_NATIVE_DOUBLE;
    field_desc.size = sizeof(double);
#endif
    descs.push_back(field_desc);

    field_desc.name = "start_sweep_rec";
    field_desc.details = "Starting record number of the sweeps comprising this event.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _start_sweep_rec);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "end_sweep_rec";
    field_desc.details = "Ending record number of the sweeps comprising this event.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(EventData, _end_sweep_rec);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);
}

void quakelib::ModelEvent::write_ascii_header(std::ostream &out_stream) {
    std::vector<FieldDesc>                  descs;
    std::vector<FieldDesc>::const_iterator  dit;

    // Write section header
    ModelEvent::get_field_descs(descs);

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_stream << "# " << dit->name << ": " << dit->details << "\n";
    }

    out_stream << "# ";

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_stream << dit->name << " ";
    }

    out_stream << "\n";
}

#ifdef HDF5_FOUND
void quakelib::ModelEvent::setup_event_hdf5(const hid_t &data_file) {
    std::vector<FieldDesc>  descs;
    size_t                  num_fields;
    unsigned int            i;
    EventData               blank_event;
    char                    **field_names, **field_details;
    size_t                  *field_offsets;
    hid_t                   *field_types;
    size_t                  *field_sizes;
    herr_t                  res;

    // Set up the section table definition
    descs.clear();
    ModelEvent::get_field_descs(descs);
    num_fields = descs.size();
    field_names = new char *[num_fields];
    field_details = new char *[num_fields];
    field_offsets = new size_t[num_fields];
    field_types = new hid_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_names[i] = new char[descs[i].name.length()+1];
        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
        field_names[i][descs[i].name.length()] = '\0';
        field_details[i] = new char[descs[i].details.length()+1];
        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
        field_details[i][descs[i].details.length()] = '\0';
        field_offsets[i] = descs[i].offset;
        field_types[i] = descs[i].type;
        field_sizes[i] = descs[i].size;
    }

    blank_event._event_number = UNDEFINED_EVENT_ID;
    blank_event._event_year = std::numeric_limits<float>::quiet_NaN();
    blank_event._event_magnitude = std::numeric_limits<float>::quiet_NaN();
    blank_event._event_trigger = UNDEFINED_ELEMENT_ID;
    blank_event._shear_stress_init = blank_event._shear_stress_final = std::numeric_limits<float>::quiet_NaN();
    blank_event._normal_stress_init = blank_event._normal_stress_final = std::numeric_limits<float>::quiet_NaN();
    blank_event._start_sweep_rec = blank_event._end_sweep_rec = UNDEFINED_EVENT_ID;

    // Create the event table
    res = H5TBmake_table("Event Table",
                         data_file,
                         ModelEvent::hdf5_table_name().c_str(),
                         num_fields,
                         0,
                         sizeof(EventData),
                         (const char **)field_names,
                         field_offsets,
                         field_types,
                         100,
                         &blank_event,
                         0,
                         NULL);

    if (res < 0) exit(-1);

    // Add the details of each field as an attribute
    for (i=0; i<num_fields; ++i) {
        std::stringstream   ss;
        ss << "FIELD_" << i << "_DETAILS";
        res = H5LTset_attribute_string(data_file,
                                       ModelEvent::hdf5_table_name().c_str(),
                                       ss.str().c_str(),
                                       field_details[i]);

        if (res < 0) exit(-1);
    }

    // Free memory for HDF5 related data
    for (i=0; i<num_fields; ++i) delete field_names[i];

    delete field_names;

    for (i=0; i<num_fields; ++i) delete field_details[i];

    delete field_details;
    delete field_offsets;
    delete field_types;
    delete field_sizes;
}

void quakelib::ModelEvent::append_event_hdf5(const hid_t &data_file) const {
    std::vector<FieldDesc>  descs;
    size_t                  num_fields;
    unsigned int            i;
    size_t                  *field_offsets;
    size_t                  *field_sizes;
    herr_t                  res;

    // Set up the section table definition
    descs.clear();
    ModelEvent::get_field_descs(descs);
    num_fields = descs.size();
    field_offsets = new size_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_offsets[i] = descs[i].offset;
        field_sizes[i] = descs[i].size;
    }

    // Append the event record
    res = H5TBappend_records(data_file,
                             ModelEvent::hdf5_table_name().c_str(),
                             1,
                             sizeof(EventData),
                             field_offsets,
                             field_sizes,
                             &(_data));

    if (res < 0) exit(-1);

    // Free memory for HDF5 related data
    delete field_offsets;
    delete field_sizes;
}
#endif

int quakelib::ModelEventSet::read_file_ascii(const std::string &event_file_name, const std::string &sweep_file_name) {
    std::ifstream   event_file, sweep_file;
    ModelSweeps     file_sweeps;

    // Try to open the event file
    event_file.open(event_file_name.c_str());

    if (!event_file.is_open()) return -1;

    // Try to open the sweeps file
    sweep_file.open(sweep_file_name.c_str());

    if (!sweep_file.is_open()) return -1;

    // Keep going until we hit the end of either file
    while (!event_file.eof() && !sweep_file.eof()) {
        ModelEvent  new_event;
        ModelSweeps new_sweeps;
        new_event.read_ascii(event_file);
        unsigned int num_rec_sweeps = new_event.getNumRecordedSweeps();
        new_sweeps.read_ascii(sweep_file, num_rec_sweeps);
        new_event.setSweeps(new_sweeps);

        if (!event_file.eof() && !sweep_file.eof()) _events.push_back(new_event);
    }

    // Close the files
    event_file.close();
    sweep_file.close();

    return 0;
}



int quakelib::ModelEventSet::read_file_hdf5(const std::string &file_name) {
#ifdef HDF5_FOUND
    hid_t       plist_id, data_file;
    herr_t      res;

    if (!H5Fis_hdf5(file_name.c_str())) return -1;

    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    if (plist_id < 0) exit(-1);

    data_file = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, plist_id);

    if (data_file < 0) exit(-1);

    read_events_hdf5(data_file);
    read_sweeps_hdf5(data_file);

    // Release HDF5 handles
    res = H5Pclose(plist_id);

    if (res < 0) exit(-1);

    res = H5Fclose(data_file);

    if (res < 0) exit(-1);

#else
    // TODO: Error out
#endif
    return 0;
}

void quakelib::ModelEventSet::read_events_hdf5(const hid_t &data_file) {
#ifdef HDF5_FOUND
    std::vector<FieldDesc>                        descs;
    std::map<UIndex, ModelEvent>::const_iterator  fit;
    hsize_t                     num_fields, num_events;
    unsigned int                i;
    EventData                   *event_data;
    size_t                      *field_offsets;
    size_t                      *field_sizes;
    herr_t                      res;

    _events.clear();
    descs.clear();
    ModelEvent::get_field_descs(descs);
    num_fields = descs.size();
    field_offsets = new size_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_offsets[i] = descs[i].offset;
        field_sizes[i] = descs[i].size;
    }

    res = H5TBget_table_info(data_file, ModelEvent::hdf5_table_name().c_str(), &num_fields, &num_events);

    if (res < 0) exit(-1);

    // TODO: check that num_fields matches the descs

    event_data = new EventData[num_events];
    res = H5TBread_records(data_file, ModelEvent::hdf5_table_name().c_str(), 0, num_events, sizeof(EventData), field_offsets, field_sizes, event_data);

    if (res < 0) exit(-1);

    // Read section data into the World
    for (i=0; i<num_events; ++i) {
        ModelEvent  new_event;
        new_event.read_data(event_data[i]);
        _events.push_back(new_event);
    }

    // Free memory for HDF5 related data
    delete event_data;
    delete field_offsets;
    delete field_sizes;
#else
    // TODO: Error out
#endif
}

void quakelib::ModelEventSet::read_sweeps_hdf5(const hid_t &data_file) {
#ifdef HDF5_FOUND
    std::vector<FieldDesc>                          descs;
    ModelEventSet::iterator                   fit;
    hsize_t                     num_fields, num_sweeps;
    unsigned int                i;
    unsigned int                start_sweep;
    unsigned int                end_sweep;
    SweepData                   *event_sweeps;
    size_t                      *field_offsets;
    size_t                      *field_sizes;
    herr_t                      res;

    descs.clear();
    ModelSweeps::get_field_descs(descs);
    num_fields = descs.size();
    field_offsets = new size_t[num_fields];
    field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_offsets[i] = descs[i].offset;
        field_sizes[i] = descs[i].size;
    }

    res = H5TBget_table_info(data_file, ModelSweeps::hdf5_table_name().c_str(), &num_fields, &num_sweeps);

    if (res < 0) exit(-1);

    event_sweeps = new SweepData[num_sweeps];
    res = H5TBread_records(data_file, ModelSweeps::hdf5_table_name().c_str(), 0, num_sweeps, sizeof(SweepData), field_offsets, field_sizes, event_sweeps);

    if (res < 0) exit(-1);

    // Read sweeps data into the ModelEventSet
    for (fit=_events.begin(); fit!=_events.end(); ++fit) {
        fit->getStartEndSweep(start_sweep, end_sweep);
        ModelSweeps new_sweeps;

        for (i=start_sweep; i<end_sweep; i++) {
            new_sweeps.read_data(event_sweeps[i]);
        }

        fit->setSweeps(new_sweeps);

    }

    delete event_sweeps;
#else
    // TODO: Error out
#endif
}

// ********************************************************************************************

void quakelib::ModelStress::read_ascii(std::istream &in_stream, const unsigned int num_records) {
    for (unsigned int i=0; i<num_records; ++i) {
        std::stringstream   ss(next_line(in_stream));
        StressData     new_stress_rec;
        ss >> new_stress_rec._element_id;
        ss >> new_stress_rec._shear_stress;
        ss >> new_stress_rec._normal_stress;
        // Put the stress record on the list
        _data.push_back(new_stress_rec);
    }
}

void quakelib::ModelStress::write_ascii(std::ostream &out_stream) const {
    std::vector<StressData>::const_iterator it;

    for (it=_data.begin(); it!=_data.end(); ++it) {
        out_stream << it->_element_id << " ";
        out_stream << it->_shear_stress << " ";
        out_stream << it->_normal_stress;

        next_line(out_stream);
    }
}

#ifdef HDF5_FOUND
void quakelib::ModelStress::setup_stress_hdf5(const hid_t &data_file) {
    std::vector<FieldDesc>  descs;
    size_t                  num_fields;
    unsigned int            i;
    StressData              blank_data;
    herr_t                  res;

    // Set up the section table definition
    descs.clear();
    ModelStress::get_field_descs(descs);
    num_fields = descs.size();
    char **field_names = new char *[num_fields];
    char **field_details = new char *[num_fields];
    size_t *field_offsets = new size_t[num_fields];
    hid_t *field_types = new hid_t[num_fields];
    size_t *field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_names[i] = new char[descs[i].name.length()+1];
        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
        field_names[i][descs[i].name.length()] = '\0';
        field_details[i] = new char[descs[i].details.length()+1];
        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
        field_details[i][descs[i].details.length()] = '\0';
        field_offsets[i] = descs[i].offset;
        field_types[i] = descs[i].type;
        field_sizes[i] = descs[i].size;
    }

    blank_data._element_id = UNDEFINED_ELEMENT_ID;
    blank_data._shear_stress = blank_data._normal_stress = std::numeric_limits<float>::quiet_NaN();

    // Create the sweep table
    res = H5TBmake_table("Stress Table",
                         data_file,
                         ModelStress::hdf5_table_name().c_str(),
                         num_fields,
                         0,
                         sizeof(StressData),
                         (const char **)field_names,
                         field_offsets,
                         field_types,
                         100,
                         &blank_data,
                         0,
                         NULL);

    if (res < 0) exit(-1);

    // Add the details of each field as an attribute
    for (i=0; i<num_fields; ++i) {
        std::stringstream   ss;
        ss << "FIELD_" << i << "_DETAILS";
        res = H5LTset_attribute_string(data_file, ModelStress::hdf5_table_name().c_str(), ss.str().c_str(), field_details[i]);

        if (res < 0) exit(-1);
    }

    // Free memory for HDF5 related data
    for (i=0; i<num_fields; ++i) delete field_names[i];

    delete field_names;

    for (i=0; i<num_fields; ++i) delete field_details[i];

    delete field_details;
    delete field_offsets;
    delete field_types;
    delete field_sizes;
}

void quakelib::ModelStress::append_stress_hdf5(const hid_t &data_file) const {
    std::vector<FieldDesc>                  descs;
    std::vector<StressData>::const_iterator it;
    herr_t                      res;
    unsigned int                i;

    // Set up the section table definition
    descs.clear();
    ModelStress::get_field_descs(descs);
    size_t num_fields = descs.size();
    size_t num_stress_recs = _data.size();
    size_t *field_offsets = new size_t[num_fields];
    size_t *field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_offsets[i] = descs[i].offset;
        field_sizes[i] = descs[i].size;
    }

    // Fill in the data for the sections
    StressData *stress_data = new StressData[num_stress_recs];

    for (i=0,it=_data.begin(); it!=_data.end(); ++i,++it) {
        memcpy(&(stress_data[i]), &(*it), sizeof(StressData));
    }

    // Create the section table
    res = H5TBappend_records(data_file,
                             ModelStress::hdf5_table_name().c_str(),
                             num_stress_recs,
                             sizeof(StressData),
                             field_offsets,
                             field_sizes,
                             stress_data);

    if (res < 0) exit(-1);

    // Free memory for HDF5 related data
    delete stress_data;

    delete field_offsets;
    delete field_sizes;
}
#endif

void quakelib::ModelStress::write_ascii_header(std::ostream &out_stream) {
    std::vector<FieldDesc>                  descs;
    std::vector<FieldDesc>::const_iterator  dit;

    // Write section header
    ModelStress::get_field_descs(descs);

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_stream << "# " << dit->name << ": " << dit->details << "\n";
    }

    out_stream << "# ";

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_stream << dit->name << " ";
    }

    out_stream << "\n";
}

void quakelib::ModelStress::get_field_descs(std::vector<quakelib::FieldDesc> &descs) {
    FieldDesc       field_desc;

    // Stress state table definition
    field_desc.name = "element_id";
    field_desc.details = "ID number of the element corresponding to the stress values.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(StressData, _element_id);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "shear_stress";
    field_desc.details = "Shear stress on the element (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(StressData, _shear_stress);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "normal_stress";
    field_desc.details = "Normal stress on the element (Pascals).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(StressData, _normal_stress);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);
}

// ********************************************************************************************

void quakelib::ModelStressState::write_ascii_header(std::ostream &out_stream) {
    std::vector<FieldDesc>                  descs;
    std::vector<FieldDesc>::const_iterator  dit;

    // Write section header
    ModelStressState::get_field_descs(descs);

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_stream << "# " << dit->name << ": " << dit->details << "\n";
    }

    out_stream << "# ";

    for (dit=descs.begin(); dit!=descs.end(); ++dit) {
        out_stream << dit->name << " ";
    }

    out_stream << "\n";
}

#ifdef HDF5_FOUND
void quakelib::ModelStressState::setup_stress_state_hdf5(const hid_t &data_file) {
    std::vector<FieldDesc>  descs;
    unsigned int            i;
    StressDataTime          blank_data;
    herr_t                  res;

    // Set up the section table definition
    descs.clear();
    ModelStressState::get_field_descs(descs);
    size_t num_fields = descs.size();
    char **field_names = new char *[num_fields];
    char **field_details = new char *[num_fields];
    size_t *field_offsets = new size_t[num_fields];
    hid_t *field_types = new hid_t[num_fields];
    size_t *field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_names[i] = new char[descs[i].name.length()+1];
        strncpy(field_names[i], descs[i].name.c_str(), descs[i].name.length());
        field_names[i][descs[i].name.length()] = '\0';
        field_details[i] = new char[descs[i].details.length()+1];
        strncpy(field_details[i], descs[i].details.c_str(), descs[i].details.length());
        field_details[i][descs[i].details.length()] = '\0';
        field_offsets[i] = descs[i].offset;
        field_types[i] = descs[i].type;
        field_sizes[i] = descs[i].size;
    }

    blank_data._year = std::numeric_limits<float>::quiet_NaN();
    blank_data._event_num = UNDEFINED_EVENT_ID;
    blank_data._sweep_num = UNDEFINED_EVENT_ID;
    blank_data._end_rec = UNDEFINED_ELEMENT_ID;
    blank_data._start_rec = UNDEFINED_ELEMENT_ID;

    // Create the sweep table
    res = H5TBmake_table("Stress State Table",
                         data_file,
                         ModelStressState::hdf5_table_name().c_str(),
                         num_fields,
                         0,
                         sizeof(StressDataTime),
                         (const char **)field_names,
                         field_offsets,
                         field_types,
                         100,
                         &blank_data,
                         0,
                         NULL);

    if (res < 0) exit(-1);

    // Add the details of each field as an attribute
    for (i=0; i<num_fields; ++i) {
        std::stringstream   ss;
        ss << "FIELD_" << i << "_DETAILS";
        res = H5LTset_attribute_string(data_file, ModelStressState::hdf5_table_name().c_str(), ss.str().c_str(), field_details[i]);

        if (res < 0) exit(-1);
    }

    // Free memory for HDF5 related data
    for (i=0; i<num_fields; ++i) delete field_names[i];

    delete field_names;

    for (i=0; i<num_fields; ++i) delete field_details[i];

    delete field_details;
    delete field_offsets;
    delete field_types;
    delete field_sizes;
}

void quakelib::ModelStressState::append_stress_state_hdf5(const hid_t &data_file) const {
    std::vector<FieldDesc>                  descs;
    std::vector<StressData>::const_iterator it;
    herr_t                      res;
    unsigned int                i;

    // Set up the section table definition
    descs.clear();
    ModelStressState::get_field_descs(descs);
    size_t num_fields = descs.size();
    size_t *field_offsets = new size_t[num_fields];
    size_t *field_sizes = new size_t[num_fields];

    for (i=0; i<num_fields; ++i) {
        field_offsets[i] = descs[i].offset;
        field_sizes[i] = descs[i].size;
    }

    // Add to the stress state table
    res = H5TBappend_records(data_file,
                             ModelStressState::hdf5_table_name().c_str(),
                             1,
                             sizeof(StressDataTime),
                             field_offsets,
                             field_sizes,
                             &_times);

    if (res < 0) exit(-1);

    // Free memory for HDF5 related data
    delete field_offsets;
    delete field_sizes;
}
#endif

void quakelib::ModelStressState::read_ascii(std::istream &in_stream) {
    std::stringstream   ss(next_line(in_stream));
    ss >> _times._year;
    ss >> _times._event_num;
    ss >> _times._sweep_num;
    ss >> _times._start_rec;
    ss >> _times._end_rec;
}

void quakelib::ModelStressState::write_ascii(std::ostream &out_stream) const {
    out_stream << _times._year << " ";
    out_stream << _times._event_num << " ";
    out_stream << _times._sweep_num << " ";
    out_stream << _times._start_rec << " ";
    out_stream << _times._end_rec;

    next_line(out_stream);
}

void quakelib::ModelStressState::get_field_descs(std::vector<quakelib::FieldDesc> &descs) {
    FieldDesc       field_desc;

    // Stress state table definition
    field_desc.name = "year";
    field_desc.details = "Simulation time this stress state was recorded (years).";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(StressDataTime, _year);
    field_desc.type = H5T_NATIVE_FLOAT;
    field_desc.size = sizeof(float);
#endif
    descs.push_back(field_desc);

    field_desc.name = "event_num";
    field_desc.details = "Event number this stress state corresponds to.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(StressDataTime, _event_num);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "sweep_num";
    field_desc.details = "Sweep number this stress state corresponds to.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(StressDataTime, _sweep_num);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "start_rec";
    field_desc.details = "Starting record of stress values for this time.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(StressDataTime, _start_rec);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);

    field_desc.name = "end_rec";
    field_desc.details = "Ending record of stress values for this time.";
#ifdef HDF5_FOUND
    field_desc.offset = HOFFSET(StressDataTime, _end_rec);
    field_desc.type = H5T_NATIVE_UINT;
    field_desc.size = sizeof(unsigned int);
#endif
    descs.push_back(field_desc);
}

int quakelib::ModelStressSet::read_file_ascii(const std::string &stress_index_file_name, const std::string &stress_file_name) {
    std::ifstream   stress_ind_file, stress_file;
    ModelSweeps     file_sweeps;

    // Try to open the stress index file
    stress_ind_file.open(stress_index_file_name.c_str());

    if (!stress_ind_file.is_open()) return -1;

    // Try to open the stress file
    stress_file.open(stress_file_name.c_str());

    if (!stress_file.is_open()) return -1;

    // Keep going until we hit the end of either file
    while (!stress_ind_file.eof() && !stress_file.eof()) {
        ModelStressState    new_stress_state;
        ModelStress         new_stresses;
        new_stress_state.read_ascii(stress_ind_file);
        unsigned int num_recs = new_stress_state.getNumStressRecords();
        new_stresses.read_ascii(stress_file, num_recs);
        new_stress_state.setStresses(new_stresses);

        if (!stress_ind_file.eof() && !stress_file.eof()) _states.push_back(new_stress_state);
    }

    // Close the files
    stress_ind_file.close();
    stress_file.close();

    return 0;
}

namespace quakelib {
    std::ostream &operator<<(std::ostream &os, const ModelSweeps &ms) {
        os << "SWEEPS(" << ms._sweeps.size() << ")";
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const ModelEvent &me) {
        os << me._data._event_number << " " << me._data._event_year;
        return os;
    }
}
