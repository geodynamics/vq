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
    field_desc.details = "Lame's parameter describing the shear modulus of the material for this element (Units?).";
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

void quakelib::ModelWorld::create_section(std::vector<unsigned int> &unused_trace_segments, const std::vector<FaultTracePoint> &trace, const LatLonDepth &base_coord, const UIndex &fault_id, const float &element_size, const std::string &section_name, const std::string &taper_method) {
    Vec<3>              cur_trace_point, next_trace_point, element_start, element_end, element_step_vec, vert_step;
    double              t, inter_t;
    std::vector<UIndex> elem_ids;
    double              x_0, x_1, x_s, y_0, y_1, y_s, l;
    double              elem_depth, elem_slip_rate, elem_aseismic;
    double              elem_rake, elem_dip;
    double              elem_lame_mu, elem_lame_lambda;
    unsigned int        num_vert_elems, ve;
    double              total_trace_length, taper_t;
    Conversion          conv(base_coord);

    if (element_size <= 0) return;

    if (trace.size() == 0 || trace.size() == 1) return;

    ModelSection &section = new_section();
    section.set_name(section_name);
    section.set_fault_id(fault_id);

    // Determine the total length of the trace
    total_trace_length = 0;
    cur_trace_point = conv.convert2xyz(trace.at(0).pos());

    for (unsigned int i=1; i<trace.size(); ++i) {
        next_trace_point = conv.convert2xyz(trace.at(i).pos());
        total_trace_length += (next_trace_point-cur_trace_point).mag();
        cur_trace_point = next_trace_point;
    }

    double dist_along_strike = 0;
    double dist_along_trace = 0;
    double taper_flow = 0;
    double taper_full = 0;
    double fault_area = 0;
    cur_trace_point = element_start = conv.convert2xyz(trace.at(0).pos());

    for (unsigned int i=1; i<trace.size(); ++i) {
        next_trace_point = conv.convert2xyz(trace.at(i).pos());
        t = 0;
        bool used_trace_segment = false;

        while (t <= 1) {
            // Treat the line between the current two trace points as a parameterized equation with t in [0,1]
            // We then want to solve for t in the equation: element_size = dist(element_start, line(t))
            // If t < 0 or t > 1, no element will fit using this trace segment so we move to the next one
            x_0 = cur_trace_point[0];
            x_1 = next_trace_point[0];
            x_s = element_start[0];
            y_0 = cur_trace_point[1];
            y_1 = next_trace_point[1];
            y_s = element_start[1];
            l = element_size;
            // Horribly ugly equation to find the next point on the parameterized line
            t = fabs((sqrt(pow(-2*x_0*x_0+2*x_0*x_1+2*x_0*x_s-2*x_1*x_s-2*y_0*y_0+2*y_0*y_1+2*y_0*y_s-2*y_1*y_s, 2)-
                           4*(x_0*x_0-2*x_0*x_1+x_1*x_1+y_0*y_0-2*y_0*y_1+y_1*y_1)*(x_0*x_0-2*x_0*x_s+x_s*x_s+y_0*y_0-2*y_0*y_s+y_s*y_s-l*l))+
                      2*x_0*x_0-2*x_0*x_1-2*x_0*x_s+2*x_1*x_s+2*y_0*y_0-2*y_0*y_1-2*y_0*y_s+2*y_1*y_s)/(2*(x_0*x_0-2*x_0*x_1+x_1*x_1+y_0*y_0-2*y_0*y_1+y_1*y_1)));

            if (t > 1) break;

            used_trace_segment = true;
            element_end = cur_trace_point*(1-t) + next_trace_point*t;
            element_step_vec = element_end-element_start;

            // Use the t value in the middle of the element for interpolation
            inter_t = t;
            elem_depth = inter_t *trace.at(i).depth_along_dip()+(1.0-inter_t)*trace.at(i-1).depth_along_dip();
            // Warn user if element depth is smaller than element size
            num_vert_elems = floor(elem_depth/element_size);

            if (num_vert_elems == 0) std::cerr << "WARNING: Depth is smaller than element size in trace segment "
                                                   << i << ". Element size may be too big." << std::endl;

            elem_slip_rate = conv.cm_per_yr2m_per_sec(inter_t *trace.at(i).slip_rate()+(1.0-inter_t)*trace.at(i-1).slip_rate());
            elem_aseismic = inter_t *trace.at(i).aseismic()+(1.0-inter_t)*trace.at(i-1).aseismic();
            elem_dip = conv.deg2rad(inter_t *trace.at(i).dip()+(1.0-inter_t)*trace.at(i-1).dip());
            elem_rake = conv.deg2rad(inter_t *trace.at(i).rake()+(1.0-inter_t)*trace.at(i-1).rake());
            elem_lame_mu = inter_t *trace.at(i).lame_mu()+(1.0-inter_t)*trace.at(i-1).lame_mu();
            elem_lame_lambda = inter_t *trace.at(i).lame_lambda()+(1.0-inter_t)*trace.at(i-1).lame_lambda();

            // Set up the vertical step to go along the dip
            vert_step = element_step_vec.rotate_around_axis(Vec<3>(0,0,-1), M_PI/2);
            vert_step = vert_step.rotate_around_axis(element_step_vec, elem_dip);

            for (ve=0; ve<num_vert_elems; ++ve) {
                // Calculate values for element based on interpolation between trace points
                taper_t = 1;

                double cur_dist = dist_along_trace+t*(next_trace_point-cur_trace_point).mag()-0.5*element_size;

                if (taper_method == "taper_full" || taper_method == "taper_renorm") {
                    double x = cur_dist/total_trace_length;
                    double z = (float(ve)+0.5)/num_vert_elems;
                    taper_t *= 4*(x-x*x)*sqrt(1-z);
                } else if (taper_method == "taper") {
                    double inside_dist = (total_trace_length/2.0 - fabs(total_trace_length/2.0-cur_dist));

                    if (inside_dist <= elem_depth) {
                        double x = inside_dist/elem_depth;
                        double z = (float(ve)+0.5)/num_vert_elems;
                        taper_t *= sqrt(x)*sqrt(1-z);
                    }
                }

                taper_flow += taper_t *elem_slip_rate*(element_size*element_size);
                taper_full += elem_slip_rate*(element_size*element_size);

                // Create the new vertices
                ModelVertex &v0 = new_vertex();
                ModelVertex &v1 = new_vertex();
                ModelVertex &v2 = new_vertex();

                // Set xyz for the vertices
                v0.set_xyz(element_start+vert_step*ve, base_coord);
                v1.set_xyz(element_start+vert_step*(ve+1), base_coord);
                v2.set_xyz(element_start+vert_step*ve+element_step_vec, base_coord);

                // Set distance along strike for each vertex
                v0.set_das(dist_along_strike);
                v1.set_das(dist_along_strike);
                v2.set_das(dist_along_strike+element_size);

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

                fault_area += element_size*element_size;
            }

            element_start = element_end;
            dist_along_strike += element_size;
        }

        dist_along_trace += (next_trace_point-cur_trace_point).mag();

        if (!used_trace_segment) unused_trace_segments.push_back(i);

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

int quakelib::ModelWorld::read_file_trace_latlon(std::vector<unsigned int> &unused_trace_segments, const std::string &file_name, const float &elem_size, const std::string &taper_method) {
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

        new_world.create_section(unused_trace_segments, trace_pts, LatLonDepth(min_lat, min_lon), fault_id, elem_size, cur_section_name, taper_method);
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
    fiterator           fit;
    UIndex              sid;
    unsigned int        i;
    double              max_alt;
    bool                element_on_trace;
    Conversion          c;

    out_file.open(file_name.c_str());

    // Write traces by section
    for (fit=begin_section(); fit!=end_section(); ++fit) {
        std::vector<FaultTracePoint>    trace_pts;

        trace_pts.clear();
        sid = fit->id();

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
        out_file << fit->fault_id() << " " << trace_pts.size() << " " << fit->name() << "\n";

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
                         "sections",
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
        res = H5LTset_attribute_string(data_file, "sections", ss.str().c_str(), field_details[i]);

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
                         "elements",
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
        res = H5LTset_attribute_string(data_file, "elements", ss.str().c_str(), field_details[i]);

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

    // TODO: factor this out?
    blank_vertex = ModelVertex().data();

    // Fill in the data for the vertices
    vertex_data = new VertexData[num_vertices];

    for (i=0,vit=_vertices.begin(); vit!=_vertices.end(); ++i,++vit) {
        vit->second.write_data(vertex_data[i]);
    }

    // Create the vertices table
    res = H5TBmake_table("Vertices",
                         data_file,
                         "vertices",
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
        res = H5LTset_attribute_string(data_file, "vertices", ss.str().c_str(), field_details[i]);

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
    quakelib::LatLonDepth           base;

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

    // Triangle elements are currently not supported
    if (geometry_data.num_triangles() > 0) {
        std::cerr << "ERROR: Currently cannot handle EQSim triangle elements. These elements will be ignored." << std::endl;
    }

    // Go through the geometry and create sections for each EQSim section
    for (sit=geometry_data.sections.begin(); sit!=geometry_data.sections.end(); ++sit) {
        eqsim_world.insert(sit->second.create_model_section());

        // Assuming aligned rectangular elements
        for (it=sit->second.rectangles.begin(); it!=sit->second.rectangles.end(); ++it) {
            quakelib::ModelElement  new_element;
            unsigned int            i;

            new_element = it->second.create_model_element();
            new_element.set_section_id(sit->second.sid());

            for (i=0; i<4; ++i) eqsim_world.insert(sit->second.vertices.find(it->second.vertex(i))->second.create_model_vertex(base));

            new_element.set_lame_lambda(friction_data.get_lame_lambda());
            new_element.set_lame_mu(friction_data.get_lame_mu());
            //new_element.set_static_strength(friction_data.get_static_strength(it->first));
            //new_element.set_dynamic_strength(friction_data.get_dynamic_strength(it->first));

            eqsim_world.insert(new_element);
        }
    }

    insert(eqsim_world);

    return 0;
}

int quakelib::ModelWorld::write_files_eqsim(const std::string &geom_file_name, const std::string &cond_file_name, const std::string &fric_file_name) {
    EQSimGeometryWriter     geometry_data;
    EQSimFrictionWriter     friction_data;
    eiterator               eit;
    fiterator               fit;
    UIndex                  fid;
    UIndex                  vind, eind;
    LatLonDepth             base = min_bound();
    Conversion              c(base);

    vind = eind = 1;
    friction_data.set_lame_lambda_mu(3.2e10, 3.0e10);

    for (fit=begin_section(); fit!=end_section(); ++fit) {
        EQSimGeometrySection &section = geometry_data.new_section();

        // Set section properties
        section.set_name(fit->name());
        section.set_fid(fit->fault_id());

        fid = fit->id();

        for (eit=begin_element(fid); eit!=end_element(fid); ++eit) {
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
                out_file << "Slip rate: " << c.m_per_sec2cm_per_yr(eit->second.slip_rate()) << " cm/year\n";
                out_file << "Rake: " << c.rad2deg(eit->second.rake()) << " degrees\n";
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
