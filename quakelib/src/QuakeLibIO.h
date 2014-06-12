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

#include "QuakeLibUtil.h"
#include "QuakeLib.h"

#ifdef QUAKELIB_HAVE_STRING_H
#include <string.h>
#endif

#ifdef HDF5_FOUND
#include "hdf5.h"
#include "hdf5_hl.h"
#endif

#include <map>

#ifndef _QUAKELIB_IO_H_
#define _QUAKELIB_IO_H_

namespace quakelib {
    typedef unsigned int UIndex;
    static const UIndex INVALID_INDEX = std::numeric_limits<unsigned int>::max();
    static const unsigned int NAME_MAX_LEN = 256;

    // Block info related definitions
#define VERTICES_TABLE_HDF5             "vertices"
#define ELEMENTS_TABLE_HDF5             "elements"

    class ModelIO {
        private:
            std::string         _comment;

        protected:
            std::string comment(void) const {
                return _comment;
            };
            void set_comment(std::string new_comment) {
                _comment = new_comment;
            };
            std::string next_line(std::istream &in_stream);
            void next_line(std::ostream &out_stream) const;
    };

    struct FieldDesc {
        std::string name;
        std::string details;
#ifdef HDF5_FOUND
        size_t      offset;
        hid_t       type;
        size_t      size;
#endif
    };

    typedef struct FieldDesc FieldDesc;

    typedef std::map<UIndex, UIndex> IndexRemap;
    class ModelRemapping {
        private:
            IndexRemap    section_remap;
            IndexRemap    element_remap;
            IndexRemap    vertex_remap;

        public:
            void remap_section(const UIndex &orig_id, const UIndex &new_id) {
                section_remap.insert(std::make_pair(orig_id, new_id));
            };
            void remap_element(const UIndex &orig_id, const UIndex &new_id) {
                element_remap.insert(std::make_pair(orig_id, new_id));
            };
            void remap_vertex(const UIndex &orig_id, const UIndex &new_id) {
                vertex_remap.insert(std::make_pair(orig_id, new_id));
            };

            // Returns index that the specified section ID will be remapped to, INVALID_INDEX if the section ID is not in this remapping
            UIndex get_section_id(const UIndex &orig_id) const;
            UIndex get_element_id(const UIndex &orig_id) const;
            UIndex get_vertex_id(const UIndex &orig_id) const;
    };

    struct VertexData {
        UIndex  _id;
        float   _lat, _lon, _alt;
        float   _das;
    };

    class ModelVertex : public ModelIO {
        private:
            VertexData      _data;
            Vec<3>          _pos;

        public:
            ModelVertex(void) {
                _data._id = INVALID_INDEX;
                _data._lat = _data._lon = _data._alt = std::numeric_limits<float>::quiet_NaN();
                _pos = Vec<3>();
                _data._das = std::numeric_limits<float>::quiet_NaN();
            };

            VertexData data(void) const {
                return _data;
            };

            UIndex id(void) const {
                return _data._id;
            };
            void set_id(const UIndex &id) {
                _data._id = id;
            };

            LatLonDepth lld(void) const {
                return LatLonDepth(_data._lat, _data._lon, _data._alt);
            };
            void set_lld(const LatLonDepth &lld, const LatLonDepth &base) {
                Conversion c(base);
                Vec<3> xyz = c.convert2xyz(lld);
                _data._lat = lld.lat();
                _data._lon = lld.lon();
                _data._alt = lld.altitude();
                _pos = xyz;
            };

            Vec<3> xyz(void) const {
                return _pos;
            };
            void set_xyz(const Vec<3> &new_xyz, const LatLonDepth &base) {
                Conversion c(base);
                LatLonDepth lld = c.convert2LatLon(new_xyz);
                _pos = new_xyz;
                _data._lat = lld.lat();
                _data._lon = lld.lon();
                _data._alt = lld.altitude();
            };

            float das(void) const {
                return _data._das;
            };
            void set_das(const float &das) {
                _data._das = das;
            };

#ifdef HDF5_FOUND
            static std::string hdf5_table_name(void) {
                return "vertices";
            };
#endif
            static void get_field_descs(std::vector<FieldDesc> &descs);

            void apply_remap(const ModelRemapping &remap);
            // Recalculate the Cartesian location of this vertex based on the new base coordinate
            //void change_base(const LatLonDepth &new_base) { Conversion c(new_base); set_xyz(c.convert2xyz(lld())); };

            void read_data(const VertexData &in_data);
            void write_data(VertexData &out_data) const;

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;
    };

    struct ElementData {
        UIndex              _id;
        UIndex              _section_id;
        UIndex              _vertices[3];
        bool                _is_quad;
        float               _slip_rate;
        float               _aseismic;
        float               _rake;
        float               _lame_mu;
        float               _lame_lambda;
        float               _max_slip;
    };

    class ModelElement : public ModelIO {
        private:
            ElementData         _data;

        public:
            ModelElement(void) {
                _data._id = _data._section_id = INVALID_INDEX;

                for (unsigned int i=0; i<3; ++i) _data._vertices[i] = INVALID_INDEX;

                _data._is_quad = false;
                _data._slip_rate = _data._aseismic = _data._rake = std::numeric_limits<float>::quiet_NaN();
                _data._lame_mu = _data._lame_lambda = _data._max_slip = std::numeric_limits<float>::quiet_NaN();
            };

            ElementData data(void) const {
                return _data;
            };

            UIndex id(void) const {
                return _data._id;
            };
            void set_id(const UIndex &id) {
                _data._id = id;
            };

            UIndex section_id(void) const {
                return _data._section_id;
            };
            void set_section_id(const UIndex &section_id) {
                _data._section_id = section_id;
            };

            UIndex vertex(const unsigned int &v) const {
                assert(v<3);
                return _data._vertices[v];
            };
            void set_vertex(const unsigned int &v, const UIndex &ind) {
                assert(v<3);
                _data._vertices[v] = ind;
            };

            bool is_quad(void) const {
                return _data._is_quad;
            };
            void set_is_quad(const bool &is_quad) {
                _data._is_quad = is_quad;
            };

            float slip_rate(void) const {
                return _data._slip_rate;
            };
            void set_slip_rate(const float &slip_rate) {
                _data._slip_rate = slip_rate;
            };

            float aseismic(void) const {
                return _data._aseismic;
            };
            void set_aseismic(const float &aseismic) {
                _data._aseismic = aseismic;
            };

            float rake(void) const {
                return _data._rake;
            };
            void set_rake(const float &rake) {
                _data._rake = rake;
            };

            float lame_mu(void) const {
                return _data._lame_mu;
            };
            void set_lame_mu(const float &lame_mu) {
                _data._lame_mu = lame_mu;
            };

            float lame_lambda(void) const {
                return _data._lame_lambda;
            };
            void set_lame_lambda(const float &lame_lambda) {
                _data._lame_lambda = lame_lambda;
            };
        
            float max_slip(void) const {
                return _data._max_slip;
            };
            void set_max_slip(const float &max_slip) {
                _data._max_slip = max_slip;
            };

#ifdef HDF5_FOUND
            static std::string hdf5_table_name(void) {
                return "elements";
            };
#endif
            static void get_field_descs(std::vector<FieldDesc> &descs);

            void apply_remap(const ModelRemapping &remap);

            void read_data(const ElementData &in_data);
            void write_data(ElementData &out_data) const;

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;
    };

    struct SectionData {
        UIndex              _id;
        UIndex              _fault_id;
        char                _name[NAME_MAX_LEN];
    };

    class ModelSection : public ModelIO {
        private:
            SectionData         _data;

        public:
            ModelSection(void) {
                _data._id = INVALID_INDEX;
                strncpy(_data._name, "(invalid)", NAME_MAX_LEN);
                _data._name[NAME_MAX_LEN-1] = '\0';
            };

            SectionData data(void) const {
                return _data;
            };

            UIndex id(void) const {
                return _data._id;
            };
            void set_id(const UIndex &id) {
                _data._id = id;
            };

            UIndex fault_id(void) const {
                return _data._fault_id;
            };
            void set_fault_id(const UIndex &fault_id) {
                _data._fault_id = fault_id;
            };

            std::string name(void) const {
                return _data._name;
            };
            void set_name(const std::string &name) {
                strncpy(_data._name, name.c_str(), NAME_MAX_LEN);
                _data._name[NAME_MAX_LEN-1] = '\0';
            };

#ifdef HDF5_FOUND
            static std::string hdf5_table_name(void) {
                return "sections";
            };
#endif
            static void get_field_descs(std::vector<FieldDesc> &descs);

            void apply_remap(const ModelRemapping &remap);

            void read_data(const SectionData &in_data);
            void write_data(SectionData &out_data) const;

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;

            // TODO: write these functions
            void write_kml_label(std::ostream &out_stream) const {};
            void write_kml_geometry(std::ostream &out_stream) const {};
    };

    class FaultTracePoint : public ModelIO {
        private:
            LatLonDepth         _pos;
            float               _depth_along_dip;
            float               _slip_rate, _aseismic, _rake, _dip;
            float               _lame_mu, _lame_lambda;
        public:
            // Create a trace point by reading an input stream
            FaultTracePoint(std::istream &in_stream) {
                read_ascii(in_stream);
            };

            // Create a trace point on a fault. By default, this will be a right-lateral
            // strike slip fault with standard earth-like material properties with long-term
            // slip of 1 cm/year
            FaultTracePoint(const LatLonDepth &pos, const float &depth_along_dip=12000,
                            const float &slip_rate=(0.01/(86400*365.25)),
                            const float &aseismic=0, const float &rake=180, const float &dip=90,
                            const float &lame_mu=3.0e10, const float &lame_lambda=3.2e10) :
                _pos(pos), _depth_along_dip(depth_along_dip), _slip_rate(slip_rate), _aseismic(aseismic),
                _rake(rake), _dip(dip), _lame_mu(lame_mu), _lame_lambda(lame_lambda) {};

            LatLonDepth pos(void) const {
                return _pos;
            };
            float depth_along_dip(void) const {
                return _depth_along_dip;
            };
            float slip_rate(void) const {
                return _slip_rate;
            };
            float aseismic(void) const {
                return _aseismic;
            };
            float rake(void) const {
                return _rake;
            };
            float dip(void) const {
                return _dip;
            };
            float lame_mu(void) const {
                return _lame_mu;
            };
            float lame_lambda(void) const {
                return _lame_lambda;
            };

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;
    };

    class fiterator {
        private:
            std::map<UIndex, ModelSection>              *_map;
            std::map<UIndex, ModelSection>::iterator    _it;

        public:
            fiterator(void) : _map(NULL) {};
            fiterator(std::map<UIndex, ModelSection> *map, std::map<UIndex, ModelSection>::iterator start) : _map(map), _it(start) {};
            fiterator &operator=(const fiterator &other) {
                _map = other._map;
                _it = other._it;
                return *this;
            };
            bool operator==(const fiterator &other) {
                return (_map == other._map && _it == other._it);
            };
            bool operator!=(const fiterator &other) {
                return (_map != other._map || _it != other._it);
            };
            fiterator &operator++(void) {
                if (_map && _it != _map->end()) _it++;

                return *this;
            };
            //fiterator &operator++(int) { fiterator tmp(*this); ++(*this); return tmp; };
            ModelSection &operator*(void) {
                return _it->second;
            };
            ModelSection *operator->(void) {
                return (&*(fiterator)*this);
            };
    };

    class eiterator {
        private:
            std::map<UIndex, ModelElement>              *_map;
            std::map<UIndex, ModelElement>::iterator    _it;
            UIndex                                      _fid;

        public:
            eiterator(void) : _map(NULL), _fid(INVALID_INDEX) {};
            eiterator(std::map<UIndex, ModelElement> *map, std::map<UIndex, ModelElement>::iterator start, const UIndex &fid) : _map(map), _it(start), _fid(fid) {
                if (_map && _fid != INVALID_INDEX) {
                    while (_it != _map->end() && _it->second.section_id() != _fid) {
                        _it++;
                    }
                }
            };
            eiterator &operator=(const eiterator &other) {
                _map = other._map;
                _it = other._it;
                _fid = other._fid;
                return *this;
            };
            bool operator==(const eiterator &other) {
                return (_map == other._map && _it == other._it && _fid == other._fid);
            };
            bool operator!=(const eiterator &other) {
                return (_map != other._map || _it != other._it || _fid != other._fid);
            };
            eiterator &operator++(void) {
                if (_map && _it != _map->end()) {
                    if (_fid == INVALID_INDEX) {
                        _it++;
                    } else {
                        do {
                            _it++;
                        } while (_it != _map->end() && _it->second.section_id() != _fid);
                    }
                }

                return *this;
            };
            //eiterator &operator++(int) { eiterator tmp(*this); ++(*this); return tmp; };
            ModelElement &operator*(void) {
                return _it->second;
            };
            ModelElement *operator->(void) {
                return (&*(eiterator)*this);
            };
    };

    class ModelWorld : public ModelIO {
        private:
            std::map<UIndex, ModelVertex>   _vertices;
            std::map<UIndex, ModelElement>  _elements;
            std::map<UIndex, ModelSection>  _sections;

#ifdef HDF5_FOUND
            void read_section_hdf5(const hid_t &data_file);
            void read_element_hdf5(const hid_t &data_file);
            void read_vertex_hdf5(const hid_t &data_file);

            void write_section_hdf5(const hid_t &data_file) const;
            void write_element_hdf5(const hid_t &data_file) const;
            void write_vertex_hdf5(const hid_t &data_file) const;
#endif

        public:
            ModelSection &new_section(void);
            ModelElement &new_element(void);
            ModelVertex &new_vertex(void);

            ModelSection &section(const UIndex &ind) throw(std::domain_error);
            ModelElement &element(const UIndex &ind) throw(std::domain_error);
            ModelVertex &vertex(const UIndex &ind) throw(std::domain_error);

            fiterator begin_section(void) {
                return fiterator(&_sections, _sections.begin());
            };
            fiterator end_section(void) {
                return fiterator(&_sections, _sections.end());
            };

            eiterator begin_element(const UIndex &fid=INVALID_INDEX) {
                return eiterator(&_elements, _elements.begin(), fid);
            };
            eiterator end_element(const UIndex &fid=INVALID_INDEX) {
                return eiterator(&_elements, _elements.end(), fid);
            };

            UIndex next_section_index(void) const {
                if (_sections.size()) return _sections.rbegin()->first+1;
                else return 0;
            };
            UIndex next_element_index(void) const {
                if (_elements.size()) return _elements.rbegin()->first+1;
                else return 0;
            };
            UIndex next_vertex_index(void) const {
                if (_vertices.size()) return _vertices.rbegin()->first+1;
                else return 0;
            };

            size_t num_sections(void) const;
            size_t num_faults(void) const;
            size_t num_elements(const UIndex &fid=INVALID_INDEX) const;
            size_t num_vertices(const UIndex &fid=INVALID_INDEX) const;

            LatLonDepth min_bound(const UIndex &fid=INVALID_INDEX) const;
            LatLonDepth max_bound(const UIndex &fid=INVALID_INDEX) const;

            void insert(const ModelWorld &other_world);
            void insert(const ModelSection &new_section);
            void insert(const ModelElement &new_element);
            void insert(const ModelVertex &new_vertex);

            void get_bounds(LatLonDepth &minimum, LatLonDepth &maximum) const;

            SimElement create_sim_element(const UIndex &element_id) const;

            bool overwrite(const ModelRemapping &remap);
            void apply_remap(const ModelRemapping &remap);
            ModelRemapping remap_indices_contiguous(const UIndex &start_section_index,
                                                    const UIndex &start_element_index,
                                                    const UIndex &start_vertex_index) const;
            ModelRemapping remove_duplicate_vertices_remap(void) const;
            void clear(void);

            void reset_base_coord(const LatLonDepth &new_base);

            void create_section(std::vector<unsigned int> &unused_trace_segments,
                                const std::vector<FaultTracePoint> &trace,
                                const LatLonDepth &base_coord,
                                const UIndex &fault_id,
                                const float &element_size,
                                const std::string &section_name,
                                const std::string &taper_method);

            int read_file_ascii(const std::string &file_name);
            int write_file_ascii(const std::string &file_name) const;

            int read_file_trace_latlon(std::vector<unsigned int> &unused_trace_segments, const std::string &file_name, const float &elem_size, const std::string &taper_method);
            int write_file_trace_latlon(const std::string &file_name, const float &depth_along_dip);

            int read_file_hdf5(const std::string &file_name);
            int write_file_hdf5(const std::string &file_name) const;

            int write_file_kml(const std::string &file_name);

            int read_files_eqsim(const std::string &geom_file_name, const std::string &cond_file_name, const std::string &fric_file_name);
            int write_files_eqsim(const std::string &geom_file_name, const std::string &cond_file_name, const std::string &fric_file_name);
    };
}

#endif
