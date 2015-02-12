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

#define UNDEFINED_ELEMENT_ID    UINT_MAX
#define UNDEFINED_FAULT_ID      UINT_MAX
#define UNDEFINED_SECTION_ID    UINT_MAX
#define UNDEFINED_EVENT_ID      UINT_MAX

namespace quakelib {
    typedef unsigned int UIndex;
    static const UIndex INVALID_INDEX = std::numeric_limits<unsigned int>::max();
    static const unsigned int NAME_MAX_LEN = 256;

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

            void read_data(const VertexData &in_data);
            void write_data(VertexData &out_data) const;

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;
    };

    struct ElementData {
        UIndex              _id;
        UIndex              _section_id;
        UIndex              _vertices[3];
        unsigned int        _is_quad;
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

    class siterator {
        private:
            std::map<UIndex, ModelSection>              *_map;
            std::map<UIndex, ModelSection>::iterator    _it;

        public:
            siterator(void) : _map(NULL) {};
            siterator(std::map<UIndex, ModelSection> *map, std::map<UIndex, ModelSection>::iterator start) : _map(map), _it(start) {};
            siterator &operator=(const siterator &other) {
                _map = other._map;
                _it = other._it;
                return *this;
            };
            bool operator==(const siterator &other) {
                return (_map == other._map && _it == other._it);
            };
            bool operator!=(const siterator &other) {
                return (_map != other._map || _it != other._it);
            };
            siterator &operator++(void) {
                if (_map && _it != _map->end()) _it++;

                return *this;
            };
            ModelSection &operator*(void) {
                return _it->second;
            };
            ModelSection *operator->(void) {
                return (&*(siterator)*this);
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

            siterator begin_section(void) {
                return siterator(&_sections, _sections.begin());
            };
            siterator end_section(void) {
                return siterator(&_sections, _sections.end());
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
                                const std::string &taper_method,
                                const bool resize_trace_elements);

            int read_file_ascii(const std::string &file_name);
            int write_file_ascii(const std::string &file_name) const;

            int read_file_trace_latlon(std::vector<unsigned int> &unused_trace_segments, const std::string &file_name, const float &elem_size, const std::string &taper_method, const bool resize_trace_elements);
            int write_file_trace_latlon(const std::string &file_name, const float &depth_along_dip);

            int read_file_hdf5(const std::string &file_name);
            int write_file_hdf5(const std::string &file_name) const;

            int write_file_kml(const std::string &file_name);

            int read_files_eqsim(const std::string &geom_file_name, const std::string &cond_file_name, const std::string &fric_file_name);
            int write_files_eqsim(const std::string &geom_file_name, const std::string &cond_file_name, const std::string &fric_file_name);
    };

    typedef std::set<UIndex> ElementIDSet;

    // Class recording data associated with a block that slipped during an event
    // mu is static, but we retain it for use in calculating magnitude.
    struct SweepData {
        UIndex      _event_number;
        UIndex      _sweep_number;
        UIndex      _element_id;
        double      _slip, _area, _mu;
        double      _shear_init, _shear_final;
        double      _normal_init, _normal_final;
    };

    /*!
     An event consists of one or more sweeps. Each sweep is a list of failed blocks and
     the slip they experienced. This gives the order that the failure propagated
     through the fault(s).
     */
    class ModelSweeps : public ModelIO {
        private:
            //! The number of this sweep and the associated event number
            UIndex      _event_number;

            //! All the sweeps associated with a single event
            std::vector<SweepData>      _sweeps;

            //! Map of sweep/element to the location in the _sweeps vector
            std::map<std::pair<UIndex, UIndex>, unsigned int>    _rel;

            // Get the index into _sweeps for the given sweep/element combination
            // If it doesn't exist, create it
            unsigned int sweepElementPos(const UIndex &sweep, const UIndex &elem) {
                std::pair<UIndex, UIndex>       se = std::make_pair(sweep, elem);

                if (_rel.count(se) == 0) {
                    SweepData new_data;
                    new_data._event_number = _event_number;
                    new_data._sweep_number = sweep;
                    new_data._element_id = elem;
                    _rel.insert(std::make_pair(se, _sweeps.size()));
                    _sweeps.push_back(new_data);
                }

                return _rel.at(se);
            }

        public:
            ModelSweeps(void) : _event_number(UNDEFINED_EVENT_ID) {};

            typedef std::vector<SweepData>::iterator       iterator;
            typedef std::vector<SweepData>::const_iterator const_iterator;

            iterator begin(void) {
                return _sweeps.begin();
            };
            iterator end(void) {
                return _sweeps.end();
            };
            const_iterator begin(void) const {
                return _sweeps.begin();
            };
            const_iterator end(void) const {
                return _sweeps.end();
            };

            void setSlipAndArea(const UIndex &sweep_number, const UIndex &element_id, const double &slip, const double &area, const double &mu) {
                unsigned int    pos = sweepElementPos(sweep_number, element_id);
                _sweeps[pos]._slip = slip;
                _sweeps[pos]._area = area;
                _sweeps[pos]._mu = mu;
            };

            void setInitStresses(const UIndex &sweep_number, const UIndex &element_id, const double &shear_init, const double &normal_init) {
                unsigned int    pos = sweepElementPos(sweep_number, element_id);
                _sweeps[pos]._shear_init = shear_init;
                _sweeps[pos]._normal_init = normal_init;
            }

            void setFinalStresses(const UIndex &sweep_number, const UIndex &element_id, const double &shear_final, const double &normal_final) {
                unsigned int    pos = sweepElementPos(sweep_number, element_id);
                _sweeps[pos]._shear_final = shear_final;
                _sweeps[pos]._normal_final = normal_final;
            }

            ElementIDSet getInvolvedElements(void) const {
                ElementIDSet element_id_set;

                for (std::vector<SweepData>::const_iterator it=_sweeps.begin(); it!=_sweeps.end(); ++it) {
                    element_id_set.insert(it->_element_id);
                }

                return element_id_set;
            }

            void clear(void) {
                _event_number = UNDEFINED_EVENT_ID;
                _sweeps.clear();
                _rel.clear();
            };
            unsigned int size(void) const {
                return _sweeps.size();
            };

#ifdef HDF5_FOUND
            static std::string hdf5_table_name(void) {
                return "sweeps";
            };
            static void setup_sweeps_hdf5(const hid_t &data_file);
            void append_sweeps_hdf5(const hid_t &data_file) const;
#endif
            static void get_field_descs(std::vector<FieldDesc> &descs);
            static void write_ascii_header(std::ostream &out_stream);

            void read_data(const SweepData &in_data);
            void write_data(SweepData &out_data) const;

            void read_ascii(std::istream &in_stream, const unsigned int num_records);
            void write_ascii(std::ostream &out_stream) const;

            SweepData &operator[](const unsigned int ind) throw(std::out_of_range) {
                if (ind >= _sweeps.size()) throw std::out_of_range("ModelSweeps[]");

                return _sweeps[ind];
            };

            const SweepData &operator[](const unsigned int ind) const throw(std::out_of_range) {
                if (ind >= _sweeps.size()) throw std::out_of_range("ModelSweeps[]");

                return _sweeps[ind];
            };

            friend std::ostream &operator<<(std::ostream &os, const ModelSweeps &ms);
    };

    typedef std::map<UIndex, SweepData> EventElementMap;

    struct EventData {
        //! The event number is used as a simple unique identifier for each event
        unsigned int    _event_number;

        //! The year the event occurred
        double          _event_year;

        //! The magnitude of the event
        double          _event_magnitude;

        //! Which element triggered the event through a static friction failure
        UIndex          _event_trigger;

        //! Initial and final sum of shear stress on all elements involved in the event
        double          _shear_stress_init, _shear_stress_final;

        //! Initial and final sum of normal stress on all elements involved in the event
        double          _normal_stress_init, _normal_stress_final;

        //! Start and end record of the sweeps comprising this event
        unsigned int    _start_sweep_rec, _end_sweep_rec;
    };

    /*!
     A ModelEvent is a rupture of one or more blocks on one or more faults in a model.
     This consists of origin information (trigger block, year) and information about
     how the event propagated through the system (event sweeps).
     */
    class ModelEvent : public ModelIO {
        private:
            EventData       _data;

            //! Whether the event trigger element is on this node or not (used in parallel simulation)
            bool            _event_trigger_on_this_node;

            //! The set of event sweeps
            ModelSweeps     _sweeps;

            //! Sum of the sweep slips, used to quickly calculate magnitude
            EventElementMap _total_slip;

            // TODO: add a flag indicating whether this is an aftershock?

        public:
            ModelEvent(void) {
                clear();
            }
            typedef EventElementMap::iterator iterator;
            typedef EventElementMap::const_iterator const_iterator;

            void setEventStresses(const double &init_shear, const double &final_shear,
                                  const double &init_normal, const double &final_normal) {
                _data._shear_stress_init = init_shear;
                _data._shear_stress_final = final_shear;
                _data._normal_stress_init = init_normal;
                _data._normal_stress_final = final_normal;
            }

            const ModelSweeps &getSweeps(void) const {
                return _sweeps;
            };

            void setStartEndSweep(const unsigned int start_sweep, const unsigned int end_sweep) {
                _data._start_sweep_rec = start_sweep;
                _data._end_sweep_rec = end_sweep;
            }
            void getStartEndSweep(unsigned int &start_sweep, unsigned int &end_sweep) const {
                start_sweep = _data._start_sweep_rec;
                end_sweep = _data._end_sweep_rec;
            }


            unsigned int getNumRecordedSweeps(void) const {
                return _data._end_sweep_rec - _data._start_sweep_rec;
            };
            double getShearStressInit(void) {
                return _data._shear_stress_init;
            };
            double getShearStressFinal(void) {
                return _data._shear_stress_final;
            };
            double getNormalStressInit(void) {
                return _data._normal_stress_init;
            };
            double getNormalStressFinal(void) {
                return _data._normal_stress_final;
            };

            //! Set the unique identifier number for this event
            void setEventNumber(const unsigned int &event_number) {
                _data._event_number = event_number;
            };
            //! Get the unique identifier number for this event
            unsigned int getEventNumber(void) const {
                return _data._event_number;
            };

            //! Set the year this event occurred
            void setEventYear(const double &event_year) {
                _data._event_year = event_year;
            };
            //! Get the year this event occurred
            double getEventYear(void) const {
                return _data._event_year;
            };

            //! Set the block that triggered this event
            void setEventTrigger(const UIndex &event_trigger) {
                _data._event_trigger = event_trigger;
            };
            //! Get the block that triggered this event
            UIndex getEventTrigger(void) const {
                return _data._event_trigger;
            };

            //! Set whether this event occurred on this node
            void setEventTriggerOnThisNode(const bool &totn) {
                _event_trigger_on_this_node = totn;
            };
            //! Get whether this event occurred on this node
            bool getEventTriggerOnThisNode(void) const {
                return _event_trigger_on_this_node;
            };

            //! Set the sweeps for this event to be those specified in the argument.
            //! Also calculate relevant values related to these sweeps.
            void setSweeps(const ModelSweeps &sweeps) {
                ModelSweeps::iterator   it;
                double                  moment = 0;

                _sweeps = sweeps;
                _total_slip.clear();

                // Sum up sweep information for total_slip records
                // and calculate magnitude at the same time
                for (it=_sweeps.begin(); it!=_sweeps.end(); ++it) {
                    it->_event_number = _data._event_number;
                    _total_slip[it->_element_id]._slip += it->_slip;
                    _total_slip[it->_element_id]._area = it->_area;
                    _total_slip[it->_element_id]._mu = it->_mu;
                    moment += it->_slip*it->_mu*it->_area;
                }

                _data._event_magnitude = (2.0/3.0)*log10(1e7*moment) - 10.7;
            }

            //! Get the total amount a given block slipped during this event
            double getEventSlip(const UIndex &element_id) const {
                return (_total_slip.count(element_id) > 0 ? _total_slip.at(element_id)._slip : 0);
            }

            //! Get the area of a given block that slipped during this event
            double getEventArea(const UIndex &element_id) const {
                return (_total_slip.count(element_id) > 0 ? _total_slip.at(element_id)._area : 0);
            }

            //! Get the value of Mu for a given block
            double getBlockMu(const UIndex &element_id) const {
                return (_total_slip.count(element_id) > 0 ? _total_slip.at(element_id)._mu : 0);
            }

            //! Get a set of block IDs for all the blocks that failed in this event.
            ElementIDSet getInvolvedElements(void) const {
                return _sweeps.getInvolvedElements();
            }

            //! Get the magnitude of the earthquake in this event based on the set of specified blocks.
            double getMagnitude(const ElementIDSet &involved_elements) const {
                ElementIDSet::const_iterator    it;
                double                          moment;

                for (moment=0,it=involved_elements.begin(); it!=involved_elements.end(); ++it) {
                    moment += getBlockMu(*it)*getEventSlip(*it)*getEventArea(*it);
                }

                return (2.0/3.0)*log10(1e7*moment) - 10.7;
            }

            //! Get the magnitude of the earthquake in this event.
            double getMagnitude(void) const {
                ElementIDSet  block_id_set = getInvolvedElements();
                return getMagnitude(block_id_set);
            }

            //! Get the total number of blocks that failed in this event.
            unsigned int size(void) const {
                return _total_slip.size();
            };

            //! Calculate the total slipped area in the event
            double calcEventRuptureArea(void) const {
                EventElementMap::const_iterator     it;
                double                              rupture_area = 0;

                for (it=_total_slip.begin(); it!=_total_slip.end(); ++it) rupture_area += it->second._area;

                return rupture_area;
            }

            //! Calculate mean slip weighted by rupture area
            double calcMeanSlip(void) const {
                EventElementMap::const_iterator     it;
                double                              sum_slip_area = 0, sum_area = 0;

                for (it=_total_slip.begin(); it!=_total_slip.end(); ++it) {
                    sum_slip_area += it->second._area*it->second._slip;
                    sum_area += it->second._area;
                }

                return sum_slip_area/sum_area;
            }

            void clear(void) {
                _event_trigger_on_this_node = false;

                _data._event_number = UNDEFINED_EVENT_ID;
                _data._event_year = _data._event_magnitude = std::numeric_limits<double>::quiet_NaN();
                _data._event_trigger = UNDEFINED_ELEMENT_ID;
                _data._shear_stress_init = _data._shear_stress_final = std::numeric_limits<double>::quiet_NaN();
                _data._normal_stress_init = _data._normal_stress_final = std::numeric_limits<double>::quiet_NaN();
                _data._start_sweep_rec = _data._end_sweep_rec = UNDEFINED_EVENT_ID;

                _sweeps.clear();
                _total_slip.clear();
            }

#ifdef HDF5_FOUND
            static std::string hdf5_table_name(void) {
                return "events";
            };
            static void setup_event_hdf5(const hid_t &data_file);
            void append_event_hdf5(const hid_t &data_file) const;
#endif
            static void get_field_descs(std::vector<FieldDesc> &descs);
            static void write_ascii_header(std::ostream &out_stream);

            void read_data(const EventData &in_data);
            void write_data(EventData &out_data) const;

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;

            friend std::ostream &operator<<(std::ostream &os, const ModelEvent &me);
    };

    class ModelEventSet {
        private:
            std::vector<ModelEvent>     _events;
            void read_events_hdf5(const hid_t &data_file);
            void read_sweeps_hdf5(const hid_t &data_file);

        public:
            typedef std::vector<ModelEvent>::iterator       iterator;
            typedef std::vector<ModelEvent>::const_iterator const_iterator;

            iterator begin(void) {
                return _events.begin();
            };
            iterator end(void) {
                return _events.end();
            };

            unsigned int size(void) const {
                return _events.size();
            }

            ModelEvent &operator[](const unsigned int ind) throw(std::out_of_range) {
                if (ind >= _events.size()) throw std::out_of_range("ModelEventSet[]");

                return _events[ind];
            };

            const ModelEvent &operator[](const unsigned int ind) const throw(std::out_of_range) {
                if (ind >= _events.size()) throw std::out_of_range("ModelEventSet[]");

                return _events[ind];
            };

            int read_file_ascii(const std::string &event_file_name, const std::string &sweep_file_name);

            int read_file_hdf5(const std::string &file_name);
    };

    /*!
     The stress state of an element in the model at a specified time in the simulation.
     32-bit floats are used to save space since there may be vast amounts of stress
     data output and high precision is not needed for analysis.
     */
    struct StressData {
        //! The ID of the element
        UIndex          _element_id;

        //! Shear and normal stress on the element at this time
        float          _shear_stress, _normal_stress;
    };

    /*!
     ModelStress is the stress state of the model at a specified point in time.
     */
    class ModelStress : public ModelIO {
        private:
            std::vector<StressData>     _data;

        public:
            ModelStress(void) {
                clear();
            }

            void clear(void) {
                _data.clear();
            }

            void add_stress_entry(UIndex element_id, float shear_stress, float normal_stress) {
                StressData new_entry;
                new_entry._element_id = element_id;
                new_entry._shear_stress = shear_stress;
                new_entry._normal_stress = normal_stress;
                _data.push_back(new_entry);
            }
#ifdef HDF5_FOUND
            static std::string hdf5_table_name(void) {
                return "stresses";
            };
            static void setup_stress_hdf5(const hid_t &data_file);
            void append_stress_hdf5(const hid_t &data_file) const;
#endif
            static void get_field_descs(std::vector<FieldDesc> &descs);
            static void write_ascii_header(std::ostream &out_stream);

            void read_ascii(std::istream &in_stream, const unsigned int num_records);
            void write_ascii(std::ostream &out_stream) const;

            unsigned int size(void) const {
                return _data.size();
            }

            StressData &operator[](const unsigned int ind) throw(std::out_of_range) {
                if (ind >= _data.size()) throw std::out_of_range("ModelStress[]");

                return _data[ind];
            };

            const StressData &operator[](const unsigned int ind) const throw(std::out_of_range) {
                if (ind >= _data.size()) throw std::out_of_range("ModelStress[]");

                return _data[ind];
            };

            friend std::ostream &operator<<(std::ostream &os, const ModelStress &ms);
    };

    struct StressDataTime {
        //! The year of the stress data
        float           _year;

        //! The event and sweep of the stress data (if applicable)
        UIndex          _event_num, _sweep_num;

        //! Starting and ending index of the actual stress entries
        UIndex          _start_rec, _end_rec;
    };

    /*!
     ModelStressState is the stress state of the model at a specified point in time, including information specifying the time.
     */
    class ModelStressState : public ModelIO {
        private:
            ModelStress         _stress;

            StressDataTime      _times;

        public:
            ModelStressState(void) {
                clear();
            }

            void clear(void) {
                _stress.clear();
                _times._year = std::numeric_limits<float>::quiet_NaN();
                _times._event_num = UNDEFINED_EVENT_ID;
                _times._sweep_num = UNDEFINED_EVENT_ID;
                _times._start_rec = UNDEFINED_EVENT_ID;
                _times._end_rec = UNDEFINED_EVENT_ID;
            }

#ifdef HDF5_FOUND
            static std::string hdf5_table_name(void) {
                return "stress_state";
            };
            static void setup_stress_state_hdf5(const hid_t &data_file);
            void append_stress_state_hdf5(const hid_t &data_file) const;
#endif
            static void get_field_descs(std::vector<FieldDesc> &descs);
            static void write_ascii_header(std::ostream &out_stream);

            void setStresses(const ModelStress &new_stresses) {
                _stress = new_stresses;
            };
            const ModelStress &stresses(void) const {
                return _stress;
            }
            ModelStress &stresses(void) {
                return _stress;
            }

            void setYear(const double year) {
                _times._year = year;
            }
            float getYear(void) const {
                return _times._year;
            }
            void setEventNum(const UIndex event_num) {
                _times._event_num = event_num;
            }
            UIndex getEventNum(void) const {
                return _times._event_num;
            }
            void setSweepNum(const UIndex sweep_num) {
                _times._sweep_num = sweep_num;
            }
            UIndex getSweepNum(void) const {
                return _times._sweep_num;
            }
            void setStartEndRecNums(const unsigned int start_rec, const unsigned int end_rec) {
                _times._start_rec = start_rec;
                _times._end_rec = end_rec;
            }
            unsigned int getNumStressRecords(void) const {
                return _times._end_rec - _times._start_rec;
            };

            void read_ascii(std::istream &in_stream);
            void write_ascii(std::ostream &out_stream) const;

            friend std::ostream &operator<<(std::ostream &os, const ModelStressState &mss);
    };

    class ModelStressSet {
        private:
            std::vector<ModelStressState>   _states;

        public:
            typedef std::vector<ModelStressState>::iterator         iterator;
            typedef std::vector<ModelStressState>::const_iterator   const_iterator;

            iterator begin(void) {
                return _states.begin();
            };
            iterator end(void) {
                return _states.end();
            };

            unsigned int size(void) const {
                return _states.size();
            }

            ModelStressState &operator[](const unsigned int ind) throw(std::out_of_range) {
                if (ind >= _states.size()) throw std::out_of_range("ModelStressSet[]");

                return _states[ind];
            };

            const ModelStressState &operator[](const unsigned int ind) const throw(std::out_of_range) {
                if (ind >= _states.size()) throw std::out_of_range("ModelStressSet[]");

                return _states[ind];
            };

            int read_file_ascii(const std::string &stress_index_file_name, const std::string &stress_file_name);
    };
}

#endif
