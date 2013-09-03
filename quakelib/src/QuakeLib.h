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

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <limits.h>
#include <float.h>

#include "QuakeLibOkada.h"

#ifndef _QUAKELIB_H_
#define _QUAKELIB_H_

// TODO: add default constructors to blank out objects

namespace quakelib {
	//! Holds a list of errors related to an EqSim model. Each error is a string and
	//! an optional associated line. The existence of errors does not affect editing
	//! or writing of a model. Erroneous models may still be written or used for a simulation.
	class EQSimErrors {
	private:
		//! A list of errors associated with an EqSim model or file, as well as
		//! the line number (< 0 if there is no associated line number).
		std::multimap<int, std::string>	error_list;
		
	public:
		EQSimErrors(void) { error_list.clear(); };
		//! Returns a count of the number of errors.
		size_t count(void) const { return error_list.size(); };
		//! Record an error and associated line number (< 0 if not from a file).
		void report(const std::string &error_msg, const int &line_num);
		//! Format errors and write them to the specified ostream.
		void write(std::ostream &os) const;
		//! Clear all errors from the list.
		void clear(void) { error_list.clear(); };
	};
	
	//! Unsigned index for EqSim objects.
	typedef unsigned int UIndex;
	//! Map of record indices with their associated line numbers.
	typedef std::map<UIndex, int> LineNumbers;
	
	//! Unsigned index for section ids
	typedef unsigned int SectionID;
	
	// Enumerations for metadata record types.
	enum RECORD_TYPE {
		META_COMMENT,
		META_SIGNATURE,
		META_INFO,
		META_TITLE,
		META_AUTHOR,
		META_DATE
	};
	
	// These classes are internal to the EqSim parsing/writing and generally can be ignored by the user.
	namespace internal {
		//! Simple wrapper class that represents an output file for the EqSim format.
		class EQSimFileWriter {
		protected:
			//! The output file stream for this file writer.
			std::ofstream		out_stream;
			
		public:
			//! Open a file with the given name to write to. If the file already exists, it will be overwritten.
			void open(const std::string &file_name) {
				out_stream.open(file_name.c_str(), std::ios_base::out|std::ios_base::trunc);
				//TODO: if (!out_stream) throw std::exception("Failed to open stream.");
			};
			//! Close the file and flush unwritten data.
			void close(void) { out_stream.close(); };
		};
		
		typedef std::string Record;
		typedef std::vector<Record> MultiRecord;
		
		//! Represents a field descriptor in the metadata section. This is comprised of a field name, data type and field description.
		class FieldDesc {
		private:
			std::string	_name;		//! A character string that gives the name of the field.
			int			_type;		//! An integer giving the type of data in the field.
			int			_field_num;	//! Field number in this record.
			std::string _desc;		//! Description of the field.
			
		public:
			FieldDesc(const std::string &name="", const int &type=0, const int &field_num=0, const std::string &desc="") : _name(name), _type(type), _field_num(field_num), _desc(desc) {};
			
			//! Write this field with the specified key to the specified output stream.
			void write(const unsigned int &key, std::ostream &out_stream) const;
		};
		
		// key is an integer in the range 1 to 100 that gives the position of the field within the data record
		typedef std::map<int, FieldDesc> MultiFieldDesc;
		
		//! Represents a record descriptor in the metadata section. This is comprised of a record name and description, and a set of field descriptors related to the record.
		class RecordDesc {
		private:
			int				line_num;	//! Line that this record was declared on
			std::string		name;		//! A character string that gives the name of the data record
			int				n_field;	//! An integer giving the number of fields in the data record
			std::string		desc;		//! Description of the data record
			MultiFieldDesc	field_descs;//! Actual field descriptors tied to this record descriptor
			
		public:
			RecordDesc(const int &_line_num=0, const std::string &_name="", const int &_n_field=0, const std::string &_desc="") :
			line_num(_line_num), name(_name), n_field(_n_field), desc(_desc) {};
			
			//! Add a field to the record description.
			bool add_field(const int &_index, const FieldDesc &_new_desc);
			
			//! Returns the line number this record was read from, -1 if not read.
			int get_line_num(void) const { return line_num; };
			//! Returns the maximum number of fields described in this record.
			int get_max_index(void) const { return n_field; };
			//! Returns the actual number of fields described in this record.
			size_t get_num_fields(void) const { return field_descs.size(); };
			
			//! Write this record with the specified key to the specified output stream.
			void write(const unsigned int &key, std::ostream &out_stream) const;
		};
		
		// key is an integer in the range 200 to 399 that gives the kind of data record being described
		typedef std::map<int, RecordDesc> MultiRecordDesc;
		
		// TODO: move metadata out of internal(?)
		//! Represents the metadata and descriptor sections of an EqSim file. This includes
		//! information such as the file title, author, date of creation, and record/field
		//! descriptors for the file data.
		class EQSimMetadata {
		private:
			// Checks whether a string has invalid characters (CR, LF, etc)
			bool has_bad_chars(const std::string &test_str);
		protected:
			//! The specification level of this data.
			//! TODO: properly support this.
			int				spec_level;
			
			//! Metadata records organized by record type.
			//! Although these are all MultiRecords, some types may only allow a single record.
			std::map<RECORD_TYPE, MultiRecord>			metadata_recs;
			
			//! Metadata record line numbers.
			std::map<RECORD_TYPE, std::vector<int> >	metadata_line_nums;
			
			//! Map of metadata record type to record number.
			std::map<RECORD_TYPE, unsigned int>			metadata_rec_nums;
			
			//! Descriptor records.
			MultiRecordDesc	record_descs;
		public:
			// TODO: fix this to be the correct signature record and spec level for a given file
			EQSimMetadata(void) {
				spec_level = 0;
				metadata_recs[META_COMMENT] = MultiRecord();
				metadata_recs[META_SIGNATURE] = MultiRecord();
				metadata_recs[META_INFO] = MultiRecord();
				metadata_recs[META_TITLE] = MultiRecord();
				metadata_recs[META_AUTHOR] = MultiRecord();
				metadata_recs[META_DATE] = MultiRecord();
				
				metadata_line_nums[META_COMMENT] = std::vector<int>();
				metadata_line_nums[META_SIGNATURE] = std::vector<int>();
				metadata_line_nums[META_INFO] = std::vector<int>();
				metadata_line_nums[META_TITLE] = std::vector<int>();
				metadata_line_nums[META_AUTHOR] = std::vector<int>();
				metadata_line_nums[META_DATE] = std::vector<int>();
				
				metadata_rec_nums[META_COMMENT] = 100;
				metadata_rec_nums[META_SIGNATURE] = 101;
				metadata_rec_nums[META_INFO] = 110;
				metadata_rec_nums[META_TITLE] = 111;
				metadata_rec_nums[META_AUTHOR] = 112;
				metadata_rec_nums[META_DATE] = 113;
			};
			
			//! Set the specification level.
			void set_spec_level(const int &new_level) { spec_level = new_level; };
			//! Get the specification level.
			int get_spec_level(void) const { return spec_level; };
			
			//! Get the number of metadata records of a given type. Throws an exception if the type is unknown.
			unsigned int meta_num_records(const RECORD_TYPE &rec_type) const
				throw(std::invalid_argument);
			
			//! Clear all metadata records of a given type. Throws an exception if the type is unknown.
			void meta_clear_record(const RECORD_TYPE &rec_type)
				throw(std::invalid_argument);
			
			//! Get a metadata record of a given type. Throws an exception if the type is unknown
			//! or the rec_num is outside of the valid range.
			std::string meta_get_record(const RECORD_TYPE &rec_type, const unsigned int &rec_num) const
				throw(std::invalid_argument, std::out_of_range);
			
			//! Append a record of a given type. If the type can only have one associated record,
			//! overwrites any previous recorded records. Throws an exception if the type is unknown.
			void meta_add_record(const RECORD_TYPE &rec_type, const std::string &new_rec)
				throw(std::invalid_argument);
			
			//! Set a metadata record of a given type. Throws an exception if the type is unknown
			//! or the rec_num is outside of the valid range.
			void meta_set_record(const RECORD_TYPE &rec_type, const unsigned int &rec_num, const std::string &new_rec)
				throw(std::invalid_argument, std::out_of_range);
			
			//! Erase the specified metadata record of the specified type. Throws an exception
			//! if the type is unknown or the rec_num is outside the valid range.
			void meta_erase_record(const RECORD_TYPE &rec_type, const unsigned int &rec_num) throw(std::invalid_argument, std::out_of_range);
			
			//! Add a record description to the metadata. An exception will be thrown if the key is not between 200 to 399.
			void add_record_desc_record(const unsigned int &rec_key, internal::RecordDesc &rec_desc) throw(std::invalid_argument);
			
			//! Perform internal correctness checks and write any errors to the specified Errors object.
			void validate(EQSimErrors &errors) const;
		};
	}
	
	//! A class of functions for reading the metadata of an EqSim style file.
	//! The EqSimFileReader works by using callback functions to process specified
	//! record numbers. This means a certain function can be specified to process
	//! record 110, and another function can process record 121, and so on. These
	//! functions can be turned on or off to enable context dependent parsing.
	//! Subclasses provide additional callbacks to handle different files.
	class EQSimMetadataReader : public internal::EQSimMetadata {
	private:
		//! Current record ID being parsed.
		int				cur_record_index;
		
		bool			parsed_sig, parsed_title, parsed_date, parse_metadata_finished, parsed_descriptors;
		
		bool parse_metadata_line(const int &rec_num, const int &line_num, std::istringstream &line_stream);
		
		virtual bool parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream) = 0;
		
		void parse_comment(const int &line_num, std::istringstream &line_stream);
		void parse_signature(const int &line_num, std::istringstream &line_stream);
		void parse_end_metadata(const int &line_num, std::istringstream &line_stream);
		void parse_end_descriptor(const int &line_num, std::istringstream &line_stream);
		void parse_info_record(const int &line_num, std::istringstream &line_stream);
		void parse_title_record(const int &line_num, std::istringstream &line_stream);
		void parse_author_record(const int &line_num, std::istringstream &line_stream);
		void parse_date_record(const int &line_num, std::istringstream &line_stream);
		void parse_record_desc_record(const int &line_num, std::istringstream &line_stream);
		void parse_field_desc_record(const int &line_num, std::istringstream &line_stream);
		void parse_eof(const int &line_num, std::istringstream &line_stream);
		
	protected:
		//! Records errors associated with parsing the metadata.
		EQSimErrors			parse_errors;
		
		//! Whether parsing is completed for the input file metadata.
		bool				finish_parsing;
		
	public:
		EQSimMetadataReader(void) : parsed_sig(false), parsed_title(false), parsed_date(false), parse_metadata_finished(false), parsed_descriptors(false) { cur_record_index = -1; };
		
		//! Parses an input stream and calls the appropriate callback functions.
		//! The callback function to use for a given line is defined by the 3 digit
		//! code at the line start. Returns true if the file exists.
		// TODO: throw exception if problem reading file (e.g. doesn't exist)
		bool parse_file(const std::string &file_name, const int &start_line_num=0);
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	//! A class for storing and writing the metadata in an EqSim style file.
	// TODO: make writer abstract like reader is
	class EQSimMetadataWriter : public internal::EQSimMetadata, public internal::EQSimFileWriter {
	private:
		bool	wrote_header;		//! Whether we have written the header lines yet
		
	public:
		EQSimMetadataWriter(void) : wrote_header(false) {};
		
		//! Write all metadata to the associated file output stream.
		void write(void);
		
		//! Close the file output stream and flush the data.
		void close(void);
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	typedef std::map<UIndex, std::pair<double, double> > StressMap;
	typedef std::map<UIndex, double> StateMap;
	
	//! Represents the records of an EqSim style file describing the initial stresses and rate state laws for a given model.
	class EQSimCondition {
	protected:
		// Record 201: Initial stress states
		//! A map of index numbers to initial stress data.
		StressMap	_stresses;		// Field 2,3: Element initial shear/normal stress (Pascal)
		
		// Record 202: Initial state
		//! A map of index numbers to rate-state data.
		StateMap	_rs_theta;		// Field 2: Element initial state for rate-state law (seconds)
		
		//! Lines of parsed initial conditions (-1 if not parsed).
		LineNumbers	_stress_lines;
		//! Lines of parsed rate-state data (-1 if not parsed).
		LineNumbers	_rst_lines;
		
	public:
		//! Returns the number of elements found in this condition data.
		size_t num_elements(void) const;
		//! Get the initial shear and normal stress (in Pascals) of the specified element.
		//! If the element is not found an exception is raised.
		double get_shear_stress(const UIndex &ind) const throw(std::out_of_range);
		double get_normal_stress(const UIndex &ind) const throw(std::out_of_range);
		//! Set the initial shear and normal stress (in Pascals) of the specified element.
		void set_stresses(const UIndex &ind, const double &shear_stress, const double &normal_stress);
		
		//! Returns the rate-state data for the specified index.
		double get_rate_state(const UIndex &ind) const throw(std::out_of_range);
		//! Sets the rate-state data for the specified index.
		void set_rate_state(const UIndex &ind, const double &rate_state);
		
		//! Write the condition data to the specified output stream.
		void write(std::ostream &out_stream) const;
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	//! Represents the records of an EqSim style file describing the initial stresses and rate state laws for a given model.
	class EQSimConditionReader : public EQSimMetadataReader, public EQSimCondition {
	private:
		// Record 200: Initial condition summary
		int				_line_num;
		unsigned int	_n_element;		// Field 1: Total number of elements in the file
		bool			_stress_flag;	// Field 2: 1 if initial stress (record 201) is included, 0 if not
		bool			_state_flag;		// Field 3: 1 if initial state (record 202) is included, 0 if not
		
		void parse_summary_record(const int &line_num, std::istringstream &line_stream);
		void parse_stress_record(const int &line_num, std::istringstream &line_stream);
		void parse_state_record(const int &line_num, std::istringstream &line_stream);
		
		bool parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream);
		
	public:
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	//! Class to output EqSim condition records to a file.
	class EQSimConditionWriter : public EQSimMetadataWriter, public EQSimCondition {
	public:
		EQSimConditionWriter(void);
		void write(void) { EQSimMetadataWriter::write(); EQSimCondition::write(out_stream); };
	};
	
	//! Represents a slip element in an EqSim style event output file.
	class EQSimEventSlipElement {		// Record 202: Slip element list entry
	private:
		int			_line_num;		// line number this record was read from (-1 if not read from a file)
		UIndex		_element_id;	// Field 1: Element ID number (positive integer), or negative of element count (zero if no element info)
		
	public:
		//! Create an EQSimEventSlipElement with the specified element ID and optional line number.
		EQSimEventSlipElement(const UIndex &element_id, const int &line=-1) : _line_num(line), _element_id(element_id) {};
		
		//! Returns element ID associated with this event slip.
		UIndex element_id(void) const { return _element_id; };
		//! Sets element ID associated with this event slip.
		void set_element_id(const UIndex &element_id) { _element_id = element_id; };
		
		//! Write the record for this event slip element.
		void write(std::ostream &out_stream) const;
	};
	
	typedef std::vector<EQSimEventSlipElement> EQSimEventSlipElementList;
	
	//! Represents a slip map for an event in an EqSim style event output file.
	class EQSimEventSlipMap {		// Record 201: Slip map entry
	private:
		int			_line_num;		//! line number this record was read from (-1 if not read from a file)
		double		_depth_lo;		//! Field 1: Lowest value of depth (meters, negative underground)
		double		_depth_hi;		//! Field 2: Highest value of depth in the rupture (meters, negative underground)
		double		_das_lo;		//! Field 3: Lowest value of distance-along-strike in the rupture (meters)
		double		_das_hi;		//! Field 4: Highest value of distance-along-strike in the rupture (meters)
		double		_area;			//! Field 5: Rupture area (square meters)
		double		_mean_slip;		//! Field 6: Average slip (meters)
		double		_moment;		//! Field 7: Seismic moment (Newton-meters)
		double		_shear_before;	//! Field 8: Shear stress before event (Pascal)
		double		_shear_after;	//! Field 9: Shear stress after event (Pascal)
		double		_normal_before;	//! Field 10: Normal stress before event (Pascal)
		double		_normal_after;	//! Field 11: Normal stress after event (Pascal)
		int			_element_id;	//! Field 12: Element ID number (positive integer), or negative of element count (zero if no element info)
		
		EQSimEventSlipElementList	slip_elements;
		
	public:
		//! Returns the lowest depth of the rupture (meters, negative underground).
		double depth_lo(void) const { return _depth_lo; };
		//! Returns the highest depth of the rupture (meters, negative underground).
		double depth_hi(void) const { return _depth_hi; };
		//! Returns the lowest distance-along-strike of the rupture (meters).
		double das_lo(void) const { return _das_lo; };
		//! Returns the highest distance-along-strike of the rupture (meters).
		double das_hi(void) const { return _das_hi; };
		//! Returns the total rupture area (in square meters).
		double area(void) const { return _area; };
		//! Mean slip (in meters).
		double mean_slip(void) const { return _mean_slip; };
		//! Rupture seismic moment (in Newton-meters).
		double moment(void) const { return _moment; };
		//! Total shear stress on ruptured elements before the event (Pascals).
		double shear_before(void) const { return _shear_before; };
		//! Total shear stress on ruptured elements after the event (Pascals).
		double shear_after(void) const { return _shear_after; };
		//! Total normal stress on ruptured elements before the event (Pascals).
		double normal_before(void) const { return _normal_before; };
		//! Total normal stress on ruptured elements after the event (Pascals).
		double normal_after(void) const { return _normal_after; };
		
		//! Set the lowest depth of the rupture (meters, negative underground).
		void set_depth_lo(const double &depth_lo) { _depth_lo = depth_lo; };
		//! Set the highest depth of the rupture (meters, negative underground).
		void set_depth_hi(const double &depth_hi) { _depth_hi = depth_hi; };
		//! Set the lowest distance-along-strike of the rupture (meters).
		void set_das_lo(const double &das_lo) { _das_lo = das_lo; };
		//! Set the highest distance-along-strike of the rupture (meters).
		void set_das_hi(const double &das_hi) { _das_hi = das_hi; };
		//! Set the total rupture area (in square meters).
		void set_area(const double &area) { _area = area; };
		//! Set the mean slip (in meters).
		void set_mean_slip(const double &mean_slip) { _mean_slip = mean_slip; };
		//! Set the rupture seismic moment (in Newton-meters).
		void set_moment(const double &moment) { _moment = moment; };
		//! Set total shear stress on ruptured elements before the event (Pascals).
		void set_shear_before(const double &shear_before) { _shear_before = shear_before; };
		//! Set total shear stress on ruptured elements after the event (Pascals).
		void set_shear_after(const double &shear_after) { _shear_after = shear_after; };
		//! Set total normal stress on ruptured elements before the event (Pascals).
		void set_normal_before(const double &normal_before) { _normal_before = normal_before; };
		//! Set total normal stress on ruptured elements after the event (Pascals).
		void set_normal_after(const double &normal_after) { _normal_after = normal_after; };
		
		//! Add an elmeent slip entry to this slip map.
		void add_slip_entry(const EQSimEventSlipElement &new_entry) { slip_elements.push_back(new_entry); };
		
		//! Parse an event slip map record from the input stream to this object.
		void parse(const int &line_num, std::istringstream &line_stream);
		//! Write the record for this event slip map.
		void write(std::ostream &out_stream) const;
	};
	
	typedef std::vector<EQSimEventSlipMap> EQSimEventSlipList;
	
	//! Represents an event summary in an EqSim style event output file.
	class EQSimEventSummary {			// Record 200: Event summary
	private:
		int				_line_num;		// line number this record was read from (-1 if not read from a file)
		UIndex			_event_id;		// Field 1: Event ID number (positive integers, in order, need not be consecutive)
		double			_magnitude;		// Field 2: Event magnitude
		double			_time;			// Field 3: Starting time of event (seconds)
		double			_duration;		// Field 4: Duration of event (seconds)
		UIndex			_sid;			// Field 5: Fault section ID number (positive integer)
		double			_depth_lo;		// Field 6: Lowest value of depth in the rupture (meters, negative underground)
		double			_depth_hi;		// Field 7: Highest value of depth in the rupture (meters, negative underground)
		double			_das_lo;		// Field 8: Lowest value of distance-along-strike in the rupture (meters)
		double			_das_hi;		// Field 9: Highest value of distance-along-strike in the rupture (meters)
		double			_hypo_depth;	// Field 10: Hypocenter depth (meters, negative underground)
		double			_hypo_das;		// Field 11: Hypocenter distance-along-strike (meters)
		double			_area;			// Field 12: Rupture area (square meters)
		double			_mean_slip;		// Field 13: Average slip (meters)
		double			_moment;		// Field 14: Seismic moment (Newton-meters)
		double			_shear_before;	// Field 15: Shear stress before event (Pascal)
		double			_shear_after;	// Field 16: Shear stress after event (Pascal)
		double			_normal_before;	// Field 17: Normal stress before event (Pascal)
		double			_normal_after;	// Field 18: Normal stress after event (Pascal)
		
	public:
		// slip maps associated with this event summary
		//! TODO: make this private?
		EQSimEventSlipList	slip_maps;
		
		//! Returns event ID.
		UIndex event_id(void) const { return _event_id; };
		//! Returns the magnitude of the event.
		double magnitude(void) const { return _magnitude; };
		//! Returns the event starting time, measured in seconds from simulation start.
		double time(void) const { return _time; };
		//! Returns the event duration (in seconds).
		double duration(void) const { return _duration; };
		//! Returns the section ID of the event. Events spanning multiple sections
		//! should be represented as separate events.
		UIndex sid(void) const { return _sid; };
		//! Returns the lowest depth of the ruptured portion of the fault (meters, negative underground).
		double depth_lo(void) const { return _depth_lo; };
		//! Returns the highest depth of the ruptured portion of the fault (meters, negative underground).
		double depth_hi(void) const { return _depth_hi; };
		//! Returns the lowest distance-along-strike of the ruptured portion of the fault (meters).
		double das_lo(void) const { return _das_lo; };
		//! Returns the highest distance-along-strike of the ruptured portion of the fault (meters).
		double das_hi(void) const { return _das_hi; };
		//! Returns the depth of the event hypocenter (meters, negative underground).
		double hypo_depth(void) const { return _hypo_depth; };
		//! Returns the distance-along-strike of the event hypocenter (meters).
		double hypo_das(void) const { return _hypo_das; };
		//! Returns the rupture area for the event (in square meters).
		double area(void) const { return _area; };
		//! Mean slip (in meters).
		double mean_slip(void) const { return _mean_slip; };
		//! Rupture seismic moment (in Newton-meters).
		double moment(void) const { return _moment; };
		//! Total shear stress on ruptured elements before the event (Pascals).
		double shear_before(void) const { return _shear_before; };
		//! Total shear stress on ruptured elements after the event (Pascals).
		double shear_after(void) const { return _shear_after; };
		//! Total normal stress on ruptured elements before the event (Pascals).
		double normal_before(void) const { return _normal_before; };
		//! Total normal stress on ruptured elements after the event (Pascals).
		double normal_after(void) const { return _normal_after; };
		
		//! Set ID of this event.
		void set_event_id(const UIndex &event_id) { _event_id = event_id; };
		//! Set the magnitude of the event.
		void set_magnitude(const double &magnitude) { _magnitude = magnitude; };
		//! Set the event starting time, measured in seconds from simulation start.
		void set_time(const double &time) { _time = time; };
		//! Set the event duration (in seconds).
		void set_duration(const double &duration) { _duration = duration; };
		//! Set the section ID of the event. Events spanning multiple sections
		//! should be represented as separate events.
		void set_sid(const UIndex &sid) { _sid = sid; };
		//! Set the lowest depth of the rupture (meters, negative underground).
		void set_depth_lo(const double &depth_lo) { _depth_lo = depth_lo; };
		//! Set the highest depth of the rupture (meters, negative underground).
		void set_depth_hi(const double &depth_hi) { _depth_hi = depth_hi; };
		//! Set the lowest distance-along-strike of the rupture (meters).
		void set_das_lo(const double &das_lo) { _das_lo = das_lo; };
		//! Set the highest distance-along-strike of the rupture (meters).
		void set_das_hi(const double &das_hi) { _das_hi = das_hi; };
		//! Set the depth of the event hypocenter (meters, negative underground).
		void set_hypo_depth(const double &hypo_depth) { _hypo_depth = hypo_depth; };
		//! Set the distance-along-strike of the event hypocenter (meters).
		void set_hypo_das(const double &hypo_das) { _hypo_das = hypo_das; };
		//! Set the total rupture area (in square meters).
		void set_area(const double &area) { _area = area; };
		//! Set the mean slip (in meters).
		void set_mean_slip(const double &mean_slip) { _mean_slip = mean_slip; };
		//! Set the rupture seismic moment (in Newton-meters).
		void set_moment(const double &moment) { _moment = moment; };
		//! Set total shear stress on ruptured elements before the event (Pascals).
		void set_shear_before(const double &shear_before) { _shear_before = shear_before; };
		//! Set total shear stress on ruptured elements after the event (Pascals).
		void set_shear_after(const double &shear_after) { _shear_after = shear_after; };
		//! Set total normal stress on ruptured elements before the event (Pascals).
		void set_normal_before(const double &normal_before) { _normal_before = normal_before; };
		//! Set total normal stress on ruptured elements after the event (Pascals).
		void set_normal_after(const double &normal_after) { _normal_after = normal_after; };
		
		//! Add a slip map to this event summary.
		void add_slip_map(const EQSimEventSlipMap &new_map) { slip_maps.push_back(new_map); };
		
		//! Parse an event summary record from the input stream to this object.
		void parse(const int &line_num, std::istringstream &line_stream);
		//! Write the record for this event summary to the specified output stream.
		void write(std::ostream &out_stream) const;
	};
	
	typedef std::vector<EQSimEventSummary> EQSimEventSummaryList;
	
	//! Represents a set of event summaries.
	class EQSimEventSet {
	public:
		//! The event summaries - note there can be multiple summaries per
		//! event number so the vector index is not necessarily the event number
		// TODO: make this private?
		EQSimEventSummaryList				event_summaries;
		
		//! Add an event summary to the event set.
		void add_event_summary(const EQSimEventSummary &new_summary) { event_summaries.push_back(new_summary); };
	};
	
	// TODO: fix write to actually write events without needing to call flush()
	//! Represents a collection of event summaries and related info to write.
	class EQSimEventWriter : public EQSimEventSet, public EQSimMetadataWriter {
	private:
		std::vector<UIndex>	entries_to_flush;	//! Ordered list of event summary indices to flush to file
		
	public:
		EQSimEventWriter(void);
		//! Add an event summary to be written for this event set.
		void add_event_summary(const EQSimEventSummary &new_summary) { entries_to_flush.push_back(event_summaries.size()); EQSimEventSet::add_event_summary(new_summary); };
		//! Flush any unwritten event summaries to the file.
		void flush(void);
	};
	
	//! Represents a collection of parse event summaries and related info.
	class EQSimEventReader : public EQSimEventSet, public EQSimMetadataReader {
	private:
		void parse_event_summary(const int &line_num, std::istringstream &line_stream);
		void parse_event_slip_map(const int &line_num, std::istringstream &line_stream);
		void parse_event_slip_element(const int &line_num, std::istringstream &line_stream);
		
		virtual bool parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream);
	};
	
	//! Represents the parameters for a rate-state equation specified by the EqSim friction file.
	// TODO: add correctness check for this
	class EQSimFrictionRateState {
	private:
		//! Field 1: Element index number (consecutive integers, starting with 1)
		double			_A;		//! Field 2: Element rate-state parameter A
		double			_B;		//! Field 3: Element rate-state parameter B
		double			_L;		//! Field 4: Element rate-state characteristic distance L (meters)
		double			_f0;	//! Field 5: Element rate-state friction coefficient f0
		double			_V0;	//! Field 6: Element rate-state reference velocity V0 (meters/second)
		
	public:
		EQSimFrictionRateState(void) : _A(0), _B(0), _L(0), _f0(0), _V0(0) {};
		EQSimFrictionRateState(const double &A, const double &B, const double &L, const double &f0, const double &V0) : _A(A), _B(B), _L(L), _f0(f0), _V0(V0) {};
		
		double A(void) const { return _A; };
		double B(void) const { return _B; };
		double L(void) const { return _L; };
		double f0(void) const { return _f0; };
		double V0(void) const { return _V0; };
		
		//! Parse rate state parameter record from the input stream to this object.
		void parse(const int &line_num, std::istringstream &line_stream);
		//! Write a record for the friction rate state parameters to the specified output stream.
		void write(std::ostream &out_stream) const;
	};
	
	typedef std::map<UIndex, std::pair<double, double> > StrengthMap;
	typedef std::map<UIndex, EQSimFrictionRateState> RateStateMap;
	
	//! Represents an EqSim style friction information file.
	class EQSimFriction {
	protected:
		// Record 201: Elastic parameters
		//! Line number of elastic parameters
		int				_line_num;
		//! Whether this friction data has specified elastic parameters.
		bool			_have_elastic_params;
		//! Field 1: Lame modulus lambda (Pascal).
		double			_lame_lambda;
		//! Field 2: Lame modulus mu, also known as the shear modulus (Pascal).
		double			_lame_mu;
		
		// Record 202: Fault strength
		//! Field 2/3: Element static yield/dynamic sliding strength (Pascal)
		StrengthMap		_strengths;
		
		//! Record 203: Rate-state parameters
		//! TODO: ensure these are implemented correctly
		RateStateMap	_rate_state;
		
		//! Lines of parsed strengths (-1 if not parsed).
		LineNumbers		_str_lines;
		//! Lines of parsed rate-state parameters (-1 if not parsed).
		LineNumbers		_rs_lines;
		
	public:
		EQSimFriction(void) : _line_num(-1), _have_elastic_params(false), _lame_lambda(0), _lame_mu(0) {
			_strengths.clear(); _rate_state.clear(); _str_lines.clear(); _rs_lines.clear();
		};
		
		//! Set the strengths for the specified element index.
		void set_strengths(const UIndex &ind, const double &static_strength, const double &dynamic_strength);
		//! Get the static strength for the specified element index.
		double get_static_strength(const UIndex &ind) const throw(std::out_of_range);
		//! Get the dynamic strength for the specified element index.
		double get_dynamic_strength(const UIndex &ind) const throw(std::out_of_range);
		
		//! Get the rate-state parameters for the specified element index.
		EQSimFrictionRateState get_rs_param(const UIndex &ind) const;
		//! Set the rate-state parameters for the specified element index.
		void set_rs_param(const UIndex &ind, const EQSimFrictionRateState &params);
		
		//! Get the constant Lame lambda parameter for this friction data (in Pascals).
		double get_lame_lambda(void) const { return _lame_lambda; };
		//! Get the constant Lame mu parameter for this friction data (in Pascals).
		double get_lame_mu(void) const { return _lame_mu; };
		//! Set the constant Lame lambda and mu parameters for this friction data (in Pascals).
		void set_lame_lambda_mu(const double &new_lambda, const double &new_mu) { _lame_lambda = new_lambda; _lame_mu = new_mu; _have_elastic_params = true; };
		
		//! Write the friction data records to the specified output stream.
		void write(std::ostream &out_stream) const;
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	//! Represents a parsed EqSim style friction information file.
	class EQSimFrictionReader : public EQSimFriction, public EQSimMetadataReader {
	private:
		// Record 200: Fault friction summary
		int				_line_num;
		unsigned int	_num_element;	//! Field 1: Total number of elements in the file
		bool			_elastic_flag;	//! Field 2: 1 if elastic parameters (record 201) are included, 0 if not
		bool			_strength_flag;	//! Field 3: 1 if fault strength (record 202) is included, 0 if not
		bool			_rate_state_flag;//! Field 4: 1 if rate-state parameters (record 203) are included, 0 if not
		
		void parse_summary_record(const int &line_num, std::istringstream &line_stream);
		void parse_elastic_params_record(const int &line_num, std::istringstream &line_stream);
		void parse_fault_strength_record(const int &line_num, std::istringstream &line_stream);
		void parse_rate_state_record(const int &line_num, std::istringstream &line_stream);
		
		virtual bool parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream);
		
	public:
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	//! Output file for EqSim friction data records.
	class EQSimFrictionWriter : public EQSimMetadataWriter, public EQSimFriction {
	public:
		EQSimFrictionWriter(void);
		
		//! Write the metadata and friction data for this model.
		void write(void) { EQSimMetadataWriter::write(); EQSimFriction::write(out_stream); };
	};
	
	//! Possible locations of a vertex on the fault trace.
	enum TraceFlag {
		UNDEFINED_TRACE_STATUS=-1,
		NOT_ON_TRACE=0,
		MIDDLE_TRACE=1,
		BEGINNING_TRACE=2,
		END_TRACE=3,
	};

	//! Represents a vertex in the EqSim file associated with a particular section.
	class EQSimGeometryVertex {
	private:
		int				_line_num;	//! line number this record was read from (-1 if not read)
		UIndex			_index;		//! Field 1: Vertex index number
		LatLonDepth		_loc;		//! Encapsulates fields 2-4
									//! Field 2: Latitude (decimal degrees, positive north)
									//! Field 3: Longitude (decimal degrees, positive east)
									//! Field 4: Depth (meters, negative underground)
		double			_das;		//! Field 5: Distance-along-strike (meters)
		TraceFlag		_trace_flag;//! Field 6: Trace flag (0 = not on trace, 1 = on trace but not initial or final,
									//! 2 = initial point on trace, 3 = final point on trace)
		
	public:
		EQSimGeometryVertex(void) : _line_num(-1), _index(0), _loc(), _das(0), _trace_flag(UNDEFINED_TRACE_STATUS) {};
		
		//! Return the line number this vertex was read on, -1 if not read.
		int line_num(void) const { return _line_num; };
		//! Return the index of this vertex.
		UIndex index(void) const { return _index; };
		//! Return a LatLonDepth object specifying the vertex location.
		LatLonDepth loc(void) const { return _loc; };
		//! Returns the distance-along-strike for this vertex (meters).
		double das(void) const { return _das; };
		//! Returns the trace flag information for this vertex.
		TraceFlag trace_flag(void) const { return _trace_flag; };
		
		//! Set the index of this vertex.
		void set_index(const UIndex &index) { _index = index; };
		//! Set the location of this vertex.
		void set_loc(const LatLonDepth &loc) { _loc = loc; };
		//! Set the distance-along-strike of this vertex.
		void set_das(const double &das) { _das = das; };
		//! Set the trace flag of this vertex.
		void set_trace_flag(const TraceFlag &trace_flag) { _trace_flag = trace_flag; };
		
		//! Parse a vertex record from the input stream to this object.
		void parse(const int &line_num, std::istringstream &line_stream);
		//! Write the record for this vertex to the specified output stream.
		void write(std::ostream &out_stream) const;
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
		
		//! Test for equality of vertices. Two vertices are equivalent if they have equal locations,
		//! distance along strike and trace flags.
		bool operator==(EQSimGeometryVertex &vertex) const { return (_loc == vertex._loc) && (_das == vertex._das) && (_trace_flag == vertex._trace_flag); };
		//! Test for inequality of vertices.
		bool operator!=(EQSimGeometryVertex &vertex) const { return (!(*this == vertex)); };
	};
	
	typedef std::map<UIndex, EQSimGeometryVertex> EQSimGeomVertexMap;
	typedef std::map<UIndex, UIndex> EQSimGeomIndexMap;
	
	//! Represents a set of remappings for vertices, triangles and rectangles.
	class EQSimGeomRemapping {
	public:
		//! Mapping of old vertex indices to new indices.
		EQSimGeomIndexMap	vert_remap;
		//! Mapping of old triangle indices to new indices.
		EQSimGeomIndexMap	tri_remap;
		//! Mapping of old rectangle indices to new indices.
		EQSimGeomIndexMap	rect_remap;
	};
	
	//! Represents a triangle in the EqSim file associated with a particular section.
	class EQSimGeometryTriangle {
	private:
		int				_line_num;		//! line number this record was read from (-1 if not read)
		UIndex			_index;			//! Field 1: Triangle index number
		UIndex			_vertex[3];		//! Fields 2-4: Vertex index number for corner #1-3
										//! (counting counterclockwise as viewed from positive side of element)
		double			_rake;			//! Field 5: Rake angle (decimal degrees, 0.0 = left lateral, 90.0 = positive side moves up)
		double			_slip_rate;		//! Field 6: Element slip rate (meters/second)
		double			_aseis_factor;	//! Field 7: Element aseismicity factor
		double			_strike;		//! Field 8: Strike angle (decimal degrees)
		double			_dip;			//! Field 9: Dip angle (decimal degrees)
		
	public:
		EQSimGeometryTriangle(void) : _line_num(-1), _index(0), _rake(0), _slip_rate(0), _aseis_factor(0), _strike(0), _dip(0)
				{ for (int i=0;i<3;++i) _vertex[i] = 0; };
		
		//! Return the line number this vertex was read on, -1 if not read.
		int line_num(void) const { return _line_num; };
		//! Return the index of this triangle.
		UIndex index(void) const { return _index; };
		
		//! Returns the index of vertices for the triangle.
		//! Throws an exception if the vertex number is outside [0,2].
		UIndex vertex(const unsigned int &vert) const throw (std::out_of_range) {
			if (vert > 2) throw std::out_of_range("quakelib::EQSimGeometryTriangle::vertex");
			return _vertex[vert];
		};
		//! Returns the rake angle of the triangle (in degrees, 0 = left lateral).
		double rake(void) const { return _rake; };
		//! Returns the slip rate of the triangle (meters per second).
		double slip_rate(void) const { return _slip_rate; };
		//! Returns the aseismic factor for this triangle.
		double aseismic(void) const { return _aseis_factor; };
		//! Returns the strike angle of the triangle (in degrees).
		double strike(void) const { return _strike; };
		//! Returns the dip angle of the triangle (in degrees).
		double dip(void) const { return _dip; };
		
		//! Set the index of this triangle.
		void set_index(const UIndex &index) { _index = index; };
		//! Set the indices for the triangle vertices.
		//! Throws an exception if the vertex number is outside [0,2].
		void set_vertex(const unsigned int &vert, const UIndex &vertex_ind) throw (std::out_of_range) {
			if (vert > 2) throw std::out_of_range("quakelib::EQSimGeometryTriangle::set_vertex");
			_vertex[vert] = vertex_ind;
		};
		//! Set the rake (in degrees) for this triangle.
		void set_rake(const double &rake) { _rake = rake; };
		//! Set the slip rate (in meters/second) for this triangle.
		void set_slip_rate(const double &slip_rate) { _slip_rate = slip_rate; };
		//! Set the aseismic factor for this triangle.
		void set_aseismic(const double &aseismic) { _aseis_factor = aseismic; };
		//! Set the strike angle (in degrees) for this triangle.
		void set_strike(const double &strike) { _strike = strike; };
		//! Set the dip angle (in degrees) for this triangle.
		void set_dip(const double &dip) { _dip = dip; };
		
		//! Apply the specified index remapping to this triangle.
		void apply_remap(const EQSimGeomRemapping &remap);
		
		//! Parse a triangle record from the input stream to this object.
		void parse(const int &line_num, std::istringstream &line_stream);
		//! Write the record for this triangle to the specified output stream.
		void write(std::ostream &out_stream) const;
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	typedef std::map<UIndex, EQSimGeometryTriangle> EQSimGeomTriangleMap;
	
	//! Represents a rectangle in the EqSim file associated with a particular section.
	class EQSimGeometryRectangle {
	private:
		int				_line_num;		//! line number this record was read from (-1 if not read)
		UIndex			_index;			//! Field 1: Rectangle index number
		UIndex			_vertex[4];		//! Fields 2-5: Vertex index number for corners #1-4
										//! (counting counterclockwise as viewed from positive side of element)
		double			_rake;			//! Field 6: Rake angle (decimal degrees, 0.0 = left lateral, 90.0 = positive side moves up)
		double			_slip_rate;		//! Field 7: Element slip rate (meters/second)
		double			_aseis_factor;	//! Field 8: Element aseismicity factor
		double			_strike;		//! Field 9: Strike angle (decimal degrees)
		double			_dip;			//! Field 10: Dip angle (decimal degrees)
		int				_perfect_flag;	//! Field 11: Perfect flag (0 = not perfect rectangle, 1 = perfect rectangle)
		
	public:
		EQSimGeometryRectangle(void) : _line_num(-1), _index(0), _rake(0), _slip_rate(0), _aseis_factor(0), _strike(0), _dip(0), _perfect_flag(-1) { for (int i=0;i<4;++i) _vertex[i] = 0; };
		
		//! Returns the line number this vertex was read on, -1 if not read.
		int line_num(void) const { return _line_num; };
		//! Returns the index of this rectangle.
		UIndex index(void) const { return _index; };
		
		// TODO: change these to 0 based index
		//! Returns the index of vertex 1 for the rectangle.
		UIndex vertex(const unsigned int &vert) const throw (std::out_of_range) {
			if (vert > 3) throw std::out_of_range("quakelib::EQSimGeometryRectangle::vertex");
			return _vertex[vert];
		};
		
		//! Returns the rake angle of the rectangle (in degrees, 0 = left lateral).
		double rake(void) const { return _rake; };
		//! Returns the slip rate of the rectangle (meters per second).
		double slip_rate(void) const { return _slip_rate; };
		//! Returns the aseismic factor for this rectangle.
		double aseismic(void) const { return _aseis_factor; };
		//! Returns the strike angle of the rectangle (in degrees).
		double strike(void) const { return _strike; };
		//! Returns the dip angle of the rectangle (in degrees).
		double dip(void) const { return _dip; };
		//! Returns whether this rectangle is flagged as perfect or not.
		int perfect_flag(void) const { return _perfect_flag; };
		
		//! Set the index of the rectangle.
		void set_index(const UIndex &index) { _index = index; };
		//! Set the indices for the rectangle vertices.
		void set_vertex(const unsigned int &vert, const UIndex &vertex_ind) throw (std::out_of_range) {
			if (vert > 3) throw std::out_of_range("quakelib::EQSimGeometryRectangle::set_vertex");
			_vertex[vert] = vertex_ind;
		};
		
		//! Set the rake (in degrees) for this rectangle.
		void set_rake(const double &rake) { _rake = rake; };
		//! Set the slip rate (in meters/second) for this rectangle.
		void set_slip_rate(const double &slip_rate) { _slip_rate = slip_rate; };
		//! Set the aseismic factor for this rectangle.
		void set_aseismic(const double &aseismic) { _aseis_factor = aseismic; };
		//! Set the strike angle (in degrees) for this rectangle.
		void set_strike(const double &strike) { _strike = strike; };
		//! Set the dip angle (in degrees) for this rectangle.
		void set_dip(const double &dip) { _dip = dip; };
		//! Set whether this rectangle is perfect.
		void set_perfect_flag(const int &perfect) { _perfect_flag = perfect; };
		
		//! Apply the specified index remapping to this rectangle.
		void apply_remap(const EQSimGeomRemapping &remap);
		
		//! Parse a rectangle record from the input stream to this object.
		void parse(const int &line_num, std::istringstream &line_stream);
		//! Write the record for this rectangle to the specified output stream.
		void write(std::ostream &out_stream) const;
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	typedef std::map<UIndex, EQSimGeometryRectangle> EQSimGeomRectangleMap;
	
	//! Class to find vertices by Cartesian location and distance.
	class VertexFinder {
	private:
		//! Map of which vertices are at which X coordinates
		std::multimap<double, UIndex>	x_vertex;
		//! Map of which vertices are at which Y coordinates
		std::multimap<double, UIndex>	y_vertex;
		//! Map of which vertices are at which Z coordinates
		std::multimap<double, UIndex>	z_vertex;
		
	public:
		//! Add a vertex index and associated Cartesian point to the finder
		void add_point(const UIndex &index, const Vec<3> &loc);
		//! Find a set of vertices within a specified distance of a Cartesian location.
		//! Indices of the vertices within this distance are put in the vertex_indices set.
		void find_vertices_in_distance(std::set<UIndex> &vertex_indices, const Vec<3> &loc, const double &distance);
	};
	
	//! Represents a complete triangular or rectangular element for use in a simulation.
	template <unsigned int nverts>
	class Element {
	private:
		//! Used to determine the intersection of convex objects in the method of separating axes.
		bool check_separation(const Vec<3> &sep_vec, const Element<nverts> &block) const;
		
	protected:
		//! Coordinates of element vertices in meters
		Vec<3>		_vert[nverts];
		//! Distance-along-strike (meters)
		double		_das[nverts];
		//! Trace flag (0 = not on trace, 1 = on trace but not initial or final,
		//! 2 = initial point on trace, 3 = final point on trace)
		TraceFlag	_trace_flag[nverts];
		//! Rake angle (radians, 0.0 = left lateral, PI/2 = positive side moves up)
		double		_rake;
		//! Element slip rate (meters/second)
		double		_slip_rate;
		//! Element aseismicity factor (in [0, 1])
		double		_aseis_factor;
		
	public:
		//! Calculate the stress tensor at a location with Lame parameters lambda and mu
		//! caused by this element moving unit_slip meters.
		Tensor<3,3> calc_stress_tensor(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);
		//! Calculate the derivatives at a location with Lame parameters lambda and mu
		//! caused by this element moving unit_slip meters.
		Vec<3> calc_dudx(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);
		Vec<3> calc_dudy(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);
		Vec<3> calc_dudz(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);
		//! Calculate the displacement vector at a location with Lame parameters lambda and mu
		//! caused by this element moving unit_slip meters.
		Vec<3> calc_displacement_vector(const Vec<3> &location, const double &unit_slip, const double &loc_lambda, const double &loc_mu) const throw(std::invalid_argument);
		
		//! Set the location of a vertex.
		void set_vert(const unsigned int &vert, const Vec<3> &new_vert) throw(std::out_of_range) { if(vert>=nverts) throw std::out_of_range("quakelib::Element::set_vert"); _vert[vert] = new_vert; };
		//! Get the location of a vertex.
		Vec<3> vert(const unsigned int &vert) const throw(std::out_of_range) { if(vert>=nverts) throw std::out_of_range("quakelib::Element::vert"); return _vert[vert]; };
		
		//! Set the distance-along-strike for a vertex.
		void set_das(const unsigned int &vert, const double &new_das) throw(std::out_of_range) { if(vert>=nverts) throw std::out_of_range("quakelib::Element::set_das"); _das[vert] = new_das; };
		//! Get the distance-along-strike for a vertex.
		double das(const unsigned int &vert) const throw(std::out_of_range) { if(vert>=nverts) throw std::out_of_range("quakelib::Element::das"); return _das[vert]; };
		
		//! Set the trace flags for a vertex.
		void set_trace_flag(const unsigned int &vert, const TraceFlag &new_flag) throw(std::out_of_range) { if(vert>=nverts) throw std::out_of_range("quakelib::Element::set_trace_flag"); _trace_flag[vert] = new_flag; };
		//! Get the trace flag for a vertex.
		TraceFlag trace_flag(const unsigned int &vert) const throw(std::out_of_range) { if(vert>=nverts) throw std::out_of_range("quakelib::Element::trace_flag"); return _trace_flag[vert]; };
		
		//! Set the slip rate in m/s for this block.
		void set_slip_rate(const double &new_slip_rate) { _slip_rate = new_slip_rate; };
		//! Get the slip rate in cm/year for this block.
		double slip_rate(void) const { return _slip_rate; };
		
		//! Set the rake angle of this block in radians.
		void set_rake(const double &new_rake) { _rake = new_rake; };
		//! Get the rake angle of this block in radians.
		double rake(void) const { return _rake; };
		
		//! Get the fraction of slip which is aseismic for this element.
		double aseismic(void) const { return _aseis_factor; };
		//! Set the fraction of slip	which is aseismic for this element.
		void set_aseismic(const double &new_aseismic) throw(std::invalid_argument) { if (new_aseismic < 0 || new_aseismic > 1) throw std::invalid_argument("quakelib::Element::set_aseismic"); _aseis_factor = new_aseismic; };
		
		//! Clear all variables for this element.
		void clear(void) { for (unsigned int i=0;i<nverts;++i) { _vert[i] = Vec<3>(); _das[i] = nan(""); _trace_flag[i] = UNDEFINED_TRACE_STATUS; } _rake = _slip_rate = _aseis_factor = nan(""); };
		
		//! Returns a unit vector along the direction of fault dip.
		double dip(void) const { Vec<3> v; v[2] = 1.0; return normal().vector_angle(v); };
		
		//! Returns unit vector along the direction of fault rake.
		Vec<3> rake_vector(void) const {
			Vec<3> vec;
			vec = (_vert[nverts-1]-_vert[0]).unit_vector();
			return vec.rotate_around_axis(normal()*(-1.0), _rake);
		};
		
		//! Get the height of this block (in m).
		double height(void) const { return max_depth()-min_depth(); };
		
		//! Get the average depth of this element (in m) as a simple average of vertices.
		double avg_depth(void) const { double d=0; for(unsigned int i=0;i<nverts;++i) d+=_vert[i][2]; return d/nverts; };
		
		//! Get the average DAS for this block.
		double avg_das(void) const { double d=0; for(unsigned int i=0;i<nverts;++i) d += _das[i]; return d/nverts; };
		
		//! Get center point of this block
		Vec<3> center() const { Vec<3> c; for (unsigned int i=0;i<nverts;++i) c += _vert[i]; return c * (1.0/nverts); };
		
		//! Returns true if any point of this block is on the trace, false otherwise
		bool on_trace() const { for (unsigned int i=0;i<nverts;++i) if (_trace_flag[i] != NOT_ON_TRACE) return true; return false; };
		
		//! Returns length along largest dimension
		double largest_dimension(void) const throw(std::domain_error) { if(nverts!=4) throw std::domain_error("quakelib::Element::largest_dimension"); return fmax((_vert[3]-_vert[0]).mag(),(_vert[1]-_vert[0]).mag()); };
		
		//! Calculates the area of this block based on vertices.
		//double get_area(void) const { assert(nverts==3 || nverts==4); Vec<3> a, b; a = _vert[nverts-2]-_vert[0]; b = _vert[nverts-1]-_vert[nverts-3]; return 0.5*a.cross(b).mag(); };
		double get_area(void) const { return (_vert[3]-_vert[0]).mag() * (_vert[1]-_vert[0]).mag(); };
		
		//! Determine minimum depth of this block.
		//! Note that depth is negative, so max depth returns the most negative point and min depth the least negative
		double min_depth(void) const { double d = -DBL_MAX; for (unsigned int i=0;i<nverts;i++) d = fmax(_vert[i][2], d); return d; };
		//! Determine maximum depth of this block.
		double max_depth(void) const { double d = DBL_MAX; for (unsigned int i=0;i<nverts;i++) d = fmin(_vert[i][2], d); return d; };
		
		//! Determine minimum distance-along-strike (DAS) of this block.
		double min_das(void) const { double d = DBL_MAX; for (unsigned int i=0;i<nverts;i++) d = fmin(_das[i], d); return d; };
		//! Determine maximum distance-along-strike (DAS) of this block.
		double max_das(void) const { double d = -DBL_MAX; for (unsigned int i=0;i<nverts;i++) d = fmax(_das[i], d); return d; };
		
		//! Calculates the Euclidean distance between the 3D midpoint of this block and another block.
		double center_distance(const Element<nverts> &a) const { return (a.center() - center()).mag(); };
		
		bool overlaps(const Element<nverts> &block) const;
		
		//! Whether this block is above another block (defined by the top of the blocks).
		bool is_above(const Element<nverts> &block) const { return (block.min_depth() <= min_depth()); };
		
		//! Get the normal unit vector to the plane of this element.
		Vec<3> normal(void) const { Vec<3> a,b; a=_vert[1]-_vert[0]; b=_vert[nverts-1]-_vert[0]; return a.cross(b).unit_vector(); };
	};
	
	typedef Element<3> ElementTri;
	typedef Element<4> ElementRect;
	
	//! Represents a geometry section (composed of vertices, triangles and rectangles) in the EqSim file.
	class EQSimGeometrySection {
	private:
		void apply_remap(const EQSimGeomRemapping &remap);
		
	protected:
		//! Index number of this section.
		UIndex			_sid;
		//! Name of this section.
		std::string		_name;
		//! Fault ID of this section. Note that multiple sections may have the same fault ID.
		UIndex			_fid;
		
	public:
		//! Map of indices to vertices.
		EQSimGeomVertexMap		vertices;
		//! Map of indices to triangles.
		EQSimGeomTriangleMap	triangles;
		//! Map of indices to rectangles.
		EQSimGeomRectangleMap	rectangles;
		
		EQSimGeometrySection(void) : _sid(0), _name("(empty)"), _fid(0) { vertices.clear(); triangles.clear(); rectangles.clear(); };
		
		//! Returns the number of vertices in this section.
		size_t num_vertices(void) const { return vertices.size(); };
		//! Returns the number of triangles in this section.
		size_t num_triangles(void) const { return triangles.size(); };
		//! Returns the number of rectangles in this section.
		size_t num_rectangles(void) const { return rectangles.size(); };
		
		//! Returns the highest (northernmost) latitude of all vertices in this section.
		double lat_hi(void) const;
		//! Returns the lowest (southernmost) latitude of all vertices in this section.
		double lat_lo(void) const;
		
		//! Returns the highest (easternmost) longitude of all vertices in this section.
		double lon_hi(void) const;
		//! Returns the lowest (westernmost) longitude of all vertices in this section.
		double lon_lo(void) const;
		
		//! Returns the highest depth of all vertices in this section (meters, negative underground).
		double depth_hi(void) const;
		//! Returns the lowest depth of all vertices in this section (meters, negative underground).
		double depth_lo(void) const;
		
		//! Returns the furthest distance-along-strike of all vertices in this section.
		double das_hi(void) const;
		//! Returns the shortest distance-along-strike of all vertices in this section.
		double das_lo(void) const;
		
		//! Returns the ID number of this section.
		UIndex sid(void) const { return _sid; };
		//! Returns the name of this section.
		std::string name(void) const { return _name; };
		//! Returns the fault ID. Note that multiple sections may have the same fault ID.
		UIndex fid(void) const { return _fid; };
		
		//! Set the ID number of this section.
		void set_sid(const UIndex &sid) { _sid = sid; };
		//! Set the name of this section. If the name contains a space or newline an exception will be raised.
		void set_name(const std::string &name) throw(std::invalid_argument) { if (name.find_first_of("\n\r ") != std::string::npos) throw std::invalid_argument("quakelib::EQSimGeometrySection::set_name."); _name = name; };
		//! Set the fault ID of this section. Note that multiple sections may have the same fault ID.
		void set_fid(const UIndex &fid) { _fid = fid; };
		
		//! Create a new vertex in this section and return a reference to it.
		//! If no more vertices can be created an exception is raised.
		EQSimGeometryVertex &new_vertex(void) throw(std::out_of_range);
		//! Create a new triangle in this section and return a reference to it.
		//! If no more triangles can be created an exception is raised.
		EQSimGeometryTriangle &new_triangle(void) throw(std::out_of_range);
		//! Create a new rectangle in this section and return a reference to it.
		//! If no more rectangles can be created an exception is raised.
		EQSimGeometryRectangle &new_rectangle(void) throw(std::out_of_range);
		
		//! Adds a rectangular element (with XYZ coordinates) into this section using the specified conversion.
		UIndex add_element(const ElementRect &elem, const Conversion &conv);
		
		//! TODO: write this
		void erase_vertex(const UIndex &ind);
		
		//! TODO: write this
		void fill_vertex_finder(const Conversion &conv, VertexFinder &finder);
		
		//! Determine if the four specified points constitute a perfect rectangle.
		//! A perfect rectangle is one where all points are coplanar and parallel
		//! sides are the same length (within the specified tolerance).
		bool is_perfect_rect(const Vec<3> v[4], const double &tolerance);
		
		//! Create a rectangle object from four input coordinates and a slip vector.
		//! The tolerance is used to determine whether the rectangle is perfect or not.
		//! Returns a reference to the newly created rectangle which has been inserted in the section.
		EQSimGeometryRectangle &rect_from_lat_lon_depth(const std::vector<LatLonDepth> &locs,
														const Vec<2> &slip,
														const double &tolerance=1.0e-2) throw (std::invalid_argument);
		//! TODO: write this
		void calculate_trace_das(void);
		//! TODO: write this
		void remove_duplicate_vertices(const double &distance);
		//! TODO: write this
		void remove_unused_vertices(void);
		//! TODO: write this
		UIndex remap_continguous(const UIndex &start_ind);
		
		//! Write all vertices, triangles and rectangles of this section to the specified output stream.
		void write(std::ostream &out_stream) const;
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	typedef std::map<UIndex, EQSimGeometrySection> EQSimGeomSectionMap;
	
	class EQSimGeometry;
	class EQSimGeometryReader;
	//! Parsed values for a geometry section specifying numbers of records and value limits in the section.
	class EQSimParsedGeometrySection {
	private:
		int				_line_num;	//! line number this record was read from (-1 if not read)
		UIndex			_sid;		//! Field 1: Section identification number (positive integer, may not be consecutive)
		std::string		_name;		//! Field 2: Section name
		unsigned int	_n_vertex;	//! Field 3: Total number of vertices in the section
		unsigned int	_n_triangle;//! Field 4: Total number of triangles in the section
		unsigned int	_n_rectangle;//! Field 5: Total number of rectangles in the section
		double			_lat_lo;	//! Field 6: Lowest value of latitude in the section (decimal degrees, positive north)
		double			_lat_hi;	//! Field 7: Highest value of latitude in the section (decimal degrees, positive north)
		double			_lon_lo;	//! Field 8: Lowest value of longitude in the section (decimal degrees, positive east)
		double			_lon_hi;	//! Field 9: Highest value of longitude in the section (decimal degrees, positive east)
		double			_depth_lo;	//! Field 10: Lowest value of depth in the section (meters, negative underground)
		double			_depth_hi;	//! Field 11: Highest value of depth in the section (meters, negative underground)
		double			_das_lo;	//! Field 12: Lowest value of distance-along-strike in the section (meters)
		double			_das_hi;	//! Field 13: Highest value of distance-along-strike in the section (meters)
		UIndex			_fid;		//! Field 14: Fault identification number (positive integer)
		
	public:
		EQSimParsedGeometrySection(void) : _line_num(-1), _sid(0), _name(""), _n_vertex(0), _n_triangle(0),
			_n_rectangle(0), _lat_lo(0), _lat_hi(0), _lon_lo(0), _lon_hi(0), _depth_lo(0), _depth_hi(0),
			_das_lo(0), _das_hi(0), _fid(0) {};
		
		//! Parse geometry section info from the input stream to this object.
		void parse(const int &line_num, std::istringstream &line_stream);
		
		friend class EQSimGeometry;
		friend class EQSimGeometryReader;
	};
	
	typedef std::map<UIndex, EQSimParsedGeometrySection> EQSimParsedGeomSectionMap;
	
	//! Represents EqSim model geometry.
	class EQSimGeometry {
	protected:
		//! Parsed geometry section data. This data is used in validate to ensure the file is consistent.
		EQSimParsedGeomSectionMap		_parsed_sections;
		
	public:
		//! Map of section IDs to sections in this model.
		EQSimGeomSectionMap				sections;
		
		//! Creates and returns a new section in this geometry with an index one greater than the highest section ID.
		EQSimGeometrySection &new_section(void) throw(std::out_of_range);
		
		//! Returns the highest (northernmost) latitude of all vertices in this model.
		double lat_hi(void) const;
		//! Returns the lowest (southernmost) latitude of all vertices in this model.
		double lat_lo(void) const;
		
		//! Returns the highest (easternmost) longitude of all vertices in this model.
		double lon_hi(void) const;
		//! Returns the lowest (westernmost) longitude of all vertices in this model.
		double lon_lo(void) const;
		
		//! Returns the highest depth of all vertices in this model (meters, negative underground).
		double depth_hi(void) const;
		//! Returns the lowest depth of all vertices in this model (meters, negative underground).
		double depth_lo(void) const;
		
		//! Returns the furthest distance-along-strike of all vertices in this model (meters).
		double das_hi(void) const;
		//! Returns the shortest distance-along-strike of all vertices in this model (meters).
		double das_lo(void) const;
		
		//! Returns the number of sections in this model.
		size_t num_sections(void) const { return sections.size(); };
		//! Returns the number of vertices in this model.
		size_t num_vertices(void) const;
		//! Returns the number of triangles in this model.
		size_t num_triangles(void) const;
		//! Returns the number of rectangles in this model.
		size_t num_rectangles(void) const;
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
		
		//! Write all sections and their associated records to the specified output stream.
		void write(std::ostream &out_stream) const;
	};
	
	//! Represents the full contents of a parsed EqSim file.
	class EQSimGeometryWriter : public EQSimGeometry, public EQSimMetadataWriter {
	public:
		EQSimGeometryWriter(void);
		
		//! Write the model metadata and geometry records to the associated file.
		void write(void) { EQSimMetadataWriter::write(); EQSimGeometry::write(out_stream); };
	};
	
	//! Represents the full contents of a parsed EqSim file.
	class EQSimGeometryReader : public EQSimGeometry, public EQSimMetadataReader {
	private:
		int					_line_num;	//! Line number this record was read from.
		unsigned int		_n_section;	//! Field 1: Total number of fault sections in the file.
		unsigned int		_n_vertex;	//! Field 2: Total number of vertices in the file.
		unsigned int		_n_triangle;//! Field 3: Total number of triangles in the file.
		unsigned int		_n_rectangle;//! Field 4: Total number of rectangles in the file.
		double				_lat_lo;	//! Field 5: Lowest value of latitude in the file (decimal degrees, positive north).
		double				_lat_hi;	//! Field 6: Highest value of latitude in the file (decimal degrees, positive north).
		double				_lon_lo;	//! Field 7: Lowest value of longitude in the file (decimal degrees, positive east).
		double				_lon_hi;	//! Field 8: Highest value of longitude in the file (decimal degrees, positive east).
		double				_depth_lo;	//! Field 9: Lowest value of depth in the file (meters, negative underground).
		double				_depth_hi;	//! Field 10: Highest value of depth in the file (meters, negative underground).
		
		//! The current section being read. This is for internal purposes only.
		UIndex				cur_section;
		
		void parse_fault_summary_record(const int &line_num, std::istringstream &line_stream);
		void parse_section_record(const int &line_num, std::istringstream &line_stream);
		void parse_vertex_record(const int &line_num, std::istringstream &line_stream);
		void parse_triangle_record(const int &line_num, std::istringstream &line_stream);
		void parse_rectangle_record(const int &line_num, std::istringstream &line_stream);
		
		virtual bool parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream);
		
	public:
		EQSimGeometryReader(void) : _line_num(-1), _n_section(0), _n_vertex(0), _n_triangle(0), _n_rectangle(0),
			_lat_lo(0), _lat_hi(0), _lon_lo(0), _lon_hi(0), _depth_lo(0), _depth_hi(0) {};
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		void validate(EQSimErrors &errors) const;
	};
	
	typedef std::map<SectionID, double> DynamicTriggeringMap;
	typedef std::map<SectionID, int> SlipScalingMap;
	
	//! Represents the records of a file defining the dynamic triggering and slip scaling on a section by section basis.
	class SectionParams {
	protected:
		// Record 201: Initial stress states
		//! A map of index numbers to initial stress data.
		DynamicTriggeringMap	_dyns;		// Field 2,3: Element initial shear/normal stress (Pascal)
		
		// Record 202: Initial state
		//! A map of index numbers to rate-state data.
		SlipScalingMap	_sts;		// Field 2: Element initial state for rate-state law (seconds)
		
		//! Lines of parsed initial conditions (-1 if not parsed).
		//LineNumbers	_stress_lines;
		//! Lines of parsed rate-state data (-1 if not parsed).
		//LineNumbers	_rst_lines;
		
	public:
		//! Returns the number of elements found in this condition data.
		//size_t num_elements(void) const;
		//! Get the initial shear and normal stress (in Pascals) of the specified element.
		//! If the element is not found an exception is raised.
		//double get_shear_stress(const UIndex &ind) const throw(std::out_of_range);
		//double get_normal_stress(const UIndex &ind) const throw(std::out_of_range);
		//! Set the initial shear and normal stress (in Pascals) of the specified element.
		//void set_stresses(const UIndex &ind, const double &shear_stress, const double &normal_stress);
		
		//! Returns the rate-state data for the specified index.
		//double get_rate_state(const UIndex &ind) const throw(std::out_of_range);
		//! Sets the rate-state data for the specified index.
		//void set_rate_state(const UIndex &ind, const double &rate_state);
		
		//! Write the condition data to the specified output stream.
		//oid write(std::ostream &out_stream) const;
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		//void validate(EQSimErrors &errors) const;
		
		//! Returns the dynamic triggering value for the specified section.
		double get_dyn_val(const SectionID &ind) const throw(std::out_of_range);
		
		//! Returns the slip scaling value for the specified section.
		double get_st_val(const SectionID &ind) const throw(std::out_of_range);
		
	};
	
	//! Represents the records of an EqSim style file describing the initial stresses and rate state laws for a given model.
	class SectionParamsReader : public SectionParams {
	private:
		// Record 200: Initial condition summary
		int				_line_num;
		//unsigned int	_n_element;		// Field 1: Total number of elements in the file
		//bool			_stress_flag;	// Field 2: 1 if initial stress (record 201) is included, 0 if not
		//bool			_state_flag;		// Field 3: 1 if initial state (record 202) is included, 0 if not
		//! Whether parsing is completed for the input file metadata.
		bool				finish_parsing;
		
		//void parse_summary_record(const int &line_num, std::istringstream &line_stream);
		//void parse_stress_record(const int &line_num, std::istringstream &line_stream);
		//void parse_state_record(const int &line_num, std::istringstream &line_stream);
		
		void parse_line(const int &line_num, std::istringstream &line_stream);
		bool parse_eof(const int &line_num, std::istringstream &line_stream);
		
	public:
		//! Parses an input stream and calls the appropriate callback functions.
		//! The callback function to use for a given line is defined by the 3 digit
		//! code at the line start. Returns true if the file exists.
		// TODO: throw exception if problem reading file (e.g. doesn't exist)
		bool parse_file(const std::string &file_name, const int &start_line_num=0);
		
		//! Perform internal correctness checks and write any errors to the specified Errors object.
		//void validate(EQSimErrors &errors) const;
	};

}

#endif
