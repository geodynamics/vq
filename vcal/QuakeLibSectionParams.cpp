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

bool quakelib::SectionParamsReader::parse_file(const std::string &file_name, const int &start_line_num) {
	std::string			next_line;
	std::istringstream	line_stream;
	int					line_num;
	std::ifstream		input_file;
	//bool				found_match;
	
	input_file.open(file_name.c_str());
	if (!input_file.good()) return false;
	
	line_num = start_line_num;
	finish_parsing = false;
	while (input_file.good() && !finish_parsing) {
		// Increment the line count
		++line_num;
		
		// Get the next line in the file
		getline(input_file, next_line);
		
		// Find the record number for the line
		line_stream.str(next_line);
		//line_stream >> section_num;
		
		parse_line(line_num, line_stream);
		
		/*
		// Call the corresponding line processor with the remainder of the line
		try {
			found_match = parse_metadata_line(record_num, line_num, line_stream);
			if (!found_match) found_match = parse_line(record_num, line_num, line_stream);
			if (!found_match) {
				std::stringstream		error_msg;
				error_msg << "Unexpected record number: " << record_num;
				parse_errors.report(error_msg.str(), line_num);
			}
		} catch (std::exception e) {
			std::stringstream		error_msg;
			error_msg << "Unexpected error: " << e.what();
			parse_errors.report(error_msg.str(), line_num);
		}
		*/
		line_stream.clear();
		//record_num = 0;
	}
	
	input_file.close();
	return true;
}


void quakelib::SectionParamsReader::parse_line(const int &line_num, std::istringstream &line_stream) {
	SectionID	sid;
	double		dyn;
	int			st;
	
	line_stream >> sid				// Field 1: Element index number (consecutive integers, starting with 1)
				>> dyn
				>> st;			// Field 2: Element initial state for rate-state law (seconds)
	
	_dyns.insert(std::make_pair(sid, dyn));
	_sts.insert(std::make_pair(sid, st));

}

/*
void quakelib::SectionParamsReader::parse_eof(const int &line_num, std::istringstream &line_stream) {
	std::string str;
	line_stream >> str;
	if (str.compare("End")) parse_errors.report("End record has incorrect format.", line_num);
	finish_parsing = true;
}
*/
/*
void quakelib::SectionParamsReader::parse_summary_record(const int &line_num, std::istringstream &line_stream) {
	_line_num = line_num;
	line_stream >> _n_element		// Field 1: Total number of elements in the file
				>> _stress_flag		// Field 2: 1 if initial stress (record 201) is included, 0 if not
				>> _state_flag;		// Field 3: 1 if initial state (record 202) is included, 0 if not
}

void quakelib::SectionParamsReader::parse_stress_record(const int &line_num, std::istringstream &line_stream) {
	UIndex		ind;
	double		shear_str, normal_str;
	
	line_stream >> ind				// Field 1: Element index number (consecutive integers, starting with 1)
				>> shear_str		// Field 2: Element initial shear stress (Pascal)
				>> normal_str;		// Field 3: Element initial normal stress (Pascal)
	
	_stresses.insert(std::make_pair(ind, std::make_pair(shear_str, normal_str)));
}

void quakelib::SectionParamsReader::parse_state_record(const int &line_num, std::istringstream &line_stream) {
	UIndex		ind;
	double		rs_t;
	
	line_stream >> ind				// Field 1: Element index number (consecutive integers, starting with 1)
				>> rs_t;			// Field 2: Element initial state for rate-state law (seconds)
	
	_rs_theta.insert(std::make_pair(ind, rs_t));
}
*/
/*
// Get the number of elements in this condition object.
size_t quakelib::SectionParams::num_elements(void) const {
	if (_stresses.size() > 0) return _stresses.size();
	else return _rs_theta.size();
}
*/

/*!
 Get the dynamic triggering and slip scaling for a given section.
 */
double quakelib::SectionParams::get_dyn_val(const SectionID &ind) const throw(std::out_of_range) {
	DynamicTriggeringMap::const_iterator it = _dyns.find(ind);
	if (it == _dyns.end()) throw std::out_of_range("quakelib::SectionParams::get_dyn_val");
	return it->second;
}

double quakelib::SectionParams::get_st_val(const SectionID &ind) const throw(std::out_of_range) {
	SlipScalingMap::const_iterator it = _sts.find(ind);
	if (it == _sts.end()) throw std::out_of_range("quakelib::SectionParams::get_st_val");
	return it->second;
}

/*!
 Set the shear and normal stresses of a given element.
 
void quakelib::SectionParams::set_stresses(const UIndex &ind, const double &shear_stress, const double &normal_stress) {
	_stress_lines.insert(std::make_pair(ind, -1));
	_stresses.insert(std::make_pair(ind, std::make_pair(shear_stress, normal_stress)));
}
*/
/*
 Get the rate-state value of a given element.
 
double quakelib::SectionParams::get_rate_state(const UIndex &ind) const throw(std::out_of_range) {
	StateMap::const_iterator it = _rs_theta.find(ind);
	if (it == _rs_theta.end()) throw std::out_of_range("quakelib::SectionParams::get_rate_state");
	return it->second;
}
*/
/*!
 Set the rate-state value of a given element.

void quakelib::SectionParams::set_rate_state(const UIndex &ind, const double &rate_state) {
	_rst_lines.insert(std::make_pair(ind, -1));
	_rs_theta.insert(std::make_pair(ind, rate_state));
}
 */
/*
// TODO: interleave elements according to specification
void quakelib::SectionParams::write(std::ostream &out_stream) const {
	StressMap::const_iterator	sit;
	StateMap::const_iterator	rit;
	
	if (_stresses.size() > 0 && _rs_theta.size() > 0 && _stresses.size() != _rs_theta.size())
		throw std::length_error("quakelib::SectionParams::write");
	
	// Write the header record
	out_stream	<< "200"
				<< " " << _stresses.size()
				<< " " << (_stresses.size() > 0)
				<< " " << (_rs_theta.size() > 0)
				<< std::endl;
	
	// Write stress information
	for (sit=_stresses.begin();sit!=_stresses.end();++sit) {
		out_stream	<< "201"
					<< " " << sit->first
					<< " " << sit->second.first
					<< " " << sit->second.second
					<< std::endl;
	}
	
	// Write rate-state information
	for (rit=_rs_theta.begin();rit!=_rs_theta.end();++rit) {
		out_stream	<< "202"
					<< " " << rit->first
					<< " " << rit->second
					<< std::endl;
	}
}

void quakelib::SectionParams::validate(EQSimErrors &errors) const {
	StressMap::const_iterator	sit;

	// Check that no stresses are negative
	for (sit=_stresses.begin();sit!=_stresses.end();++sit) {
		if (sit->second.first < 0 || sit->second.second < 0) {
			std::stringstream		error_msg;
			error_msg << "Element " << sit->first << " has negative initial stress (shear="
			<< sit->second.first << ", normal=" << sit->second.second << ").";
			errors.report(error_msg.str(), _stress_lines.find(sit->first)->second);
		}
	}
}
*/
/*!
 Confirms that the condition file represented by this class is self consistent.
 To do so, it ensures that the number of shear and normal stress elements
 actually stored is the same as the number specified in the file. It also
 does the same check for the number of rate-state entries.

void quakelib::SectionParamsReader::validate(EQSimErrors &errors) const {
	EQSimMetadataReader::validate(errors);
	SectionParams::validate(errors);
	
	if (_stress_flag) {
		// If stress flag was true, make sure the number of entries was as specified
		if (_stresses.size() != _n_element) {
			std::stringstream		error_msg;
			error_msg << "Specified number of stress entries (" << _n_element
			<< ") does not match actual number in file (" << _stresses.size() << ").";
			errors.report(error_msg.str(), _line_num);
		}
	} else {
		// Check the shear and normal stress entries
		if (_stresses.size() > 0) {
			std::stringstream		error_msg;
			error_msg << "Stress flag indicates no stress entries, but found "
			<< _stresses.size() << " in the file.";
			errors.report(error_msg.str(), _line_num);
		}
	}

	if (_state_flag) {
		// If state flag was true, make sure the number of entries was as specified
		if (_rs_theta.size() != _n_element) {
			std::stringstream		error_msg;
			error_msg << "Specified number of rate-state entries (" << _n_element
			<< ") does not match actual number in file (" << _rs_theta.size() << ").";
			errors.report(error_msg.str(), _line_num);
		}
	} else {
		// Check the rate-state entries
		if (_rs_theta.size() > 0) {
			std::stringstream		error_msg;
			error_msg << "State flag indicates no rate-state entries, but found "
			<< _rs_theta.size() << " in the file.";
			errors.report(error_msg.str(), _line_num);
		}
	}
}
*/
