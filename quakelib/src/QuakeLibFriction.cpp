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

quakelib::EQSimFrictionWriter::EQSimFrictionWriter(void) {
	internal::RecordDesc	fric_summary_rec, elastic_params_rec, fault_strength_rec, rate_state_rec;
	
	set_spec_level(1);
	meta_add_record(META_SIGNATURE, "EQSim_Input_Friction_2");
	
	fric_summary_rec = internal::RecordDesc(0, "summary", 4, "Record 200: Fault friction summary");
	fric_summary_rec.add_field(1, internal::FieldDesc("n_element", 1, 1, "Total number of elements in the file"));
	fric_summary_rec.add_field(2, internal::FieldDesc("elastic_flag", 1, 2, "1 if elastic parameters (record 201) are included, 0 if not"));
	fric_summary_rec.add_field(3, internal::FieldDesc("strength_flag", 1, 3, "1 if fault strength (record 202) is included, 0 if not"));
	fric_summary_rec.add_field(4, internal::FieldDesc("rate_state_flag", 1, 4, "1 if rate-state parameters (record 203) are included, 0 if not"));
	elastic_params_rec = internal::RecordDesc(0, "elastic_param", 2, "Record 201: Elastic parameters");
	elastic_params_rec.add_field(1, internal::FieldDesc("lame_lambda", 2, 1, "Lame modulus lambda (Pascal)"));
	elastic_params_rec.add_field(2, internal::FieldDesc("lame_mu", 2, 2, "Lame modulus mu, also known as the shear modulus (Pascal)"));
	fault_strength_rec = internal::RecordDesc(0, "fault_strength", 3, "Record 202: Fault strength");
	fault_strength_rec.add_field(1, internal::FieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
	fault_strength_rec.add_field(2, internal::FieldDesc("static_strength", 2, 2, "Element static yield strength (Pascal)"));
	fault_strength_rec.add_field(3, internal::FieldDesc("dynamic_strength", 2, 3, "Element dynamic sliding strength (Pascal)"));
	rate_state_rec = internal::RecordDesc(0, "rate_state", 6, "Record 203: Rate-state parameters");
	rate_state_rec.add_field(1, internal::FieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
	rate_state_rec.add_field(2, internal::FieldDesc("A", 2, 2, "Element rate-state parameter A"));
	rate_state_rec.add_field(3, internal::FieldDesc("B", 2, 3, "Element rate-state parameter B"));
	rate_state_rec.add_field(4, internal::FieldDesc("L", 2, 4, "Element rate-state characteristic distance L (meters)"));
	rate_state_rec.add_field(5, internal::FieldDesc("f0", 2, 5, "Element rate-state friction coefficient f0"));
	rate_state_rec.add_field(6, internal::FieldDesc("V0", 2, 6, "Element rate-state reference velocity V0 (meters/second)"));
	
	add_record_desc_record(200, fric_summary_rec);
	add_record_desc_record(201, elastic_params_rec);
	add_record_desc_record(202, fault_strength_rec);
	add_record_desc_record(203, rate_state_rec);
}

void quakelib::EQSimFriction::write(std::ostream &out_stream) const {
	StrengthMap::const_iterator		sit;
	RateStateMap::const_iterator	rit;
	
	if (_strengths.size() > 0 && _rate_state.size() > 0 && _strengths.size() != _rate_state.size())
		throw std::length_error("stresses");
	
	// Write the summary
	out_stream	<< "200"
				<< " " << _strengths.size()
				<< " " << _have_elastic_params
				<< " " << (_strengths.size() > 0)
				<< " " << (_rate_state.size() > 0)
				<< std::endl;
	
	// Write the elastic parameters
	if (_have_elastic_params) {
		out_stream	<< "201"
					<< " " << _lame_lambda
					<< " " << _lame_mu
					<< std::endl;
	}
	
	// Write strength information
	for (sit=_strengths.begin();sit!=_strengths.end();++sit) {
		out_stream	<< "202"
					<< " " << sit->first
					<< " " << sit->second.first
					<< " " << sit->second.second
					<< std::endl;
	}
	
	// Write rate-state information
	for (rit=_rate_state.begin();rit!=_rate_state.end();++rit) {
		out_stream<< "203 " << rit->first << " ";
		rit->second.write(out_stream);
		out_stream << std::endl;
	}
}

void quakelib::EQSimFrictionRateState::write(std::ostream &out_stream) const {
	out_stream	<< _A
				<< " " << _B
				<< " " << _L
				<< " " << _f0
				<< " " << _V0;
}

bool quakelib::EQSimFrictionReader::parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream) {
	switch (rec_num) {
		case 200: parse_summary_record(line_num, line_stream); return true;
		case 201: parse_elastic_params_record(line_num, line_stream); return true;
		case 202: parse_fault_strength_record(line_num, line_stream); return true;
		case 203: parse_rate_state_record(line_num, line_stream); return true;
		default: return false;
	}
}

void quakelib::EQSimFrictionReader::parse_summary_record(const int &line_num, std::istringstream &line_stream) {
	_line_num = line_num;
	line_stream >> _num_element		// Field 1: Total number of elements in the file
				>> _elastic_flag	// Field 2: 1 if elastic parameters (record 201) are included, 0 if not
				>> _strength_flag	// Field 3: 1 if fault strength (record 202) is included, 0 if not
				>> _rate_state_flag;// Field 4: 1 if rate-state parameters (record 203) are included, 0 if not
}

void quakelib::EQSimFrictionReader::parse_elastic_params_record(const int &line_num, std::istringstream &line_stream) {
	_line_num = line_num;
	_have_elastic_params = true;
	line_stream >> _lame_lambda		// Field 1: Lame modulus lambda (Pascal)
				>> _lame_mu;		// Field 2: Lame modulus mu, also known as the shear modulus (Pascal)
}

void quakelib::EQSimFrictionReader::parse_fault_strength_record(const int &line_num, std::istringstream &line_stream) {
	UIndex		ind;
	double		stat_str, dyn_str;
	
	line_stream >> ind				// Field 1: Element index number (consecutive integers, starting with 1)
				>> stat_str			// Field 2: Element static yield strength (Pascal)
				>> dyn_str;			// Field 3: Element dynamic sliding strength (Pascal)
	
	_str_lines.insert(std::make_pair(ind, _line_num));
	_strengths.insert(std::make_pair(ind, std::make_pair(stat_str, dyn_str)));
}

void quakelib::EQSimFrictionRateState::parse(const int &line_num, std::istringstream &line_stream) {
	line_stream >> _A			// Field 2: Element rate-state parameter A
				>> _B			// Field 3: Element rate-state parameter B
				>> _L			// Field 4: Element rate-state characteristic distance L (meters)
				>> _f0			// Field 5: Element rate-state friction coefficient f0
				>> _V0;			// Field 6: Element rate-state reference velocity V0 (meters/second)
}

void quakelib::EQSimFrictionReader::parse_rate_state_record(const int &line_num, std::istringstream &line_stream) {
	EQSimFrictionRateState	new_rate_state;
	UIndex					_index;
	
	line_stream >> _index;		// Field 1: Element index number (consecutive integers, starting with 1)
	new_rate_state.parse(line_num, line_stream);
	
	_rs_lines.insert(std::make_pair(_index, _line_num));
	_rate_state.insert(std::make_pair(_index, new_rate_state));
}

double quakelib::EQSimFriction::get_static_strength(const UIndex &ind) const throw(std::out_of_range) {
	StrengthMap::const_iterator it = _strengths.find(ind);
	if (it == _strengths.end()) throw std::out_of_range("quakelib::EQSimFriction::get_static_strength");
	return it->second.first;
}

double quakelib::EQSimFriction::get_dynamic_strength(const UIndex &ind) const throw(std::out_of_range) {
	StrengthMap::const_iterator it = _strengths.find(ind);
	if (it == _strengths.end()) throw std::out_of_range("quakelib::EQSimFriction::get_dynamic_strength");
	return it->second.second;
}

void quakelib::EQSimFriction::set_strengths(const UIndex &ind, const double &static_strength, const double &dynamic_strength) {
	_strengths.insert(std::make_pair(ind, std::make_pair(static_strength, dynamic_strength)));
}

quakelib::EQSimFrictionRateState quakelib::EQSimFriction::get_rs_param(const UIndex &ind) const {
	RateStateMap::const_iterator it = _rate_state.find(ind);
	if (it == _rate_state.end()) throw std::out_of_range("quakelib::EQSimFriction::rate_state_params");
	return it->second;
}

void quakelib::EQSimFriction::set_rs_param(const UIndex &ind, const EQSimFrictionRateState &params) {
	_rate_state.insert(std::make_pair(ind, params));
}

void quakelib::EQSimFriction::validate(EQSimErrors &errors) const {
	StrengthMap::const_iterator		it;
	
	// Ensure correct parameter values
	if (_lame_mu <= 0) {
		std::stringstream		error_msg;
		error_msg << "Lame mu parameter (" << _lame_mu << ") should be greater than 0.";
		errors.report(error_msg.str(), _line_num);
	}
	
	// Ensure strengths have valid values
	for (it=_strengths.begin();it!=_strengths.end();++it) {
		if (it->second.first < 0) {
			std::stringstream		error_msg;
			error_msg << "Static yield strength (" << it->second.first << ") should not be below 0.";
			errors.report(error_msg.str(), _str_lines.find(it->first)->second);
		}
		if (it->second.second < 0) {
			std::stringstream		error_msg;
			error_msg << "Dynamic yield strength (" << it->second.second << ") should not be below 0.";
			errors.report(error_msg.str(), _str_lines.find(it->first)->second);
		}
	}
	
	// TODO: check rate-state params
}

void quakelib::EQSimFrictionReader::validate(EQSimErrors &errors) const {
	EQSimMetadataReader::validate(errors);
	EQSimFriction::validate(errors);
	
	// If strength_flag is false, make sure there were no strength entries
	if (_strength_flag) {
		// If strength flag was true, make sure the number of entries was as specified
		if (_strengths.size() != _num_element) {
			std::stringstream		error_msg;
			error_msg << "Specified number of strength entries (" << _num_element
			<< ") does not match actual number in file (" << _strengths.size() << ").";
			errors.report(error_msg.str(), _line_num);
		}
	} else {
		// Check the static and dynamic strength entries
		if (_strengths.size() > 0) {
			std::stringstream		error_msg;
			error_msg << "Strength flag indicates no strength entries, but found "
			<< _strengths.size() << " in the file.";
			errors.report(error_msg.str(), _line_num);
		}
	}

	// If rate_state_flag is false, make sure there were no rate_state entries
	if (_rate_state_flag) {
		// If strength flag was true, make sure the number of entries was as specified
		if (_rate_state.size() != _num_element) {
			std::stringstream		error_msg;
			error_msg << "Specified number of rate-state entries (" << _num_element
			<< ") does not match actual number in file (" << _rate_state.size() << ").";
			errors.report(error_msg.str(), _line_num);
		}
	} else {
		// Check the rate-state entries
		if (_rate_state.size() > 0) {
			std::stringstream		error_msg;
			error_msg << "Rate-state flag indicates no rate-state entries, but found "
			<< _rate_state.size() << " in the file.";
			errors.report(error_msg.str(), _line_num);
		}
	}
}
