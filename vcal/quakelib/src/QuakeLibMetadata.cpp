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

bool quakelib::EQSimMetadataReader::parse_metadata_line(const int &rec_num, const int &line_num, std::istringstream &line_stream) {
	// If we haven't parsed the signature yet, we must do that first
	if (!parsed_sig) {
		if (rec_num == 101) {
			parse_signature(line_num, line_stream);
			return true;
		} else {
			return false;
		}
	}
	
	// If we're past the signature but not yet past the metadata, parse the appropriate lines
	if (!parse_metadata_finished) {
		switch (rec_num) {
			case 100: parse_comment(line_num, line_stream); return true;
			case 102: parse_end_metadata(line_num, line_stream); return true;
			case 110: parse_info_record(line_num, line_stream); return true;
			case 111: parse_title_record(line_num, line_stream); return true;
			case 112: parse_author_record(line_num, line_stream); return true;
			case 113: parse_date_record(line_num, line_stream); return true;
			case 999: parse_eof(line_num, line_stream); return true;
			default: return false;
		}
	}
	
	// Otherwise, assuming we're not done parsing descriptors, we read them in
	if (!parsed_descriptors) {
		switch (rec_num) {
			case 100: parse_comment(line_num, line_stream); return true;
			case 103: parse_end_descriptor(line_num, line_stream); return true;
			case 120: parse_record_desc_record(line_num, line_stream); return true;
			case 121: parse_field_desc_record(line_num, line_stream); return true;
			case 999: parse_eof(line_num, line_stream); return true;
			default: return false;
		}
	}
	
	// Once we're past the descriptors, we can only accept EOF and comments
	switch (rec_num) {
		case 100: parse_comment(line_num, line_stream); return true;
		case 999: parse_eof(line_num, line_stream); return true;
		default: return false;
	}
}

bool quakelib::EQSimMetadataReader::parse_file(const std::string &file_name, const int &start_line_num) {
	std::string			next_line;
	std::istringstream	line_stream;
	int					record_num, line_num;
	std::ifstream		input_file;
	bool				found_match;
	
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
		line_stream >> record_num;
		
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
		line_stream.clear();
		record_num = 0;
	}
	
	input_file.close();
	return true;
}

bool quakelib::internal::EQSimMetadata::has_bad_chars(const std::string &test_str) {
	if (test_str.find_first_of("\n\r") != std::string::npos) return true;
	return false;
}

unsigned int quakelib::internal::EQSimMetadata::meta_num_records(const RECORD_TYPE &rec_type) const throw(std::invalid_argument) {
	switch (rec_type) {
		case META_COMMENT:
		case META_INFO:
		case META_AUTHOR:
		case META_TITLE:
		case META_SIGNATURE:
		case META_DATE:
			return metadata_recs.find(rec_type)->second.size();
		default:
			throw std::invalid_argument("quakelib::EQSimMetadata::meta_num_records");
	}
}

void quakelib::internal::EQSimMetadata::meta_clear_record(const RECORD_TYPE &rec_type) throw(std::invalid_argument) {
	switch (rec_type) {
		case META_COMMENT:
		case META_INFO:
		case META_AUTHOR:
		case META_TITLE:
		case META_SIGNATURE:
		case META_DATE:
			metadata_recs[rec_type].clear();
			metadata_line_nums[rec_type].clear();
			break;
		default:
			throw std::invalid_argument("quakelib::EQSimMetadata::meta_clear_record");
	}
}

std::string quakelib::internal::EQSimMetadata::meta_get_record(const RECORD_TYPE &rec_type, const unsigned int &rec_num) const throw(std::invalid_argument, std::out_of_range) {
	switch (rec_type) {
		case META_COMMENT:
		case META_INFO:
		case META_AUTHOR:
		case META_TITLE:
		case META_SIGNATURE:
		case META_DATE:
			if (rec_num >= metadata_recs.find(rec_type)->second.size()) throw std::out_of_range("quakelib::EQSimMetadata::meta_get_record");
			return metadata_recs.find(rec_type)->second.at(rec_num);
		default:
			throw std::invalid_argument("quakelib::EQSimMetadata::meta_get_record");
	}
}

void quakelib::internal::EQSimMetadata::meta_add_record(const RECORD_TYPE &rec_type, const std::string &new_rec) throw(std::invalid_argument) {
	if (has_bad_chars(new_rec)) throw std::invalid_argument("quakelib::EQSimMetadata::meta_add_record");
	switch (rec_type) {
		case META_COMMENT:
		case META_INFO:
		case META_AUTHOR:
			metadata_recs.find(rec_type)->second.push_back(new_rec);
			metadata_line_nums.find(rec_type)->second.push_back(-1);
			break;
		case META_TITLE:
		case META_SIGNATURE:
		case META_DATE:
			metadata_recs.find(rec_type)->second.clear();
			metadata_line_nums.find(rec_type)->second.clear();
			metadata_recs.find(rec_type)->second.push_back(new_rec);
			metadata_line_nums.find(rec_type)->second.push_back(-1);
			break;
		default:
			throw std::invalid_argument("quakelib::EQSimMetadata::meta_add_record");
	}
}

void quakelib::internal::EQSimMetadata::meta_set_record(const RECORD_TYPE &rec_type, const unsigned int &rec_num, const std::string &new_rec) throw(std::invalid_argument, std::out_of_range) {
	switch (rec_type) {
		case META_COMMENT:
		case META_INFO:
		case META_AUTHOR:
			if (rec_num >= metadata_recs.find(rec_type)->second.size()) throw std::out_of_range("quakelib::EQSimMetadata::meta_set_record");
			metadata_recs.find(rec_type)->second.at(rec_num) = new_rec;
			metadata_line_nums.find(rec_type)->second.at(rec_num) = -1;
			break;
		case META_TITLE:
		case META_SIGNATURE:
		case META_DATE:
			if (rec_num >= metadata_recs.find(rec_type)->second.size()) throw std::out_of_range("quakelib::EQSimMetadata::meta_set_record");
			metadata_recs.find(rec_type)->second.clear();
			metadata_line_nums.find(rec_type)->second.clear();
			metadata_recs.find(rec_type)->second.push_back(new_rec);
			metadata_line_nums.find(rec_type)->second.push_back(-1);
			break;
		default:
			throw std::invalid_argument("quakelib::EQSimMetadata::meta_set_record");
	}
}

void quakelib::internal::EQSimMetadata::meta_erase_record(const RECORD_TYPE &rec_type, const unsigned int &rec_num) throw(std::invalid_argument, std::out_of_range) {
	switch (rec_type) {
		case META_COMMENT:
		case META_INFO:
		case META_AUTHOR:
		case META_TITLE:
		case META_SIGNATURE:
		case META_DATE:
			if (rec_num >= metadata_recs.find(rec_type)->second.size()) throw std::out_of_range("quakelib::EQSimMetadata::meta_erase_record");
			metadata_recs.find(rec_type)->second.erase(metadata_recs.find(rec_type)->second.begin()+rec_num);
			metadata_line_nums.find(rec_type)->second.erase(metadata_line_nums.find(rec_type)->second.begin()+rec_num);
			break;
		default:
			throw std::invalid_argument("quakelib::EQSimMetadata::meta_erase_record");
	}
}

void quakelib::internal::EQSimMetadata::add_record_desc_record(const unsigned int &rec_key, internal::RecordDesc &rec_desc) throw(std::invalid_argument) {
	if(rec_key<200 || rec_key>=400) throw std::invalid_argument("quakelib::EQSimMetadata::add_record_desc_record");
	record_descs.insert(std::make_pair(rec_key, rec_desc));
}

void quakelib::EQSimMetadataReader::parse_comment(const int &line_num, std::istringstream &line_stream) {
	std::string		comment;
	getline(line_stream, comment);
	comment.erase(0,1);		// erase the extra space getline adds to the beginning
	metadata_recs.find(META_COMMENT)->second.push_back(comment);
	metadata_line_nums.find(META_COMMENT)->second.push_back(line_num);
}

void quakelib::EQSimMetadataReader::parse_signature(const int &line_num, std::istringstream &line_stream) {
	std::string		sig_rec;
	
	line_stream >> sig_rec
				>> spec_level;
	
	metadata_recs.find(META_SIGNATURE)->second.clear();
	metadata_line_nums.find(META_SIGNATURE)->second.clear();
	metadata_recs.find(META_SIGNATURE)->second.push_back(sig_rec);
	metadata_line_nums.find(META_SIGNATURE)->second.push_back(line_num);
	
	parsed_sig = true;
}

void quakelib::EQSimMetadataReader::parse_end_metadata(const int &line_num, std::istringstream &line_stream) {
	std::string str;
	line_stream >> str;
	if (str.compare("End_Metadata")) parse_errors.report("Metadata end record has incorrect format.", line_num);
	if (parse_metadata_finished) parse_errors.report("Metadata end record previously found.", line_num);
	parse_metadata_finished = true;
}

void quakelib::EQSimMetadataReader::parse_end_descriptor(const int &line_num, std::istringstream &line_stream) {
	std::string str;
	line_stream >> str;
	if (str.compare("End_Descriptor")) parse_errors.report("Descriptor end record has incorrect format.", line_num);
	parsed_descriptors = true;
}

void quakelib::EQSimMetadataReader::parse_info_record(const int &line_num, std::istringstream &line_stream) {
	std::string new_info;
	getline(line_stream, new_info);
	new_info.erase(0,1);		// erase the extra space getline adds to the beginning
	metadata_recs.find(META_INFO)->second.push_back(new_info);
	metadata_line_nums.find(META_INFO)->second.push_back(line_num);
}

void quakelib::EQSimMetadataReader::parse_title_record(const int &line_num, std::istringstream &line_stream) {
	std::string		title_rec;
	
	getline(line_stream, title_rec);
	title_rec.erase(0,1);		// erase the extra space getline adds to the beginning
	metadata_recs.find(META_TITLE)->second.clear();
	metadata_line_nums.find(META_TITLE)->second.clear();
	metadata_recs.find(META_TITLE)->second.push_back(title_rec);
	metadata_line_nums.find(META_TITLE)->second.push_back(line_num);
	if (parsed_title) parse_errors.report("Title record previously found.", line_num);
	parsed_title = true;
}

void quakelib::EQSimMetadataReader::parse_author_record(const int &line_num, std::istringstream &line_stream) {
	std::string	new_author;
	getline(line_stream, new_author);
	new_author.erase(0,1);		// erase the extra space getline adds to the beginning
	metadata_recs.find(META_AUTHOR)->second.push_back(new_author);
	metadata_line_nums.find(META_AUTHOR)->second.push_back(line_num);
}

void quakelib::EQSimMetadataReader::parse_date_record(const int &line_num, std::istringstream &line_stream) {
	std::string		date_rec;
	
	getline(line_stream, date_rec);
	date_rec.erase(0,1);		// erase the extra space getline adds to the beginning
	metadata_recs.find(META_DATE)->second.clear();
	metadata_line_nums.find(META_DATE)->second.clear();
	metadata_recs.find(META_DATE)->second.push_back(date_rec);
	metadata_line_nums.find(META_DATE)->second.push_back(line_num);
	
	if (parsed_date) parse_errors.report("Date record previously found.", line_num);
	parsed_date = true;
}

void quakelib::EQSimMetadataReader::parse_record_desc_record(const int &line_num, std::istringstream &line_stream) {
	int				rec_kind;
	std::string		rec_name;
	int				rec_n_field;
	
	line_stream >> rec_kind
				>> rec_name
				>> rec_n_field;
	
	if (record_descs.count(rec_kind)) {
		parse_errors.report("Record descriptor with this ID already exists.", line_num);
	} else if (rec_kind < 200 || rec_kind > 399) {
		parse_errors.report("Record kind must be between 200 and 399.", line_num);
	} else if (rec_n_field < 1 || rec_n_field > 100) {
		parse_errors.report("Record descriptor must have between 1 and 100 fields.", line_num);
	} else {
		record_descs.insert(std::make_pair(rec_kind, internal::RecordDesc(line_num, rec_name, rec_n_field)));
		cur_record_index = rec_kind;
	}
}

void quakelib::EQSimMetadataReader::parse_field_desc_record(const int &line_num, std::istringstream &line_stream) {
	int				field_index;
	std::string		field_name;
	int				field_type;
	
	line_stream >> field_index
				>> field_name
				>> field_type;
	
	if (!record_descs.count(cur_record_index)) {
		parse_errors.report("Field descriptor may not be declared without corresponding record descriptor.", line_num);
	} else if (field_index < 1 || field_index > 100) {
		parse_errors.report("Field index must be between 1 and 100.", line_num);
	} else if (field_index > record_descs[cur_record_index].get_num_fields()) {
		parse_errors.report("Field index must be less than record maximum index.", line_num);
	} else if (field_type < 1 || field_type > 3) {
		parse_errors.report("Field type must be between 1 and 3.", line_num);
	} else {
		bool		added_field;
		added_field = record_descs[cur_record_index].add_field(field_index, internal::FieldDesc(field_name, field_type));
		if (!added_field) {
			parse_errors.report("Field index must be unique within record.", line_num);
		}
	}
}

void quakelib::EQSimMetadataReader::parse_eof(const int &line_num, std::istringstream &line_stream) {
	std::string str;
	line_stream >> str;
	if (str.compare("End")) parse_errors.report("End record has incorrect format.", line_num);
	finish_parsing = true;
}

void quakelib::EQSimMetadataReader::validate(EQSimErrors &errors) const {
	internal::MultiRecordDesc::const_iterator		it;
	std::stringstream					error_msg;
	
	internal::EQSimMetadata::validate(errors);
	
	for (it=record_descs.begin();it!=record_descs.end();++it) {
		if (it->second.get_num_fields() != it->second.get_num_fields()) {
			error_msg << "Number of field records (" << it->second.get_num_fields()
					<< ") does not match specified number for record (" << it->second.get_num_fields() << ").";
			errors.report(error_msg.str(), it->second.get_max_index());
		}
	}
}

void quakelib::EQSimMetadataWriter::validate(EQSimErrors &errors) const {
	internal::EQSimMetadata::validate(errors);
}

void quakelib::internal::EQSimMetadata::validate(EQSimErrors &errors) const {
	unsigned int				i;
	std::stringstream			error_msg;
	std::vector<RECORD_TYPE>	rec_types;
	std::vector<RECORD_TYPE>::const_iterator	it;
	
	rec_types.push_back(META_COMMENT);
	rec_types.push_back(META_SIGNATURE);
	rec_types.push_back(META_INFO);
	rec_types.push_back(META_TITLE);
	rec_types.push_back(META_AUTHOR);
	rec_types.push_back(META_DATE);
	
	for (it=rec_types.begin();it!=rec_types.end();++it) {
		for (i=0;i<metadata_recs.find(*it)->second.size();++i) {
			if (metadata_recs.find(*it)->second.at(i).empty()) {
				error_msg << "Metadata records may not be blank (record "
					<< metadata_rec_nums.find(*it)->second << " number "
					<< i << ").";
				errors.report(error_msg.str(), metadata_line_nums.find(*it)->second.at(i));
			}
		}
	}
}

void quakelib::EQSimMetadataWriter::write(void) {
	internal::MultiRecordDesc::const_iterator		it;
	internal::MultiRecord::const_iterator			rit;
	
	if (wrote_header) return;
	
	for (rit=metadata_recs.find(META_SIGNATURE)->second.begin();rit!=metadata_recs.find(META_SIGNATURE)->second.end();++rit)
		out_stream << metadata_rec_nums.find(META_SIGNATURE)->second << " " << *rit << " " << spec_level << "\n";
	
	for (rit=metadata_recs.find(META_COMMENT)->second.begin();rit!=metadata_recs.find(META_COMMENT)->second.end();++rit)
		out_stream << metadata_rec_nums.find(META_COMMENT)->second << " " << *rit << "\n";
	for (rit=metadata_recs.find(META_INFO)->second.begin();rit!=metadata_recs.find(META_INFO)->second.end();++rit)
		out_stream << metadata_rec_nums.find(META_INFO)->second << " " << *rit << "\n";
	for (rit=metadata_recs.find(META_TITLE)->second.begin();rit!=metadata_recs.find(META_TITLE)->second.end();++rit)
		out_stream << metadata_rec_nums.find(META_TITLE)->second << " " << *rit << "\n";
	for (rit=metadata_recs.find(META_AUTHOR)->second.begin();rit!=metadata_recs.find(META_AUTHOR)->second.end();++rit)
		out_stream << metadata_rec_nums.find(META_AUTHOR)->second << " " << *rit << "\n";
	for (rit=metadata_recs.find(META_DATE)->second.begin();rit!=metadata_recs.find(META_DATE)->second.end();++rit)
		out_stream << metadata_rec_nums.find(META_DATE)->second << " " << *rit << "\n";
	
	out_stream << "102 End_Metadata\n";
	
	for (it=record_descs.begin();it!=record_descs.end();++it) it->second.write(it->first, out_stream);
	
	out_stream << "103 End_Descriptor\n";
	out_stream.flush();
	
	wrote_header = true;
}

void quakelib::EQSimMetadataWriter::close(void) {
	out_stream << "999 End\n";
	out_stream.flush();
}

void quakelib::internal::FieldDesc::write(const unsigned int &key, std::ostream &out_stream) const {
	out_stream << "121 " << key << " " << _name << " " << _type << "    " << "Field " << _field_num << ": " << _desc << "\n";
}

void quakelib::internal::RecordDesc::write(const unsigned int &key, std::ostream &out_stream) const {
	MultiFieldDesc::const_iterator		it;
	out_stream << "120 " << key << " " << name << " " << field_descs.size() << "    " << desc << "\n";
	for (it=field_descs.begin();it!=field_descs.end();++it) it->second.write(it->first, out_stream);
}

bool quakelib::internal::RecordDesc::add_field(const int &_index, const FieldDesc &_new_desc) {
	if (field_descs.count(_index)) return false;
	else field_descs[_index] = _new_desc;
	return true;
}

void quakelib::EQSimErrors::report(const std::string &error_msg, const int &line_num=-1) {
	error_list.insert(std::make_pair(line_num, error_msg));
}

void quakelib::EQSimErrors::write(std::ostream &os) const {
	std::multimap<int, std::string>::const_iterator		it;
	
	for (it=error_list.begin();it!=error_list.end();++it) {
		if (it->first >= 0) os << "LINE " << it->first << ": ";
		os << it->second << std::endl;
	}
}
