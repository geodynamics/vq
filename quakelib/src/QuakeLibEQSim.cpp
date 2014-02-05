// Copyright (c) 2012-2013 Eric M. Heien
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

#include "QuakeLibEQSim.h"
#include <iomanip>

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
            case 100:
                parse_comment(line_num, line_stream);
                return true;

            case 102:
                parse_end_metadata(line_num, line_stream);
                return true;

            case 110:
                parse_info_record(line_num, line_stream);
                return true;

            case 111:
                parse_title_record(line_num, line_stream);
                return true;

            case 112:
                parse_author_record(line_num, line_stream);
                return true;

            case 113:
                parse_date_record(line_num, line_stream);
                return true;

            case 999:
                parse_eof(line_num, line_stream);
                return true;

            default:
                return false;
        }
    }

    // Otherwise, assuming we're not done parsing descriptors, we read them in
    if (!parsed_descriptors) {
        switch (rec_num) {
            case 100:
                parse_comment(line_num, line_stream);
                return true;

            case 103:
                parse_end_descriptor(line_num, line_stream);
                return true;

            case 120:
                parse_record_desc_record(line_num, line_stream);
                return true;

            case 121:
                parse_field_desc_record(line_num, line_stream);
                return true;

            case 999:
                parse_eof(line_num, line_stream);
                return true;

            default:
                return false;
        }
    }

    // Once we're past the descriptors, we can only accept EOF and comments
    switch (rec_num) {
        case 100:
            parse_comment(line_num, line_stream);
            return true;

        case 999:
            parse_eof(line_num, line_stream);
            return true;

        default:
            return false;
    }
}

bool quakelib::EQSimMetadataReader::parse_file(const std::string &file_name, const int &start_line_num) {
    std::string         next_line;
    std::istringstream  line_stream;
    int                 record_num, line_num;
    std::ifstream       input_file;
    bool                found_match;

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
                std::stringstream       error_msg;
                error_msg << "Unexpected record number: " << record_num;
                parse_errors.report(error_msg.str(), line_num);
            }
        } catch (std::exception e) {
            std::stringstream       error_msg;
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
    if (rec_key<200 || rec_key>=400) throw std::invalid_argument("quakelib::EQSimMetadata::add_record_desc_record");

    record_descs.insert(std::make_pair(rec_key, rec_desc));
}

void quakelib::EQSimMetadataReader::parse_comment(const int &line_num, std::istringstream &line_stream) {
    std::string     comment;
    getline(line_stream, comment);
    comment.erase(0,1);     // erase the extra space getline adds to the beginning
    metadata_recs.find(META_COMMENT)->second.push_back(comment);
    metadata_line_nums.find(META_COMMENT)->second.push_back(line_num);
}

void quakelib::EQSimMetadataReader::parse_signature(const int &line_num, std::istringstream &line_stream) {
    std::string     sig_rec;

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
    new_info.erase(0,1);        // erase the extra space getline adds to the beginning
    metadata_recs.find(META_INFO)->second.push_back(new_info);
    metadata_line_nums.find(META_INFO)->second.push_back(line_num);
}

void quakelib::EQSimMetadataReader::parse_title_record(const int &line_num, std::istringstream &line_stream) {
    std::string     title_rec;

    getline(line_stream, title_rec);
    title_rec.erase(0,1);       // erase the extra space getline adds to the beginning
    metadata_recs.find(META_TITLE)->second.clear();
    metadata_line_nums.find(META_TITLE)->second.clear();
    metadata_recs.find(META_TITLE)->second.push_back(title_rec);
    metadata_line_nums.find(META_TITLE)->second.push_back(line_num);

    if (parsed_title) parse_errors.report("Title record previously found.", line_num);

    parsed_title = true;
}

void quakelib::EQSimMetadataReader::parse_author_record(const int &line_num, std::istringstream &line_stream) {
    std::string new_author;
    getline(line_stream, new_author);
    new_author.erase(0,1);      // erase the extra space getline adds to the beginning
    metadata_recs.find(META_AUTHOR)->second.push_back(new_author);
    metadata_line_nums.find(META_AUTHOR)->second.push_back(line_num);
}

void quakelib::EQSimMetadataReader::parse_date_record(const int &line_num, std::istringstream &line_stream) {
    std::string     date_rec;

    getline(line_stream, date_rec);
    date_rec.erase(0,1);        // erase the extra space getline adds to the beginning
    metadata_recs.find(META_DATE)->second.clear();
    metadata_line_nums.find(META_DATE)->second.clear();
    metadata_recs.find(META_DATE)->second.push_back(date_rec);
    metadata_line_nums.find(META_DATE)->second.push_back(line_num);

    if (parsed_date) parse_errors.report("Date record previously found.", line_num);

    parsed_date = true;
}

void quakelib::EQSimMetadataReader::parse_record_desc_record(const int &line_num, std::istringstream &line_stream) {
    int             rec_kind;
    std::string     rec_name;
    int             rec_n_field;

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
    int             field_index;
    std::string     field_name;
    int             field_type;

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
        bool        added_field;
        added_field = record_descs[cur_record_index].add_field(field_index, internal::EQSimFieldDesc(field_name, field_type));

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
    internal::MultiRecordDesc::const_iterator       it;
    std::stringstream                   error_msg;

    internal::EQSimMetadata::validate(errors);

    for (it=record_descs.begin(); it!=record_descs.end(); ++it) {
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
    unsigned int                i;
    std::stringstream           error_msg;
    std::vector<RECORD_TYPE>    rec_types;
    std::vector<RECORD_TYPE>::const_iterator    it;

    rec_types.push_back(META_COMMENT);
    rec_types.push_back(META_SIGNATURE);
    rec_types.push_back(META_INFO);
    rec_types.push_back(META_TITLE);
    rec_types.push_back(META_AUTHOR);
    rec_types.push_back(META_DATE);

    for (it=rec_types.begin(); it!=rec_types.end(); ++it) {
        for (i=0; i<metadata_recs.find(*it)->second.size(); ++i) {
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
    internal::MultiRecordDesc::const_iterator       it;
    internal::MultiRecord::const_iterator           rit;

    if (wrote_header) return;

    for (rit=metadata_recs.find(META_SIGNATURE)->second.begin(); rit!=metadata_recs.find(META_SIGNATURE)->second.end(); ++rit)
        out_stream << metadata_rec_nums.find(META_SIGNATURE)->second << " " << *rit << " " << spec_level << "\n";

    for (rit=metadata_recs.find(META_COMMENT)->second.begin(); rit!=metadata_recs.find(META_COMMENT)->second.end(); ++rit)
        out_stream << metadata_rec_nums.find(META_COMMENT)->second << " " << *rit << "\n";

    for (rit=metadata_recs.find(META_INFO)->second.begin(); rit!=metadata_recs.find(META_INFO)->second.end(); ++rit)
        out_stream << metadata_rec_nums.find(META_INFO)->second << " " << *rit << "\n";

    for (rit=metadata_recs.find(META_TITLE)->second.begin(); rit!=metadata_recs.find(META_TITLE)->second.end(); ++rit)
        out_stream << metadata_rec_nums.find(META_TITLE)->second << " " << *rit << "\n";

    for (rit=metadata_recs.find(META_AUTHOR)->second.begin(); rit!=metadata_recs.find(META_AUTHOR)->second.end(); ++rit)
        out_stream << metadata_rec_nums.find(META_AUTHOR)->second << " " << *rit << "\n";

    for (rit=metadata_recs.find(META_DATE)->second.begin(); rit!=metadata_recs.find(META_DATE)->second.end(); ++rit)
        out_stream << metadata_rec_nums.find(META_DATE)->second << " " << *rit << "\n";

    out_stream << "102 End_Metadata\n";

    for (it=record_descs.begin(); it!=record_descs.end(); ++it) it->second.write(it->first, out_stream);

    out_stream << "103 End_Descriptor\n";
    out_stream.flush();

    wrote_header = true;
}

void quakelib::EQSimMetadataWriter::close(void) {
    out_stream << "999 End\n";
    out_stream.flush();
}

void quakelib::internal::EQSimFieldDesc::write(const unsigned int &key, std::ostream &out_stream) const {
    out_stream << "121 " << key << " " << _name << " " << _type << "    " << "Field " << _field_num << ": " << _desc << "\n";
}

void quakelib::internal::RecordDesc::write(const unsigned int &key, std::ostream &out_stream) const {
    MultiFieldDesc::const_iterator      it;
    out_stream << "120 " << key << " " << name << " " << field_descs.size() << "    " << desc << "\n";

    for (it=field_descs.begin(); it!=field_descs.end(); ++it) it->second.write(it->first, out_stream);
}

bool quakelib::internal::RecordDesc::add_field(const int &_index, const EQSimFieldDesc &_new_desc) {
    if (field_descs.count(_index)) return false;
    else field_descs[_index] = _new_desc;

    return true;
}

void quakelib::EQSimErrors::report(const std::string &error_msg, const int &line_num=-1) {
    error_list.insert(std::make_pair(line_num, error_msg));
}

void quakelib::EQSimErrors::write(std::ostream &os) const {
    std::multimap<int, std::string>::const_iterator     it;

    for (it=error_list.begin(); it!=error_list.end(); ++it) {
        if (it->first >= 0) os << "LINE " << it->first << ": ";

        os << it->second << std::endl;
    }
}

quakelib::EQSimConditionWriter::EQSimConditionWriter(void) {
    internal::RecordDesc    summary_rec, init_stress_rec, init_state_rec;

    meta_add_record(META_SIGNATURE, "EQSim_Input_Friction_2");
    set_spec_level(1);

    summary_rec = internal::RecordDesc(0, "summary", 18, "Record 200: Initial condition summary");
    summary_rec.add_field(1, internal::EQSimFieldDesc("n_element", 1, 1, "Total number of elements in the file"));
    summary_rec.add_field(2, internal::EQSimFieldDesc("stress_flag", 1, 2, "1 if initial stress (record 201) is included, 0 if not"));
    summary_rec.add_field(3, internal::EQSimFieldDesc("state_flag", 1, 3, "1 if initial state (record 202) is included, 0 if not"));

    init_stress_rec = internal::RecordDesc(0, "initial_stress", 3, "Record 201: Slip map entry");
    init_stress_rec.add_field(1, internal::EQSimFieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
    init_stress_rec.add_field(2, internal::EQSimFieldDesc("shear_stress", 2, 2, "Element initial shear stress (Pascal)"));
    init_stress_rec.add_field(3, internal::EQSimFieldDesc("normal_stress", 2, 3, "Element initial normal stress (Pascal)"));

    init_state_rec = internal::RecordDesc(0, "initial_state", 2, "Record 202: Initial state");
    init_state_rec.add_field(1, internal::EQSimFieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
    init_state_rec.add_field(2, internal::EQSimFieldDesc("rs_theta", 2, 2, "Element initial state for rate-state law (seconds)"));

    add_record_desc_record(200, summary_rec);
    add_record_desc_record(201, init_stress_rec);
    add_record_desc_record(202, init_state_rec);
}

bool quakelib::EQSimConditionReader::parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream) {
    switch (rec_num) {
        case 200:
            parse_summary_record(line_num, line_stream);
            return true;

        case 201:
            parse_stress_record(line_num, line_stream);
            return true;

        case 202:
            parse_state_record(line_num, line_stream);
            return true;

        default:
            return false;
    }
}

void quakelib::EQSimConditionReader::parse_summary_record(const int &line_num, std::istringstream &line_stream) {
    _line_num = line_num;
    line_stream >> _n_element       // Field 1: Total number of elements in the file
                >> _stress_flag     // Field 2: 1 if initial stress (record 201) is included, 0 if not
                >> _state_flag;     // Field 3: 1 if initial state (record 202) is included, 0 if not
}

void quakelib::EQSimConditionReader::parse_stress_record(const int &line_num, std::istringstream &line_stream) {
    UIndex      ind;
    double      shear_str, normal_str;

    line_stream >> ind              // Field 1: Element index number (consecutive integers, starting with 1)
                >> shear_str        // Field 2: Element initial shear stress (Pascal)
                >> normal_str;      // Field 3: Element initial normal stress (Pascal)

    _stresses.insert(std::make_pair(ind, std::make_pair(shear_str, normal_str)));
}

void quakelib::EQSimConditionReader::parse_state_record(const int &line_num, std::istringstream &line_stream) {
    UIndex      ind;
    double      rs_t;

    line_stream >> ind              // Field 1: Element index number (consecutive integers, starting with 1)
                >> rs_t;            // Field 2: Element initial state for rate-state law (seconds)

    _rs_theta.insert(std::make_pair(ind, rs_t));
}

// Get the number of elements in this condition object.
size_t quakelib::EQSimCondition::num_elements(void) const {
    if (_stresses.size() > 0) return _stresses.size();
    else return _rs_theta.size();
}

/*!
 Get the shear and normal stresses of a given element.
 */
double quakelib::EQSimCondition::get_shear_stress(const UIndex &ind) const throw(std::out_of_range) {
    StressMap::const_iterator it = _stresses.find(ind);

    if (it == _stresses.end()) throw std::out_of_range("quakelib::EQSimCondition::get_shear_stress");

    return it->second.first;
}

double quakelib::EQSimCondition::get_normal_stress(const UIndex &ind) const throw(std::out_of_range) {
    StressMap::const_iterator it = _stresses.find(ind);

    if (it == _stresses.end()) throw std::out_of_range("quakelib::EQSimCondition::get_normal_stress");

    return it->second.second;
}

/*!
 Set the shear and normal stresses of a given element.
 */
void quakelib::EQSimCondition::set_stresses(const UIndex &ind, const double &shear_stress, const double &normal_stress) {
    _stress_lines.insert(std::make_pair(ind, -1));
    _stresses.insert(std::make_pair(ind, std::make_pair(shear_stress, normal_stress)));
}

/*
 Get the rate-state value of a given element.
 */
double quakelib::EQSimCondition::get_rate_state(const UIndex &ind) const throw(std::out_of_range) {
    StateMap::const_iterator it = _rs_theta.find(ind);

    if (it == _rs_theta.end()) throw std::out_of_range("quakelib::EQSimCondition::get_rate_state");

    return it->second;
}

/*!
 Set the rate-state value of a given element.
 */
void quakelib::EQSimCondition::set_rate_state(const UIndex &ind, const double &rate_state) {
    _rst_lines.insert(std::make_pair(ind, -1));
    _rs_theta.insert(std::make_pair(ind, rate_state));
}

// TODO: interleave elements according to specification
void quakelib::EQSimCondition::write(std::ostream &out_stream) const {
    StressMap::const_iterator   sit;
    StateMap::const_iterator    rit;

    if (_stresses.size() > 0 && _rs_theta.size() > 0 && _stresses.size() != _rs_theta.size())
        throw std::length_error("quakelib::EQSimCondition::write");

    // Write the header record
    out_stream  << "200"
                << " " << _stresses.size()
                << " " << (_stresses.size() > 0)
                << " " << (_rs_theta.size() > 0)
                << std::endl;

    // Write stress information
    for (sit=_stresses.begin(); sit!=_stresses.end(); ++sit) {
        out_stream  << "201"
                    << " " << sit->first
                    << " " << sit->second.first
                    << " " << sit->second.second
                    << std::endl;
    }

    // Write rate-state information
    for (rit=_rs_theta.begin(); rit!=_rs_theta.end(); ++rit) {
        out_stream  << "202"
                    << " " << rit->first
                    << " " << rit->second
                    << std::endl;
    }
}

void quakelib::EQSimCondition::validate(EQSimErrors &errors) const {
    StressMap::const_iterator   sit;

    // Check that no stresses are negative
    for (sit=_stresses.begin(); sit!=_stresses.end(); ++sit) {
        if (sit->second.first < 0 || sit->second.second < 0) {
            std::stringstream       error_msg;
            error_msg << "Element " << sit->first << " has negative initial stress (shear="
                      << sit->second.first << ", normal=" << sit->second.second << ").";
            errors.report(error_msg.str(), _stress_lines.find(sit->first)->second);
        }
    }
}

/*!
 Confirms that the condition file represented by this class is self consistent.
 To do so, it ensures that the number of shear and normal stress elements
 actually stored is the same as the number specified in the file. It also
 does the same check for the number of rate-state entries.
 */
void quakelib::EQSimConditionReader::validate(EQSimErrors &errors) const {
    EQSimMetadataReader::validate(errors);
    EQSimCondition::validate(errors);

    if (_stress_flag) {
        // If stress flag was true, make sure the number of entries was as specified
        if (_stresses.size() != _n_element) {
            std::stringstream       error_msg;
            error_msg << "Specified number of stress entries (" << _n_element
                      << ") does not match actual number in file (" << _stresses.size() << ").";
            errors.report(error_msg.str(), _line_num);
        }
    } else {
        // Check the shear and normal stress entries
        if (_stresses.size() > 0) {
            std::stringstream       error_msg;
            error_msg << "Stress flag indicates no stress entries, but found "
                      << _stresses.size() << " in the file.";
            errors.report(error_msg.str(), _line_num);
        }
    }

    if (_state_flag) {
        // If state flag was true, make sure the number of entries was as specified
        if (_rs_theta.size() != _n_element) {
            std::stringstream       error_msg;
            error_msg << "Specified number of rate-state entries (" << _n_element
                      << ") does not match actual number in file (" << _rs_theta.size() << ").";
            errors.report(error_msg.str(), _line_num);
        }
    } else {
        // Check the rate-state entries
        if (_rs_theta.size() > 0) {
            std::stringstream       error_msg;
            error_msg << "State flag indicates no rate-state entries, but found "
                      << _rs_theta.size() << " in the file.";
            errors.report(error_msg.str(), _line_num);
        }
    }
}

void quakelib::EQSimEventWriter::flush(void) {
    std::vector<unsigned int>::const_iterator   it;

    quakelib::EQSimMetadataWriter::write();

    for (it=entries_to_flush.begin(); it!=entries_to_flush.end(); ++it) {
        event_summaries.at(*it).write(out_stream);
    }

    entries_to_flush.clear();
    out_stream.flush();
}

void quakelib::EQSimEventSummary::parse(const int &line_num, std::istringstream &line_stream) {
    _line_num = line_num;

    line_stream >> _event_id    // Field 1: Event ID number (positive integers, in order, need not be consecutive)
                >> _magnitude       // Field 2: Event magnitude
                >> _time            // Field 3: Starting time of event (seconds)
                >> _duration        // Field 4: Duration of event (seconds)
                >> _sid             // Field 5: Fault section ID number (positive integer)
                >> _depth_lo        // Field 6: Lowest value of depth in the rupture (meters, negative underground)
                >> _depth_hi        // Field 7: Highest value of depth in the rupture (meters, negative underground)
                >> _das_lo          // Field 8: Lowest value of distance-along-strike in the rupture (meters)
                >> _das_hi          // Field 9: Highest value of distance-along-strike in the rupture (meters)
                >> _hypo_depth      // Field 10: Hypocenter depth (meters, negative underground)
                >> _hypo_das        // Field 11: Hypocenter distance-along-strike (meters)
                >> _area            // Field 12: Rupture area (square meters)
                >> _mean_slip       // Field 13: Average slip (meters)
                >> _moment          // Field 14: Seismic moment (Newton-meters)
                >> _shear_before    // Field 15: Shear stress before event (Pascal)
                >> _shear_after     // Field 16: Shear stress after event (Pascal)
                >> _normal_before   // Field 17: Normal stress before event (Pascal)
                >> _normal_after;   // Field 18: Normal stress after event (Pascal)
}

void quakelib::EQSimEventReader::parse_event_summary(const int &line_num, std::istringstream &line_stream) {
    quakelib::EQSimEventSummary es;

    es.parse(line_num, line_stream);

    event_summaries.push_back(es);
}

void quakelib::EQSimEventSlipMap::parse(const int &line_num, std::istringstream &line_stream) {
    _line_num = line_num;

    line_stream >> _depth_lo        // Field 1: Lowest value of depth (meters, negative underground)
                >> _depth_hi        // Field 2: Highest value of depth in the rupture (meters, negative underground)
                >> _das_lo          // Field 3: Lowest value of distance-along-strike in the rupture (meters)
                >> _das_hi          // Field 4: Highest value of distance-along-strike in the rupture (meters)
                >> _area            // Field 5: Rupture area (square meters)
                >> _mean_slip       // Field 6: Average slip (meters)
                >> _moment          // Field 7: Seismic moment (Newton-meters)
                >> _shear_before    // Field 8: Shear stress before event (Pascal)
                >> _shear_after     // Field 9: Shear stress after event (Pascal)
                >> _normal_before   // Field 10: Normal stress before event (Pascal)
                >> _normal_after    // Field 11: Normal stress after event (Pascal)
                >> _element_id;     // Field 12: Element ID number (positive integer), or negative of element count (zero if no element info)
}

void quakelib::EQSimEventReader::parse_event_slip_map(const int &line_num, std::istringstream &line_stream) {
    quakelib::EQSimEventSlipMap sm;

    if (event_summaries.size() == 0) throw std::exception();

    sm.parse(line_num, line_stream);

    event_summaries.back().add_slip_map(sm);
}

void quakelib::EQSimEventReader::parse_event_slip_element(const int &line_num, std::istringstream &line_stream) {
    unsigned int        elem_id;

    if (event_summaries.size() == 0 || event_summaries.back().slip_maps.size() == 0) throw std::exception();

    line_stream >> elem_id;

    event_summaries.back().slip_maps.back().add_slip_entry(EQSimEventSlipElement(elem_id, line_num));
}

void quakelib::EQSimEventSummary::write(std::ostream &out_stream) const {
    std::vector<EQSimEventSlipMap>::const_iterator      it;
    std::streamsize                                     old_prec;

    old_prec = out_stream.precision();
    out_stream  << "200"
                << " " << _event_id         // Field 1: Event ID number (positive integers, in order, need not be consecutive)
                << " " << _magnitude            // Field 2: Event magnitude
                << " " << std::setprecision(16)
                << _time                        // Field 3: Starting time of event (seconds)
                << " " << std::setprecision(old_prec)
                << _duration                    // Field 4: Duration of event (seconds)
                << " " << _sid              // Field 5: Fault section ID number (positive integer)
                << " " << _depth_lo         // Field 6: Lowest value of depth in the rupture (meters, negative underground)
                << " " << _depth_hi         // Field 7: Highest value of depth in the rupture (meters, negative underground)
                << " " << _das_lo           // Field 8: Lowest value of distance-along-strike in the rupture (meters)
                << " " << _das_hi           // Field 9: Highest value of distance-along-strike in the rupture (meters)
                << " " << _hypo_depth       // Field 10: Hypocenter depth (meters, negative underground)
                << " " << _hypo_das         // Field 11: Hypocenter distance-along-strike (meters)
                << " " << _area             // Field 12: Rupture area (square meters)
                << " " << _mean_slip            // Field 13: Average slip (meters)
                << " " << _moment           // Field 14: Seismic moment (Newton-meters)
                << " " << _shear_before     // Field 15: Shear stress before event (Pascal)
                << " " << _shear_after      // Field 16: Shear stress after event (Pascal)
                << " " << _normal_before        // Field 17: Normal stress before event (Pascal)
                << " " << _normal_after     // Field 18: Normal stress after event (Pascal)
                << "\n";

    for (it=slip_maps.begin(); it!=slip_maps.end(); ++it) {
        it->write(out_stream);
    }
}

void quakelib::EQSimEventSlipMap::write(std::ostream &out_stream) const {
    std::vector<EQSimEventSlipElement>::const_iterator      it;

    out_stream  << "201"
                << " " << _depth_lo         // Field 1: Lowest value of depth (meters, negative underground)
                << " " << _depth_hi         // Field 2: Highest value of depth in the rupture (meters, negative underground)
                << " " << _das_lo           // Field 3: Lowest value of distance-along-strike in the rupture (meters)
                << " " << _das_hi           // Field 4: Highest value of distance-along-strike in the rupture (meters)
                << " " << _area             // Field 5: Rupture area (square meters)
                << " " << _mean_slip        // Field 6: Average slip (meters)
                << " " << _moment           // Field 7: Seismic moment (Newton-meters)
                << " " << _shear_before     // Field 8: Shear stress before event (Pascal)
                << " " << _shear_after      // Field 9: Shear stress after event (Pascal)
                << " " << _normal_before    // Field 10: Normal stress before event (Pascal)
                << " " << _normal_after;    // Field 11: Normal stress after event (Pascal)

    // Field 12: Element ID number (positive integer), or negative of element count (zero if no element info)
    if (slip_elements.size() == 0) {
        out_stream << " " << 0 << "\n";
    } else if (slip_elements.size() == 1) {
        out_stream << " " << slip_elements.back().element_id() << "\n";
    } else {
        out_stream << " -" << slip_elements.size() << "\n";

        for (it=slip_elements.begin(); it!=slip_elements.end(); ++it) it->write(out_stream);
    }
}

void quakelib::EQSimEventSlipElement::write(std::ostream &out_stream) const {
    out_stream  << "202"
                << " " << _element_id
                << "\n";
}

bool quakelib::EQSimEventReader::parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream) {
    switch (rec_num) {
        case 200:
            parse_event_summary(line_num, line_stream);
            return true;

        case 201:
            parse_event_slip_map(line_num, line_stream);
            return true;

        case 202:
            parse_event_slip_element(line_num, line_stream);
            return true;

        default:
            return false;
    }
}

quakelib::EQSimEventWriter::EQSimEventWriter(void) {
    quakelib::internal::RecordDesc  event_summary_rec, slip_map_rec, slip_elem_rec;

    set_spec_level(2);
    meta_add_record(META_SIGNATURE, "EQSim_Output_Event_2");

    event_summary_rec = internal::RecordDesc(0, "event_summary", 18, "Record 200: Event summary");
    event_summary_rec.add_field(1, internal::EQSimFieldDesc("event_id", 1, 1, "Event ID number (positive integers, in order, need not be consecutive)"));
    event_summary_rec.add_field(2, internal::EQSimFieldDesc("magnitude", 2, 2, "Event magnitude"));
    event_summary_rec.add_field(3, internal::EQSimFieldDesc("time", 2, 3, "Starting time of event (seconds)"));
    event_summary_rec.add_field(4, internal::EQSimFieldDesc("duration", 2, 4, "Duration of event (seconds)"));
    event_summary_rec.add_field(5, internal::EQSimFieldDesc("sid", 1, 5, "Fault section ID number (positive integer)"));
    event_summary_rec.add_field(6, internal::EQSimFieldDesc("depth_lo", 2, 6, "Lowest value of depth in the rupture (meters, negative underground)"));
    event_summary_rec.add_field(7, internal::EQSimFieldDesc("depth_hi", 2, 7, "Highest value of depth in the rupture (meters, negative underground)"));
    event_summary_rec.add_field(8, internal::EQSimFieldDesc("das_lo", 2, 8, "Lowest value of distance-along-strike in the rupture (meters)"));
    event_summary_rec.add_field(9, internal::EQSimFieldDesc("das_hi", 2, 9, "Highest value of distance-along-strike in the rupture (meters)"));
    event_summary_rec.add_field(10, internal::EQSimFieldDesc("hypo_depth", 2, 10, "Hypocenter depth (meters, negative underground)"));
    event_summary_rec.add_field(11, internal::EQSimFieldDesc("hypo_das", 2, 11, "Hypocenter distance-along-strike (meters)"));
    event_summary_rec.add_field(12, internal::EQSimFieldDesc("area", 2, 12, "Rupture area (square meters)"));
    event_summary_rec.add_field(13, internal::EQSimFieldDesc("mean_slip", 2, 13, "Average slip (meters)"));
    event_summary_rec.add_field(14, internal::EQSimFieldDesc("moment", 2, 14, "Seismic moment (Newton-meters)"));
    event_summary_rec.add_field(15, internal::EQSimFieldDesc("shear_before", 2, 15, "Shear stress before event (Pascal)"));
    event_summary_rec.add_field(16, internal::EQSimFieldDesc("shear_after", 2, 16, "Shear stress after event (Pascal)"));
    event_summary_rec.add_field(17, internal::EQSimFieldDesc("normal_before", 2, 17, "Normal stress before event (Pascal)"));
    event_summary_rec.add_field(18, internal::EQSimFieldDesc("normal_after", 2, 18, "Normal stress after event (Pascal)"));

    slip_map_rec = internal::RecordDesc(0, "slip_map", 12, "Record 201: Slip map entry");
    slip_map_rec.add_field(1, internal::EQSimFieldDesc("depth_lo", 2, 1, "Lowest value of depth (meters, negative underground)"));
    slip_map_rec.add_field(2, internal::EQSimFieldDesc("depth_hi", 2, 2, "Highest value of depth (meters, negative underground)"));
    slip_map_rec.add_field(3, internal::EQSimFieldDesc("das_lo", 2, 3, "Lowest value of distance-along-strike (meters)"));
    slip_map_rec.add_field(4, internal::EQSimFieldDesc("das_hi", 2, 4, "Highest value of distance-along-strike (meters)"));
    slip_map_rec.add_field(5, internal::EQSimFieldDesc("area", 2, 5, "Area (square meters)"));
    slip_map_rec.add_field(6, internal::EQSimFieldDesc("mean_slip", 2, 6, "Average slip (meters)"));
    slip_map_rec.add_field(7, internal::EQSimFieldDesc("moment", 2, 7, "Seismic moment (Newton-meters)"));
    slip_map_rec.add_field(8, internal::EQSimFieldDesc("shear_before", 2, 8, "Shear stress before event (Pascal)"));
    slip_map_rec.add_field(9, internal::EQSimFieldDesc("shear_after", 2, 9, "Shear stress after event (Pascal)"));
    slip_map_rec.add_field(10, internal::EQSimFieldDesc("normal_before", 2, 10, "Normal stress before event (Pascal)"));
    slip_map_rec.add_field(11, internal::EQSimFieldDesc("normal_after", 2, 11, "Normal stress after event (Pascal)"));
    slip_map_rec.add_field(12, internal::EQSimFieldDesc("element_id", 1, 12, "Element ID number (positive integer), or negative of element count (zero if no element info)"));

    slip_elem_rec = internal::RecordDesc(0, "slip_element", 1, "Record 202: Slip element list entry");
    slip_elem_rec.add_field(1, internal::EQSimFieldDesc("element_id", 1, 1, "Element ID number (positive integer)"));

    add_record_desc_record(200, event_summary_rec);
    add_record_desc_record(201, slip_map_rec);
    add_record_desc_record(202, slip_elem_rec);
}

quakelib::EQSimFrictionWriter::EQSimFrictionWriter(void) {
    internal::RecordDesc    fric_summary_rec, elastic_params_rec, fault_strength_rec, rate_state_rec;

    set_spec_level(1);
    meta_add_record(META_SIGNATURE, "EQSim_Input_Friction_2");

    fric_summary_rec = internal::RecordDesc(0, "summary", 4, "Record 200: Fault friction summary");
    fric_summary_rec.add_field(1, internal::EQSimFieldDesc("n_element", 1, 1, "Total number of elements in the file"));
    fric_summary_rec.add_field(2, internal::EQSimFieldDesc("elastic_flag", 1, 2, "1 if elastic parameters (record 201) are included, 0 if not"));
    fric_summary_rec.add_field(3, internal::EQSimFieldDesc("strength_flag", 1, 3, "1 if fault strength (record 202) is included, 0 if not"));
    fric_summary_rec.add_field(4, internal::EQSimFieldDesc("rate_state_flag", 1, 4, "1 if rate-state parameters (record 203) are included, 0 if not"));
    elastic_params_rec = internal::RecordDesc(0, "elastic_param", 2, "Record 201: Elastic parameters");
    elastic_params_rec.add_field(1, internal::EQSimFieldDesc("lame_lambda", 2, 1, "Lame modulus lambda (Pascal)"));
    elastic_params_rec.add_field(2, internal::EQSimFieldDesc("lame_mu", 2, 2, "Lame modulus mu, also known as the shear modulus (Pascal)"));
    fault_strength_rec = internal::RecordDesc(0, "fault_strength", 3, "Record 202: Fault strength");
    fault_strength_rec.add_field(1, internal::EQSimFieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
    fault_strength_rec.add_field(2, internal::EQSimFieldDesc("static_strength", 2, 2, "Element static yield strength (Pascal)"));
    fault_strength_rec.add_field(3, internal::EQSimFieldDesc("dynamic_strength", 2, 3, "Element dynamic sliding strength (Pascal)"));
    rate_state_rec = internal::RecordDesc(0, "rate_state", 6, "Record 203: Rate-state parameters");
    rate_state_rec.add_field(1, internal::EQSimFieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
    rate_state_rec.add_field(2, internal::EQSimFieldDesc("A", 2, 2, "Element rate-state parameter A"));
    rate_state_rec.add_field(3, internal::EQSimFieldDesc("B", 2, 3, "Element rate-state parameter B"));
    rate_state_rec.add_field(4, internal::EQSimFieldDesc("L", 2, 4, "Element rate-state characteristic distance L (meters)"));
    rate_state_rec.add_field(5, internal::EQSimFieldDesc("f0", 2, 5, "Element rate-state friction coefficient f0"));
    rate_state_rec.add_field(6, internal::EQSimFieldDesc("V0", 2, 6, "Element rate-state reference velocity V0 (meters/second)"));

    add_record_desc_record(200, fric_summary_rec);
    add_record_desc_record(201, elastic_params_rec);
    add_record_desc_record(202, fault_strength_rec);
    add_record_desc_record(203, rate_state_rec);
}

void quakelib::EQSimFriction::write(std::ostream &out_stream) const {
    StrengthMap::const_iterator     sit;
    RateStateMap::const_iterator    rit;

    if (_strengths.size() > 0 && _rate_state.size() > 0 && _strengths.size() != _rate_state.size())
        throw std::length_error("stresses");

    // Write the summary
    out_stream  << "200"
                << " " << _strengths.size()
                << " " << _have_elastic_params
                << " " << (_strengths.size() > 0)
                << " " << (_rate_state.size() > 0)
                << std::endl;

    // Write the elastic parameters
    if (_have_elastic_params) {
        out_stream  << "201"
                    << " " << _lame_lambda
                    << " " << _lame_mu
                    << std::endl;
    }

    // Write strength information
    for (sit=_strengths.begin(); sit!=_strengths.end(); ++sit) {
        out_stream  << "202"
                    << " " << sit->first
                    << " " << sit->second.first
                    << " " << sit->second.second
                    << std::endl;
    }

    // Write rate-state information
    for (rit=_rate_state.begin(); rit!=_rate_state.end(); ++rit) {
        out_stream<< "203 " << rit->first << " ";
        rit->second.write(out_stream);
        out_stream << std::endl;
    }
}

void quakelib::EQSimFrictionRateState::write(std::ostream &out_stream) const {
    out_stream  << _A
                << " " << _B
                << " " << _L
                << " " << _f0
                << " " << _V0;
}

bool quakelib::EQSimFrictionReader::parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream) {
    switch (rec_num) {
        case 200:
            parse_summary_record(line_num, line_stream);
            return true;

        case 201:
            parse_elastic_params_record(line_num, line_stream);
            return true;

        case 202:
            parse_fault_strength_record(line_num, line_stream);
            return true;

        case 203:
            parse_rate_state_record(line_num, line_stream);
            return true;

        default:
            return false;
    }
}

void quakelib::EQSimFrictionReader::parse_summary_record(const int &line_num, std::istringstream &line_stream) {
    _line_num = line_num;
    line_stream >> _num_element     // Field 1: Total number of elements in the file
                >> _elastic_flag    // Field 2: 1 if elastic parameters (record 201) are included, 0 if not
                >> _strength_flag   // Field 3: 1 if fault strength (record 202) is included, 0 if not
                >> _rate_state_flag;// Field 4: 1 if rate-state parameters (record 203) are included, 0 if not
}

void quakelib::EQSimFrictionReader::parse_elastic_params_record(const int &line_num, std::istringstream &line_stream) {
    _line_num = line_num;
    _have_elastic_params = true;
    line_stream >> _lame_lambda     // Field 1: Lame modulus lambda (Pascal)
                >> _lame_mu;        // Field 2: Lame modulus mu, also known as the shear modulus (Pascal)
}

void quakelib::EQSimFrictionReader::parse_fault_strength_record(const int &line_num, std::istringstream &line_stream) {
    UIndex      ind;
    double      stat_str, dyn_str;

    line_stream >> ind              // Field 1: Element index number (consecutive integers, starting with 1)
                >> stat_str         // Field 2: Element static yield strength (Pascal)
                >> dyn_str;         // Field 3: Element dynamic sliding strength (Pascal)

    _str_lines.insert(std::make_pair(ind, _line_num));
    _strengths.insert(std::make_pair(ind, std::make_pair(stat_str, dyn_str)));
}

void quakelib::EQSimFrictionRateState::parse(const int &line_num, std::istringstream &line_stream) {
    line_stream >> _A           // Field 2: Element rate-state parameter A
                >> _B           // Field 3: Element rate-state parameter B
                >> _L           // Field 4: Element rate-state characteristic distance L (meters)
                >> _f0          // Field 5: Element rate-state friction coefficient f0
                >> _V0;         // Field 6: Element rate-state reference velocity V0 (meters/second)
}

void quakelib::EQSimFrictionReader::parse_rate_state_record(const int &line_num, std::istringstream &line_stream) {
    EQSimFrictionRateState  new_rate_state;
    UIndex                  _index;

    line_stream >> _index;      // Field 1: Element index number (consecutive integers, starting with 1)
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
    StrengthMap::const_iterator     it;

    // Ensure correct parameter values
    if (_lame_mu <= 0) {
        std::stringstream       error_msg;
        error_msg << "Lame mu parameter (" << _lame_mu << ") should be greater than 0.";
        errors.report(error_msg.str(), _line_num);
    }

    // Ensure strengths have valid values
    for (it=_strengths.begin(); it!=_strengths.end(); ++it) {
        if (it->second.first < 0) {
            std::stringstream       error_msg;
            error_msg << "Static yield strength (" << it->second.first << ") should not be below 0.";
            errors.report(error_msg.str(), _str_lines.find(it->first)->second);
        }

        if (it->second.second < 0) {
            std::stringstream       error_msg;
            error_msg << "Dynamic sliding strength (" << it->second.second << ") should not be below 0.";
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
            std::stringstream       error_msg;
            error_msg << "Specified number of strength entries (" << _num_element
                      << ") does not match actual number in file (" << _strengths.size() << ").";
            errors.report(error_msg.str(), _line_num);
        }
    } else {
        // Check the static and dynamic strength entries
        if (_strengths.size() > 0) {
            std::stringstream       error_msg;
            error_msg << "Strength flag indicates no strength entries, but found "
                      << _strengths.size() << " in the file.";
            errors.report(error_msg.str(), _line_num);
        }
    }

    // If rate_state_flag is false, make sure there were no rate_state entries
    if (_rate_state_flag) {
        // If strength flag was true, make sure the number of entries was as specified
        if (_rate_state.size() != _num_element) {
            std::stringstream       error_msg;
            error_msg << "Specified number of rate-state entries (" << _num_element
                      << ") does not match actual number in file (" << _rate_state.size() << ").";
            errors.report(error_msg.str(), _line_num);
        }
    } else {
        // Check the rate-state entries
        if (_rate_state.size() > 0) {
            std::stringstream       error_msg;
            error_msg << "Rate-state flag indicates no rate-state entries, but found "
                      << _rate_state.size() << " in the file.";
            errors.report(error_msg.str(), _line_num);
        }
    }
}

void quakelib::EQSimGeometry::write(std::ostream &out_stream) const {
    quakelib::EQSimGeomSectionMap::const_iterator       sit;

    // TODO: write the metadata
    out_stream  << "200"
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
    for (sit=sections.begin(); sit!=sections.end(); ++sit) {
        sit->second.write(out_stream);
    }
}

void quakelib::EQSimGeometrySection::write(std::ostream &out_stream) const {
    quakelib::EQSimGeomVertexMap::const_iterator            vit;
    quakelib::EQSimGeomTriangleMap::const_iterator      tit;
    quakelib::EQSimGeomRectangleMap::const_iterator     rit;

    out_stream  << "201"
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
    for (vit=vertices.begin(); vit!=vertices.end(); ++vit) {
        vit->second.write(out_stream);
    }

    // Write the triangles
    for (tit=triangles.begin(); tit!=triangles.end(); ++tit) {
        tit->second.write(out_stream);
    }

    // Write the rectangles
    for (rit=rectangles.begin(); rit!=rectangles.end(); ++rit) {
        rit->second.write(out_stream);
    }
}

void quakelib::EQSimGeometryVertex::write(std::ostream &out_stream) const {
    out_stream  << "202"
                << " " << _index
                << " " << _loc.lat()
                << " " << _loc.lon()
                << " " << _loc.altitude()
                << " " << _das
                << " " << _trace_flag
                << std::endl;
}

void quakelib::EQSimGeometryTriangle::write(std::ostream &out_stream) const {
    out_stream  << "203"
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
    out_stream  << "204"
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

quakelib::EQSimGeometryVertex &quakelib::EQSimGeometrySection::new_vertex(void) throw(std::out_of_range) {
    UIndex                                          new_index;
    EQSimGeometryVertex                             new_vertex;
    std::pair<EQSimGeomVertexMap::iterator, bool>   res;

    if (vertices.size() == 0) new_index = 1;
    else new_index = vertices.rbegin()->first + 1;

    new_vertex.set_index(new_index);
    res = vertices.insert(std::make_pair(new_index, new_vertex));

    if (!res.second) throw std::out_of_range("new_vertex");

    return res.first->second;
}

quakelib::EQSimGeometryTriangle &quakelib::EQSimGeometrySection::new_triangle(void) throw(std::out_of_range) {
    UIndex                                          new_index;
    EQSimGeometryTriangle                           new_triangle;
    std::pair<EQSimGeomTriangleMap::iterator, bool> res;

    if (triangles.size() == 0) new_index = 1;
    else new_index = triangles.rbegin()->first + 1;

    new_triangle.set_index(new_index);
    res = triangles.insert(std::make_pair(new_index, new_triangle));

    if (!res.second) throw std::out_of_range("new_triangle");

    return res.first->second;
}

quakelib::EQSimGeometryRectangle &quakelib::EQSimGeometrySection::new_rectangle(void) throw(std::out_of_range) {
    UIndex                                      new_index;
    EQSimGeometryRectangle                      new_rectangle;
    std::pair<EQSimGeomRectangleMap::iterator, bool>    res;

    if (rectangles.size() == 0) new_index = 1;
    else new_index = rectangles.rbegin()->first + 1;

    new_rectangle.set_index(new_index);
    res = rectangles.insert(std::make_pair(new_index, new_rectangle));

    if (!res.second) throw std::out_of_range("new_rectangle");

    return res.first->second;
}

quakelib::EQSimGeometrySection &quakelib::EQSimGeometry::new_section(void) throw(std::out_of_range) {
    UIndex                                  new_sid;
    EQSimGeometrySection                    new_section;
    std::pair<EQSimGeomSectionMap::iterator, bool>  res;

    if (sections.size() == 0) new_sid = 1;
    else new_sid = sections.rbegin()->first + 1;

    new_section.set_sid(new_sid);
    res = sections.insert(std::make_pair(new_sid, new_section));

    if (!res.second) throw std::out_of_range("new_section");

    return res.first->second;
}

bool quakelib::EQSimGeometrySection::is_perfect_rect(const Vec<3> v[4], const double &tolerance) {
    Vec<3>      v01, v02, v03, v12, v23, v30;
    double      vol, rel_diff1, rel_diff2, rel_vol;

    // Check if the points are coplanar by getting the volume
    // of the parallelepiped formed by vectors between them
    v01 = v[1] - v[0];
    v02 = v[2] - v[0];
    v03 = v[3] - v[0];
    vol = v01.dot_product(v02.cross(v03));
    // Check that the volume is relatively small compared to the edge size
    rel_vol = fabs(vol / (v01.mag() * v02.mag() * v03.mag()));

    if (rel_vol < tolerance) {
        // Ensure that parallel sides are of the same length (within a specified relative tolerance)
        v12 = v[2] - v[1];
        v23 = v[3] - v[2];
        v30 = v[0] - v[3];
        rel_diff1 = fabs((v01.mag()-v23.mag())/v01.mag());
        rel_diff2 = fabs((v12.mag()-v30.mag())/v12.mag());

        if (rel_diff1 < tolerance && rel_diff2 < tolerance) return true;
    }

    return false;
}

size_t quakelib::EQSimGeometry::num_vertices(void) const {
    EQSimGeomSectionMap::const_iterator it;
    unsigned int        total = 0;

    for (it=sections.begin(); it!=sections.end(); ++it) total += it->second.num_vertices();

    return total;
}

size_t quakelib::EQSimGeometry::num_triangles(void) const {
    EQSimGeomSectionMap::const_iterator it;
    unsigned int        total = 0;

    for (it=sections.begin(); it!=sections.end(); ++it) total += it->second.num_triangles();

    return total;
}

size_t quakelib::EQSimGeometry::num_rectangles(void) const {
    EQSimGeomSectionMap::const_iterator it;
    unsigned int        total = 0;

    for (it=sections.begin(); it!=sections.end(); ++it) total += it->second.num_rectangles();

    return total;
}

quakelib::ModelSection quakelib::EQSimGeometrySection::create_model_section(void) const {
    ModelSection  new_section;

    new_section.set_id(_sid);
    new_section.set_fault_id(_fid);
    new_section.set_name(_name);

    return new_section;
}

quakelib::ModelElement quakelib::EQSimGeometryRectangle::create_model_element(void) const {
    ModelElement    new_element;

    new_element.set_id(_index);
    new_element.set_vertex(0, _vertex[0]);
    new_element.set_vertex(1, _vertex[1]);
    new_element.set_vertex(2, _vertex[3]);
    new_element.set_is_quad(true);
    new_element.set_rake(Conversion::deg2rad(_rake));
    new_element.set_slip_rate(_slip_rate);
    new_element.set_aseismic(_aseis_factor);

    return new_element;
}

quakelib::ModelVertex quakelib::EQSimGeometryVertex::create_model_vertex(const quakelib::LatLonDepth &base) const {
    ModelVertex     new_vertex;

    new_vertex.set_id(_index);
    new_vertex.set_lld(_loc, base);
    new_vertex.set_das(_das);

    return new_vertex;
}

double quakelib::EQSimGeometry::lat_hi(void) const {
    EQSimGeomSectionMap::const_iterator it;
    double      cur_high = -INFINITY;

    for (it=sections.begin(); it!=sections.end(); ++it) cur_high = fmax(cur_high, it->second.lat_hi());

    return cur_high;
}

double quakelib::EQSimGeometrySection::lat_hi(void) const {
    EQSimGeomVertexMap::const_iterator  it;
    double      cur_high = -INFINITY;

    for (it=vertices.begin(); it!=vertices.end(); ++it) cur_high = fmax(cur_high, it->second.loc().lat());

    return cur_high;
}

double quakelib::EQSimGeometry::lat_lo(void) const {
    EQSimGeomSectionMap::const_iterator it;
    double      cur_low = INFINITY;

    for (it=sections.begin(); it!=sections.end(); ++it) cur_low = fmin(cur_low, it->second.lat_lo());

    return cur_low;
}

double quakelib::EQSimGeometrySection::lat_lo(void) const {
    EQSimGeomVertexMap::const_iterator  it;
    double      cur_low = INFINITY;

    for (it=vertices.begin(); it!=vertices.end(); ++it) cur_low = fmin(cur_low, it->second.loc().lat());

    return cur_low;
}

double quakelib::EQSimGeometry::lon_hi(void) const {
    EQSimGeomSectionMap::const_iterator it;
    double      cur_high = -INFINITY;

    for (it=sections.begin(); it!=sections.end(); ++it) cur_high = fmax(cur_high, it->second.lon_hi());

    return cur_high;
}

double quakelib::EQSimGeometrySection::lon_hi(void) const {
    EQSimGeomVertexMap::const_iterator  it;
    double      cur_high = -INFINITY;

    for (it=vertices.begin(); it!=vertices.end(); ++it) cur_high = fmax(cur_high, it->second.loc().lon());

    return cur_high;
}

double quakelib::EQSimGeometry::lon_lo(void) const {
    EQSimGeomSectionMap::const_iterator it;
    double      cur_low = INFINITY;

    for (it=sections.begin(); it!=sections.end(); ++it) cur_low = fmin(cur_low, it->second.lon_lo());

    return cur_low;
}

double quakelib::EQSimGeometrySection::lon_lo(void) const {
    EQSimGeomVertexMap::const_iterator  it;
    double      cur_low = INFINITY;

    for (it=vertices.begin(); it!=vertices.end(); ++it) cur_low = fmin(cur_low, it->second.loc().lon());

    return cur_low;
}

double quakelib::EQSimGeometry::depth_hi(void) const {
    EQSimGeomSectionMap::const_iterator it;
    double      cur_high = -INFINITY;

    for (it=sections.begin(); it!=sections.end(); ++it) cur_high = fmax(cur_high, it->second.depth_hi());

    return cur_high;
}

double quakelib::EQSimGeometrySection::depth_hi(void) const {
    EQSimGeomVertexMap::const_iterator  it;
    double      cur_high = -INFINITY;

    for (it=vertices.begin(); it!=vertices.end(); ++it) cur_high = fmax(cur_high, it->second.loc().altitude());

    return cur_high;
}

double quakelib::EQSimGeometry::depth_lo(void) const {
    EQSimGeomSectionMap::const_iterator it;
    double      cur_low = INFINITY;

    for (it=sections.begin(); it!=sections.end(); ++it) cur_low = fmin(cur_low, it->second.depth_lo());

    return cur_low;
}

double quakelib::EQSimGeometrySection::depth_lo(void) const {
    EQSimGeomVertexMap::const_iterator  it;
    double      cur_low = INFINITY;

    for (it=vertices.begin(); it!=vertices.end(); ++it) cur_low = fmin(cur_low, it->second.loc().altitude());

    return cur_low;
}

double quakelib::EQSimGeometry::das_hi(void) const {
    EQSimGeomSectionMap::const_iterator it;
    double      cur_high = -INFINITY;

    for (it=sections.begin(); it!=sections.end(); ++it) cur_high = fmax(cur_high, it->second.das_hi());

    return cur_high;
}

double quakelib::EQSimGeometrySection::das_hi(void) const {
    EQSimGeomVertexMap::const_iterator  it;
    double      cur_high = -INFINITY;

    for (it=vertices.begin(); it!=vertices.end(); ++it) cur_high = fmax(cur_high, it->second.das());

    return cur_high;
}

double quakelib::EQSimGeometry::das_lo(void) const {
    EQSimGeomSectionMap::const_iterator it;
    double      cur_low = INFINITY;

    for (it=sections.begin(); it!=sections.end(); ++it) cur_low = fmin(cur_low, it->second.das_lo());

    return cur_low;
}

double quakelib::EQSimGeometrySection::das_lo(void) const {
    EQSimGeomVertexMap::const_iterator  it;
    double      cur_low = INFINITY;

    for (it=vertices.begin(); it!=vertices.end(); ++it) cur_low = fmin(cur_low, it->second.das());

    return cur_low;
}

bool quakelib::EQSimGeometryReader::parse_line(const int &rec_num, const int &line_num, std::istringstream &line_stream) {
    switch (rec_num) {
        case 200:
            parse_fault_summary_record(line_num, line_stream);
            return true;

        case 201:
            parse_section_record(line_num, line_stream);
            return true;

        case 202:
            parse_vertex_record(line_num, line_stream);
            return true;

        case 203:
            parse_triangle_record(line_num, line_stream);
            return true;

        case 204:
            parse_rectangle_record(line_num, line_stream);
            return true;

        default:
            return false;
    }
}

void quakelib::EQSimGeometryReader::parse_fault_summary_record(const int &line_num, std::istringstream &line_stream) {
    _line_num = line_num;
    line_stream >> _n_section   // Field 1: Total number of fault sections in the file
                >> _n_vertex    // Field 2: Total number of vertices in the file
                >> _n_triangle  // Field 3: Total number of triangles in the file
                >> _n_rectangle // Field 4: Total number of rectangles in the file
                >> _lat_lo      // Field 5: Lowest value of latitude in the file (decimal degrees, positive north)
                >> _lat_hi      // Field 6: Highest value of latitude in the file (decimal degrees, positive north)
                >> _lon_lo      // Field 7: Lowest value of longitude in the file (decimal degrees, positive east)
                >> _lon_hi      // Field 8: Highest value of longitude in the file (decimal degrees, positive east)
                >> _depth_lo    // Field 9: Lowest value of depth in the file (meters, negative underground)
                >> _depth_hi;   // Field 10: Highest value of depth in the file (meters, negative underground)
}

void quakelib::EQSimParsedGeometrySection::parse(const int &line_num, std::istringstream &line_stream) {
    _line_num = line_num;
    line_stream >> _sid         // Field 1: Section identification number (positive integer, may not be consecutive)
                >> _name        // Field 2: Section name
                >> _n_vertex    // Field 3: Total number of vertices in the section
                >> _n_triangle  // Field 4: Total number of triangles in the section
                >> _n_rectangle // Field 5: Total number of rectangles in the section
                >> _lat_lo      // Field 6: Lowest value of latitude in the section (decimal degrees, positive north)
                >> _lat_hi      // Field 7: Highest value of latitude in the section (decimal degrees, positive north)
                >> _lon_lo      // Field 8: Lowest value of longitude in the section (decimal degrees, positive east)
                >> _lon_hi      // Field 9: Highest value of longitude in the section (decimal degrees, positive east)
                >> _depth_lo    // Field 10: Lowest value of depth in the section (meters, negative underground)
                >> _depth_hi    // Field 11: Highest value of depth in the section (meters, negative underground)
                >> _das_lo      // Field 12: Lowest value of distance-along-strike in the section (meters)
                >> _das_hi      // Field 13: Highest value of distance-along-strike in the section (meters)
                >> _fid;        // Field 14: Fault identification number (positive integer)
}

void quakelib::EQSimGeometryReader::parse_section_record(const int &line_num, std::istringstream &line_stream) {
    EQSimParsedGeometrySection  pgs;
    EQSimGeometrySection        gs;

    pgs.parse(line_num, line_stream);
    gs.set_sid(pgs._sid);
    gs.set_name(pgs._name);
    gs.set_fid(pgs._fid);

    cur_section = gs.sid();

    _parsed_sections.insert(std::make_pair(cur_section, pgs));
    sections.insert(std::make_pair(cur_section, gs));
}

void quakelib::EQSimGeometryVertex::parse(const int &line_num, std::istringstream &line_stream) {
    double      new_lat, new_lon, new_depth;
    int         trace_flag;

    _line_num = line_num;
    line_stream >> _index       // Field 1: Vertex index number
                >> new_lat      // Field 2: Latitude (decimal degrees, positive north)
                >> new_lon      // Field 3: Longitude (decimal degrees, positive east)
                >> new_depth    // Field 4: Depth (meters, negative underground)
                >> _das         // Field 5: Distance-along-strike (meters)
                >> trace_flag;  // Field 6: Trace flag (0 = not on trace, 1 = on trace but not initial or final, 2 = initial point on trace, 3 = final point on trace)

    switch (trace_flag) {
        case 0:
            this->set_trace_flag(NOT_ON_TRACE);
            break;

        case 1:
            this->set_trace_flag(MIDDLE_TRACE);
            break;

        case 2:
            this->set_trace_flag(BEGINNING_TRACE);
            break;

        case 3:
            this->set_trace_flag(END_TRACE);
            break;

        default:
            this->set_trace_flag(UNDEFINED_TRACE_STATUS);
            break;
    }

    this->set_loc(LatLonDepth(new_lat, new_lon, new_depth));
}

void quakelib::EQSimGeometryReader::parse_vertex_record(const int &line_num, std::istringstream &line_stream) {
    quakelib::EQSimGeometryVertex       gv;

    gv.parse(line_num, line_stream);

    sections[cur_section].vertices.insert(std::make_pair(gv.index(), gv));
}

void quakelib::EQSimGeometryTriangle::parse(const int &line_num, std::istringstream &line_stream) {
    _line_num = line_num;
    line_stream >> _index           // Field 1: Vertex index number
                >> _vertex[0]       // Field 2: Vertex index number for corner #1 (counting counterclockwise as viewed from positive side of element)
                >> _vertex[1]       // Field 3: Vertex index number for corner #2 (counting counterclockwise as viewed from positive side of element)
                >> _vertex[2]       // Field 4: Vertex index number for corner #3 (counting counterclockwise as viewed from positive side of element)
                >> _rake            // Field 5: Rake angle (decimal degrees, 0.0 = left lateral, 90.0 = positive side moves up)
                >> _slip_rate       // Field 6: Element slip rate (meters/second)
                >> _aseis_factor;   // Field 7: Element aseismicity factor
}

void quakelib::EQSimGeometryReader::parse_triangle_record(const int &line_num, std::istringstream &line_stream) {
    quakelib::EQSimGeometryTriangle     gt;

    gt.parse(line_num, line_stream);

    sections[cur_section].triangles.insert(std::make_pair(gt.index(), gt));
}

void quakelib::EQSimGeometryRectangle::parse(const int &line_num, std::istringstream &line_stream) {
    _line_num = line_num;

    line_stream >> _index       // Field 1: Rectangle index number
                >> _vertex[0]   // Field 2: Vertex index number for corner #1 (counting counterclockwise as viewed from positive side of element)
                >> _vertex[1]   // Field 3: Vertex index number for corner #2 (counting counterclockwise as viewed from positive side of element)
                >> _vertex[2]   // Field 4: Vertex index number for corner #3 (counting counterclockwise as viewed from positive side of element)
                >> _vertex[3]   // Field 5: Vertex index number for corner #4 (counting counterclockwise as viewed from positive side of element)
                >> _rake        // Field 6: Rake angle (decimal degrees, 0.0 = left lateral, 90.0 = positive side moves up)
                >> _slip_rate   // Field 7: Element slip rate (meters/second)
                >> _aseis_factor// Field 8: Element aseismicity factor
                >> _strike      // Field 9: Strike angle (decimal degrees)
                >> _dip         // Field 10: Dip angle (decimal degrees)
                >> _perfect_flag;// Field 11: Perfect flag (0 = not perfect rectangle, 1 = perfect rectangle)
}

void quakelib::EQSimGeometryReader::parse_rectangle_record(const int &line_num, std::istringstream &line_stream) {
    quakelib::EQSimGeometryRectangle        gr;

    gr.parse(line_num, line_stream);

    sections[cur_section].rectangles.insert(std::make_pair(gr.index(), gr));
}

void quakelib::EQSimGeometryVertex::validate(EQSimErrors &errors) const {
    // Check that the trace flag is from 0 to 3
    if (_trace_flag < 0 || _trace_flag > 3) {
        std::stringstream       error_msg;
        error_msg << "Trace flag of vertex (" << _index
                  << ") must be 0, 1, 2 or 3 (currently " << _trace_flag << ").";
        errors.report(error_msg.str(), _line_num);
    }
}

void quakelib::EQSimGeometryTriangle::validate(EQSimErrors &errors) const {
    // TODO: assert vertices are in correct order

    // Check that the rectangle aseismicity factor is in acceptable bounds
    if (_aseis_factor < 0 || _aseis_factor > 1) {
        std::stringstream       error_msg;
        error_msg << "Aseismicity factor " << _aseis_factor << " of triangle (" << _index
                  << ") is outside acceptable bounds.";
        errors.report(error_msg.str(), _line_num);
    }
}

void quakelib::EQSimGeometryRectangle::validate(EQSimErrors &errors) const {
    // TODO: assert vertices are in correct order

    // Check that the rectangle aseismicity factor is in acceptable bounds
    if (_aseis_factor < 0 || _aseis_factor > 1) {
        std::stringstream       error_msg;
        error_msg << "Aseismicity factor " << _aseis_factor << " of rectangle (" << _index
                  << ") is outside acceptable bounds.";
        errors.report(error_msg.str(), _line_num);
    }

    // Check that the rectangle perfect flag is either 0 or 1
    if (_perfect_flag != 0 && _perfect_flag != 1) {
        std::stringstream       error_msg;
        error_msg << "Perfect flag of rectangle (" << _index
                  << ") must be either 0 or 1 (currently " << _perfect_flag << ").";
        errors.report(error_msg.str(), _line_num);
    }
}

void quakelib::EQSimGeometrySection::validate(EQSimErrors &errors) const {
    EQSimGeomTriangleMap::const_iterator        tit;
    EQSimGeomRectangleMap::const_iterator       rit;
    EQSimGeomVertexMap::const_iterator          vit;
    double                                  das1, das2, das3, das4, dep1, dep2, dep3, dep4;
    double                                  min_das, max_das, min_dep, max_dep;

    for (vit=vertices.begin(); vit!=vertices.end(); ++vit) {
        // Perform internal vertex correctness check
        vit->second.validate(errors);

        // Check that vertex index matches the map index
        if (vit->second.index() != vit->first) {
            std::stringstream       error_msg;
            error_msg << "Index mismatch for vertex " << vit->first << " (internal index is " << vit->second.index() << ").";
            errors.report(error_msg.str(), vit->second.line_num());
        }
    }

    for (tit=triangles.begin(); tit!=triangles.end(); ++tit) {
        // Perform internal triangle correctness check
        tit->second.validate(errors);

        // Check that triangle index matches the map index
        if (tit->second.index() != tit->first) {
            std::stringstream       error_msg;
            error_msg << "Index mismatch for triangle " << tit->first << " (internal index is " << tit->second.index() << ").";
            errors.report(error_msg.str(), tit->second.line_num());
        }

        // Check that the vertices used by triangles are within their own sections
        // Check vertex 1
        if (vertices.count(tit->second.vertex(0)) == 0) {
            std::stringstream       error_msg;
            error_msg << "Vertex 1 of triangle " << tit->first << " not found in section " << _sid << ".";
            errors.report(error_msg.str(), tit->second.line_num());
        }

        // Check vertex 2
        if (vertices.count(tit->second.vertex(1)) == 0) {
            std::stringstream       error_msg;
            error_msg << "Vertex 2 of triangle " << tit->first << " not found in section " << _sid << ".";
            errors.report(error_msg.str(), tit->second.line_num());
        }

        // Check vertex 3
        if (vertices.count(tit->second.vertex(2)) == 0) {
            std::stringstream       error_msg;
            error_msg << "Vertex 3 of triangle " << tit->first << " not found in section " << _sid << ".";
            errors.report(error_msg.str(), tit->second.line_num());
        }
    }

    for (rit=rectangles.begin(); rit!=rectangles.end(); ++rit) {
        // Perform internal rectangle correctness check
        rit->second.validate(errors);

        // Check that rectangle index matches the map index
        if (rit->second.index() != rit->first) {
            std::stringstream       error_msg;
            error_msg << "Index mismatch for rectangle " << rit->first << " (internal index is " << rit->second.index() << ").";
            errors.report(error_msg.str(), rit->second.line_num());
        }

        // Check that the vertices used by rectangles are within their own sections
        // Check vertex 1
        if (vertices.count(rit->second.vertex(0)) == 0) {
            std::stringstream       error_msg;
            error_msg << "Vertex 1 of rectangle " << rit->first << " not found in section " << _sid << ".";
            errors.report(error_msg.str(), rit->second.line_num());
        }

        // Check vertex 2
        if (vertices.count(rit->second.vertex(1)) == 0) {
            std::stringstream       error_msg;
            error_msg << "Vertex 2 of rectangle " << rit->first << " not found in section " << _sid << ".";
            errors.report(error_msg.str(), rit->second.line_num());
        }

        // Check vertex 3
        if (vertices.count(rit->second.vertex(2)) == 0) {
            std::stringstream       error_msg;
            error_msg << "Vertex 3 of rectangle " << rit->first << " not found in section " << _sid << ".";
            errors.report(error_msg.str(), rit->second.line_num());
        }

        // Check vertex 4
        if (vertices.count(rit->second.vertex(3)) == 0) {
            std::stringstream       error_msg;
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
            std::stringstream       error_msg;
            error_msg << "Rectangle " << rit->first << " is not rectangular in DAS.";
            errors.report(error_msg.str(), rit->second.line_num());
        }

        if ((dep1 != min_dep && dep1 != max_dep) || (dep2 != min_dep && dep2 != max_dep) ||
            (dep3 != min_dep && dep3 != max_dep) || (dep4 != min_dep && dep4 != max_dep)) {
            std::stringstream       error_msg;
            error_msg << "Rectangle " << rit->first << " is not rectangular in depth.";
            errors.report(error_msg.str(), rit->second.line_num());
        }
    }
}

// TODO: add check for consecutively increasing indices
void quakelib::EQSimGeometry::validate(EQSimErrors &errors) const {
    EQSimGeomSectionMap::const_iterator             sit;
    EQSimParsedGeomSectionMap::const_iterator       pit;
    EQSimGeomVertexMap::const_iterator              vit;

    // Assert correctness of each parsed section
    for (pit=_parsed_sections.begin(); pit!=_parsed_sections.end(); ++pit) {
        sit = sections.find(pit->second._sid);

        // Check that the section has the correct number of vertices
        if (sit->second.num_vertices() != pit->second._n_vertex) {
            std::stringstream       error_msg;
            error_msg << "Specified number of vertices (" << pit->second._n_vertex << ") for section "
                      << pit->second._sid << " does not match actual count in file (" << sit->second.num_vertices() << ").";
            errors.report(error_msg.str(), pit->second._line_num);
        }

        // Check that the section has the correct number of triangles
        if (sit->second.num_triangles() != pit->second._n_triangle) {
            std::stringstream       error_msg;
            error_msg << "Specified number of triangles (" << pit->second._n_triangle << ") for section "
                      << pit->second._sid << " does not match actual count in file (" << sit->second.num_triangles() << ").";
            errors.report(error_msg.str(), pit->second._line_num);
        }

        // Check that the section has the correct number of rectangles
        if (sit->second.num_rectangles() != pit->second._n_rectangle) {
            std::stringstream       error_msg;
            error_msg << "Specified number of rectangles (" << pit->second._n_rectangle << ") for section "
                      << pit->second._sid << " does not match actual count in file (" << sit->second.num_rectangles() << ").";
            errors.report(error_msg.str(), pit->second._line_num);
        }

        // Check that the bounds of each vertex is within the section bounds
        for (vit=sit->second.vertices.begin(); vit!=sit->second.vertices.end(); ++vit) {
            // Check that the vertex latitude is within the section bounds
            if (vit->second.loc().lat() < pit->second._lat_lo || vit->second.loc().lat() > pit->second._lat_hi) {
                std::stringstream       error_msg;
                error_msg << "Vertex " << vit->first << " latitude (" << vit->second.loc().lat() << ") is outside section bounds.";
                errors.report(error_msg.str(), vit->second.line_num());
            }

            // Check that the vertex longitude is within the section bounds
            if (vit->second.loc().lon() < pit->second._lon_lo || vit->second.loc().lon() > pit->second._lon_hi) {
                std::stringstream       error_msg;
                error_msg << "Vertex " << vit->first << " longitude (" << vit->second.loc().lon() << ") is outside section bounds.";
                errors.report(error_msg.str(), vit->second.line_num());
            }

            // Check that the vertex depth is within the section bounds
            if (vit->second.loc().altitude() < pit->second._depth_lo || vit->second.loc().altitude() > pit->second._depth_hi) {
                std::stringstream       error_msg;
                error_msg << "Vertex " << vit->first << " depth (" << vit->second.loc().altitude() << ") is outside section bounds.";
                errors.report(error_msg.str(), vit->second.line_num());
            }

            // Check that the vertex DAS is within the section bounds
            if (vit->second.das() < pit->second._das_lo || vit->second.das() > pit->second._das_hi) {
                std::stringstream       error_msg;
                error_msg << "Vertex " << vit->first << " DAS (" << vit->second.das() << ") is outside section bounds.";
                errors.report(error_msg.str(), vit->second.line_num());
            }
        }
    }

    // Assert correctness of each section
    for (sit=sections.begin(); sit!=sections.end(); ++sit) {
        sit->second.validate(errors);
    }
}

void quakelib::EQSimGeometryReader::validate(EQSimErrors &errors) const {
    quakelib::EQSimMetadataReader::validate(errors);

    // Check that the number of sections was the same as specified
    if (sections.size() != _n_section) {
        std::stringstream       error_msg;
        error_msg << "Specified number of fault sections (" << _n_section
                  << ") does not match actual count in file (" << sections.size() << ").";
        errors.report(error_msg.str(), _line_num);
    }

    // Check that the number of vertices was the same as specified
    if (num_vertices() != _n_vertex) {
        std::stringstream       error_msg;
        error_msg << "Specified number of vertices (" << _n_vertex
                  << ") does not match actual count in file (" << num_vertices() << ").";
        errors.report(error_msg.str(), _line_num);
    }

    // Check that the number of triangles was the same as specified
    if (num_triangles() != _n_triangle) {
        std::stringstream       error_msg;
        error_msg << "Specified number of triangles (" << _n_triangle
                  << ") does not match actual count in file (" << num_triangles() << ").";
        errors.report(error_msg.str(), _line_num);
    }

    // Check that the number of triangles was the same as specified
    if (num_rectangles() != _n_rectangle) {
        std::stringstream       error_msg;
        error_msg << "Specified number of rectangles (" << _n_rectangle
                  << ") does not match actual count in file (" << num_rectangles() << ").";
        errors.report(error_msg.str(), _line_num);
    }

    // Check that the latitude, longitude and depths are within the specified limits
    if (lat_hi() > _lat_hi || lat_lo() < _lat_lo) {
        std::stringstream       error_msg;
        error_msg << "Vertex latitude limits in file (" << lat_lo() << "," << lat_hi()
                  << ") are outside specified latitude limits (" << _lat_lo << "," << _lat_hi << ").";
        errors.report(error_msg.str(), _line_num);
    }

    if (lon_hi() > _lon_hi || lon_lo() < _lon_lo) {
        std::stringstream       error_msg;
        error_msg << "Vertex longitude limits in file (" << lon_lo() << "," << lon_hi()
                  << ") are outside specified longitude limits (" << _lon_lo << "," << _lon_hi << ").";
        errors.report(error_msg.str(), _line_num);
    }

    if (depth_hi() > _depth_hi || depth_lo() < _depth_lo) {
        std::stringstream       error_msg;
        error_msg << "Vertex depth limits in file (" << depth_lo() << "," << depth_hi()
                  << ") are outside specified depth limits (" << _depth_lo << "," << _depth_hi << ").";
        errors.report(error_msg.str(), _line_num);
    }

    // TODO: assert section parsed bounds are within model parse bounds

    EQSimGeometry::validate(errors);
}

quakelib::EQSimGeometryWriter::EQSimGeometryWriter(void) : EQSimGeometry(), EQSimMetadataWriter() {
    internal::RecordDesc    geom_summary_rec, geom_section_rec, geom_vertex_rec, geom_triangle_rec, geom_rectangle_rec;

    set_spec_level(2);
    meta_add_record(META_SIGNATURE, "EQSim_Input_Geometry_2");

    geom_summary_rec = internal::RecordDesc(0, "summary", 10, "Record 200: Fault system summary");
    geom_summary_rec.add_field(1, internal::EQSimFieldDesc("n_section", 1, 1, "Total number of fault sections in the file"));
    geom_summary_rec.add_field(2, internal::EQSimFieldDesc("n_vertex", 1, 2, "Total number of vertices in the file"));
    geom_summary_rec.add_field(3, internal::EQSimFieldDesc("n_triangle", 1, 3, "Total number of triangles in the file"));
    geom_summary_rec.add_field(4, internal::EQSimFieldDesc("n_rectangle", 1, 4, "Total number of rectangles in the file"));
    geom_summary_rec.add_field(5, internal::EQSimFieldDesc("lat_lo", 2, 5, "Lowest value of latitude in the file (decimal degrees, positive north)"));
    geom_summary_rec.add_field(6, internal::EQSimFieldDesc("lat_hi", 2, 6, "Highest value of latitude in the file (decimal degrees, positive north)"));
    geom_summary_rec.add_field(7, internal::EQSimFieldDesc("lon_lo", 2, 7, "Lowest value of longitude in the file (decimal degrees, positive east)"));
    geom_summary_rec.add_field(8, internal::EQSimFieldDesc("lon_hi", 2, 8, "Highest value of longitude in the file (decimal degrees, positive east)"));
    geom_summary_rec.add_field(9, internal::EQSimFieldDesc("depth_lo", 2, 9, "Lowest value of depth in the file (meters, negative underground)"));
    geom_summary_rec.add_field(10, internal::EQSimFieldDesc("depth_hi", 2, 10, "Highest value of depth in the file (meters, negative underground)"));

    geom_section_rec = internal::RecordDesc(0, "section", 14, "Record 201: Fault section information");
    geom_section_rec.add_field(1, internal::EQSimFieldDesc("sid", 1, 1, "Section identification number (positive integer, may not be consecutive)"));
    geom_section_rec.add_field(2, internal::EQSimFieldDesc("name", 3, 2, "Section name"));
    geom_section_rec.add_field(3, internal::EQSimFieldDesc("n_vertex", 1, 3, "Total number of vertices in the section"));
    geom_section_rec.add_field(4, internal::EQSimFieldDesc("n_triangle", 1, 4, "Total number of triangles in the section"));
    geom_section_rec.add_field(5, internal::EQSimFieldDesc("n_rectangle", 1, 5, "Total number of rectangles in the section"));
    geom_section_rec.add_field(6, internal::EQSimFieldDesc("lat_lo", 2, 6, "Lowest value of latitude in the section (decimal degrees, positive north)"));
    geom_section_rec.add_field(7, internal::EQSimFieldDesc("lat_hi", 2, 7, "Highest value of latitude in the section (decimal degrees, positive north)"));
    geom_section_rec.add_field(8, internal::EQSimFieldDesc("lon_lo", 2, 8, "Lowest value of longitude in the section (decimal degrees, positive east)"));
    geom_section_rec.add_field(9, internal::EQSimFieldDesc("lon_hi", 2, 9, "Highest value of longitude in the section (decimal degrees, positive east)"));
    geom_section_rec.add_field(10, internal::EQSimFieldDesc("depth_lo", 2, 10, "Lowest value of depth in the section (meters, negative underground)"));
    geom_section_rec.add_field(11, internal::EQSimFieldDesc("depth_hi", 2, 11, "Highest value of depth in the section (meters, negative underground)"));
    geom_section_rec.add_field(12, internal::EQSimFieldDesc("das_lo", 2, 12, "Lowest value of distance-along-strike in the section (meters)"));
    geom_section_rec.add_field(13, internal::EQSimFieldDesc("das_hi", 2, 13, "Highest value of distance-along-strike in the section (meters)"));
    geom_section_rec.add_field(14, internal::EQSimFieldDesc("fault_id", 1, 14, "Fault identification number (positive integer)"));

    geom_vertex_rec = internal::RecordDesc(0, "vertex", 6, "Record 202: Vertex");
    geom_vertex_rec.add_field(1, internal::EQSimFieldDesc("index", 1, 1, "Vertex index number (consecutive integers, starting with 1)"));
    geom_vertex_rec.add_field(2, internal::EQSimFieldDesc("lat", 2, 2, "Latitude (decimal degrees, positive north)"));
    geom_vertex_rec.add_field(3, internal::EQSimFieldDesc("lon", 2, 3, "Longitude (decimal degrees, positive east)"));
    geom_vertex_rec.add_field(4, internal::EQSimFieldDesc("depth", 2, 4, "Depth (meters, negative underground)"));
    geom_vertex_rec.add_field(5, internal::EQSimFieldDesc("das", 2, 5, "Distance-along-strike (meters)"));
    geom_vertex_rec.add_field(6, internal::EQSimFieldDesc("trace_flag", 1, 6, "Trace flag (0 = not on trace, 1 = on trace but not initial or final, 2 = initial point on trace, 3 = final point on trace)"));

    geom_triangle_rec = internal::RecordDesc(0, "triangle", 9, "Record 203: Triangle");
    geom_triangle_rec.add_field(1, internal::EQSimFieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
    geom_triangle_rec.add_field(2, internal::EQSimFieldDesc("vertex_1", 1, 2, "Vertex index number for corner #1 (counting counterclockwise as viewed from positive side of element)"));
    geom_triangle_rec.add_field(3, internal::EQSimFieldDesc("vertex_2", 1, 3, "Vertex index number for corner #2 (counting counterclockwise as viewed from positive side of element)"));
    geom_triangle_rec.add_field(4, internal::EQSimFieldDesc("vertex_3", 1, 4, "Vertex index number for corner #3 (counting counterclockwise as viewed from positive side of element)"));
    geom_triangle_rec.add_field(5, internal::EQSimFieldDesc("rake", 2, 5, "Rake angle (decimal degrees)"));
    geom_triangle_rec.add_field(6, internal::EQSimFieldDesc("slip_rate", 2, 6, "Element slip rate (meters/second)"));
    geom_triangle_rec.add_field(7, internal::EQSimFieldDesc("aseis_factor", 2, 7, "Element aseismicity factor"));
    geom_triangle_rec.add_field(8, internal::EQSimFieldDesc("strike", 2, 8, "Strike angle (decimal degrees)"));
    geom_triangle_rec.add_field(9, internal::EQSimFieldDesc("dip", 2, 9, "Dip angle (decimal degrees)"));

    geom_rectangle_rec = internal::RecordDesc(0, "rectangle", 11, "Record 204: Rectangle");
    geom_rectangle_rec.add_field(1, internal::EQSimFieldDesc("index", 1, 1, "Element index number (consecutive integers, starting with 1)"));
    geom_rectangle_rec.add_field(2, internal::EQSimFieldDesc("vertex_1", 1, 2, "Vertex index number for corner #1 (counting counterclockwise as viewed from positive side of element)"));
    geom_rectangle_rec.add_field(3, internal::EQSimFieldDesc("vertex_2", 1, 3, "Vertex index number for corner #2 (counting counterclockwise as viewed from positive side of element)"));
    geom_rectangle_rec.add_field(4, internal::EQSimFieldDesc("vertex_3", 1, 4, "Vertex index number for corner #3 (counting counterclockwise as viewed from positive side of element)"));
    geom_rectangle_rec.add_field(5, internal::EQSimFieldDesc("vertex_4", 1, 5, "Vertex index number for corner #4 (counting counterclockwise as viewed from positive side of element)"));
    geom_rectangle_rec.add_field(6, internal::EQSimFieldDesc("rake", 2, 6, "Rake angle (decimal degrees)"));
    geom_rectangle_rec.add_field(7, internal::EQSimFieldDesc("slip_rate", 2, 7, "Element slip rate (meters/second)"));
    geom_rectangle_rec.add_field(8, internal::EQSimFieldDesc("aseis_factor", 2, 8, "Element aseismicity factor"));
    geom_rectangle_rec.add_field(9, internal::EQSimFieldDesc("strike", 2, 9, "Strike angle (decimal degrees)"));
    geom_rectangle_rec.add_field(10, internal::EQSimFieldDesc("dip", 2, 10, "Dip angle (decimal degrees)"));
    geom_rectangle_rec.add_field(11, internal::EQSimFieldDesc("perfect_flag", 1, 11, "Perfect flag (0 = not perfect rectangle, 1 = perfect rectangle)"));

    add_record_desc_record(200, geom_summary_rec);
    add_record_desc_record(201, geom_section_rec);
    add_record_desc_record(202, geom_vertex_rec);
    add_record_desc_record(203, geom_triangle_rec);
    add_record_desc_record(204, geom_rectangle_rec);
}
