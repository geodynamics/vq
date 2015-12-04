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

#include "config.h"
#include "QuakeLibIO.h"
#include "QuakeLibEQSim.h"

#include <iomanip>
#include <getopt.h>

/* options descriptor */
static struct option longopts[] = {
    { "export_file",                required_argument,      NULL,           'e' },
    { "export_file_type",           required_argument,      NULL,           'f' },
    { "import_file",                required_argument,      NULL,           'i' },
    { "import_file_type",           required_argument,      NULL,           'j' },
    { "export_eqsim_condition",     required_argument,      NULL,           'D' },
    { "export_eqsim_friction",      required_argument,      NULL,           'R' },
    { "export_eqsim_geometry",      required_argument,      NULL,           'M' },
    { "import_eqsim_condition",     required_argument,      NULL,           'C' },
    { "import_eqsim_friction",      required_argument,      NULL,           'F' },
    { "import_eqsim_geometry",      required_argument,      NULL,           'G' },
    { "import_trace_element_size",  required_argument,      NULL,           'l' },
    { "merge_duplicate_verts",      no_argument,            NULL,           'm' },
    { "delete_unused",              no_argument,            NULL,           'd' },
    { "taper_fault_method",         required_argument,      NULL,           't' },
    { "print_statistics",           required_argument,      NULL,           's' },
    { NULL,                         0,                      NULL,           0 }
};

//    { "do_not_compute_stress_drops",no_argument,            NULL,           'x' },
//    { "stress_drop_factor",         required_argument,      NULL,           'q' },

std::string mem_string(const double &num_bytes) {
    std::stringstream       ss;
    int unit_place = (num_bytes > 0 ? (int)(log(num_bytes)/log(1024)) : 0);
    ss << num_bytes/pow(1024, unit_place) << " " << std::string(" KMGTPE").at(unit_place) << "B";
    return ss.str();
}

void print_statistics(quakelib::ModelWorld &world, const std::string &file_name) {
    size_t                  num_elements, num_vertices;
    quakelib::UIndex        sid;
    quakelib::siterator     sit;
    quakelib::eiterator     eit;
    std::vector<double>     rake_vals, slip_rate_vals;
    std::ofstream           out_file;
    int                     section_field_width, section_name_width;
    int                     elem_field_width, vert_field_width, val_field_width;
    double                  mem_req;
    quakelib::Conversion    c;

    out_file.open(file_name.c_str());

    out_file << "STATISTICS\n";
    out_file << "Number of Sections: " << world.num_sections() << "\n";
    out_file << "Number of Faults: " << world.num_faults() << "\n";
    out_file << "Number of Elements: " << world.num_elements() << "\n";
    out_file << "Number of Vertices: " << world.num_vertices() << "\n";
    mem_req = (world.num_elements()*world.num_elements())*2*8;
    out_file << "Expected memory requirement: " << mem_string(mem_req) << "\n";

    section_field_width = fmax(5, log10(world.num_sections())+1);
    section_name_width = 10;
    elem_field_width = fmax(5, log10(world.num_elements())+1);
    vert_field_width = fmax(5, log10(world.num_vertices())+1);
    val_field_width = 5;
    out_file << std::setw(section_field_width) << "Section\t";
    out_file << std::setw(section_name_width) << "Name\t";
    out_file << std::setw(elem_field_width) << "Elems\t";
    out_file << std::setw(vert_field_width) << "Verts\t";
    out_file << std::setw(val_field_width) << "\t";
    out_file << std::setw(val_field_width) << "Rake\t";
    out_file << std::setw(val_field_width) << "\t";
    out_file << std::setw(val_field_width) << "\t";
    out_file << std::setw(val_field_width) << "Rate\t";
    out_file << std::setw(val_field_width) << "\t";
    out_file << "\n";

    out_file << std::setw(section_field_width) << "\t";
    out_file << std::setw(section_name_width) << "\t";
    out_file << std::setw(elem_field_width) << "\t";
    out_file << std::setw(vert_field_width) << "\t";
    out_file << std::setw(val_field_width) << "Min\t";
    out_file << std::setw(val_field_width) << "Median\t";
    out_file << std::setw(val_field_width) << "Max\t";
    out_file << std::setw(val_field_width) << "Min\t";
    out_file << std::setw(val_field_width) << "Median\t";
    out_file << std::setw(val_field_width) << "Max\t";
    out_file << "\n";

    for (sit=world.begin_section(); sit!=world.end_section(); ++sit) {
        rake_vals.clear();
        slip_rate_vals.clear();
        sid = sit->id();
        num_elements = world.num_elements(sid);
        num_vertices = world.num_vertices(sid);

        for (eit=world.begin_element(sid); eit!=world.end_element(sid); ++eit) {
            rake_vals.push_back(c.rad2deg(eit->rake()));
            slip_rate_vals.push_back(c.m_per_sec2cm_per_yr(eit->slip_rate()));
        }

        if (num_elements == 0) {
            rake_vals.push_back(std::numeric_limits<double>::quiet_NaN());
            slip_rate_vals.push_back(std::numeric_limits<double>::quiet_NaN());
        }

        std::sort(rake_vals.begin(), rake_vals.end());
        std::sort(slip_rate_vals.begin(), slip_rate_vals.end());

        out_file << std::setw(section_field_width) << sit->id() << "\t";
        out_file << std::setw(section_name_width) << sit->name().substr(0, section_name_width-1) << "\t";
        out_file << std::setw(elem_field_width) << num_elements << "\t";
        out_file << std::setw(vert_field_width) << num_vertices << "\t";

        out_file << std::setw(val_field_width) << rake_vals[0] << "\t";
        out_file << std::setw(val_field_width) << rake_vals[rake_vals.size()/2] << "\t";
        out_file << std::setw(val_field_width) << rake_vals[rake_vals.size()-1] << "\t";

        out_file << std::setw(val_field_width) << slip_rate_vals[0] << "\t";
        out_file << std::setw(val_field_width) << slip_rate_vals[slip_rate_vals.size()/2] << "\t";
        out_file << std::setw(val_field_width) << slip_rate_vals[slip_rate_vals.size()-1] << "\t";
        out_file << "\n";
    }

    out_file.close();
}

void print_usage(int argc, char **argv) {
    std::cerr << "Used for model file manipulation with VC. Imports one or more files, manipulates them based on user arguments and exports one or more files." << std::endl;
    std::cerr << std::endl;

    std::cerr << argv[0] << " [options]" << std::endl;
    std::cerr << "-s FILE, --print_statistics=FILE" << std::endl;
    std::cerr << "\tPrint statistics regarding final model to the specified file." << std::endl;
    std::cerr << "-m, --merge_duplicate_verts" << std::endl;
    std::cerr << "\tMerge duplicate vertices after importing files." << std::endl;
    std::cerr << "-d, --delete_unused" << std::endl;
    std::cerr << "\tDelete unused vertices after importing files." << std::endl;
    std::cerr << "-r, --resize_trace_elements" << std::endl;
    std::cerr << "\tResize elements generated on traces to better match fault length." << std::endl;
    std::cerr << "\tThis will only decrease and at most halve the element size." << std::endl;
    std::cerr << "-x, --do_not_compute_stress_drops" << std::endl;
    std::cerr << "\tDo not compute stress drops, must specify from EQSim friction instead." << std::endl;
    std::cerr << "-q FACTOR, --stress_drop_factor=FACTOR" << std::endl;
    std::cerr << "\tSpecify the stress drop factor (usually 0.2-0.6). It's a multiplier for the computed stress drops (it's logarithmic, 0.1 increase means multiplying stress drops by 1.4)." << std::endl;

    std::cerr << std::endl;
    std::cerr << "FILE IMPORT" << std::endl;
    std::cerr << "-i FILE, --import_file=FILE" << std::endl;
    std::cerr << "\tSpecify a model file to import and merge. Must have a paired import_file_type." << std::endl;
    std::cerr << "-j TYPE, --import_file_type=TYPE" << std::endl;
    std::cerr << "\tSpecify a model file type for importing. Must have a paired import_file." << std::endl;
    std::cerr << "-l SIZE, --import_trace_element_size=SIZE" << std::endl;
    std::cerr << "\tSpecify the element size (in meters) to use for trace file meshing. Must have a paired trace type file import." << std::endl;
    std::cerr << "-t METHOD, --taper_trace_method=METHOD" << std::endl;
    std::cerr << "\tSpecify the how to taper the imported trace when meshing. Must have a paired trace type file import." << std::endl;
    std::cerr << "-C FILE, --import_eqsim_condition=FILE" << std::endl;
    std::cerr << "\tSpecify an EQSim condition file to import for the model." << std::endl;
    std::cerr << "-F FILE, --import_eqsim_friction=FILE" << std::endl;
    std::cerr << "\tSpecify an EQSim friction file to import for the model." << std::endl;
    std::cerr << "-G FILE, --import_eqsim_geometry=FILE" << std::endl;
    std::cerr << "\tSpecify an EQSim geometry file to import for the model." << std::endl;

    std::cerr << std::endl;
    std::cerr << "FILE EXPORT" << std::endl;
    std::cerr << "-e FILE, --export_file=FILE" << std::endl;
    std::cerr << "\tSpecify a file to export the completed model to. Must have a paired export_file_type." << std::endl;
    std::cerr << "-f TYPE, --export_file_type=TYPE" << std::endl;
    std::cerr << "\tSpecify a file type to export the completed model. Must have a paired export_file." << std::endl;
    std::cerr << "-D FILE, --export_eqsim_condition=FILE" << std::endl;
    std::cerr << "\tSpecify an EQSim condition file to export for the model." << std::endl;
    std::cerr << "-R FILE, --export_eqsim_friction=FILE" << std::endl;
    std::cerr << "\tSpecify an EQSim friction file to export for the model." << std::endl;
    std::cerr << "-M FILE, --export_eqsim_geometry=FILE" << std::endl;
    std::cerr << "\tSpecify an EQSim geometry file to export for the model." << std::endl;
}

int main (int argc, char **argv) {
    quakelib::ModelWorld        world;
    bool                        delete_unused, merge_duplicate_vertices, resize_trace_elements, compute_stress_drops=true;
    bool                        arg_error, failed;
    std::string                 names[2] = {"import", "export"};
    std::string                 eqsim_geom_in_file, eqsim_fric_in_file, eqsim_cond_in_file;
    std::string                 eqsim_geom_out_file, eqsim_fric_out_file, eqsim_cond_out_file;
    std::string                 stat_out_file;
    std::vector<std::string>    files[2], types[2];
    std::vector<std::string>    taper_fault_methods;
    std::vector<double>         trace_element_sizes;
    double                      stress_drop_factor = 0.3;  // default value is a reasonable 0.3
    int                         ch, res;
    unsigned int                i, n, j, num_trace_files;

    arg_error = delete_unused = merge_duplicate_vertices = resize_trace_elements = false;
    eqsim_geom_in_file = eqsim_fric_in_file = eqsim_cond_in_file = "";
    eqsim_geom_out_file = eqsim_fric_out_file = eqsim_cond_out_file = "";

    while ((ch = getopt_long(argc, argv, "mdrs:D:R:M:C:F:G:i:j:e:f:l:t:", longopts, NULL)) != -1) {
        switch (ch) {
            case 'd':
                delete_unused = true;
                break;

            case 'm':
                merge_duplicate_vertices = true;
                break;

            case 'r':
                resize_trace_elements = true;
                break;
                
//            case 'x':
//                compute_stress_drops = false;
//                break;

            case 's':
                stat_out_file = optarg;
                break;

//            case 'q':
//                stress_drop_factor = atof(optarg);
//                break;

            case 'D':
                eqsim_cond_out_file = optarg;
                break;

            case 'R':
                eqsim_fric_out_file = optarg;
                break;

            case 'M':
                eqsim_geom_out_file = optarg;
                break;

            case 'C':
                eqsim_cond_in_file = optarg;
                break;

            case 'F':
                eqsim_fric_in_file = optarg;
                break;

            case 'G':
                eqsim_geom_in_file = optarg;
                break;

            case 'i':
                files[0].push_back(optarg);
                break;

            case 'j':
                types[0].push_back(optarg);
                break;

            case 'e':
                files[1].push_back(optarg);
                break;

            case 'f':
                types[1].push_back(optarg);
                break;

            case 'l':
                trace_element_sizes.push_back(atof(optarg));
                break;

            case 't':
                taper_fault_methods.push_back(optarg);
                break;

            default:
                std::cerr << "Unknown argument " << argv[optind] << std::endl;
                print_usage(argc, argv);
                exit(1);
        }
    }

    // Make sure there's at least some input files
    if (files[0].size() == 0 && eqsim_geom_in_file.empty()) {
        std::cerr << "ERROR: Must specify an input file." << std::endl;
        arg_error = true;
    }

    // Check validity of input parameters
    num_trace_files = 0;

    for (i=0; i<2; ++i) {
        if (files[i].size() != types[i].size()) {
            std::cerr << "ERROR: must declare a type for each " << names[i] << "ed file. (" << files[i].size() << " files and " << types[i].size() << " types)" << std::endl;
            arg_error = true;
        }

        for (n=0; n<types[i].size(); ++n) {
            if (types[i][n] != "text" && types[i][n] != "hdf5" && types[i][n] != "trace" && types[i][n] != "kml") {
                std::cerr << "ERROR: " << names[i] << " type must be one of: text, hdf5, trace, kml." << std::endl;
                arg_error = true;
            }

            if (i == 0 && types[i][n] == "trace") num_trace_files++;
        }
    }

    // Check that there are the same number of element sizes as input trace files, unless specifying taper and input from eqsim not trace
    if ((num_trace_files != trace_element_sizes.size()) && eqsim_geom_in_file.empty()) {
        std::cerr << "ERROR: Incorrect number of element sizes (" << trace_element_sizes.size()
                  << ") to input trace files (" << num_trace_files << ").  EQSIM_FILE empty? " << eqsim_geom_in_file.empty() << std::endl;
        arg_error = true;
    }

    // Check that the taper methods are valid
    for (i=0; i<taper_fault_methods.size(); ++i) {
        if (taper_fault_methods[i] != "none" && taper_fault_methods[i] != "taper"
            && taper_fault_methods[i] != "taper_full" && taper_fault_methods[i] != "taper_renorm") {
            std::cerr << "ERROR: Taper method " << taper_fault_methods[i] << " must be one of: none, taper, taper_full, taper_renorm." << std::endl;
            arg_error = true;
        }
    }

    if (!(compute_stress_drops) && (eqsim_fric_in_file.empty() || eqsim_geom_in_file.empty() )) {
        std::cerr << "ERROR: If not computing stress drops, must specify EQSim geometry and EQSim friction file. EQSim friction should contain stress drops." << std::endl;
        arg_error = true;
    }


    // If any argument was malformed, print the correct usage
    if (arg_error) {
        print_usage(argc, argv);
        exit(-1);
    }

    // Print an informational summary for the user
    std::cout << "# *** VQ Mesher Version " << VQ_VERSION_STR << " ***" << std::endl;
#ifdef VQ_GIT_SHA1
    std::cout << "# *** Git revision ID " << VQ_GIT_SHA1 << " ***" << std::endl;
#endif
    std::cout << "# *** " << quakelib::quakelib_info() << " ***" << std::endl;;

    std::cout << "*** Summary of edits ***" << std::endl;

    // *** INPUT SECTION ***
    // Do the actual editing. First step is to import files
    failed = false;

    for (j=0,n=0; n<files[0].size(); ++n) {
        quakelib::ModelWorld    new_world;
        std::vector<unsigned int>   unused_trace_segments;

        std::cout << "File " << names[0] << " " << files[0][n] << " with type " << types[0][n] << "... ";
        unused_trace_segments.clear();

        // TODO: change these to fail if file is not in expected format
        if (types[0][n] == "text") {
            res = new_world.read_file_ascii(files[0][n]);
        } else if (types[0][n] == "hdf5") {
            res = new_world.read_file_hdf5(files[0][n]);
        } else if (types[0][n] == "trace") {
            res = new_world.read_file_trace_latlon(unused_trace_segments, files[0][n], trace_element_sizes.at(j), taper_fault_methods.at(j), resize_trace_elements);
            ++j;
        }

        if (res) {
            std::cout << "*** Error reading file " << files[0][n] << std::endl;
            failed = true;
        } else {
            world.insert(new_world);
            std::cout << "done." << std::endl;

            if (unused_trace_segments.size() > 0) {
                std::cout << "WARNING: Unused trace segment" << (unused_trace_segments.size() > 1 ? "s" : "");

                for (i=0; i<unused_trace_segments.size(); ++i) std::cout << (i == 0 ? " " : ", ") << unused_trace_segments[i];

                std::cout << ". Element size may be too big." << std::endl;
            }
        }
    }

    // Import any EQSim files as well
    // TODO: add console output here
    if (!eqsim_geom_in_file.empty()) {
        quakelib::ModelWorld        new_world;

        new_world.read_files_eqsim(eqsim_geom_in_file, eqsim_cond_in_file, eqsim_fric_in_file, taper_fault_methods.at(0));
        world.insert(new_world);
    }

    if (failed) {
        std::cout << "*** Quitting due to errors." << std::endl;
        return -1;
    }

    // *** PROCESSING/MANIPULATION SECTION ***
    // If requested, merge duplicate vertices
    if (merge_duplicate_vertices) {
        std::cout << "Merge duplicate vertices" << std::endl;
        quakelib::ModelRemapping    remap;
        remap = world.remove_duplicate_vertices_remap();
        world.apply_remap(remap);
    }

    // TODO: If requested, delete unused elements/vertices
    if (delete_unused) {
        std::cout << "Delete unused elements and vertices" << std::endl;
    }

    // ***CREATE ModelFault OBJECTS***
    // Here we take the final fault model and create ModelFault objects, which correctly rewrites DAS for
    // each vertex to be relative to fault, applies horizontal tapering wrt length of fault
    world.create_faults(taper_fault_methods.at(0));


    // Schultz: Moving stress drop computation here (used to be in UpdateBlockStress.cpp.
    // ------------------------------------------------------------------------------------------
    // Compute the stress drops. Default behavior is to compute them. To prescribe stress drops,
    //   one must specify the values in an EQSim friction file.
    if (compute_stress_drops) {
        std::cout << "Computing stress drops with stress_drop_factor=" << stress_drop_factor << std::endl;
        world.compute_stress_drops(stress_drop_factor);
    }


    // *** OUTPUT SECTION ***
    // Finally, export the appropriate file types
    for (n=0; n<files[1].size(); ++n) {
        std::cout << "File " << names[1] << " " << files[1][n] << " with type " << types[1][n] << "... ";

        if (types[1][n] == "text") res = world.write_file_ascii(files[1][n]);
        else if (types[1][n] == "hdf5") res = world.write_file_hdf5(files[1][n]);
        else if (types[1][n] == "kml") res = world.write_file_kml(files[1][n]);
        else if (types[1][n] == "trace") res = world.write_file_trace_latlon();

        if (res) std::cout << "error." << std::endl;
        else std::cout << "done." << std::endl;
    }

    // Export EQSim model if needed
    if (!eqsim_geom_out_file.empty()) {
        world.write_files_eqsim(eqsim_geom_out_file, "", eqsim_fric_out_file);
    }

    // Print final statistics
    if (!stat_out_file.empty()) {
        std::cout << "Print statistics to " << stat_out_file << std::endl;
        print_statistics(world, stat_out_file);
    }

    return 0;
}
