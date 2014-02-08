// Copyright (c) 2012-2014 Eric M. Heien, Michael K. Sachs, John B. Rundle
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

#include "VCParams.h"
#include "ConfigFile.h"

/*!
 Parse a configuration file for each of the possible parameters, assigning default values
 for parameters that are not explicitly specified.
 */
void VCParams::read_params(const std::string &param_file_name) {
    std::string     greens_method, friction_law_str;
    ConfigFile param_file(param_file_name);

    valid = true;

    version = param_file.read<string>("sim.version", "");

    year = param_file.read<double>("sim.start_year", 0.0);
    sim_end_year = param_file.read<double>("sim.end_year", 1000.0);

    noise_event = param_file.read<double>("sim.noise.event", 0.0);
    noise_slip_deficit = param_file.read<double>("sim.noise.slip_deficit", 0.0);

    fault_kill_cff = param_file.read<double>("sim.fault.kill_cff", 0.0);

    checkpoint_period = param_file.read<int>("sim.checkpoint.period", 0);
    checkpoint_save_prefix = param_file.read<string>("sim.checkpoint.prefix", "sim_state_");

    progress_period = param_file.read<unsigned int>("sim.progress.period", 0);

    dynamic = param_file.read<double>("sim.dynamic", -1.0);
    friction_law_str = param_file.read<string>("sim.friction_law", "original");
    slip_scaling_threshold = param_file.read<unsigned int>("sim.slip_scaling_threshold", 10);

    greens_kill_distance = param_file.read<double>("sim.greens.kill_distance", 0.0);
    greens_method = param_file.read<string>("sim.greens.method", "original");
    barnes_hut_theta = param_file.read<double>("sim.greens.bh_theta", 0.0);
    greens_infile = param_file.read<string>("sim.greens.infile", "");

    event_outfile = param_file.read<string>("sim.file.output_event", "");
    sweep_outfile = param_file.read<string>("sim.file.output_sweep", "");
    event_outfile_type = param_file.read<string>("sim.file.output_event_type", "");

    bass_max_generations = param_file.read<unsigned int>("sim.bass.max_generations", 0);
    bass_min_magnitude_mm = param_file.read<double>("sim.bass.min_magnitude_mm", 4.0);
    bass_aftershock_strength_dm = param_file.read<double>("sim.bass.aftershock_strength_dm", 1.25);
    bass_frequency_scale_b = param_file.read<double>("sim.bass.frequency_scale_b", 1.0);
    bass_aftershock_start_c = param_file.read<double>("sim.bass.aftershock_start_c", 0.1);
    bass_time_decay_p = param_file.read<double>("sim.bass.time_decay_p", 1.25);
    bass_distance_d = param_file.read<double>("sim.bass.distance_d", 300);
    bass_distance_decay_q = param_file.read<double>("sim.bass.distance_decay_q", 1.35);

    sanity_check = param_file.read<bool>("sim.sanity_check", false);
    do_normal_stress = param_file.read<bool>("sim.do_normal_stress", true);
    use_transpose_matrix = param_file.read<bool>("sim.use_transpose_matrix", true);

    greens_outfile = param_file.read<string>("sim.greens.outfile", "");
    events_file = param_file.read<string>("sim.events.file", "");

    section_params_file = param_file.read<string>("sim.section_params.file", "");

    hdf5_file = param_file.read<string>("sim.hdf5.file", "vc_result.h5");

    eqsim_condition_file = param_file.read<string>("sim.eqsim.file.condition", "");
    eqsim_friction_file = param_file.read<string>("sim.eqsim.file.friction", "");
    eqsim_geometry_file = param_file.read<string>("sim.eqsim.file.geometry", "");
    eqsim_output_file = param_file.read<string>("sim.eqsim.file.output", "");
    eqsim_slipmap_mag = param_file.read<double>("sim.eqsim.file.slipmap_mag", 7.5);

    // Parse the Greens calculation method string
    if (!greens_method.compare("file")) {
        greens_calc_method = GREENS_FILE_PARSE;
    } else if (!greens_method.compare("bh")) {
        greens_calc_method = GREENS_CALC_BARNES_HUT;
    } else if (!greens_method.compare("standard")) {
        greens_calc_method = GREENS_CALC_STANDARD;
    } else {
        greens_calc_method = GREENS_CALC_UNDEFINED;
    }

    // Parse the friction law method string
    if (!friction_law_str.compare("original")) {
        friction_law_method = FRIC_LAW_ORIG;
    } else if (!friction_law_str.compare("stepped")) {
        friction_law_method = FRIC_LAW_STEPPED;
    } else {
        friction_law_method = FRIC_LAW_UNDEFINED;
    }

    //std::cout << *this;
}

std::ostream &operator<<(std::ostream &os, const GreensCalcMethod &calc_method) {
    switch (calc_method) {
        case GREENS_FILE_PARSE:
            os << "file";
            break;

        case GREENS_CALC_NONE:
            os << "none";
            break;

        case GREENS_CALC_BARNES_HUT:
            os << "bh";
            break;

        default:
            os << "undefined";
            break;
    }

    return os;
}

std::ostream &operator<<(std::ostream &os, const FrictionLawMethod &calc_method) {
    switch (calc_method) {
        case FRIC_LAW_ORIG:
            os << "original";
            break;

        case FRIC_LAW_STEPPED:
            os << "stepped";
            break;

        default:
            os << "undefined";
            break;
    }

    return os;
}

std::ostream &operator<<(std::ostream &os, const VCParams &params) {
    os << "sim.version\t\t\t\t\t\t\t= " << params.version << "\n";

    os << "sim.start_year\t\t\t\t\t\t= " << params.year << "\n";
    os << "sim.end_year\t\t\t\t\t\t\t= " << params.sim_end_year << "\n";

    os << "sim.noise.event\t\t\t\t\t\t= " << params.noise_event << "\n";
    os << "sim.noise.slip_deficit\t\t\t\t= " << params.noise_slip_deficit << "\n";

    os << "sim.fault.kill_cff\t\t\t\t\t= " << params.fault_kill_cff << "\n";

    os << "sim.checkpoint.period\t\t\t\t\t\t= " << params.checkpoint_period << "\n";
    os << "sim.checkpoint.prefix\t\t\t\t\t\t= " << params.checkpoint_save_prefix << "\n";

    os << "sim.progress.period\t\t\t\t\t= " << params.progress_period << "\n";

    os << "sim.dynamic\t\t\t\t\t\t\t= " << params.dynamic << "\n";
    os << "sim.friction_law\t\t\t\t\t= " << params.friction_law_method << "\n";

    os << "sim.greens.kill_distance\t\t\t= " << params.greens_kill_distance << "\n";
    os << "sim.greens.method\t\t\t\t\t\t= " << params.greens_calc_method << "\n";
    os << "sim.greens.bh_theta\t\t\t\t= " << params.barnes_hut_theta << "\n";
    os << "sim.greens.infile\t\t\t\t\t= " << params.greens_infile << "\n";

    os << "sim.file.output_event\t\t\t= " << params.event_outfile << "\n";
    os << "sim.file.output_sweep\t\t\t= " << params.sweep_outfile << "\n";
    os << "sim.file.output_event_type\t\t= " << params.event_outfile_type << "\n";

    os << "sim.bass.max_generations\t\t\t= " << params.bass_max_generations << "\n";
    os << "sim.bass.min_magnitude_mm\t\t\t= " << params.bass_min_magnitude_mm << "\n";
    os << "sim.bass.aftershock_strength_dm\t= " << params.bass_aftershock_strength_dm << "\n";
    os << "sim.bass.frequency_scale_b\t\t\t= " << params.bass_frequency_scale_b << "\n";
    os << "sim.bass.aftershock_start_c\t\t\t= " << params.bass_aftershock_start_c << "\n";
    os << "sim.bass.time_decay_p\t\t\t\t= " << params.bass_time_decay_p << "\n";
    os << "sim.bass.distance_d\t\t\t\t\t= " << params.bass_distance_d << "\n";
    os << "sim.bass.distance_decay_q\t\t\t= " << params.bass_distance_decay_q << "\n";

    os << "sim.sanity_check\t\t\t\t\t\t= " << params.sanity_check << "\n";
    os << "sim.do_normal_stress\t\t\t\t\t= " << params.do_normal_stress << "\n";
    os << "sim.use_transpose_matrix\t\t\t\t= " << params.use_transpose_matrix << "\n";

    os << "sim.greens.outfile\t\t\t\t\t= " << params.greens_outfile << "\n";
    os << "sim.events.file\t\t\t\t\t\t= " << params.events_file << "\n";

    os << "sim.section_params.file\t\t\t\t\t\t= " << params.section_params_file << "\n";

    os << "sim.hdf5.file\t\t\t\t\t= " << params.hdf5_file << "\n";

    os << "sim.eqsim.file.condition\t\t\t= " << params.eqsim_condition_file << "\n";
    os << "sim.eqsim.file.friction\t\t\t\t= " << params.eqsim_friction_file << "\n";
    os << "sim.eqsim.file.geometry\t\t\t\t= " << params.eqsim_geometry_file << "\n";
    os << "sim.eqsim.file.output\t\t\t\t= " << params.eqsim_output_file << "\n";
    os << "sim.eqsim.file.output.slipmap_mag\t= " << params.eqsim_slipmap_mag << "\n";

    return os;
}
