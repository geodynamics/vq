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
    std::string     greens_method;
    ConfigFile param_file(param_file_name);

    valid = true;

    version = param_file.read<string>("sim.version", "");

    year = param_file.read<double>("sim.time.start_year", 0.0);
    sim_end_year = param_file.read<double>("sim.time.end_year", 1000.0);

    sanity_check = param_file.read<bool>("sim.system.sanity_check", false);
    use_transpose_matrix = param_file.read<bool>("sim.system.transpose_matrix", true);
    progress_period = param_file.read<unsigned int>("sim.system.progress_period", 0);
    checkpoint_period = param_file.read<int>("sim.system.checkpoint_period", 0);
    checkpoint_save_prefix = param_file.read<string>("sim.system.checkpoint_prefix", "sim_state_");

    dynamic = param_file.read<double>("sim.friction.dynamic", INFINITY);

    greens_method = param_file.read<string>("sim.greens.method", "standard");
    greens_infile = param_file.read<string>("sim.greens.input", "");
    greens_outfile = param_file.read<string>("sim.greens.output", "");
    do_normal_stress = param_file.read<bool>("sim.greens.use_normal", true);
    greens_kill_distance = param_file.read<double>("sim.greens.kill_distance", 0.0);
    greens_sample_distance = param_file.read<double>("sim.greens.sample_distance", 1000.0);

    if (greens_sample_distance <= 0) greens_sample_distance = 1000;

    input_model_file = param_file.read<string>("sim.file.input", "");
    input_model_file_type = param_file.read<string>("sim.file.input_type", "");
    event_outfile = param_file.read<string>("sim.file.output_event", "");
    sweep_outfile = param_file.read<string>("sim.file.output_sweep", "");
    event_outfile_type = param_file.read<string>("sim.file.output_event_type", "");

    bass_max_generations = param_file.read<unsigned int>("sim.bass.max_generations", 0);
    bass_min_magnitude_mm = param_file.read<double>("sim.bass.mm", 4.0);
    bass_aftershock_strength_dm = param_file.read<double>("sim.bass.dm", 1.25);
    bass_frequency_scale_b = param_file.read<double>("sim.bass.b", 1.0);
    bass_aftershock_start_c = param_file.read<double>("sim.bass.c", 0.1);
    bass_time_decay_p = param_file.read<double>("sim.bass.p", 1.25);
    bass_distance_d = param_file.read<double>("sim.bass.d", 300);
    bass_distance_decay_q = param_file.read<double>("sim.bass.q", 1.35);

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

std::ostream &operator<<(std::ostream &os, const VCParams &params) {
    os << "sim.version\t\t\t\t= " << params.version << "\n";

    os << "sim.time.start_year\t\t\t= " << params.year << "\n";
    os << "sim.time.end_year\t\t\t= " << params.sim_end_year << "\n";

    os << "sim.system.sanity_check\t\t\t= " << params.sanity_check << "\n";
    os << "sim.system.transpose_matrix\t\t= " << params.use_transpose_matrix << "\n";
    os << "sim.system.progress_period\t\t= " << params.progress_period << "\n";
    os << "sim.system.checkpoint_period\t\t= " << params.checkpoint_period << "\n";
    os << "sim.system.checkpoint_prefix\t\t= " << params.checkpoint_save_prefix << "\n";

    os << "sim.friction.dynamic\t\t\t= " << params.dynamic << "\n";

    os << "sim.greens.method\t\t\t= " << params.greens_calc_method << "\n";
    os << "sim.greens.input\t\t\t= " << params.greens_infile << "\n";
    os << "sim.greens.output\t\t\t= " << params.greens_outfile << "\n";
    os << "sim.greens.sample_distance\t\t= " << params.greens_sample_distance << "\n";
    os << "sim.greens.use_normal\t\t\t= " << params.do_normal_stress << "\n";
    os << "sim.greens.kill_distance\t\t= " << params.greens_kill_distance << "\n";

    os << "sim.file.input\t\t\t\t=" << params.input_model_file << "\n";
    os << "sim.file.input_type\t\t\t=" << params.input_model_file_type << "\n";
    os << "sim.file.output_event\t\t\t= " << params.event_outfile << "\n";
    os << "sim.file.output_sweep\t\t\t= " << params.sweep_outfile << "\n";
    os << "sim.file.output_event_type\t\t= " << params.event_outfile_type << "\n";

    os << "sim.bass.max_generations\t\t= " << params.bass_max_generations << "\n";
    os << "sim.bass.mm\t\t\t\t= " << params.bass_min_magnitude_mm << "\n";
    os << "sim.bass.dm\t\t\t\t= " << params.bass_aftershock_strength_dm << "\n";
    os << "sim.bass.b\t\t\t\t= " << params.bass_frequency_scale_b << "\n";
    os << "sim.bass.c\t\t\t\t= " << params.bass_aftershock_start_c << "\n";
    os << "sim.bass.p\t\t\t\t= " << params.bass_time_decay_p << "\n";
    os << "sim.bass.d\t\t\t\t= " << params.bass_distance_d << "\n";
    os << "sim.bass.q\t\t\t\t= " << params.bass_distance_decay_q << "\n";

    return os;
}
