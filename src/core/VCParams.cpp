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

/*!
 Parse a configuration file for each of the possible parameters, assigning default values
 for parameters that are not explicitly specified.
 */
void VCParams::read_params(const std::string &param_file_name) {
    params = ConfigFile(param_file_name);

    // Call each of the parameter checks to set the default value if it's not already set
    params.readSet<string>("sim.version", "");

    params.readSet<double>("sim.time.start_year", 0.0);
    params.readSet<double>("sim.time.end_year", 1000.0);

    // in terms of # of events between state saves
    params.readSet<int>("sim.system.checkpoint_period", 0);
    params.readSet<string>("sim.system.checkpoint_prefix", "sim_state_");

    params.readSet<unsigned int>("sim.system.progress_period", 0);

    params.readSet<double>("sim.friction.dynamic", INFINITY);

    params.readSet<double>("sim.greens.kill_distance", 0.0);
    double dist = params.readSet<double>("sim.greens.sample_distance", 1000.0);

    if (dist <= 0) params.add("sim.greens.sample_distance", "1000.0");

    std::string greens_method = params.readSet<string>("sim.greens.method", "standard");

    // Parse the Greens calculation method string
    if (greens_method.compare("file") && greens_method.compare("bh") && greens_method.compare("standard")) {
        params.add("sim.greens.method", "undefined");
    }

    // controls how much smoothing occurs in Barnes-Hut approximation
    params.readSet<double>("sim.greens.bh_theta", 0.0);
    params.readSet<string>("sim.greens.input", "");

    params.readSet<unsigned int>("sim.bass.max_generations", 0);
    params.readSet<double>("sim.bass.mm", 4.0);
    params.readSet<double>("sim.bass.dm", 1.25);
    params.readSet<double>("sim.bass.b", 1.0);
    params.readSet<double>("sim.bass.c", 0.1);
    params.readSet<double>("sim.bass.p", 1.25);
    params.readSet<double>("sim.bass.d", 300);
    params.readSet<double>("sim.bass.q", 1.35);

    params.readSet<bool>("sim.system.sanity_check", false);
    params.readSet<bool>("sim.greens.use_normal", true);
    params.readSet<bool>("sim.system.transpose_matrix", true);

    params.readSet<string>("sim.file.input", "");
    params.readSet<string>("sim.file.input_type", "");

    params.readSet<string>("sim.greens.output", "");

    params.readSet<string>("sim.file.output_event", "");
    params.readSet<string>("sim.file.output_sweep", "");
    params.readSet<string>("sim.file.output_event_type", "");

}

void VCParams::write_params(const std::string &param_file_name) {
    std::ofstream   param_out_file(param_file_name.c_str());
    param_out_file << params;

    param_out_file.close();
}
