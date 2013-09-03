// Copyright (c) 2010, John B. Rundle <rundle@cse.ucdavis.edu>, 
// All rights reserved.
// 
// Redistribution and use of this code or any derivative works are
// permitted provided that the following conditions are met:
// 
// * Redistributions may not be sold, nor may they be used in a
// commercial product or activity.
// 
// * Redistributions that are modified from the original source must
// include the complete source code, including the source code for all
// components used by a binary built from the modified
// sources. However, as a special exception, the source code
// distributed need not include anything that is normally distributed
// (in either source or binary form) with the major components
// (compiler, kernel, and so on) of the operating system on which the
// executable runs, unless that component itself accompanies the
// executable.
// 
// * Redistributions must reproduce the above copyright notice, this list
// of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "VCParams.h"
#include "ConfigFile.h"

/*!
 Parse a configuration file for each of the possible parameters, assigning default values
 for parameters that are not explicitly specified.
 */
void VCParams::read_params(const std::string &param_file_name) {
	std::string		greens_method, spec_exec_str, friction_law_str;
	ConfigFile param_file(param_file_name);
	
	valid = true;
	
	version = param_file.read<string>("sim.version", "");
	
	spec_exec_str = param_file.read<string>("sim.spec_exec", "none");
	spec_exec_dist = param_file.read<double>("sim.spec_exec_distance", 0.0);
	
	year = param_file.read<double>("sim.start_year", 0.0);
	sim_end_year = param_file.read<double>("sim.end_year", 1000.0);
	event_start_year = param_file.read<double>("sim.events.start_year", 0);
	
	noise_event = param_file.read<double>("sim.noise.event", 0.0);
	noise_slip_deficit = param_file.read<double>("sim.noise.slip_deficit", 0.0);
	noise_stress = param_file.read<double>("sim.noise.stress", 0.0);
	noise_stress_resolution = param_file.read<double>("sim.noise.stress.resolution", 10.0);
	
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
	
	asperity_num = param_file.read<int>("sim.asperity.num", 0);
	asperity_width = param_file.read<int>("sim.asperity.width", 1);
	asperity_mag = param_file.read<double>("sim.asperity.mag", 2.0);
	
	bass_max_generations = param_file.read<unsigned int>("sim.bass.max_generations", 0);
	bass_min_magnitude_mm = param_file.read<double>("sim.bass.min_magnitude_mm", 4.0);
	bass_aftershock_strength_dm = param_file.read<double>("sim.bass.aftershock_strength_dm", 1.25);
	bass_frequency_scale_b = param_file.read<double>("sim.bass.frequency_scale_b", 1.0);
	bass_aftershock_start_c = param_file.read<double>("sim.bass.aftershock_start_c", 0.1);
	bass_time_decay_p = param_file.read<double>("sim.bass.time_decay_p", 1.25);
	bass_distance_d = param_file.read<double>("sim.bass.distance_d", 300);
	bass_distance_decay_q = param_file.read<double>("sim.bass.distance_decay_q", 1.35);
	
	bg_event_mean_interevent = param_file.read<double>("sim.bg_event.mean_intervent", -1);
	bg_event_min_mag = param_file.read<double>("sim.bg_event.min_magnitude", 4.0);
	bg_event_max_mag = param_file.read<double>("sim.bg_event.max_magnitude", 7.0);
	bg_event_distance = param_file.read<double>("sim.bg_event.distance", 300);
	bg_event_distance_decay = param_file.read<double>("sim.bg_event.distance_decay", 1.35);
	
	depth_dependent_velocity = param_file.read<bool>("sim.depth_dependent_velocity", false);
	sanity_check = param_file.read<bool>("sim.sanity_check", false);
	do_normal_stress = param_file.read<bool>("sim.do_normal_stress", true);
	use_transpose_matrix = param_file.read<bool>("sim.use_transpose_matrix", true);
	
	model_lat0 = param_file.read<double>("sim.system.rawfile.lat0", 31.5);
	model_lon0 = param_file.read<double>("sim.system.rawfile.lon0", -126.0);
	
	state_begin_file = param_file.read<string>("sim.state.file.begin", "");
	system_outfile = param_file.read<string>("sim.system.outfile", "");
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
    } else if (!greens_method.compare("2011")) {
		greens_calc_method = GREENS_CALC_2011;
	} else {
		greens_calc_method = GREENS_CALC_UNDEFINED;
	}
	
	// Parse the speculative execution method string
	if (!spec_exec_str.compare("none")) {
		spec_exec_method = SPEC_EXEC_NONE;
	} else if (!spec_exec_str.compare("fixed")) {
		spec_exec_method = SPEC_EXEC_FIXED_DIST;
	} else if (!spec_exec_str.compare("adaptive")) {
		spec_exec_method = SPEC_EXEC_ADAPTIVE;
	} else {
		spec_exec_method = SPEC_EXEC_UNDEFINED;
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

std::ostream& operator<<(std::ostream& os, const GreensCalcMethod& calc_method) {
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

std::ostream& operator<<(std::ostream& os, const SpecExecMethod& calc_method) {
	switch (calc_method) {
		case SPEC_EXEC_NONE:
			os << "none";
			break;
		case SPEC_EXEC_FIXED_DIST:
			os << "fixed";
			break;
		case SPEC_EXEC_ADAPTIVE:
			os << "adaptive";
			break;
		default:
			os << "undefined";
			break;
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, const FrictionLawMethod& calc_method) {
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

std::ostream& operator<<(std::ostream& os, const VCParams& params) {
	os << "sim.version\t\t\t\t\t\t\t= " << params.version << "\n";
	
	os << "sim.spec_exec\t\t\t\t\t\t= " << params.spec_exec_method << "\n";
	os << "sim.spec_exec_distance\t\t\t\t= " << params.spec_exec_dist << "\n";
	
	os << "sim.start_year\t\t\t\t\t\t= " << params.year << "\n";
	os << "sim.end_year\t\t\t\t\t\t\t= " << params.sim_end_year << "\n";
	os << "sim.events.start_year\t\t\t\t= " << params.event_start_year << "\n";
	
	os << "sim.noise.event\t\t\t\t\t\t= " << params.noise_event << "\n";
	os << "sim.noise.slip_deficit\t\t\t\t= " << params.noise_slip_deficit << "\n";
	os << "sim.noise.stress\t\t\t\t\t\t= " << params.noise_stress << "\n";
	os << "sim.noise.stress.resolution\t\t\t= " << params.noise_stress_resolution << "\n";
	
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
	
	os << "sim.asperity.num\t\t\t\t\t\t= " << params.asperity_num << "\n";
	os << "sim.asperity.width\t\t\t\t\t= " << params.asperity_width << "\n";
	os << "sim.asperity.mag\t\t\t\t\t\t= " << params.asperity_mag << "\n";
	
	os << "sim.bass.max_generations\t\t\t= " << params.bass_max_generations << "\n";
	os << "sim.bass.min_magnitude_mm\t\t\t= " << params.bass_min_magnitude_mm << "\n";
	os << "sim.bass.aftershock_strength_dm\t= " << params.bass_aftershock_strength_dm << "\n";
	os << "sim.bass.frequency_scale_b\t\t\t= " << params.bass_frequency_scale_b << "\n";
	os << "sim.bass.aftershock_start_c\t\t\t= " << params.bass_aftershock_start_c << "\n";
	os << "sim.bass.time_decay_p\t\t\t\t= " << params.bass_time_decay_p << "\n";
	os << "sim.bass.distance_d\t\t\t\t\t= " << params.bass_distance_d << "\n";
	os << "sim.bass.distance_decay_q\t\t\t= " << params.bass_distance_decay_q << "\n";
	
	os << "sim.bg_event.mean_intervent\t\t\t= " << params.bg_event_mean_interevent << "\n";
	os << "sim.bg_event.min_magnitude\t\t\t= " << params.bg_event_min_mag << "\n";
	os << "sim.bg_event.max_magnitude\t\t\t= " << params.bg_event_max_mag << "\n";
	os << "sim.bg_event.distance\t\t\t\t= " << params.bg_event_distance << "\n";
	os << "sim.bg_event.distance_decay\t\t\t= " << params.bg_event_distance_decay << "\n";
	
	os << "sim.depth_dependent_velocity\t\t= " << params.depth_dependent_velocity << "\n";
	os << "sim.sanity_check\t\t\t\t\t\t= " << params.sanity_check << "\n";
	os << "sim.do_normal_stress\t\t\t\t\t= " << params.do_normal_stress << "\n";
	os << "sim.use_transpose_matrix\t\t\t\t= " << params.use_transpose_matrix << "\n";
	
	os << "sim.system.rawfile.lat0\t\t\t\t= " << params.model_lat0 << "\n";
	os << "sim.system.rawfile.lon0\t\t\t\t= " << params.model_lon0 << "\n";
	
	os << "sim.state.file.begin\t\t\t\t\t= " << params.state_begin_file << "\n";
	os << "sim.system.outfile\t\t\t\t\t= " << params.system_outfile << "\n";
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
