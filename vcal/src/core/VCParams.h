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

//#ifdef HAVE_CONFIG_H
#include "config.h"
//#endif

#include <string>
#include "QuakeLibUtil.h"

#ifndef _VCPARAMS_H_
#define _VCPARAMS_H_

enum GreensCalcMethod {
	GREENS_CALC_UNDEFINED,		// undefined Greens function behavior
	GREENS_CALC_NONE,			// do not calculate Greens function
	GREENS_FILE_PARSE,			// use file parsing method
	GREENS_CALC_BARNES_HUT,		// use Barnes Hut method to calculate Greens function
    GREENS_CALC_2011            // use the new Okada class to calculate Greens functions
};

enum SpecExecMethod {
	SPEC_EXEC_UNDEFINED,		// undefined behavior
	SPEC_EXEC_NONE,				// do not use speculative execution
	SPEC_EXEC_FIXED_DIST,		// use speculative execution with fixed boundary distance
	SPEC_EXEC_ADAPTIVE			// use speculative execution with adaptive prediction
};

enum FrictionLawMethod {
	FRIC_LAW_UNDEFINED,			// undefined behavior
	FRIC_LAW_ORIG,				// Original friction law where blocks slip full amount each time
	FRIC_LAW_STEPPED			// Updated friction law where slip is proportional to number of ruptured blocks
};

/*!
 The set of possible parameters for a VC simulation.
 These are described in detail in the example/sample_params.d file.
 */
class VCParams {
private:
	bool				valid;
	
	std::string			version;
	
	SpecExecMethod		spec_exec_method;
	double				spec_exec_dist;
	
	double				year;
	double				sim_end_year;
	double				event_start_year;
	
	double				noise_event;
	double				noise_slip_deficit;
	double				noise_stress;
	double				noise_stress_resolution;
	
	double				fault_kill_cff;
	
	int					checkpoint_period;	// in terms of # of events between state saves
	std::string			checkpoint_save_prefix;
	
	unsigned int		progress_period;
	
	double				dynamic;
	FrictionLawMethod	friction_law_method;
	unsigned int		slip_scaling_threshold;
	
	double				greens_kill_distance;
	GreensCalcMethod	greens_calc_method;
	double				barnes_hut_theta;		// controls how much smoothing occurs in Barnes-Hutt approximation
	std::string			greens_infile;
	
	int					asperity_num;
	int					asperity_width;
	double				asperity_mag;
	
	unsigned int		bass_max_generations;
	double				bass_min_magnitude_mm;
	double				bass_aftershock_strength_dm;
	double				bass_frequency_scale_b;
	double				bass_aftershock_start_c;
	double				bass_time_decay_p;
	double				bass_distance_d;
	double				bass_distance_decay_q;
	
	double				bg_event_mean_interevent;
	double				bg_event_min_mag;
	double				bg_event_max_mag;
	double				bg_event_distance;
	double				bg_event_distance_decay;
	
	bool				depth_dependent_velocity;
	bool				sanity_check;
	bool				do_normal_stress;
	bool				use_transpose_matrix;
	
	double				model_lat0;
	double				model_lon0;
	
	std::string			state_begin_file;
	std::string			system_outfile;
	std::string			greens_outfile;
	std::string			events_file;
	std::string			hdf5_file;
	std::string			section_params_file;
    
	std::string			eqsim_condition_file;
	std::string			eqsim_friction_file;
	std::string			eqsim_geometry_file;
	std::string			eqsim_output_file;
	double				eqsim_slipmap_mag;
	
public:
	VCParams(void) : valid(false) {};
	void read_params(const std::string &param_file_name);
	
	std::string getVersion(void) const { return version; };
	
	SpecExecMethod getSpecExecMethod(void) const { return spec_exec_method; };
	double getSpecExecDistance(void) const { return spec_exec_dist; };
	
	double getSimStart(void) const { return year; };
	double getSimDuration(void) const { return sim_end_year; };
	double getRecordingStart(void) const { return event_start_year; };
	
	double getEventNoise(void) const { return noise_event; };
	double getSlipDeficitNoise(void) const { return noise_slip_deficit; };
	double getStressNoise(void) const { return noise_stress; };
	double getStressNoiseResolution(void) const { return noise_stress_resolution; };
	
	double getFaultKillCFF(void) const { return fault_kill_cff; };
	
	int getCheckpointPeriod(void) const { return checkpoint_period; };
	std::string getCheckpointPrefix(void) const { return checkpoint_save_prefix; };
	
	unsigned int getProgressPeriod(void) const { return progress_period; };
	
	double getDynamic(void) const { return dynamic; };
	FrictionLawMethod getFrictionLaw(void) const { return friction_law_method; };
	unsigned int getSlipScalingThreshold(void) const { return slip_scaling_threshold; };
	
	double getGreensKillDistance(void) const { return greens_kill_distance; };
	GreensCalcMethod getGreensCalcMethod(void) const { return greens_calc_method; };
	double getBarnesHutTheta(void) const { return barnes_hut_theta; };
	std::string getGreensInputfile(void) const { return greens_infile; };
	
	int getNumAsperities(void) const { return asperity_num; };
	int getAsperityWidth(void) const { return asperity_width; };
	double getAsperityMagnitude(void) const { return asperity_mag; };
	
	unsigned int getBASSMaxGenerations(void) const { return bass_max_generations; };
	double getBASSMinMagnitude(void) const { return bass_min_magnitude_mm; };
	double getBASSAftershockStrength(void) const { return bass_aftershock_strength_dm; };
	double getBASSFrequencyScale(void) const { return bass_frequency_scale_b; };
	double getBASSAftershockStart(void) const { return bass_aftershock_start_c; };
	double getBASSTimeDecay(void) const { return bass_time_decay_p; };
	double getBASSDistance(void) const { return bass_distance_d; };
	double getBASSDistanceDecay(void) const { return bass_distance_decay_q; };
	
	double getBGEventMeanInterevent(void) const { return bg_event_mean_interevent; };
	double getBGEventMinMagnitude(void) const { return bg_event_min_mag; };
	double getBGEventMaxMagnitude(void) const { return bg_event_max_mag; };
	double getBGEventDistance(void) const { return bg_event_distance; };
	double getBGEventDistanceDecay(void) const { return bg_event_distance_decay; };
	
	bool doDepthDependentVelocity(void) const { return depth_dependent_velocity; };
	bool doSanityCheck(void) const { return sanity_check; };
	bool doNormalStress(void) const { return do_normal_stress; };
	bool useTransposedMatrix(void) const { return use_transpose_matrix; };
	
	quakelib::LatLonDepth getBaseLatLon(void) const { return quakelib::LatLonDepth(model_lat0, model_lon0, 0); };
	
	std::string getStateBeginFile(void) const { return state_begin_file; };
	std::string getSystemOutfile(void) const { return system_outfile; };
	std::string getGreensOutfile(void) const { return greens_outfile; };
	std::string getEventsFile(void) const { return events_file; };
	std::string getHDF5File(void) const { return hdf5_file; };
    std::string getSectionParamsFile(void) const { return section_params_file; };
	
	std::string getEqSimConditionFile(void) const { return eqsim_condition_file; };
	std::string getEqSimFrictionFile(void) const { return eqsim_friction_file; };
	std::string getEqSimGeometryFile(void) const { return eqsim_geometry_file; };
	std::string getEqSimOutputFile(void) const { return eqsim_output_file; };
	double getEqSimOutputSlipMapMag(void) const { return eqsim_slipmap_mag; };
	
    friend std::ostream& operator<<(std::ostream& os, const VCParams& params);
};

std::ostream& operator<<(std::ostream& os, const VCParams& params);
std::ostream& operator<<(std::ostream& os, const GreensCalcMethod& calc_method);
std::ostream& operator<<(std::ostream& os, const FrictionLawMethod& calc_method);

#endif
