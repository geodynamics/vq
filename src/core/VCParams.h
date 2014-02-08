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

#include "config.h"

#include <string>
#include "QuakeLibUtil.h"

#ifndef _VCPARAMS_H_
#define _VCPARAMS_H_

enum GreensCalcMethod {
    GREENS_CALC_UNDEFINED,      // undefined Greens function behavior
    GREENS_CALC_NONE,           // do not calculate Greens function
    GREENS_FILE_PARSE,          // use file parsing method
    GREENS_CALC_BARNES_HUT,     // use Barnes Hut method to calculate Greens function
    GREENS_CALC_2011            // use the new Okada class to calculate Greens functions
};

enum FrictionLawMethod {
    FRIC_LAW_UNDEFINED,         // undefined behavior
    FRIC_LAW_ORIG,              // Original friction law where blocks slip full amount each time
    FRIC_LAW_STEPPED            // Updated friction law where slip is proportional to number of ruptured blocks
};

/*!
 The set of possible parameters for a VC simulation.
 These are described in detail in the example/sample_params.d file.
 */
class VCParams {
    private:
        bool                valid;

        std::string         version;

        double              year;
        double              sim_end_year;
        double              event_start_year;

        double              noise_event;
        double              noise_slip_deficit;

        double              fault_kill_cff;

        int                 checkpoint_period;  // in terms of # of events between state saves
        std::string         checkpoint_save_prefix;

        unsigned int        progress_period;

        double              dynamic;
        FrictionLawMethod   friction_law_method;
        unsigned int        slip_scaling_threshold;

        double              greens_kill_distance;
        GreensCalcMethod    greens_calc_method;
        double              barnes_hut_theta;       // controls how much smoothing occurs in Barnes-Hutt approximation
        std::string         greens_infile;

        unsigned int        bass_max_generations;
        double              bass_min_magnitude_mm;
        double              bass_aftershock_strength_dm;
        double              bass_frequency_scale_b;
        double              bass_aftershock_start_c;
        double              bass_time_decay_p;
        double              bass_distance_d;
        double              bass_distance_decay_q;

        bool                sanity_check;
        bool                do_normal_stress;
        bool                use_transpose_matrix;

        double              model_lat0;
        double              model_lon0;

        std::string         greens_outfile;
        std::string         events_file;
        std::string         hdf5_file;
        std::string         section_params_file;

        std::string         eqsim_condition_file;
        std::string         eqsim_friction_file;
        std::string         eqsim_geometry_file;
        std::string         eqsim_output_file;
        double              eqsim_slipmap_mag;

        std::string         event_outfile, sweep_outfile;
        std::string         event_outfile_type;

    public:
        VCParams(void) : valid(false) {};
        void read_params(const std::string &param_file_name);

        std::string getVersion(void) const {
            return version;
        };

        double getSimStart(void) const {
            return year;
        };
        double getSimDuration(void) const {
            return sim_end_year;
        };
        double getRecordingStart(void) const {
            return event_start_year;
        };

        double getEventNoise(void) const {
            return noise_event;
        };
        double getSlipDeficitNoise(void) const {
            return noise_slip_deficit;
        };

        double getFaultKillCFF(void) const {
            return fault_kill_cff;
        };

        int getCheckpointPeriod(void) const {
            return checkpoint_period;
        };
        std::string getCheckpointPrefix(void) const {
            return checkpoint_save_prefix;
        };

        unsigned int getProgressPeriod(void) const {
            return progress_period;
        };

        double getDynamic(void) const {
            return dynamic;
        };
        FrictionLawMethod getFrictionLaw(void) const {
            return friction_law_method;
        };
        unsigned int getSlipScalingThreshold(void) const {
            return slip_scaling_threshold;
        };

        double getGreensKillDistance(void) const {
            return greens_kill_distance;
        };
        GreensCalcMethod getGreensCalcMethod(void) const {
            return greens_calc_method;
        };
        double getBarnesHutTheta(void) const {
            return barnes_hut_theta;
        };
        std::string getGreensInputfile(void) const {
            return greens_infile;
        };

        unsigned int getBASSMaxGenerations(void) const {
            return bass_max_generations;
        };
        double getBASSMinMagnitude(void) const {
            return bass_min_magnitude_mm;
        };
        double getBASSAftershockStrength(void) const {
            return bass_aftershock_strength_dm;
        };
        double getBASSFrequencyScale(void) const {
            return bass_frequency_scale_b;
        };
        double getBASSAftershockStart(void) const {
            return bass_aftershock_start_c;
        };
        double getBASSTimeDecay(void) const {
            return bass_time_decay_p;
        };
        double getBASSDistance(void) const {
            return bass_distance_d;
        };
        double getBASSDistanceDecay(void) const {
            return bass_distance_decay_q;
        };

        bool doSanityCheck(void) const {
            return sanity_check;
        };
        bool doNormalStress(void) const {
            return do_normal_stress;
        };
        bool useTransposedMatrix(void) const {
            return use_transpose_matrix;
        };

        quakelib::LatLonDepth getBaseLatLon(void) const {
            return quakelib::LatLonDepth(model_lat0, model_lon0, 0);
        };
        std::string getGreensOutfile(void) const {
            return greens_outfile;
        };
        std::string getEventsFile(void) const {
            return events_file;
        };
        std::string getHDF5File(void) const {
            return hdf5_file;
        };
        std::string getSectionParamsFile(void) const {
            return section_params_file;
        };

        std::string getEqSimConditionFile(void) const {
            return eqsim_condition_file;
        };
        std::string getEqSimFrictionFile(void) const {
            return eqsim_friction_file;
        };
        std::string getEqSimGeometryFile(void) const {
            return eqsim_geometry_file;
        };
        std::string getEqSimOutputFile(void) const {
            return eqsim_output_file;
        };
        double getEqSimOutputSlipMapMag(void) const {
            return eqsim_slipmap_mag;
        };

        std::string getEventOutfile(void) const {
            return event_outfile;
        };
        std::string getSweepOutfile(void) const {
            return sweep_outfile;
        };
        std::string getEventOutfileType(void) const {
            return event_outfile_type;
        };

        friend std::ostream &operator<<(std::ostream &os, const VCParams &params);
};

std::ostream &operator<<(std::ostream &os, const VCParams &params);
std::ostream &operator<<(std::ostream &os, const GreensCalcMethod &calc_method);
std::ostream &operator<<(std::ostream &os, const FrictionLawMethod &calc_method);

#endif
