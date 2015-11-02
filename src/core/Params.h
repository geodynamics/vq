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
#include "ConfigFile.h"

#ifndef _PARAMS_H_
#define _PARAMS_H_

enum GreensCalcMethod {
    GREENS_CALC_UNDEFINED,      // undefined Greens function behavior
    GREENS_CALC_NONE,           // do not calculate Greens function
    GREENS_FILE_PARSE,          // use file parsing method
    GREENS_CALC_BARNES_HUT,     // use Barnes Hut method to calculate Greens function
    GREENS_CALC_STANDARD        // use the new Okada class to calculate Greens functions
};

/*!
 The set of possible parameters for a VC simulation.
 These are described in detail in the example/sample_params.d file.
 */
class VCParams {
    private:
        ConfigFile          params;

    public:
        VCParams(void) {};
        void read_params(const std::string &param_file_name);
        void write_params(const std::string &param_file_name);

        std::string getVersion(void) const {
            return params.read<string>("sim.version");
        };

        double getSimStart(void) const {
            return params.read<double>("sim.time.start_year");
        };
        double getSimDuration(void) const {
            return params.read<double>("sim.time.end_year");
        };

        // in terms of # of events between state saves
        int getCheckpointPeriod(void) const {
            return params.read<int>("sim.system.checkpoint_period");
        };
        std::string getCheckpointPrefix(void) const {
            return params.read<string>("sim.system.checkpoint_prefix");
        };

        unsigned int getProgressPeriod(void) const {
            return params.read<unsigned int>("sim.system.progress_period");
        };

        double getDynamic(void) const {
            return params.read<double>("sim.friction.dynamic");
        };

        double getGreensKillDistance(void) const {
            return params.read<double>("sim.greens.kill_distance");
        };
        double getGreensSampleDistance(void) const {
            return params.read<double>("sim.greens.sample_distance");
        };
        GreensCalcMethod getGreensCalcMethod(void) const {
            std::string greens_method = params.read<string>("sim.greens.method");

            // Parse the Greens calculation method string
            if (!greens_method.compare("file")) return GREENS_FILE_PARSE;

            if (!greens_method.compare("bh")) return GREENS_CALC_BARNES_HUT;

            if (!greens_method.compare("standard")) return GREENS_CALC_STANDARD;

            return GREENS_CALC_UNDEFINED;
        };
        // controls how much smoothing occurs in Barnes-Hut approximation
        double getBarnesHutTheta(void) const {
            return params.read<double>("sim.greens.bh_theta");
        };
        std::string getGreensInputfile(void) const {
            return params.read<string>("sim.greens.input");
        };

        unsigned int getBASSMaxGenerations(void) const {
            return params.read<unsigned int>("sim.bass.max_generations");
        };
        double getBASSMinMagnitude(void) const {
            return params.read<double>("sim.bass.mm");
        };
        double getBASSAftershockStrength(void) const {
            return params.read<double>("sim.bass.dm");
        };
        double getBASSFrequencyScale(void) const {
            return params.read<double>("sim.bass.b");
        };
        double getBASSAftershockStart(void) const {
            return params.read<double>("sim.bass.c");
        };
        double getBASSTimeDecay(void) const {
            return params.read<double>("sim.bass.p");
        };
        double getBASSDistance(void) const {
            return params.read<double>("sim.bass.d");
        };
        double getBASSDistanceDecay(void) const {
            return params.read<double>("sim.bass.q");
        };

        bool doSanityCheck(void) const {
            return params.read<bool>("sim.system.sanity_check");
        };
        bool doNormalStress(void) const {
            return params.read<bool>("sim.greens.use_normal");
        };
        bool useTransposedMatrix(void) const {
            return params.read<bool>("sim.system.transpose_matrix");
        };

        std::string getModelFile(void) const {
            return params.read<string>("sim.file.input");
        };
        std::string getModelFileType(void) const {
            return params.read<string>("sim.file.input_type");
        };

        std::string getGreensOutfile(void) const {
            return params.read<string>("sim.greens.output");
        };

        std::string getEventOutfile(void) const {
            return params.read<string>("sim.file.output_event");
        };
        std::string getSweepOutfile(void) const {
            return params.read<string>("sim.file.output_sweep");
        };
        std::string getEventOutfileType(void) const {
            return params.read<string>("sim.file.output_event_type");
        };

        std::string getStressOutfile(void) const {
            return params.read<string>("sim.file.output_stress");
        };
        std::string getStressIndexOutfile(void) const {
            return params.read<string>("sim.file.output_stress_index");
        };
        std::string getStressOutfileType(void) const {
            return params.read<string>("sim.file.output_stress_type");
        };

        int getStressOutInterval(void) const {
            return params.read<unsigned int>("sim.file.output_stress_num_events");
        };

        std::string getStressInfile(void) const {
            return params.read<string>("sim.file.input_stress");
        };
        std::string getStressIndexInfile(void) const {
            return params.read<string>("sim.file.input_stress_index");
        };
        std::string getStressInfileType(void) const {
            return params.read<string>("sim.file.input_stress_type");
        };

        //
        // yoder: helper functions for greens function max/min values:
        // can we make a (set of) generic pramGetter(str.pram_name) function(s)? then, we can add params at will...
        // retiring the universal max/min values. BUT, leave these in place and default to the off-diagonal elements.
        double getGreenShearMax(void) const {
            //return params.read<double>("sim.greens.shear_max");
            return params.read<double>("sim.greens.shear_offdiag_max");
        };
        double getGreenShearMin(void) const {
            //return params.read<double>("sim.greens.shear_min");
            return params.read<double>("sim.greens.shear_offdiag_min");
        };
        double getGreenNormalMax(void) const {
            //return params.read<double>("sim.greens.normal_max");
            return params.read<double>("sim.greens.normal_offdiag_max");
        };
        double getGreenNormalMin(void) const {
            //return params.read<double>("sim.greens.normal_min");
            return params.read<double>("sim.greens.normal_offdiag_min");
        };

        // Kasey: new parameter for loading specified stress drops (False) or
        // computing them from Green's functions (True)
        bool computeStressDrops(void) const {
            return params.read<bool>("sim.friction.compute_stress_drops");
        };

        // Kasey: new parameter. Compute dynamic stress drops (True) or
        // computing from standard method (False)
        bool doDynamicStressDrops(void) const {
            return params.read<bool>("sim.friction.dynamic_stress_drops");
        };
        
        // Schultz: Prescribe a constant coefficient of friction for all blocks
        double getFrictionCoefficient(void) const {
            return params.read<double>("sim.friction.coefficient");
        };

        // The constant adjustment to the stress drops, should be between 0.3 and 0.5 globally
        //     with tuning the larger faults may need 0.6 to 0.8. The larger the number, the larger
        //     the stress drops. Currently it's a global constant
        double stressDropFactor(void) const {
            return params.read<double>("sim.friction.stress_drop_factor");
        };

        //
        // yoder: and overload so we can easily distinguish (off)diagonal elements (defined in GreensInit.cpp)
        double getGreenShearMax(const int block_1, const int block_2) const;
        double getGreenShearMin(const int block_1, const int block_2) const;
        double getGreenNormalMax(const int block_1, const int block_2) const;
        double getGreenNormalMin(const int block_1, const int block_2) const;

        //
        // yoder: ... and for (non)diag specific...
        double getGreenShearDiagMax(void) const {
            return params.read<double>("sim.greens.shear_diag_max");
        };
        double getGreenShearDiagMin(void) const {
            return params.read<double>("sim.greens.shear_diag_min");
        };
        double getGreenNormalDiagMax(void) const {
            return params.read<double>("sim.greens.normal_diag_max");
        };
        double getGreenNormalDiagMin(void) const {
            return params.read<double>("sim.greens.normal_diag_min");
        };
        //
        double getGreenShearOffDiagMax(void) const {
            return params.read<double>("sim.greens.shear_offdiag_max");
        };
        double getGreenShearOffDiagMin(void) const {
            return params.read<double>("sim.greens.shear_offdiag_min");
        };
        double getGreenNormalOffDiagMax(void) const {
            return params.read<double>("sim.greens.normal_offdiag_max");
        };
        double getGreenNormalOffDiagMin(void) const {
            return params.read<double>("sim.greens.normal_offdiag_min");
        };

};

#endif
