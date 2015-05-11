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

#include "SimFramework.h"
#include "Params.h"
#include "Block.h"
#include "Event.h"
#include "SimData.h"
#include "Comm.h"
#include "CommPartition.h"
#include "HDF5Data.h"

#ifdef VQ_HAVE_LIMITS_H
#include <limits.h>
#endif

#ifdef VQ_HAVE_FLOAT_H
#include <float.h>
#endif

#ifndef _Simulation_H_
#define _Simulation_H_

#include <string>

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iomanip>

enum PartitionMethod {
    PARTITION_UNDEFINED,
    PARTITION_BLOCK,
    PARTITION_DISTANCE
};

/*!
 Simulation is an instantiation of a Virtual California simulation using the SimFramework.
 It contains functions related to the simulation, parameter retrieval functions and basic
 functions for earthquake analysis.
 */
class Simulation : public SimFramework, public VCParams, public VCSimData, public VCCommPartition, public VCComm {
    public:
        Simulation(int argc, char **argv);
        ~Simulation(void);

        void init(void);

        double getYear(void) const {
            return year;
        };
        void setYear(const double &new_year) {
            year = new_year;
        };
        void incrementYear(const double &year_step) {
            year += year_step;
        };

        int numFaults(void) const;
        void getInitialFinalStresses(const quakelib::ElementIDSet &block_set, double &shear_init, double &shear_final, double &normal_init, double &normal_final) const;

        void sumStresses(const quakelib::ElementIDSet &block_set, double &shear_stress_sum, double &shear_stress0_sum, double &normal_stress_sum, double &normal_stress0_sum) const;

        double sumGreenShear(const BlockID &r) const {
            double sum = 0;

            for (int i=0; i<numGlobalBlocks(); ++i) sum += greenShear()->val(i, r);

            return sum;
        };
        double getGreenShear(const BlockID &r, const BlockID &c) const {
            return greenShear()->val(getLocalInd(r), c);
        };
        // yoder: move content to Simulation.cpp (enforcing min/max values for greens values).
        void setGreens(const BlockID &r, const BlockID &c, const double &new_green_shear, const double &new_green_normal);
        // yoder:
        void debug_out(std::string str_in);

        double getGreenNormal(const BlockID &r, const BlockID &c) const {
            return greenNormal()->val(getLocalInd(r), c);
        };
        bool compressShearRow(const unsigned int &row, const float &ratio) {
            return greenShear()->compressRow(row, ratio);
        };
        bool compressNormalRow(const unsigned int &row, const float &ratio) {
            return greenNormal()->compressRow(row, ratio);
        };
        bool decompressShearRow(const unsigned int &row) {
            return greenShear()->decompressRow(row);
        };
        bool decompressNormalRow(const unsigned int &row) {
            return greenNormal()->decompressRow(row);
        };
        void allocateShearNormalRows(const unsigned int &row) {
            greenNormal()->allocateRow(row);
            greenShear()->allocateRow(row);
        };

        void determineBlockNeighbors(void);
        void computeCFFs(void);
        void calcCFF(const BlockID gid);
        void matrixVectorMultiplyAccum(double *c, const quakelib::DenseMatrix<GREEN_VAL> *a, const double *b, const bool dense);
        void multiplySumRow(double *c, const double *b, const GREEN_VAL *a, const int n, const bool dense);
        void multiplyRow(double *c, const double *b, const GREEN_VAL *a, const int n);
        void distributeUpdateField(void);
        void broadcastUpdateField(void);
        void distributeBlocks(const quakelib::ElementIDSet &local_id_list, BlockIDProcMapping &global_id_list);
        void collectEventSweep(quakelib::ModelSweeps &sweeps);

        std::pair<quakelib::ElementIDSet::const_iterator, quakelib::ElementIDSet::const_iterator> getNeighbors(const BlockID &bid) const;
        void printTimers(void);

        bool isLocalBlockID(const BlockID &block_id) const {
            return (block_node_map.at(block_id) == node_rank);
        };

        double getSlipDeficit(const BlockID gid) {
            return slip_deficit[gid];
        };
        void setSlipDeficit(const BlockID gid, const double new_slip_deficit) {
            slip_deficit[gid] = new_slip_deficit;
        };

        double getCFF(const BlockID gid) {
            return cff[gid];
        };
        void setCFF(const BlockID gid, const double new_cff) {
            cff[gid] = new_cff;
        };
        double getCFF0(const BlockID gid) {
            return cff0[gid];
        };

        //! Calculates the static friction of this block.
        void calcFriction(const BlockID gid) {
            friction[gid] = (fabs(stress_drop[gid])/rhogd[gid]);
        }
        //! Get the static friction of this block.
        double getFriction(const BlockID gid) {
            return friction[gid];
        }

        //! Whether the block experienced static friction failure.
        //! This occurs if the Coulomb failure function goes over 0.
        bool cffFailure(const BlockID gid) const {
            return (cff[gid] >= 0);
        }

        //! Whether the block experienced dynamic friction failure.
        //! This occurs if the block is on the same fault as the trigger
        //! and it undergoes a significant CFF shift during this event
        //! (where significant is defined in terms of dynamic_val).
        bool dynamicFailure(const BlockID gid, const FaultID event_fault) const {
            return (getBlock(gid).getFaultID() == event_fault &&
                    cff[gid] > cff0[gid] &&
                    //state.cff > getStressDrop() &&
                    fabs((cff0[gid]-cff[gid])/cff0[gid]) > dynamic_val[gid]);
        }

        //! Calculate the expected recurrence of this block in years.
        // TODO: fix this, it is currently incorrect for multiple block models
        double getRecurrence(const BlockID gid) const {
            return stress_drop[gid]/(getBlock(gid).slip_rate()*(self_shear[gid]-self_normal[gid]*friction[gid]));
        };

        //! Set the Greens function self shear of this block.
        void setSelfStresses(const BlockID gid, const double new_self_shear, const double new_self_normal) {
            self_shear[gid] = new_self_shear;
            self_normal[gid] = new_self_normal;
        };
        //! Get the Greens function self shear of this block.
        double getSelfStresses(const BlockID gid) const {
            return self_shear[gid] - friction[gid]*self_normal[gid];
        };

        //! Get the stress drop for this block in bars.
        double getStressDrop(const BlockID gid) const {
            return stress_drop[gid];
        };
        //! Set the stress drop for this block in bars.
        void setStressDrop(const BlockID gid, const double new_stress_drop) {
            stress_drop[gid] = new_stress_drop;
            calcFriction(gid);
        };

        double getRhogd(const BlockID gid) {
            return rhogd[gid];
        };
        void setRhogd(const BlockID gid, const double new_rhogd) {
            rhogd[gid] = new_rhogd;
            calcFriction(gid);
        };

        void setDynamicVal(const BlockID gid, const double new_dynamic_val) {
            dynamic_val[gid] = new_dynamic_val;
        };

        void setInitShearNormalStress(const BlockID gid, const double new_init_shear, const double new_init_normal) {
            init_shear_stress[gid] = new_init_shear;
            init_normal_stress[gid] = new_init_normal;
        };

        double getInitShearStress(const BlockID gid) const {
            if (init_shear_stress.count(gid)) return init_shear_stress.find(gid)->second;

            return std::numeric_limits<double>::quiet_NaN();
        };
        double getInitNormalStress(const BlockID gid) const {
            if (init_normal_stress.count(gid)) return init_normal_stress.find(gid)->second;

            return std::numeric_limits<double>::quiet_NaN();
        };

        bool getFailed(const BlockID bid) const {
            return failed[bid];
        };
        void setFailed(const BlockID bid, const bool in_failed) {
            failed[bid] = in_failed;
        };

        double getUpdateField(const BlockID &b) const {
            return update_field[b];
        };
        void setUpdateField(const BlockID &b, const double new_val) {
            update_field[b] = new_val;
        };

        double getNormalStress(const BlockID &b) const {
            return normal_stress[b];
        };
        void setNormalStress(const BlockID &b, const double &new_val) {
            normal_stress[b] = new_val;
        };

        double getShearStress(const BlockID &b) const {
            return shear_stress[b];
        };
        void setShearStress(const BlockID &b, const double &new_val) {
            shear_stress[b] = new_val;
        };

        void saveStresses(const BlockID bid) {
            shear_stress0[bid] = shear_stress[bid];
            normal_stress0[bid] = normal_stress[bid];
            cff0[bid] = cff[bid];
        }

        void partitionBlocks(void);

        void output_stress(quakelib::UIndex event_num, quakelib::UIndex sweep_num);

    private:
#ifdef DEBUG
        int                         mult_timer;

        //! Counter for number of matrix multiplies performed
        int                         num_mults;
#endif
        //! Current simulation year
        double                      year;

        //! Temporary buffer used to speed up calculations
        double                      *mult_buffer;
        GREEN_VAL                   *decompress_buf;

        //! Files to write stress records to
        std::ofstream       stress_index_outfile, stress_outfile;

#ifdef HDF5_FOUND
        // HDF5 handle to stress data file
        hid_t               stress_data_file;
        void open_stress_hdf5_file(const std::string &hdf5_file_name);
#endif

        //! Number of stress records written to files, used for keeping track of indices
        unsigned int        num_stress_recs;

        //! Map of which blocks have which neighbors
        std::map<BlockID, quakelib::ElementIDSet>   neighbor_map;
};

#endif
