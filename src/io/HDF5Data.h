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

#ifdef MPI_C_FOUND
#include <mpi.h>
#endif

#ifdef HDF5_FOUND
#include "hdf5.h"
#include "hdf5_hl.h"
#endif

#include "VCBlock.h"
#include "VCEvent.h"
#include <map>
#include <set>
#include <fstream>

#ifndef _HDF5_DATA_H_
#define _HDF5_DATA_H_

// HDF5 greens file names
#define GREEN_SHEAR_HDF5            "greens_shear"
#define GREEN_NORMAL_HDF5           "greens_normal"

// HDF5 file data definitions
#define SIM_YEARS_HDF5              "sim_years"
#define BASE_LAT_LON_HDF5           "base_lat_lon"

// Event info related definitions
#define EVENT_TABLE_HDF5            "event_table"
#define EVENT_NUM_ENTRIES_HDF5      10

// Event sweeps table definitions
#define SWEEP_TABLE_HDF5            "event_sweep_table"
#define SWEEP_NUM_ENTRIES_HDF5      10

// State checkpoint table definitions
#define CHECKPOINT_STATE_HDF5       "checkpoint_state"
#define CHECKPOINT_YEAR_HDF5        "checkpoint_year"
#define CHECKPOINT_EVENT_HDF5       "checkpoint_event"
#define CHECKPOINT_NUM_ENTRIES_HDF5 6
#define CHECKPOINT_SLIP_DFCT_HDF5   "slipDeficit"
#define CHECKPOINT_CFF_HDF5         "cff"
#define CHECKPOINT_STRESS_S_HDF5    "stressS"
#define CHECKPOINT_STRESS_N_HDF5    "stressN"
#define CHECKPOINT_UPDATE_FLD_HDF5  "updateField"

struct EventInfo {
    unsigned int    event_number;
    double          event_year;
    BlockID         event_trigger;
    double          event_magnitude;
    unsigned int    start_sweep_rec, end_sweep_rec;
    double          init_shear, final_shear, init_normal, final_normal;
};

typedef struct EventInfo EventInfo;

struct EventSweepInfo {
    unsigned int    event_number;
    unsigned int    sweep_number;
    BlockID         block_id;
    double          block_slip;
    double          block_area;
    double          block_mu;
    double          shear_init, shear_final;
    double          normal_init, normal_final;
};

typedef struct EventSweepInfo EventSweepInfo;

#ifdef HDF5_FOUND

// Classes representing a file containing checkpoint data
class HDF5Checkpoint {
    protected:
        // HDF5 handle to checkpoint data file
        hid_t               data_file;

        // Dataspace of checkpoint state
        hid_t               state_dataspace, year_event_dataspace;

        // Handle to checkpoint state dataset
        hid_t               state_dataset, year_dataset, event_dataset;

        // Access control handle
        hid_t               plist_id;
};

class HDF5CheckpointReader : public HDF5Checkpoint {
    public:
        HDF5CheckpointReader(const std::string &ckpt_file_name,
                             double &checkpoint_year,
                             unsigned int &checkpoint_event,
                             CheckpointSet &checkpoints);
};

class HDF5CheckpointWriter : public HDF5Checkpoint {
    public:
        HDF5CheckpointWriter(const std::string &ckpt_file_name,
                             const unsigned int &nblocks,
                             const double &cur_year,
                             const unsigned int &cur_event,
                             const CheckpointSet &checkpoints);
};

// Classes representing a file containing Greens function calculation data
class HDF5GreensData {
    protected:
        // HDF5 handle to data file
        hid_t               data_file;

        // Handles to Greens data in the file
        hid_t               green_norm_set, green_shear_set;

        // Handles to data space specifications
        hid_t               green_dataspace;

        // Dimension of Greens matrix
        unsigned int        greens_dim;

    public:
        HDF5GreensData(void) {};
        ~HDF5GreensData(void);
};


class HDF5GreensDataReader : public HDF5GreensData {
    public:
        HDF5GreensDataReader(const std::string &hdf5_file_name);

        void getGreensVals(const int &bid, double *shear_vals, double *norm_vals);
};

class HDF5GreensDataWriter : public HDF5GreensData {
    public:
        HDF5GreensDataWriter(const std::string &hdf5_file_name, const unsigned int &nblocks);
        void setGreensVals(const int &bid, const double *shear_vals, const double *norm_vals);
};

class HDF5Data {
    protected:
        // HDF5 handle to data file
        hid_t               data_file;

        // Handles to data in the file
        hid_t               sim_years_set;

        // Handles to data space specifications
        hid_t               pair_val_dataspace;

        // Values read from shared memory (pointers are set to within shared memory segment)
        unsigned int    num_blocks;

        void createH5Handles(void);

    public:
        HDF5Data(void) {};
        ~HDF5Data(void);
        unsigned int modelDim(void) const {
            return num_blocks;
        };
};

class HDF5DataWriter : public HDF5Data {
    public:
        HDF5DataWriter(const std::string &hdf5_file_name);
        void setStartEndYears(const double &new_start_year, const double &new_end_year);
        void flush(void);

        void writeEvent(VCEvent &event);
};

#endif

#endif
