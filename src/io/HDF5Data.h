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

#include "Block.h"
#include "Event.h"
#include <map>
#include <set>
#include <fstream>

#ifndef _HDF5_DATA_H_
#define _HDF5_DATA_H_

// HDF5 greens file names
#define GREEN_SHEAR_HDF5            "greens_shear"
#define GREEN_NORMAL_HDF5           "greens_normal"

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
        unsigned int getGreensDim(void) const {
            return greens_dim;
        };
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

#endif

#endif
