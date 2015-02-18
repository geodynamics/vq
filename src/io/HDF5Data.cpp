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

#include "HDF5Data.h"
#include "QuakeLibIO.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string.h>
#include <sstream>

#ifdef HDF5_FOUND

HDF5CheckpointReader::HDF5CheckpointReader(const std::string &ckpt_file_name,
                                           double &checkpoint_year,
                                           unsigned int &checkpoint_event,
                                           CheckpointSet &checkpoints) : HDF5Checkpoint() {

    if (!H5Fis_hdf5(ckpt_file_name.c_str())) exit(-1);

    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    if (plist_id < 0) exit(-1);

#ifdef HDF5_IS_PARALLEL
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif

    data_file = H5Fopen(ckpt_file_name.c_str(), H5F_ACC_RDONLY, plist_id);

    if (data_file < 0) exit(-1);
}

// TODO: check status in this function
HDF5CheckpointWriter::HDF5CheckpointWriter(const std::string &ckpt_file_name,
                                           const unsigned int &nblocks,
                                           const double &cur_year,
                                           const unsigned int &cur_event,
                                           const CheckpointSet &checkpoints) : HDF5Checkpoint() {
    /*hsize_t                         dims[2] = {nblocks, CHECKPOINT_NUM_ENTRIES_HDF5};
    hsize_t                         single_val[1] = {1};
    herr_t                          status;
    CheckpointSet::const_iterator   it;
    StateCheckpointData             *mem_state;
    unsigned int                    i, num_local;
    hsize_t                         start[2], count[2];
    hid_t                           file_select, mem_select, xfer_plist_id;
    herr_t                          res;

    // Create access properties
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    if (plist_id < 0) exit(-1);

    #ifdef HDF5_IS_PARALLEL
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    #endif

    // Create the data file, overwriting any old files
    data_file = H5Fcreate(ckpt_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    if (data_file < 0) exit(-1);

    // Create data spaces with the proper dimensions
    state_dataspace = H5Screate_simple(2, dims, dims);
    year_event_dataspace = H5Screate_simple(1, single_val, single_val);

    // Create the checkpoint data sets for state information and current year/event
    state_dataset = H5Dcreate2(data_file, CHECKPOINT_STATE_HDF5, H5T_NATIVE_DOUBLE,
                               state_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (state_dataset < 0) exit(-1);

    year_dataset = H5Dcreate2(data_file, CHECKPOINT_YEAR_HDF5, H5T_NATIVE_DOUBLE,
                              year_event_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (year_dataset < 0) exit(-1);

    event_dataset = H5Dcreate2(data_file, CHECKPOINT_EVENT_HDF5, H5T_NATIVE_UINT,
                               year_event_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (event_dataset < 0) exit(-1);

    // Write the current simulation year
    status = H5Dwrite(year_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cur_year);

    // Write the current simulation event number
    status = H5Dwrite(event_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cur_event);

    // Start with a blank file and memory hyperslab selection
    num_local = checkpoints.size();

    file_select = H5Scopy(state_dataspace);
    mem_select = H5Scopy(state_dataspace);

    status = H5Sselect_none(file_select);
    status = H5Sselect_none(mem_select);

    // Add the correct state dimensions to the memory hyperslab
    start[0] = start[1] = 0;
    count[0] = num_local;
    count[1] = CHECKPOINT_NUM_ENTRIES_HDF5;
    status = H5Sselect_hyperslab(mem_select, H5S_SELECT_OR, start, NULL, count, NULL);

    // Copy block state into memory and assign the proper file hyperslab elements
    mem_state = new StateCheckpointData[num_local];
    count[0] = 1;

    for (i=0,it=checkpoints.begin(); it!=checkpoints.end(); ++i,++it) {
        start[0] = it->first;       // assign this block to the correct BlockID position
        memcpy(&(mem_state[i]), &(it->second), sizeof(StateCheckpointData));
        status = H5Sselect_hyperslab(file_select, H5S_SELECT_OR, start, NULL, count, NULL);
    }

    // Write all block state data in parallel
    xfer_plist_id = H5Pcreate(H5P_DATASET_XFER);
    #ifdef HDF5_IS_PARALLEL
    H5Pset_dxpl_mpio(xfer_plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    status = H5Dwrite(state_dataset, H5T_NATIVE_DOUBLE, mem_select, file_select, xfer_plist_id, mem_state);

    res = H5Sclose(file_select);
    res = H5Sclose(mem_select);
    res = H5Sclose(state_dataspace);
    res = H5Sclose(year_event_dataspace);

    res = H5Dclose(state_dataset);
    res = H5Dclose(year_dataset);
    res = H5Dclose(event_dataset);

    res = H5Pclose(xfer_plist_id);
    res = H5Pclose(plist_id);

    res = H5Fclose(data_file);

    delete mem_state;*/
}

/*!
 Initialize the HDF5 Greens value reader.
 */
HDF5GreensDataReader::HDF5GreensDataReader(const std::string &hdf5_file_name) : HDF5GreensData() {
    int         ndims;
    hsize_t     *dims;

    if (!H5Fis_hdf5(hdf5_file_name.c_str())) exit(-1);

    // Open the data file in read only mode
    data_file = H5Fopen(hdf5_file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    if (data_file < 0) exit(-1);

    // Open the data sets
    green_shear_set = H5Dopen2(data_file, GREEN_SHEAR_HDF5, H5P_DEFAULT);

    if (green_shear_set < 0) exit(-1);

    green_norm_set = H5Dopen2(data_file, GREEN_NORMAL_HDF5, H5P_DEFAULT);

    if (green_norm_set < 0) exit(-1);

    green_dataspace = H5Dget_space(green_shear_set);

    if (green_dataspace < 0) exit(-1);

    ndims = H5Sget_simple_extent_ndims(green_dataspace);

    if (ndims < 0) exit(-1);

    dims = new hsize_t[ndims];
    ndims = H5Sget_simple_extent_dims(green_dataspace, dims, NULL);

    if (ndims < 0) exit(-1);

    greens_dim = dims[0];
    delete dims;
}

HDF5GreensData::~HDF5GreensData(void) {
    herr_t      res;

    // Close the handles we've used
    res = H5Dclose(green_norm_set);

    if (res < 0) exit(-1);

    res = H5Dclose(green_shear_set);

    if (res < 0) exit(-1);

    res = H5Sclose(green_dataspace);

    if (res < 0) exit(-1);

    res = H5Fclose(data_file);

    if (res < 0) exit(-1);
}

/*!
 Initialize the HDF5 writer using the specified model dimensions.
 */
HDF5GreensDataWriter::HDF5GreensDataWriter(const std::string &hdf5_file_name, const unsigned int &nblocks) : HDF5GreensData() {
    hsize_t             dimsf[2];
    hid_t               plist_id;

    greens_dim = nblocks;

    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    if (plist_id < 0) exit(-1);

#ifdef HDF5_IS_PARALLEL
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
    // Create the data file, overwriting any old files
    data_file = H5Fcreate(hdf5_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    if (data_file < 0) exit(-1);

    // Specify the dimensions of the Green's matrix
    dimsf[0] = greens_dim;
    dimsf[1] = greens_dim;

    // Create the Green's matrix dataspace
    green_dataspace = H5Screate_simple(2, dimsf, NULL);

    if (green_dataspace < 0) exit(-1);

    // Create the Greens value data sets
    green_shear_set = H5Dcreate2(data_file, GREEN_SHEAR_HDF5, H5T_NATIVE_DOUBLE, green_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (green_shear_set < 0) exit(-1);

    green_norm_set = H5Dcreate2(data_file, GREEN_NORMAL_HDF5, H5T_NATIVE_DOUBLE, green_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (green_norm_set < 0) exit(-1);

    H5Pclose(plist_id);
}

/*!
 Set the Greens function normal and shear values for a specified block.
 */
void HDF5GreensDataWriter::setGreensVals(const int &bid, const double *shear_vals, const double *norm_vals) {
    herr_t      status;
    hsize_t     file_start[2], mem_start[2], count[2];
    hid_t       file_select, mem_select, plist_id;

    // Copy the selector for the entire dataspace
    file_select = H5Scopy(green_dataspace);
    mem_select = H5Scopy(green_dataspace);

    if (bid == UNDEFINED_ELEMENT_ID) {
        file_start[0] = file_start[1] = 0;
        count[0] = count[1] = 0;
    } else {
        file_start[0] = bid;                // start at xth block
        file_start[1] = 0;
        count[0] = 1;                       // 1xN set of values
        count[1] = greens_dim;
    }

    mem_start[0] = mem_start[1] = 0;    // start at element 0 in memory array

    // Select the hyperslabs for the memory and file dataspace
    status = H5Sselect_hyperslab(file_select, H5S_SELECT_SET, file_start, NULL, count, NULL);
    status = H5Sselect_hyperslab(mem_select, H5S_SELECT_SET, mem_start, NULL, count, NULL);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef HDF5_IS_PARALLEL
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif
    status = H5Dwrite(green_shear_set, H5T_NATIVE_DOUBLE, mem_select, file_select, plist_id, shear_vals);
    status = H5Dwrite(green_norm_set, H5T_NATIVE_DOUBLE, mem_select, file_select, plist_id, norm_vals);

    H5Pclose(plist_id);
    H5Sclose(file_select);
    H5Sclose(mem_select);
}

/*!
 Set the Greens function normal and shear values for a specified block.
 */
void HDF5GreensDataReader::getGreensVals(const int &bid, double *shear_vals, double *norm_vals) {
    herr_t      status;
    hsize_t     file_start[2], mem_start[2], count[2];
    hid_t       file_select, mem_select, plist_id;

    // Copy the selector for the entire dataspace
    file_select = H5Scopy(green_dataspace);

    if (file_select < 0) exit(-1);

    mem_select = H5Scopy(green_dataspace);

    if (mem_select < 0) exit(-1);

    file_start[0] = bid;                // start at xth block
    file_start[1] = 0;
    mem_start[0] = mem_start[1] = 0;    // start at element 0 in memory array
    count[0] = 1;                       // 1xN set of values
    count[1] = greens_dim;

    // Select the hyperslabs for the memory and file dataspace
    status = H5Sselect_hyperslab(file_select, H5S_SELECT_SET, file_start, NULL, count, NULL);

    if (status < 0) exit(-1);

    status = H5Sselect_hyperslab(mem_select, H5S_SELECT_SET, mem_start, NULL, count, NULL);

    if (status < 0) exit(-1);

    plist_id = H5Pcreate(H5P_DATASET_XFER);

    if (plist_id < 0) exit(-1);

    status = H5Dread(green_shear_set, H5T_NATIVE_DOUBLE, mem_select, file_select, plist_id, shear_vals);

    if (status < 0) exit(-1);

    status = H5Dread(green_norm_set, H5T_NATIVE_DOUBLE, mem_select, file_select, plist_id, norm_vals);

    if (status < 0) exit(-1);

    H5Pclose(plist_id);
    H5Sclose(file_select);
    H5Sclose(mem_select);
}

#endif
