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

#include "HDF5Data.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string.h>
#include <sstream>

HDF5CheckpointReader::HDF5CheckpointReader(const std::string &ckpt_file_name,
										   double &checkpoint_year,
										   unsigned int &checkpoint_event,
										   CheckpointSet &checkpoints) : HDF5Checkpoint() {
#ifdef HAVE_HDF5
	if (!H5Fis_hdf5(ckpt_file_name.c_str())) exit(-1);
	
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	if (plist_id < 0) exit(-1);
#ifdef HAVE_MPI
#ifdef H5_HAVE_PARALLEL
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
#endif
	
	data_file = H5Fopen(ckpt_file_name.c_str(), H5F_ACC_RDONLY, plist_id);
	if (data_file < 0) exit(-1);
	
	
#endif
}

// TODO: check status in this function
HDF5CheckpointWriter::HDF5CheckpointWriter(const std::string &ckpt_file_name,
										   const unsigned int &nblocks,
										   const double &cur_year,
										   const unsigned int &cur_event,
										   const CheckpointSet &checkpoints) : HDF5Checkpoint() {
#ifdef HAVE_HDF5
	hsize_t							dims[2] = {nblocks, CHECKPOINT_NUM_ENTRIES_HDF5};
	hsize_t							single_val[1] = {1};
	herr_t							status;
	CheckpointSet::const_iterator	it;
	StateCheckpointData				*mem_state;
	unsigned int					i, num_local;
	hsize_t							start[2], count[2];
	hid_t							file_select, mem_select, xfer_plist_id;
	herr_t							res;
	
	// Create access properties
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	if (plist_id < 0) exit(-1);
#ifdef HAVE_MPI
#ifdef H5_HAVE_PARALLEL
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
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
	for (i=0,it=checkpoints.begin();it!=checkpoints.end();++i,++it) {
		start[0] = it->first;		// assign this block to the correct BlockID position
		memcpy(&(mem_state[i]), &(it->second), sizeof(StateCheckpointData));
		status = H5Sselect_hyperslab(file_select, H5S_SELECT_OR, start, NULL, count, NULL);
	}

	// Write all block state data in parallel
	xfer_plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
#ifdef H5_HAVE_PARALLEL
	H5Pset_dxpl_mpio(xfer_plist_id, H5FD_MPIO_COLLECTIVE);
#endif
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
	
	delete mem_state;
#endif
}

/*!
 Initialize the HDF5 Greens value reader.
 */
HDF5GreensDataReader::HDF5GreensDataReader(const std::string &hdf5_file_name) : HDF5GreensData() {
#ifdef HAVE_HDF5
	int			ndims;
	hsize_t		*dims;
	
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
#endif
}

HDF5GreensData::~HDF5GreensData(void) {
#ifdef HAVE_HDF5
	herr_t		res;
	
	// Close the handles we've used
	res = H5Dclose(green_norm_set);
	if (res < 0) exit(-1);
	res = H5Dclose(green_shear_set);
	if (res < 0) exit(-1);
	
	res = H5Sclose(green_dataspace);
	if (res < 0) exit(-1);
	
	res = H5Fclose(data_file);
	if (res < 0) exit(-1);
#endif
}

/*!
 Initialize the HDF5 writer using the specified model dimensions.
 */
HDF5GreensDataWriter::HDF5GreensDataWriter(const std::string &hdf5_file_name, const unsigned int &nblocks) : HDF5GreensData() {
#ifdef HAVE_HDF5
	hsize_t				dimsf[2];
	hid_t				plist_id;
	
	greens_dim = nblocks;
	
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	if (plist_id < 0) exit(-1);
#ifdef HAVE_MPI
#ifdef H5_HAVE_PARALLEL
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
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
#endif
}

/*!
 Set the Greens function normal and shear values for a specified block.
 */
void HDF5GreensDataWriter::setGreensVals(const int &bid, const double *shear_vals, const double *norm_vals) {
#ifdef HAVE_HDF5
	herr_t		status;
	hsize_t		file_start[2], mem_start[2], count[2];
	hid_t		file_select, mem_select, plist_id;
	
	// Copy the selector for the entire dataspace
	file_select = H5Scopy(green_dataspace);
	mem_select = H5Scopy(green_dataspace);
	
	file_start[0] = bid;				// start at xth block
	file_start[1] = 0;
	mem_start[0] = mem_start[1] = 0;	// start at element 0 in memory array
	count[0] = 1;						// 1xN set of values
	count[1] = greens_dim;
	
	// Select the hyperslabs for the memory and file dataspace
	status = H5Sselect_hyperslab(file_select, H5S_SELECT_SET, file_start, NULL, count, NULL);
	status = H5Sselect_hyperslab(mem_select, H5S_SELECT_SET, mem_start, NULL, count, NULL);
	
	plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
#ifdef H5_HAVE_PARALLEL
	//H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif
#endif	
	status = H5Dwrite(green_shear_set, H5T_NATIVE_DOUBLE, mem_select, file_select, plist_id, shear_vals);
	status = H5Dwrite(green_norm_set, H5T_NATIVE_DOUBLE, mem_select, file_select, plist_id, norm_vals);
	
	H5Pclose(plist_id);
	H5Sclose(file_select);
	H5Sclose(mem_select);
#endif
}

/*!
 Set the Greens function normal and shear values for a specified block.
 */
void HDF5GreensDataReader::getGreensVals(const int &bid, double *shear_vals, double *norm_vals) {
#ifdef HAVE_HDF5
	herr_t		status;
	hsize_t		file_start[2], mem_start[2], count[2];
	hid_t		file_select, mem_select, plist_id;
	
	// Copy the selector for the entire dataspace
	file_select = H5Scopy(green_dataspace);
	mem_select = H5Scopy(green_dataspace);
	
	file_start[0] = bid;				// start at xth block
	file_start[1] = 0;
	mem_start[0] = mem_start[1] = 0;	// start at element 0 in memory array
	count[0] = 1;						// 1xN set of values
	count[1] = greens_dim;
	
	// Select the hyperslabs for the memory and file dataspace
	status = H5Sselect_hyperslab(file_select, H5S_SELECT_SET, file_start, NULL, count, NULL);
	status = H5Sselect_hyperslab(mem_select, H5S_SELECT_SET, mem_start, NULL, count, NULL);
	
	plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
#ifdef H5_HAVE_PARALLEL
	//H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif
#endif	
	status = H5Dread(green_shear_set, H5T_NATIVE_DOUBLE, mem_select, file_select, plist_id, shear_vals);
	status = H5Dread(green_norm_set, H5T_NATIVE_DOUBLE, mem_select, file_select, plist_id, norm_vals);
	
	H5Pclose(plist_id);
	H5Sclose(file_select);
	H5Sclose(mem_select);
#endif
}

/*!
 Initialize the HDF5 reader.
 This opens the memory mapped file and sets the class pointers
 to the correct positions in the file.
 */
HDF5DataReader::HDF5DataReader(const std::string &hdf5_file_name) : HDF5Data() {
#ifdef HAVE_HDF5
	hsize_t		*dims;
	
	last_event_read = 0;
	
	if (!H5Fis_hdf5(hdf5_file_name.c_str())) exit(-1);
	
	// Open the data file in read only mode
	data_file = H5Fopen(hdf5_file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (data_file < 0) exit(-1);
	
	createH5Handles();
	
	sim_years_set = H5Dopen2(data_file, SIM_YEARS_HDF5, H5P_DEFAULT);
	if (sim_years_set < 0) exit(-1);
	base_lat_lon_set = H5Dopen2(data_file, BASE_LAT_LON_HDF5, H5P_DEFAULT);
	if (base_lat_lon_set < 0) exit(-1);
	
	// Open the data sets
	num_blocks = dims[0];
	delete dims;
#endif
}

HDF5Data::~HDF5Data(void) {
#ifdef HAVE_HDF5
	herr_t		res;
	
	// Close the handles we've used
	res = H5Dclose(sim_years_set);
	if (res < 0) exit(-1);
	res = H5Dclose(base_lat_lon_set);
	if (res < 0) exit(-1);
	
	res = H5Sclose(block_info_dataspace);
	if (res < 0) exit(-1);
	
	res = H5Tclose(fault_name_datatype);
	if (res < 0) exit(-1);
	
	res = H5Fclose(data_file);
	if (res < 0) exit(-1);
#endif
}

/*!
 Initialize the HDF5 writer using the specified model dimensions.
 */
HDF5DataWriter::HDF5DataWriter(const std::string &hdf5_file_name, const int &nblocks) : HDF5Data() {
#ifdef HAVE_HDF5
	hsize_t				dimsf[2];
	herr_t				status;
	hid_t				plist_id;
	BlockInfo			empty_binfo;
	//std::cout << "there " << HAVE_HDF5 << std::endl;
	
	num_blocks = nblocks;
	
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	if (plist_id < 0) exit(-1);
#ifdef HAVE_MPI
#ifdef H5_HAVE_PARALLEL
	//H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
#endif
	// Create the data file, overwriting any old files
	data_file = H5Fcreate(hdf5_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	if (data_file < 0) exit(-1);
	
	createH5Handles();
	
	// Specify the dimensions of the Green's matrix
	dimsf[0] = num_blocks;
	dimsf[1] = num_blocks;
	
	// Create the block info dataspace
	block_info_dataspace = H5Screate_simple(1, dimsf, NULL);
	if (block_info_dataspace < 0) exit(-1);
	
	// Create entries for the simulation start/stop years and base longitude/latitude
	sim_years_set = H5Dcreate2(data_file, SIM_YEARS_HDF5, H5T_NATIVE_DOUBLE, pair_val_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (sim_years_set < 0) exit(-1);
	base_lat_lon_set = H5Dcreate2(data_file, BASE_LAT_LON_HDF5, H5T_NATIVE_DOUBLE, pair_val_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (base_lat_lon_set < 0) exit(-1);
	
	// Create the block info table
	status = H5TBmake_table("Block Info Table",
							data_file,
							B_INFO_TABLE_HDF5,
							B_INFO_NUM_ENTRIES_HDF5,
							num_blocks,
							sizeof(BlockInfo),
							binfo_field_names,
							binfo_field_offsets,
							binfo_field_types,
							num_blocks,
							&empty_binfo,
							0,
							NULL);
	if (status < 0) exit(-1);
	
	// Create the event table
	status = H5TBmake_table("Event Table",
							data_file,
							EVENT_TABLE_HDF5,
							EVENT_NUM_ENTRIES_HDF5,
							0,
							sizeof(EventInfo),
							event_field_names,
							event_field_offsets,
							event_field_types,
							100,
							NULL,
							0,
							NULL);
	if (status < 0) exit(-1);
	
	// Create the event sweeps table
	status = H5TBmake_table("Sweeps Table",
							data_file,
							SWEEP_TABLE_HDF5,
							SWEEP_NUM_ENTRIES_HDF5,
							0,
							sizeof(EventSweepInfo),
							sweep_field_names,
							sweep_field_offsets,
							sweep_field_types,
							1000,
							NULL,
							0,
							NULL);
	if (status < 0) exit(-1);
	
	// Create the aftershock table
	status = H5TBmake_table("Aftershock Table",
							data_file,
							AFTERSHOCK_TABLE_HDF5,
							AFTERSHOCK_NUM_ENTRIES_HDF5,
							0,
							sizeof(AftershockBGInfo),
							aftershock_field_names,
							aftershock_field_offsets,
							aftershock_field_types,
							1000,
							NULL,
							0,
							NULL);
	if (status < 0) exit(-1);
	
	// Create the background event table (incomplete)
	status = H5TBmake_table("Background Event Table",
							data_file,
							BG_EVENT_TABLE_HDF5,
							AFTERSHOCK_NUM_ENTRIES_HDF5,
							0,
							sizeof(AftershockBGInfo),
							aftershock_field_names,
							aftershock_field_offsets,
							aftershock_field_types,
							1000,
							NULL,
							0,
							NULL);
	if (status < 0) exit(-1);
	
	H5Pclose(plist_id);
#endif
}

void HDF5Data::createH5Handles(void) {
#ifdef HAVE_HDF5
	hsize_t		dimsf[2];
	
	num_layers = 4;
	
	// Create dataspace for fault name strings
	dimsf[0] = FAULT_NAME_MAX_LEN;
	fault_name_dataspace = H5Screate_simple(1, dimsf, NULL);
	
	// Create the datatype for the fault name strings
	fault_name_datatype = H5Tcopy(H5T_C_S1);
	H5Tset_size(fault_name_datatype, (size_t)FAULT_NAME_MAX_LEN);
	
	binfo_field_names[0] = B_INFO_BLOCK_ID_HDF5;
	binfo_field_names[1] = B_INFO_FAULT_ID_HDF5;
	binfo_field_names[2] = B_INFO_SECTION_ID_HDF5;
	binfo_field_names[3] = B_INFO_M_X_PT1_HDF5;
	binfo_field_names[4] = B_INFO_M_Y_PT1_HDF5;
	binfo_field_names[5] = B_INFO_M_Z_PT1_HDF5;
	binfo_field_names[6] = B_INFO_M_DAS_PT1_HDF5;
	binfo_field_names[7] = B_INFO_M_TRACE_FLAG_PT1_HDF5;
	binfo_field_names[8] = B_INFO_M_X_PT2_HDF5;
	binfo_field_names[9] = B_INFO_M_Y_PT2_HDF5;
	binfo_field_names[10] = B_INFO_M_Z_PT2_HDF5;
	binfo_field_names[11] = B_INFO_M_DAS_PT2_HDF5;
	binfo_field_names[12] = B_INFO_M_TRACE_FLAG_PT2_HDF5;
	binfo_field_names[13] = B_INFO_M_X_PT3_HDF5;
	binfo_field_names[14] = B_INFO_M_Y_PT3_HDF5;
	binfo_field_names[15] = B_INFO_M_Z_PT3_HDF5;
	binfo_field_names[16] = B_INFO_M_DAS_PT3_HDF5;
	binfo_field_names[17] = B_INFO_M_TRACE_FLAG_PT3_HDF5;
	binfo_field_names[18] = B_INFO_M_X_PT4_HDF5;
	binfo_field_names[19] = B_INFO_M_Y_PT4_HDF5;
	binfo_field_names[20] = B_INFO_M_Z_PT4_HDF5;
	binfo_field_names[21] = B_INFO_M_DAS_PT4_HDF5;
	binfo_field_names[22] = B_INFO_M_TRACE_FLAG_PT4_HDF5;
	binfo_field_names[23] = B_INFO_SLIP_VELOCITY_HDF5;
	binfo_field_names[24] = B_INFO_ASEISMICITY_HDF5;
	binfo_field_names[25] = B_INFO_RAKE_HDF5;
	binfo_field_names[26] = B_INFO_DIP_HDF5;
	binfo_field_names[27] = B_INFO_DYNAMIC_STRENGTH_HDF5;
	binfo_field_names[28] = B_INFO_STATIC_STRENGTH_HDF5;
	binfo_field_names[29] = B_INFO_LAME_MU_HDF5;
	binfo_field_names[30] = B_INFO_LAME_LAMBDA_HDF5;
	binfo_field_names[31] = B_INFO_FAULT_NAME_HDF5;
	
	binfo_field_offsets[0] = HOFFSET(BlockInfo, bid);
	binfo_field_offsets[1] = HOFFSET(BlockInfo, fid);
	binfo_field_offsets[2] = HOFFSET(BlockInfo, sid);
	binfo_field_offsets[3] = HOFFSET(BlockInfo, x_pt[0]);
	binfo_field_offsets[4] = HOFFSET(BlockInfo, y_pt[0]);
	binfo_field_offsets[5] = HOFFSET(BlockInfo, z_pt[0]);
	binfo_field_offsets[6] = HOFFSET(BlockInfo, das_pt[0]);
	binfo_field_offsets[7] = HOFFSET(BlockInfo, trace_flag_pt[0]);
	binfo_field_offsets[8] = HOFFSET(BlockInfo, x_pt[1]);
	binfo_field_offsets[9] = HOFFSET(BlockInfo, y_pt[1]);
	binfo_field_offsets[10] = HOFFSET(BlockInfo, z_pt[1]);
	binfo_field_offsets[11] = HOFFSET(BlockInfo, das_pt[1]);
	binfo_field_offsets[12] = HOFFSET(BlockInfo, trace_flag_pt[1]);
	binfo_field_offsets[13] = HOFFSET(BlockInfo, x_pt[2]);
	binfo_field_offsets[14] = HOFFSET(BlockInfo, y_pt[2]);
	binfo_field_offsets[15] = HOFFSET(BlockInfo, z_pt[2]);
	binfo_field_offsets[16] = HOFFSET(BlockInfo, das_pt[2]);
	binfo_field_offsets[17] = HOFFSET(BlockInfo, trace_flag_pt[2]);
	binfo_field_offsets[18] = HOFFSET(BlockInfo, x_pt[3]);
	binfo_field_offsets[19] = HOFFSET(BlockInfo, y_pt[3]);
	binfo_field_offsets[20] = HOFFSET(BlockInfo, z_pt[3]);
	binfo_field_offsets[21] = HOFFSET(BlockInfo, das_pt[3]);
	binfo_field_offsets[22] = HOFFSET(BlockInfo, trace_flag_pt[3]);
	binfo_field_offsets[23] = HOFFSET(BlockInfo, slip_velocity);
	binfo_field_offsets[24] = HOFFSET(BlockInfo, aseismicity);
	binfo_field_offsets[25] = HOFFSET(BlockInfo, rake);
	binfo_field_offsets[26] = HOFFSET(BlockInfo, dip);
	binfo_field_offsets[27] = HOFFSET(BlockInfo, dynamic_strength);
	binfo_field_offsets[28] = HOFFSET(BlockInfo, static_strength);
	binfo_field_offsets[29] = HOFFSET(BlockInfo, lame_mu);
	binfo_field_offsets[30] = HOFFSET(BlockInfo, lame_lambda);
	binfo_field_offsets[31] = HOFFSET(BlockInfo, fault_name);
	
	binfo_field_types[0] = H5T_NATIVE_UINT;
	binfo_field_types[1] = H5T_NATIVE_UINT;
	binfo_field_types[2] = H5T_NATIVE_UINT;
	binfo_field_types[3] = H5T_NATIVE_DOUBLE;
	binfo_field_types[4] = H5T_NATIVE_DOUBLE;
	binfo_field_types[5] = H5T_NATIVE_DOUBLE;
	binfo_field_types[6] = H5T_NATIVE_DOUBLE;
	binfo_field_types[7] = H5T_NATIVE_UINT;
	binfo_field_types[8] = H5T_NATIVE_DOUBLE;
	binfo_field_types[9] = H5T_NATIVE_DOUBLE;
	binfo_field_types[10] = H5T_NATIVE_DOUBLE;
	binfo_field_types[11] = H5T_NATIVE_DOUBLE;
	binfo_field_types[12] = H5T_NATIVE_UINT;
	binfo_field_types[13] = H5T_NATIVE_DOUBLE;
	binfo_field_types[14] = H5T_NATIVE_DOUBLE;
	binfo_field_types[15] = H5T_NATIVE_DOUBLE;
	binfo_field_types[16] = H5T_NATIVE_DOUBLE;
	binfo_field_types[17] = H5T_NATIVE_UINT;
	binfo_field_types[18] = H5T_NATIVE_DOUBLE;
	binfo_field_types[19] = H5T_NATIVE_DOUBLE;
	binfo_field_types[20] = H5T_NATIVE_DOUBLE;
	binfo_field_types[21] = H5T_NATIVE_DOUBLE;
	binfo_field_types[22] = H5T_NATIVE_UINT;
	binfo_field_types[23] = H5T_NATIVE_DOUBLE;
	binfo_field_types[24] = H5T_NATIVE_DOUBLE;
	binfo_field_types[25] = H5T_NATIVE_DOUBLE;
	binfo_field_types[26] = H5T_NATIVE_DOUBLE;
	binfo_field_types[27] = H5T_NATIVE_DOUBLE;
	binfo_field_types[28] = H5T_NATIVE_DOUBLE;
	binfo_field_types[29] = H5T_NATIVE_DOUBLE;
	binfo_field_types[30] = H5T_NATIVE_DOUBLE;
	binfo_field_types[31] = fault_name_datatype;
	
	binfo_field_sizes[0] = sizeof(BlockID);
	binfo_field_sizes[1] = sizeof(FaultID);
	binfo_field_sizes[2] = sizeof(quakelib::SectionID);
	binfo_field_sizes[3] = sizeof(double);
	binfo_field_sizes[4] = sizeof(double);
	binfo_field_sizes[5] = sizeof(double);
	binfo_field_sizes[6] = sizeof(double);
	binfo_field_sizes[7] = sizeof(quakelib::TraceFlag);
	binfo_field_sizes[8] = sizeof(double);
	binfo_field_sizes[9] = sizeof(double);
	binfo_field_sizes[10] = sizeof(double);
	binfo_field_sizes[11] = sizeof(double);
	binfo_field_sizes[12] = sizeof(quakelib::TraceFlag);
	binfo_field_sizes[13] = sizeof(double);
	binfo_field_sizes[14] = sizeof(double);
	binfo_field_sizes[15] = sizeof(double);
	binfo_field_sizes[16] = sizeof(double);
	binfo_field_sizes[17] = sizeof(quakelib::TraceFlag);
	binfo_field_sizes[18] = sizeof(double);
	binfo_field_sizes[19] = sizeof(double);
	binfo_field_sizes[20] = sizeof(double);
	binfo_field_sizes[21] = sizeof(double);
	binfo_field_sizes[22] = sizeof(quakelib::TraceFlag);
	binfo_field_sizes[23] = sizeof(double);
	binfo_field_sizes[24] = sizeof(double);
	binfo_field_sizes[25] = sizeof(double);
	binfo_field_sizes[26] = sizeof(double);
	binfo_field_sizes[27] = sizeof(double);
	binfo_field_sizes[28] = sizeof(double);
	binfo_field_sizes[29] = sizeof(double);
	binfo_field_sizes[30] = sizeof(double);
	binfo_field_sizes[31] = sizeof(char)*FAULT_NAME_MAX_LEN;
	
	event_field_names[0] = EVENT_NUM_HDF5;
	event_field_names[1] = EVENT_YEAR_HDF5;
	event_field_names[2] = EVENT_TRIGGER_HDF5;
	event_field_names[3] = EVENT_MAGNITUDE_HDF5;
	event_field_names[4] = EVENT_SHEAR_INIT_HDF5;
	event_field_names[5] = EVENT_NORMAL_INIT_HDF5;
	event_field_names[6] = EVENT_SHEAR_FINAL_HDF5;
	event_field_names[7] = EVENT_NORMAL_FINAL_HDF5;
	event_field_names[8] = EVENT_START_SWEEP_HDF5;
	event_field_names[9] = EVENT_END_SWEEP_HDF5;
	event_field_names[10] = EVENT_START_AS_HDF5;
	event_field_names[11] = EVENT_END_AS_HDF5;
	
	event_field_offsets[0] = HOFFSET(EventInfo, event_number);
	event_field_offsets[1] = HOFFSET(EventInfo, event_year);
	event_field_offsets[2] = HOFFSET(EventInfo, event_trigger);
	event_field_offsets[3] = HOFFSET(EventInfo, event_magnitude);
	event_field_offsets[4] = HOFFSET(EventInfo, init_shear);
	event_field_offsets[5] = HOFFSET(EventInfo, init_normal);
	event_field_offsets[6] = HOFFSET(EventInfo, final_shear);
	event_field_offsets[7] = HOFFSET(EventInfo, final_normal);
	event_field_offsets[8] = HOFFSET(EventInfo, start_sweep_rec);
	event_field_offsets[9] = HOFFSET(EventInfo, end_sweep_rec);
	event_field_offsets[10] = HOFFSET(EventInfo, start_aftershock_rec);
	event_field_offsets[11] = HOFFSET(EventInfo, end_aftershock_rec);
	
	event_field_types[0] = H5T_NATIVE_UINT;
	event_field_types[1] = H5T_NATIVE_DOUBLE;
	event_field_types[2] = H5T_NATIVE_UINT;
	event_field_types[3] = H5T_NATIVE_DOUBLE;
	event_field_types[4] = H5T_NATIVE_DOUBLE;
	event_field_types[5] = H5T_NATIVE_DOUBLE;
	event_field_types[6] = H5T_NATIVE_DOUBLE;
	event_field_types[7] = H5T_NATIVE_DOUBLE;
	event_field_types[8] = H5T_NATIVE_UINT;
	event_field_types[9] = H5T_NATIVE_UINT;
	event_field_types[10] = H5T_NATIVE_UINT;
	event_field_types[11] = H5T_NATIVE_UINT;
	
	event_field_sizes[0] = sizeof(unsigned int);
	event_field_sizes[1] = sizeof(double);
	event_field_sizes[2] = sizeof(BlockID);
	event_field_sizes[3] = sizeof(double);
	event_field_sizes[4] = sizeof(double);
	event_field_sizes[5] = sizeof(double);
	event_field_sizes[6] = sizeof(double);
	event_field_sizes[7] = sizeof(double);
	event_field_sizes[8] = sizeof(unsigned int);
	event_field_sizes[9] = sizeof(unsigned int);
	event_field_sizes[10] = sizeof(unsigned int);
	event_field_sizes[11] = sizeof(unsigned int);
	
	sweep_field_names[0] = SWEEP_EVENT_NUM_HDF5;
	sweep_field_names[1] = SWEEP_NUM_HDF5;
	sweep_field_names[2] = SWEEP_BLOCK_ID_HDF5;
	sweep_field_names[3] = SWEEP_SLIP_HDF5;
	sweep_field_names[4] = SWEEP_AREA_HDF5;
	sweep_field_names[5] = SWEEP_MU_HDF5;
	sweep_field_names[6] = SWEEP_SHEAR_INIT_HDF5;
	sweep_field_names[7] = SWEEP_NORMAL_INIT_HDF5;
	sweep_field_names[8] = SWEEP_SHEAR_FINAL_HDF5;
	sweep_field_names[9] = SWEEP_NORMAL_FINAL_HDF5;
	sweep_field_offsets[0] = HOFFSET(EventSweepInfo, event_number);
	sweep_field_offsets[1] = HOFFSET(EventSweepInfo, sweep_num);
	sweep_field_offsets[2] = HOFFSET(EventSweepInfo, block_id);
	sweep_field_offsets[3] = HOFFSET(EventSweepInfo, block_slip);
	sweep_field_offsets[4] = HOFFSET(EventSweepInfo, block_area);
	sweep_field_offsets[5] = HOFFSET(EventSweepInfo, block_mu);
	sweep_field_offsets[6] = HOFFSET(EventSweepInfo, shear_init);
	sweep_field_offsets[7] = HOFFSET(EventSweepInfo, normal_init);
	sweep_field_offsets[8] = HOFFSET(EventSweepInfo, shear_final);
	sweep_field_offsets[9] = HOFFSET(EventSweepInfo, normal_final);
	sweep_field_types[0] = H5T_NATIVE_UINT;
	sweep_field_types[1] = H5T_NATIVE_UINT;
	sweep_field_types[2] = H5T_NATIVE_UINT;
	sweep_field_types[3] = H5T_NATIVE_DOUBLE;
	sweep_field_types[4] = H5T_NATIVE_DOUBLE;
	sweep_field_types[5] = H5T_NATIVE_DOUBLE;
	sweep_field_types[6] = H5T_NATIVE_DOUBLE;
	sweep_field_types[7] = H5T_NATIVE_DOUBLE;
	sweep_field_types[8] = H5T_NATIVE_DOUBLE;
	sweep_field_types[9] = H5T_NATIVE_DOUBLE;
	sweep_field_sizes[0] = sizeof(unsigned int);
	sweep_field_sizes[1] = sizeof(unsigned int);
	sweep_field_sizes[2] = sizeof(BlockID);
	sweep_field_sizes[3] = sizeof(double);
	sweep_field_sizes[4] = sizeof(double);
	sweep_field_sizes[5] = sizeof(double);
	sweep_field_sizes[6] = sizeof(double);
	sweep_field_sizes[7] = sizeof(double);
	sweep_field_sizes[8] = sizeof(double);
	sweep_field_sizes[9] = sizeof(double);
	
	aftershock_field_names[0] = AFTERSHOCK_EVT_NUM_HDF5;
	aftershock_field_names[1] = AFTERSHOCK_GEN_HDF5;
	aftershock_field_names[2] = AFTERSHOCK_MAG_HDF5;
	aftershock_field_names[3] = AFTERSHOCK_TIME_HDF5;
	aftershock_field_names[4] = AFTERSHOCK_X_HDF5;
	aftershock_field_names[5] = AFTERSHOCK_Y_HDF5;
	aftershock_field_offsets[0] = HOFFSET(AftershockBGInfo, event_number);
	aftershock_field_offsets[1] = HOFFSET(AftershockBGInfo, gen);
	aftershock_field_offsets[2] = HOFFSET(AftershockBGInfo, mag);
	aftershock_field_offsets[3] = HOFFSET(AftershockBGInfo, time);
	aftershock_field_offsets[4] = HOFFSET(AftershockBGInfo, x);
	aftershock_field_offsets[5] = HOFFSET(AftershockBGInfo, y);
	aftershock_field_types[0] = H5T_NATIVE_UINT;
	aftershock_field_types[1] = H5T_NATIVE_UINT;
	aftershock_field_types[2] = H5T_NATIVE_FLOAT;
	aftershock_field_types[3] = H5T_NATIVE_FLOAT;
	aftershock_field_types[4] = H5T_NATIVE_FLOAT;
	aftershock_field_types[5] = H5T_NATIVE_FLOAT;
	aftershock_field_sizes[0] = sizeof(unsigned int);
	aftershock_field_sizes[1] = sizeof(unsigned int);
	aftershock_field_sizes[2] = sizeof(float);
	aftershock_field_sizes[3] = sizeof(float);
	aftershock_field_sizes[4] = sizeof(float);
	aftershock_field_sizes[5] = sizeof(float);
	
	// Create dataspace for pairs of values
	dimsf[0] = 2;
	pair_val_dataspace = H5Screate_simple(1, dimsf, NULL);
#endif
}

void HDF5DataReader::getFilterDims(int &dim_x, int &dim_y) const {
	//dim_x = filter_dim_x;
	//dim_y = filter_dim_y;
}

/*!
 Adds a layer to be filtered from any data reads.
 */
/*void HDF5DataReader::addLayerFilter(const int &add_layer) {
	filter.layer_filter.insert(add_layer);
}*/

/*!
 Adds a fault to be filtered from any data reads.
 */
void HDF5DataReader::addFaultFilter(const FaultID &add_fault) {
	filter.fault_filter.insert(add_fault);
}

/*!
 Removes a set of layers from the data read filtering
 */
/*void HDF5DataReader::removeLayerFilter(const int &remove_layer) {
	filter.layer_filter.erase(remove_layer);
}*/

/*!
 Removes a set of faults from the data read filtering
 */
void HDF5DataReader::removeFaultFilter(const FaultID &remove_fault) {
	filter.fault_filter.erase(remove_fault);
}

/*!
 Set the range of years (minimum and maximum) to filter events.
 */
void HDF5DataReader::setYearFilterRange(const double &min_year, const double &max_year) {
	filter.year_range = std::make_pair(min_year, max_year);
}

/*!
 Get the range of years (minimum and maximum) to filter events.
 */
void HDF5DataReader::getYearFilterRange(double &min_year, double &max_year) const {
	min_year = filter.year_range.first;
	max_year = filter.year_range.second;
}

/*!
 For a set of blocks, returns a mapping of which fault each block is associated with.
 */
void HDF5DataReader::getFaultBlockMapping(FaultBlockMapping &fault_block_mapping, const BlockIDSet &event_blocks) const {
	/*BlockIDSet::const_iterator		it;
	FaultID							fault_id;
	
	for (it=event_blocks.begin();it!=event_blocks.end();++it) {
		fault_id = block_info[*it].fid;
		fault_block_mapping[fault_id].insert(*it);
	}*/
}

/*!
 Given the current fault and layer filters, creates a set of filtered data.
 */
void HDF5DataReader::filterData(const unsigned int &max_dim_x, const unsigned int &max_dim_y, const bool &green_log) {
#ifdef HAVE_HDF5
	/*unsigned int	bid, ibid, i, n;
	BlockInfo	binfo;
	herr_t		status;
	hsize_t		green_start[2], mem_start[2], count[2], green_stride[2], mem_stride[2];
	hid_t		green_select, mem_select;
	
	// Copy the selector for the entire dataspace
	green_select = H5Scopy(green_dataspace);
	mem_select = H5Scopy(green_dataspace);
	
	green_start[0] = green_start[1] = 0;		// start at element 0 in the file
	mem_start[0] = mem_start[1] = 0;			// start at element 0 in memory array
	count[0] = max_dim_x;						// Read max_dim elements
	count[1] = max_dim_y;
	green_stride[0] = 1;		// TODO: change this
	green_stride[1] = 1;
	mem_stride[0] = mem_stride[1] = 1;
	
	// Select the hyperslabs for the memory and file dataspace
	status = H5Sselect_hyperslab(green_select, H5S_SELECT_SET, green_start, green_stride, count, NULL);
	status = H5Sselect_hyperslab(mem_select, H5S_SELECT_SET, mem_start, mem_stride, count, NULL);
	
	if (filter_shear) delete filter_shear;
	if (filter_norm) delete filter_norm;
	
	filter_shear = new float[max_dim_x*max_dim_y];
	filter_norm = new float[max_dim_x*max_dim_y];
	filter_dim_x = max_dim_x;
	filter_dim_y = max_dim_y;
	
	status = H5Dread(green_shear_set, float_datatype, mem_select, green_select, H5P_DEFAULT, filter_shear);
	status = H5Dread(green_norm_set, float_datatype, mem_select, green_select, H5P_DEFAULT, filter_norm);
	
	H5Sclose(green_select);
	H5Sclose(mem_select);
	
	if (green_log) {
		for (i=0;i<max_dim_x*max_dim_y;++i) {
			if (filter_shear[i] != 0) filter_shear[i] = log(fabs(filter_shear[i]));
			else filter_shear[i] = 0;
			if (filter_norm[i] != 0) filter_norm[i] = log(fabs(filter_norm[i]));
			else filter_norm[i] = 0;
		}
	}*/
	/*
	filtered_dim = 0;
	for (bid=0;bid<num_blocks;++bid) {
		if (!isBlockInFaultFilter(bid) && !isBlockInLayerFilter(bid)) {
			filtered_dim++;
		}
	}
	
	if (filtered_green_shear) delete filtered_green_shear;
	if (filtered_green_normal) delete filtered_green_normal;
	filtered_green_shear = new float[filtered_dim*filtered_dim];
	filtered_green_normal = new float[filtered_dim*filtered_dim];

	// Filter the greens functions (shear and normal) with the fault and layer filters
	for (n=0,bid=0;bid<num_blocks;++bid) {
		for (ibid=0;ibid<num_blocks;++ibid) {
			if (!isBlockInFaultFilter(bid) && !isBlockInLayerFilter(bid) &&
				!isBlockInFaultFilter(ibid) && !isBlockInLayerFilter(ibid)) {
				filtered_green_shear[n] = green_shear[bid*num_blocks+ibid];
				filtered_green_normal[n] = green_normal[bid*num_blocks+ibid];
				n++;
			}
		}
	}
	*/
	// TODO: Filter the event set with the max and min range
#endif
}

/*!
 Returns true if the specified block should be displayed given the current fault filter.
 */
bool HDF5DataReader::isBlockInFaultFilter(const BlockID &block_id) const {
	BlockInfo	binfo;
	
	binfo = getBlockInfo(block_id);
	return (filter.fault_filter.count(binfo.fid) > 0);
}

/*!
 Returns true if the specified block should be displayed given the current layer filter.
 Assumes layers are of equal size and contiguous.
 */
bool HDF5DataReader::isBlockInLayerFilter(const BlockID &block_id) const {
	int			layer_size, layer_num;
	
	layer_size = num_blocks/num_layers;
	layer_num = int(block_id/layer_size);
	return (filter.layer_filter.count(layer_num) > 0);
}

/*!
 Returns the block information for the specified block ID.
 */
BlockInfo HDF5DataReader::getBlockInfo(const BlockID &block_id) const {
#ifdef HAVE_HDF5
	BlockInfo	binfo;
	herr_t		status;
	
	status = H5TBread_records(data_file, B_INFO_TABLE_HDF5, block_id, 1, sizeof(BlockInfo), binfo_field_offsets, binfo_field_sizes, &binfo);
	if (status < 0) exit(-1);
	
	return binfo;
#endif
}

/*!
 Sets the block information for the specified block.
 */
void HDF5DataWriter::setBlockInfo(const Block &block) {
#ifdef HAVE_HDF5
	BlockID			id;
	BlockInfo		binfo;
	herr_t			status;
	
	id = block.getBlockID();
	binfo = block.getBlockInfo();
	
	status = H5TBwrite_records(data_file, B_INFO_TABLE_HDF5, id, 1, sizeof(BlockInfo), binfo_field_offsets, binfo_field_sizes, &binfo);
	if (status < 0) exit(-1);
#endif
}

/*!
 Reads from shared memory a map of fault IDs with their associated name.
 */
void HDF5DataReader::getFaultNames(std::set<std::pair<FaultID, std::string> > &fault_name_map) const {
	unsigned int	i;
	std::string		fault_name;
	FaultID			fid;
	BlockInfo		binfo;
	
	for (i=0;i<num_blocks;++i) {
		binfo = getBlockInfo(i);
		fid = binfo.fid;
		fault_name = std::string(binfo.fault_name, FAULT_NAME_MAX_LEN);
		fault_name_map.insert(std::make_pair(fid, fault_name));
	}
}

/*!
 Set the start and end years of the simulation.
 */
void HDF5DataWriter::setStartEndYears(const double &new_start_year, const double &new_end_year) {
#ifdef HAVE_HDF5
	herr_t		status;
	double		tmp[2];
	
	tmp[0] = new_start_year;
	tmp[1] = new_end_year;
	status = H5Dwrite(sim_years_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
#endif
}

/*!
 Write the base latitude and longitude
 */
void HDF5DataWriter::setLatLon0(const quakelib::LatLonDepth &new_lat_lon) {
#ifdef HAVE_HDF5
	herr_t		status;
	double		tmp[2];
	
	tmp[0] = new_lat_lon.lat();
	tmp[1] = new_lat_lon.lon();
	status = H5Dwrite(base_lat_lon_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
#endif
}

void HDF5DataReader::getStartEndYears(double &start_year, double &end_year) const {
#ifdef HAVE_HDF5
	herr_t		status;
	double		tmp[2];
	
	status = H5Dread(sim_years_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
	start_year = tmp[0];
	end_year = tmp[1];
#endif
}

void HDF5DataReader::getLatLon0(double &lat0, double &lon0) const {
#ifdef HAVE_HDF5
	herr_t		status;
	double		tmp[2];
	
	status = H5Dread(base_lat_lon_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
	lat0 = tmp[0];
	lon0 = tmp[1];
#endif
}

double HDF5DataReader::getSimStart(void) const {
	double		start_year, end_year;
	
	getStartEndYears(start_year, end_year);
	
	return start_year;
}

double HDF5DataReader::getSimLength(void) const {
	double		start_year, end_year;
	
	getStartEndYears(start_year, end_year);
	
	return end_year-start_year;
}

double HDF5DataReader::getSimEnd(void) const {
	double		start_year, end_year;
	
	getStartEndYears(start_year, end_year);
	
	return end_year;
}

double HDF5DataReader::getLat0(void) const {
	double		lat0, lon0;
	
	getLatLon0(lat0, lon0);
	
	return lat0;
}

double HDF5DataReader::getLon0(void) const {
	double		lat0, lon0;
	
	getLatLon0(lat0, lon0);
	
	return lon0;
}

/*!
 Reads all events written to the event file so far.
 */
void HDF5DataReader::readAllAvailableEvents(void) {
#ifdef HAVE_HDF5
	VCEvent				*new_event;
	VCEventSweep		new_sweep;
	VCEventAftershock	new_aftershock;
	VCGeneralEvent		new_bg_event;
	EventSweeps			sweep_list;
	EventInfo			*einfo_array;
	EventSweepInfo		*swinfo_array;
	AftershockBGInfo	*asinfo_array;
	hsize_t				total_fields, total_recs;
	herr_t				status;
	unsigned int		i, n, events_to_read;
	unsigned int		sweep_num, start_sweep_rec, end_sweep_rec, num_sweeps_to_read;
	unsigned int		start_as_rec, end_as_rec, num_as_to_read;
	
	// Get the total length of the event table
	status = H5TBget_table_info(data_file, EVENT_TABLE_HDF5, &total_fields, &total_recs);
	if (status < 0) exit(-1);
	
	// Figure out how many records remain unread
	events_to_read = total_recs - last_event_read;
	
	// Read the remaining event records
	einfo_array = new EventInfo[events_to_read];
	status = H5TBread_records(data_file,
							  EVENT_TABLE_HDF5,
							  last_event_read,
							  events_to_read,
							  sizeof(EventInfo),
							  event_field_offsets,
							  event_field_sizes,
							  einfo_array);
	if (status < 0) exit(-1);
	
	// Determine which sweep records to import and read them
	start_sweep_rec = einfo_array[0].start_sweep_rec;
	end_sweep_rec = einfo_array[events_to_read-1].end_sweep_rec;
	num_sweeps_to_read = end_sweep_rec - start_sweep_rec;
	swinfo_array = new EventSweepInfo[num_sweeps_to_read];
	status = H5TBread_records(data_file,
							  SWEEP_TABLE_HDF5,
							  start_sweep_rec,
							  num_sweeps_to_read,
							  sizeof(EventSweepInfo),
							  sweep_field_offsets,
							  sweep_field_sizes,
							  swinfo_array);
	if (status < 0) exit(-1);
	
	// Determine which aftershock records to import and read them
	start_as_rec = einfo_array[0].start_aftershock_rec;
	end_as_rec = einfo_array[events_to_read-1].end_aftershock_rec;
	num_as_to_read = end_as_rec - start_as_rec;
	if (num_as_to_read > 0) {
		asinfo_array = new AftershockBGInfo[num_as_to_read];
		status = H5TBread_records(data_file,
								  AFTERSHOCK_TABLE_HDF5,
								  start_as_rec,
								  num_as_to_read,
								  sizeof(AftershockBGInfo),
								  aftershock_field_offsets,
								  aftershock_field_sizes,
								  asinfo_array);
		if (status < 0) exit(-1);
	}
	
	for (i=0;i<events_to_read;++i) {
		new_event = new VCEvent;
		new_event->setEventNumber(einfo_array[i].event_number);
		new_event->setEventYear(einfo_array[i].event_year);
		new_event->setEventTrigger(einfo_array[i].event_trigger);
		
		// Read sweeps and associate them with the event
		sweep_list.clear();
		for (n=einfo_array[i].start_sweep_rec,sweep_num=0;n<einfo_array[i].end_sweep_rec;++n) {
			new_sweep.clear();
			EventSweepInfo &cur_sweep = swinfo_array[n-start_sweep_rec];
			if (sweep_num != cur_sweep.sweep_num) sweep_num++;
			new_sweep.setSlipAndArea(cur_sweep.block_id, cur_sweep.block_slip, cur_sweep.block_area, cur_sweep.block_mu);
			new_sweep.setInitStresses(cur_sweep.block_id, cur_sweep.shear_init, cur_sweep.normal_init);
			new_sweep.setFinalStresses(cur_sweep.block_id, cur_sweep.shear_final, cur_sweep.normal_final);
			sweep_list.push_back(new_sweep);
		}
		new_event->addSweeps(sweep_list);
		
		// Read aftershocks and associate them with the event
		for (n=einfo_array[i].start_aftershock_rec;n<einfo_array[i].end_aftershock_rec;++n) {
			AftershockBGInfo &cur_as = asinfo_array[n-start_as_rec];
			VCEventAftershock new_as(cur_as.mag, cur_as.time, cur_as.x, cur_as.y, cur_as.gen);
			new_event->addAftershock(new_as);
		}
		
		event_set.insert(std::make_pair(einfo_array[i].event_year, new_event));
	}
	delete swinfo_array;
	delete einfo_array;
	
	/*	
	
	// Read the background events (one per line)
	for (i=0;i<num_bg_events;++i) {
		new_bg_event.clear();
		event_file.read(reinterpret_cast<char *>(&as_mag),sizeof(float));
		event_file.read(reinterpret_cast<char *>(&as_time),sizeof(float));
		event_file.read(reinterpret_cast<char *>(&as_x),sizeof(float));
		event_file.read(reinterpret_cast<char *>(&as_y),sizeof(float));
		
		new_bg_event = VCGeneralEvent(as_mag, as_time, as_x, as_y);
		bg_events.insert(std::make_pair(as_time, new_bg_event));
	}
	*/
#endif
}

/*!
 Write the information for an event to the HDF5 file.
 */
void HDF5DataWriter::writeEvent(VCEvent &event, VCGeneralEventSet &bg_events) {
#ifdef HAVE_HDF5
	EventSweeps::iterator		it;
	VCEventSweep::iterator		eit;
	AftershockSet::iterator		ait;
	VCGeneralEventSet::iterator	git;
	BlockIDSet					involved_blocks;
	EventInfo					e_info;
	EventSweepInfo				*s_info_array;
	AftershockBGInfo			*a_info_array;
	herr_t						status;
	unsigned int				i, sweep_num, rec_num, event_num;
	hsize_t						start_fields, end_fields, start_recs, end_recs;
	
	event_num = event.getEventNumber();
	
	// Count the number of sweep records for preallocation
	for (it=event.sweepBegin(),rec_num=0;it!=event.sweepEnd();++it) {
		for (eit=it->begin();eit!=it->end();++eit) {
			rec_num++;
		}
	}
	
	// Check the number of records before appending to the table
	status = H5TBget_table_info(data_file, SWEEP_TABLE_HDF5, &start_fields, &start_recs);
	if (status < 0) exit(-1);
	if (rec_num > 0) {
		s_info_array = new EventSweepInfo[rec_num];
		
		// Write the event sweep information
		for (it=event.sweepBegin(),i=0,sweep_num=0;it!=event.sweepEnd();++it,++sweep_num) {
			for (eit=it->begin();eit!=it->end();++eit,++i) {
				s_info_array[i].event_number = event_num;
				s_info_array[i].sweep_num = sweep_num;
				s_info_array[i].block_id = eit->first;
				s_info_array[i].block_slip = eit->second.slip;
				s_info_array[i].block_area = eit->second.area;
				s_info_array[i].block_mu = eit->second.mu;
				s_info_array[i].shear_init = eit->second.shear_init;
				s_info_array[i].shear_final = eit->second.shear_final;
				s_info_array[i].normal_init = eit->second.normal_init;
				s_info_array[i].normal_final = eit->second.normal_final;
			}
		}
		
		H5TBappend_records(data_file,
						   SWEEP_TABLE_HDF5,
						   rec_num,
						   sizeof(EventSweepInfo),
						   sweep_field_offsets,
						   sweep_field_sizes,
						   s_info_array);
		delete s_info_array;
	}
	// Get the number of records after appending to the table
	status = H5TBget_table_info(data_file, SWEEP_TABLE_HDF5, &end_fields, &end_recs);
	if (status < 0) exit(-1);
	e_info.start_sweep_rec = start_recs;
	e_info.end_sweep_rec = end_recs;
	
	// Write the aftershocks to the file
	for (ait=event.aftershockBegin(),rec_num=0;ait!=event.aftershockEnd();++ait) { rec_num++; }
	
	// Check the number of records before appending to the table
	status = H5TBget_table_info(data_file, AFTERSHOCK_TABLE_HDF5, &start_fields, &start_recs);
	if (status < 0) exit(-1);
	
	if (rec_num > 0) {
		a_info_array = new AftershockBGInfo[rec_num];
		
		for (ait=event.aftershockBegin(),i=0;ait!=event.aftershockEnd();++ait,++i) {
			a_info_array[i].event_number = event_num;
			a_info_array[i].mag = ait->mag;
			a_info_array[i].time = ait->t;
			a_info_array[i].x = ait->x;
			a_info_array[i].y = ait->y;
			a_info_array[i].gen = ait->gen;
		}
		
		H5TBappend_records(data_file,
						   AFTERSHOCK_TABLE_HDF5,
						   rec_num,
						   sizeof(AftershockBGInfo),
						   aftershock_field_offsets,
						   aftershock_field_sizes,
						   a_info_array);
		delete a_info_array;
	}
	// Get the number of records after appending to the table
	status = H5TBget_table_info(data_file, AFTERSHOCK_TABLE_HDF5, &end_fields, &end_recs);
	if (status < 0) exit(-1);
	e_info.start_aftershock_rec = start_recs;
	e_info.end_aftershock_rec = end_recs;
	
	// Write background events
	for (git=bg_events.begin(),rec_num=0;git!=bg_events.end();++git) { rec_num++; }
	
	if (rec_num > 0) {
		a_info_array = new AftershockBGInfo[rec_num];
		
		for (git=bg_events.begin(),i=0;git!=bg_events.end();++git,++i) {
			a_info_array[i].event_number = 0;
			a_info_array[i].gen = 0;
			a_info_array[i].mag = git->mag;
			a_info_array[i].time = git->t;
			a_info_array[i].x = git->x;
			a_info_array[i].y = git->y;
		}
		
		H5TBappend_records(data_file,
						   BG_EVENT_TABLE_HDF5,
						   rec_num,
						   sizeof(AftershockBGInfo),
						   aftershock_field_offsets,
						   aftershock_field_sizes,
						   a_info_array);
		delete a_info_array;
	}
	
	event.getInvolvedBlocks(involved_blocks);
	
	e_info.event_number = event.getEventNumber();
	e_info.event_year = event.getEventYear();
	e_info.event_trigger = event.getEventTrigger();
	e_info.event_magnitude = event.getMagnitude(involved_blocks);
	
	e_info.init_shear = event.getShearStressInit();
	e_info.init_normal = event.getNormalStressInit();
	e_info.final_shear = event.getShearStressFinal();
	e_info.final_normal = event.getNormalStressFinal();
	
	H5TBappend_records(data_file, EVENT_TABLE_HDF5, 1, sizeof(EventInfo), event_field_offsets, event_field_sizes, &e_info);
#endif
}

void HDF5DataWriter::flush(void) {
#ifdef HAVE_HDF5
	H5Fflush(data_file, H5F_SCOPE_GLOBAL);
#endif
}
