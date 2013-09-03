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
#define GREEN_SHEAR_HDF5			"greens_shear"
#define GREEN_NORMAL_HDF5			"greens_normal"

// HDF5 file data definitions
#define SIM_YEARS_HDF5				"sim_years"
#define BASE_LAT_LON_HDF5			"base_lat_lon"

// Block info related definitions
#define B_INFO_TABLE_HDF5               "block_info_table"
#define B_INFO_NUM_ENTRIES_HDF5         32
#define B_INFO_BLOCK_ID_HDF5            "block_id"
#define B_INFO_FAULT_ID_HDF5            "fault_id"
#define B_INFO_SECTION_ID_HDF5          "section_id"
#define B_INFO_M_X_PT1_HDF5             "m_x_pt1"
#define B_INFO_M_Y_PT1_HDF5             "m_y_pt1"
#define B_INFO_M_Z_PT1_HDF5             "m_z_pt1"
#define B_INFO_M_DAS_PT1_HDF5           "m_das_pt1"
#define B_INFO_M_TRACE_FLAG_PT1_HDF5    "m_trace_flag_pt1"
#define B_INFO_M_X_PT2_HDF5             "m_x_pt2"
#define B_INFO_M_Y_PT2_HDF5             "m_y_pt2"
#define B_INFO_M_Z_PT2_HDF5             "m_z_pt2"
#define B_INFO_M_DAS_PT2_HDF5           "m_das_pt2"
#define B_INFO_M_TRACE_FLAG_PT2_HDF5    "m_trace_flag_pt2"
#define B_INFO_M_X_PT3_HDF5             "m_x_pt3"
#define B_INFO_M_Y_PT3_HDF5             "m_y_pt3"
#define B_INFO_M_Z_PT3_HDF5             "m_z_pt3"
#define B_INFO_M_DAS_PT3_HDF5           "m_das_pt3"
#define B_INFO_M_TRACE_FLAG_PT3_HDF5    "m_trace_flag_pt3"
#define B_INFO_M_X_PT4_HDF5             "m_x_pt4"
#define B_INFO_M_Y_PT4_HDF5             "m_y_pt4"
#define B_INFO_M_Z_PT4_HDF5             "m_z_pt4"
#define B_INFO_M_DAS_PT4_HDF5           "m_das_pt4"
#define B_INFO_M_TRACE_FLAG_PT4_HDF5    "m_trace_flag_pt4"
#define B_INFO_SLIP_VELOCITY_HDF5       "slip_velocity"
#define B_INFO_ASEISMICITY_HDF5         "aseismicity"
#define B_INFO_RAKE_HDF5                "rake_rad"
#define B_INFO_DIP_HDF5                 "dip_rad"
#define B_INFO_DYNAMIC_STRENGTH_HDF5	"dynamic_strength"
#define B_INFO_STATIC_STRENGTH_HDF5		"static_strength"
#define B_INFO_LAME_MU_HDF5				"lame_mu"
#define B_INFO_LAME_LAMBDA_HDF5			"lame_lambda"
#define B_INFO_FAULT_NAME_HDF5          "fault_name"

// Event info related definitions
#define EVENT_TABLE_HDF5			"event_table"
#define EVENT_NUM_ENTRIES_HDF5		12
#define EVENT_NUM_HDF5				"event_number"
#define EVENT_YEAR_HDF5				"event_year"
#define EVENT_TRIGGER_HDF5			"event_trigger"
#define EVENT_MAGNITUDE_HDF5		"event_magnitude"
#define EVENT_SHEAR_INIT_HDF5		"event_shear_init"
#define EVENT_NORMAL_INIT_HDF5		"event_normal_init"
#define EVENT_SHEAR_FINAL_HDF5		"event_shear_final"
#define EVENT_NORMAL_FINAL_HDF5		"event_normal_final"
#define EVENT_START_SWEEP_HDF5		"start_sweep_rec"
#define EVENT_END_SWEEP_HDF5		"end_sweep_rec"
#define EVENT_START_AS_HDF5			"start_aftershock_rec"
#define EVENT_END_AS_HDF5			"end_aftershock_rec"

// Event sweeps table definitions
#define SWEEP_TABLE_HDF5			"event_sweep_table"
#define SWEEP_NUM_ENTRIES_HDF5		10
#define SWEEP_EVENT_NUM_HDF5		"event_number"
#define SWEEP_NUM_HDF5				"sweep_num"
#define SWEEP_BLOCK_ID_HDF5			"block_id"
#define SWEEP_SLIP_HDF5				"slip"
#define SWEEP_AREA_HDF5				"area"
#define SWEEP_MU_HDF5				"mu"
#define SWEEP_SHEAR_INIT_HDF5       "shear_init"
#define SWEEP_SHEAR_FINAL_HDF5      "shear_final"
#define SWEEP_NORMAL_INIT_HDF5      "normal_init"
#define SWEEP_NORMAL_FINAL_HDF5     "normal_final"

// Aftershock/background table definitions
#define AFTERSHOCK_TABLE_HDF5		"aftershock_table"
#define AFTERSHOCK_NUM_ENTRIES_HDF5	6
#define AFTERSHOCK_EVT_NUM_HDF5		"event_number"
#define AFTERSHOCK_GEN_HDF5			"generation"
#define AFTERSHOCK_MAG_HDF5			"magnitude"
#define AFTERSHOCK_TIME_HDF5		"time"
#define AFTERSHOCK_X_HDF5			"x"
#define AFTERSHOCK_Y_HDF5			"y"

// State checkpoint table definitions
#define CHECKPOINT_STATE_HDF5		"checkpoint_state"
#define CHECKPOINT_YEAR_HDF5		"checkpoint_year"
#define CHECKPOINT_EVENT_HDF5		"checkpoint_event"
#define CHECKPOINT_NUM_ENTRIES_HDF5	6
#define CHECKPOINT_SLIP_DFCT_HDF5	"slipDeficit"
#define CHECKPOINT_SLIP_CUM_HDF5	"slipCumulative"
#define CHECKPOINT_CFF_HDF5			"cff"
#define CHECKPOINT_STRESS_S_HDF5	"stressS"
#define CHECKPOINT_STRESS_N_HDF5	"stressN"
#define CHECKPOINT_UPDATE_FLD_HDF5	"updateField"

#define BG_EVENT_TABLE_HDF5			"background_event_table"

#define FILE_NAME_MAX_LEN			256

struct EventInfo {
	unsigned int	event_number;
	double			event_year;
	BlockID			event_trigger;
	double			event_magnitude;
	unsigned int	start_sweep_rec, end_sweep_rec;
	unsigned int	start_aftershock_rec, end_aftershock_rec;
    double          init_shear, final_shear, init_normal, final_normal;
};

typedef struct EventInfo EventInfo;

struct EventSweepInfo {
	unsigned int	event_number;
	unsigned int	sweep_num;
	BlockID			block_id;
	double			block_slip;
	double			block_area;
	double			block_mu;
    double          shear_init, shear_final;
	double			normal_init, normal_final;
};

typedef struct EventSweepInfo EventSweepInfo;

struct AftershockBGInfo {
	unsigned int	event_number;
	unsigned int	gen;
	float			mag;
	float			time;
	float			x;
	float			y;
};

typedef struct AftershockBGInfo AftershockBGInfo;

// Classes representing a file containing checkpoint data
class HDF5Checkpoint {
protected:
#ifdef HDF5_FOUND
	// HDF5 handle to checkpoint data file
	hid_t				data_file;
	
	// Dataspace of checkpoint state
	hid_t				state_dataspace, year_event_dataspace;
	
	// Handle to checkpoint state dataset
	hid_t				state_dataset, year_dataset, event_dataset;
	
	// Access control handle
	hid_t				plist_id;
#endif
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
#ifdef HDF5_FOUND
	// HDF5 handle to data file
	hid_t				data_file;
	
	// Handles to Greens data in the file
	hid_t				green_norm_set, green_shear_set;
	
	// Handles to data space specifications
	hid_t				green_dataspace;
	
	// Dimension of Greens matrix
	unsigned int		greens_dim;
#endif
	
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
#ifdef HDF5_FOUND
	// HDF5 handle to data file
	hid_t				data_file;
	
	// Handles to data in the file
	hid_t				sim_years_set, base_lat_lon_set;
	
	// Handles to data types
	hid_t				fault_name_datatype;
	
	// Handles to data space specifications
	hid_t				block_info_dataspace, pair_val_dataspace, fault_name_dataspace;
	
	// Names, types, offsets and sizes for block info table
	const char *binfo_field_names[B_INFO_NUM_ENTRIES_HDF5];
	size_t binfo_field_offsets[B_INFO_NUM_ENTRIES_HDF5];
	hid_t binfo_field_types[B_INFO_NUM_ENTRIES_HDF5];
	size_t binfo_field_sizes[B_INFO_NUM_ENTRIES_HDF5];
	
	// Names, types, offsets and sizes for event table
	const char *event_field_names[EVENT_NUM_ENTRIES_HDF5];
	size_t event_field_offsets[EVENT_NUM_ENTRIES_HDF5];
	hid_t event_field_types[EVENT_NUM_ENTRIES_HDF5];
	size_t event_field_sizes[EVENT_NUM_ENTRIES_HDF5];
	
	// Names, types, offsets and sizes for event sweep table
	const char *sweep_field_names[SWEEP_NUM_ENTRIES_HDF5];
	size_t sweep_field_offsets[SWEEP_NUM_ENTRIES_HDF5];
	hid_t sweep_field_types[SWEEP_NUM_ENTRIES_HDF5];
	size_t sweep_field_sizes[SWEEP_NUM_ENTRIES_HDF5];
	
	// Names, types, offsets and sizes for aftershock/background event table
	const char *aftershock_field_names[AFTERSHOCK_NUM_ENTRIES_HDF5];
	size_t aftershock_field_offsets[AFTERSHOCK_NUM_ENTRIES_HDF5];
	hid_t aftershock_field_types[AFTERSHOCK_NUM_ENTRIES_HDF5];
	size_t aftershock_field_sizes[AFTERSHOCK_NUM_ENTRIES_HDF5];
#endif
	
	// Values read from shared memory (pointers are set to within shared memory segment)
	unsigned int	num_blocks, num_layers;
	
	void createH5Handles(void);
	
public:
	HDF5Data(void) {};
	~HDF5Data(void);
	unsigned int modelDim(void) const { return num_blocks; };
};

typedef std::pair<const double, VCEvent*> EventYear;
typedef std::multimap<double, VCEvent*> EventYearMap;
typedef std::map<double, VCGeneralEvent*> BGEventMap;

/*!
 Represents a filter for GUI information based on faults, layers and time.
 */
class GUIFilter {
public:
	// Sets of which fault IDs and layers to filter out of the Greens functions and events
	FaultIDSet					fault_filter;
	std::set<int>				layer_filter;
	
	// Range of years to return events for
	std::pair<double,double>	year_range;
};

class HDF5DataReader : public HDF5Data {
private:
	//! Current offset into the events file, so new events can immediately be read when they are appended
	unsigned int				last_event_read;
	
	//! Set of events indexed by year
	EventYearMap				event_set;
	
	//! Set of background events indexed by year
	BGEventMap					bg_events;
	
	//! Filter for events/Greens function
	GUIFilter					filter;
	
	void getStartEndYears(double &start_year, double &end_year) const;
	void getLatLon0(double &lat0, double &lon0) const;
	
public:
	typedef EventYearMap::iterator			iterator;
	typedef EventYearMap::const_iterator	const_iterator;
	
	iterator begin(const double &lower_bound=-1) {
		if (lower_bound >= 0) return event_set.lower_bound(lower_bound);
		else return event_set.lower_bound(filter.year_range.first);
	};
	iterator end(const double &upper_bound=-1) {
		if (upper_bound >= 0) return event_set.upper_bound(upper_bound);
		else return event_set.upper_bound(filter.year_range.second);
	};
	
	BGEventMap::iterator bgbegin(const double &lower_bound=-1) {
		if (lower_bound >= 0) return bg_events.lower_bound(lower_bound);
		else return bg_events.lower_bound(filter.year_range.first);
	};
	BGEventMap::iterator bgend(const double &upper_bound=-1) {
		if (upper_bound >= 0) return bg_events.upper_bound(upper_bound);
		else return bg_events.upper_bound(filter.year_range.second);
	};
	
	HDF5DataReader(const std::string &hdf5_file_name);
	
	void getFilterDims(int &dim_x, int &dim_y) const;
	
	void filterData(const unsigned int &max_dim_x, const unsigned int &max_dim_y, const bool &green_log);
	//void addLayerFilter(const int &add_layer);
	//void removeLayerFilter(const int &remove_layer);
	void addFaultFilter(const FaultID &add_fault);
	void removeFaultFilter(const FaultID &remove_fault);
	void setYearFilterRange(const double &min_year, const double &max_year);
	void getYearFilterRange(double &min_year, double &max_year) const;
	
	bool isBlockInFaultFilter(const BlockID &block_id) const;
	bool isBlockInLayerFilter(const BlockID &block_id) const;
	BlockInfo getBlockInfo(const BlockID &block_id) const;
	
	void getFaultNames(std::set<std::pair<FaultID, std::string> > &fault_name_map) const;
	void getFaultBlockMapping(FaultBlockMapping &fault_block_mapping, const BlockIDSet &event_blocks) const;
	unsigned int getNumLayers(void) const { return num_layers; };
	unsigned int getNumEvents(void) const { return event_set.size(); };
	double getSimStart(void) const;
	double getSimLength(void) const;
	double getSimEnd(void) const;
	double getLat0(void) const;
	double getLon0(void) const;
	
	void readAllAvailableEvents(void);
};

class HDF5DataWriter : public HDF5Data {
public:
	HDF5DataWriter(const std::string &hdf5_file_name, const int &nblocks);
	void setBlockInfo(const Block &block);
	void setStartEndYears(const double &new_start_year, const double &new_end_year);
	void setLatLon0(const quakelib::LatLonDepth &new_lat_lon);
	void flush(void);
	
	void writeEvent(VCEvent &event, VCGeneralEventSet &bg_events);
};

#endif
