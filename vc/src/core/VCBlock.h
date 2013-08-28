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

#include "QuakeLib.h"

#ifndef _VCBLOCK_H_
#define _VCBLOCK_H_

// Whether to store the matrix in transpose format
// This enables much faster calculation for sparse vector multiplication with no change in accuracy
#define GREEN_MATRIX_TRANSPOSE

#define GREEN_VAL		double

class VCRand {
public:
    enum {NTAB=32, INT_SIZE = NTAB+3};
    struct State {
        int idum1;
        int idum2;
        int iy;
        int iv[NTAB];
    };
	
    VCRand(int seed = 17) { init(seed); }
    VCRand(const VCRand::State& state);
    void init(int seed);
	
    double nextDouble();
    int nextInt(int max) { return ((int)(max*nextDouble()))%max; }
	
    void save(VCRand::State& state);
    void save(int* stream);
    void init(int* stream);
private:
    int idum1;
    int idum2;
    int iy;
    int iv[NTAB];
};

typedef	unsigned int BlockID;
typedef unsigned int FaultID;

#define UNDEFINED_BLOCK_ID		UINT_MAX
#define UNDEFINED_FAULT_ID		UINT_MAX
#define UNDEFINED_SECTION_ID	UINT_MAX
#define ALL_FAULTS_ID			UINT_MAX-1

#define FAULT_NAME_MAX_LEN		256

struct BlockInfo {
	BlockID		bid;
	FaultID		fid;
    quakelib::SectionID	sid;
	double		x_pt[4], y_pt[4], z_pt[4], das_pt[4];
	quakelib::TraceFlag   trace_flag_pt[4];
	double		slip_velocity;
	double		aseismicity;
    double		rake;
    double		dip;
	double		dynamic_strength, static_strength;
	double		lame_mu, lame_lambda;
	char		fault_name[FAULT_NAME_MAX_LEN];
};

typedef struct BlockInfo BlockInfo;

// TODO: add VCRand?
struct StateCheckpointData {
	double		slipDeficit;
	double		slipCumulative;
	double		cff;
	double		stressS;
	double		stressN;
	double		updateField;
};

typedef struct StateCheckpointData StateCheckpointData;

typedef std::map<BlockID, StateCheckpointData> CheckpointSet;

struct BlockVal {
	double		val;
	BlockID		block_id;
};

typedef struct BlockVal BlockVal;

enum BlockValOp {
	BLOCK_VAL_UNDEFINED,
	BLOCK_VAL_MIN,			// Get the minimum value and associated block
	BLOCK_VAL_MAX,			// Get the maximum value and associated block
	BLOCK_VAL_SUM,			// Get the sum over all blocks
};

struct BlockSweepVals {
	double		slip, init_shear, init_normal, final_shear, final_normal;
	BlockID		block_id;
};

typedef struct BlockSweepVals BlockSweepVals;

std::ostream& operator<<(std::ostream& os, const BlockVal& bv);

/*!
 Class to keep dynamic block attributes. Not all of them are needed to save between simulations,
 only slipDeficit, but they are saved anyway to provide debugging and other quickly accessible info.
 */
class State {
public:
	State(void) : slipDeficit(0), slipCumulative(0), cff(0), Fcff(0), Fcff0(0),
		 cff0(0), stressS0(0), stressN0(0), FstressS0(0), FstressN0(0), stressS(NULL), stressN(NULL), FstressS(NULL), FstressN(NULL), updateField(NULL) {};
	//! slip - Vt
    double slipDeficit;
	//! cumulative slip over all times.
    double slipCumulative;
    //! coulomb failure function for failed block
	double Fcff;
    //! coulomb failure function for failed block after a failure
	double Fcff0;
	//! coulomb failure function
	double cff;
	//! cff before an event
	double cff0;
	//! shear stress before an event
	double stressS0;
	//! normal stress before an event
	double stressN0;
	//! slip deficit before an event
	double slipDeficit0;
    //! shear stress before an event
	double FstressS0;
	//! normal stress before an event
	double FstressN0;
	//! ptr to shear stress in simulation
	double *stressS;
	//! ptr to normal stress in simulation
	double *stressN;
    double *FstressS;
	//! ptr to normal stress in simulation
	double *FstressN;
	//! ptr to update_field in simulation
	double *updateField;
	
	//! To have identical mpi and serial results
    VCRand			rand;
	
	StateCheckpointData readCheckpointData(void) const;
	void storeCheckpointData(const StateCheckpointData &ckpt_data);
	
    friend std::ostream& operator<<(std::ostream& os, const State& b);
};

/*!
 Class of block attributes. These attributes are static during the simulation except rand.
*/
class Block : public quakelib::Element<4> {
private:
    BlockID			id;					// block id
    FaultID			fid;				// fault id
    quakelib::SectionID		sid;				// section id
    
    // constants during simulation
	double			self_shear;			// Greens function shear by block on itself
	double			self_normal;		// Greens function normal by block on itself
	double			dynamic_val;		// dynamic failure value for this block
	double			active_dynamic_val;	// active dynamic failure value for this block. This is 1 until a neighbor fails. then
                                        // it takes on dynamic val
	int				slip_scaling_threshold; // the slip scaling threshold for this block
    
	double			stress_drop;		// predetermined stress drop (in pascals)
	
	double			init_shear_stress;	// initial shear stress
	double			init_normal_stress;	// initial normal stress
	
	double			rhogd;				// normal stress
	
	double			friction_val;		// friction value for this block
    
    double          section_layers;
	
	double			dynamic_strength;
	double			static_strength;
	
	std::string		fault_name;
	bool			is_valid;
    bool            failed;
    bool            failed_this_sweep;
    
    int             fault_size;
    
    /*!
     New properties defined for the updated okada functions
     */
    
    double			lambda;
    double			mu;
    double			unit_slip;
	
public:
    
    double getUnitSlip() const {return unit_slip;}
    void setUnitSlip(double inunit_slip) {unit_slip = inunit_slip;}
    
    double getLambda() const {return lambda;}
    double getMu() const {return mu;}
    void setLambda(double in_lambda) {lambda = in_lambda;}
    void setMu(double in_mu) {mu = in_mu;}
    
    void get_rake_and_normal_stress_due_to_block(double stresses[2], const Block &source_block) const;
	
	double getDynamicStrength(void) const { return dynamic_strength; };
	double getStaticStrength(void) const { return static_strength; };
	void setStrengths(double new_dynamic, double new_static) { dynamic_strength = new_dynamic; static_strength = new_static; };
	
    void dynamicOn(void) {active_dynamic_val = dynamic_val;};
    void dynamicOff(void) {active_dynamic_val = 1.0;};
	
    State			state;
	
	//! Resets the values of this block.
	void clear(void);
	//! Whether the values of this block are valid
	bool valid(void) const { return is_valid; };
	
    //! Calculates the static friction of this block.
	void calcFriction(void);
	//! Get the static friction of this block.
	double friction(void) const { return friction_val; };
    
    void setSectionLayers(int in_layers) { section_layers = (double)in_layers; };
	double getSectionLayers(void) const { return section_layers; };
    
	//! Whether the block experienced static friction failure.
	//! This occurs if the Coulomb failure function goes over 0.
	bool cffFailure(void) const;
	//! Whether the block experienced dynamic friction failure.
	//! This occurs if the block is on the same fault as the trigger
	//! and it undergoes a significant CFF shift during this event
	//! (where significant is defined in terms of dynamic_val).
	bool dynamicFailure(const FaultID &event_fault) const;
    bool additionalFailure(void) const;
    bool getFailed(void) const { return failed; };
    void setFailed(bool in_failed) { failed = in_failed; };
    bool getFailedThisSweep(void) const { return failed_this_sweep; };
    void setFailedThisSweep(bool in_failed_this_sweep) { failed_this_sweep = in_failed_this_sweep; };
	
	//! Returns the precalculated CFF of this block.
	double getCFF(void) const { return state.cff; };
    //! Returns the precalculated FCFF of this block.
	double getFCFF(void) const { return state.Fcff; };
    //! Returns the slip deficit of this block.
	double getSlipDeficit(void) const { return state.slipDeficit; };
    //! Returns the cumulative slip of this block.
	double getSlipCumulative(void) const { return state.slipCumulative; };
    //! Returns the slip deficit of this block.
	double getShearStress(void) const { return state.stressS[0]; };
    //! Returns the slip deficit of this block.
	double getFShearStress(void) const { return state.FstressS[0]; };
    //! Returns the slip deficit of this block.
	double getNormalStress(void) const { return state.stressN[0]; };
	//! Calculates and stores the CFF of this block.
	void calcCFF(bool in_event);
	//! Record the current stresses and CFF.
	void saveStresses(void) { state.stressS0 = state.stressS[0]; state.stressN0 = state.stressN[0]; state.cff0 = state.cff; state.slipDeficit0 = state.slipDeficit; };
	//! Restore the saved stresses and CFF.
	void restoreStresses(void) { *(state.stressS) = state.stressS0; *(state.stressN) = state.stressN0; state.cff = state.cff0; state.slipDeficit = state.slipDeficit0; };
	//! Record the current Fstresses and FCFF.
	void saveFStresses(void) { state.FstressS0 = state.FstressS[0]; state.FstressN0 = state.FstressN[0]; state.Fcff0 = state.Fcff; };
	//! Get the recorded shear stress.
	double getStressS0(void) const { return state.stressS0; };
	//! Get the recorded normal stress.
	double getStressN0(void) const { return state.stressN0; };
	//! Get the recorded CFF.
	double getCFF0(void) const { return state.cff0; };
	//! Set pointers into arrays for the stress values of this block
	void setStatePtrs(double *stressS, double *stressN, double *FstressS, double *FstressN, double *update_field);
	
	// Functions for manipulation of block parameters
	//! Set the block ID of this block.
	void setBlockID(const BlockID &new_id) { id = new_id; };
	//! Return the block ID of this block.
	BlockID getBlockID(void) const { return id; };
	
	//! Set the fault ID of this block.
	void setFaultID(const FaultID &new_fid) { fid = new_fid; };
	//! Return the fault ID of this block.
	FaultID getFaultID(void) const { return fid; };
    //! returns the number of elements in the block's fault
    int getFaultSize(void) const { return fault_size; };
    void setFaultSize(const int new_fault_size) { fault_size = new_fault_size; };
    
    //! Set the section ID of this block.
	void setSectionID(const quakelib::SectionID &new_sid) { sid = new_sid; };
	//! Return the section ID of this block.
	quakelib::SectionID getSectionID(void) const { return sid; };
    
	//! Calculate the expected recurrence of this block in years.
	double getRecurrence(void) const { return stress_drop/(_slip_rate*(self_shear-self_normal*friction_val)); };
	//! Set the Greens function self shear of this block.
	void setSelfStresses(const double &new_self_shear, const double &new_self_normal) {
		self_shear = new_self_shear;
		self_normal = new_self_normal;
	};
	//! Get the Greens function self shear of this block.
	double getSelfStresses(void) const { return self_shear - friction()*self_normal; };
	
	//! Get the dynamic triggering value for this block.
	double getDynamicVal(void) const { return dynamic_val; };
	//! Set the dynamic triggering value for this block.
	void setDynamicVal(const double &new_dyn_val) { dynamic_val = new_dyn_val; };
    //! Get the active dynamic triggering value for this block.
	double getActiveDynamicVal(void) const { return active_dynamic_val; };
	//! Set the active dynamic triggering value for this block.
	void setActiveDynamicVal(const double &new_dyn_val) { active_dynamic_val = new_dyn_val; };
	
	//! Get the slip scaling value for this block.
	int getSlipScalingThreshold(void) const { return slip_scaling_threshold; };
	//! Set the slip scaling value for this block.
	void setSlipScalingThreshold(const double &new_st_val) { slip_scaling_threshold = new_st_val; };
	
	//! Set the initial shear stress of this block in bars.
	void setInitShearStress(const double &new_stress) { init_shear_stress = new_stress; };
	//! Get the initial shear stress of this block in bars.
	double getInitShearStress(void) const { return init_shear_stress; };
	//! Set the initial normal stress of this block in bars.
	void setInitNormalStress(const double &new_stress) { init_normal_stress = new_stress; };
	//! Get the initial normal stress of this block in bars.
	double getInitNormalStress(void) const { return init_normal_stress; };
	
	//! Get the stress drop for this block in bars.
	double getStressDrop(void) const { return stress_drop; };
	//! Set the stress drop for this block in bars.
	void setStressDrop(const double &new_stress_drop) { stress_drop = new_stress_drop; calcFriction(); };
    double getRhogd(void) const { return rhogd; };
	void setRhogd(const double &new_rhogd) { rhogd = new_rhogd; calcFriction(); };
	
	//! Get the fault name for this block.
	std::string getFaultName(void) const { return fault_name; };
	//! Set the fault name for this block.
	void setFaultName(const std::string &new_fault_name) { fault_name = new_fault_name; };
	
	//! Whether the values of this block are valid or not.
	bool isValid(void) const { return is_valid; };
	//! Set whether the block values are valid.
	void setValid(const bool &new_valid) { is_valid = new_valid; };
	
	//! Get the string header for state files.
	static std::string getStateHeader(void) { return "       cff updateField   stressS   stressN"; };
	std::string header(void) const;
	std::string	headerUnits(void) const;
	
	BlockInfo getBlockInfo(void) const;
};

typedef std::vector<Block> BlockList;
typedef std::vector<BlockID> BlockIDList;
typedef std::map<BlockID, int> BlockIDMap;

std::ostream& operator<<(std::ostream& os, const Block& b);

#endif
