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

#include "QuakeLib.h"

#include <map>

#ifndef _VCBLOCK_H_
#define _VCBLOCK_H_

// Whether to store the matrix in transpose format
// This enables much faster calculation for sparse vector multiplication with no change in accuracy
#define GREEN_MATRIX_TRANSPOSE

#define GREEN_VAL       double

class VCRand {
    public:
        enum {NTAB=32, INT_SIZE = NTAB+3};
        struct State {
            int idum1;
            int idum2;
            int iy;
            int iv[NTAB];
        };

        VCRand(int seed = 17) {
            init(seed);
        }
        VCRand(const VCRand::State &state);
        void init(int seed);

        double nextDouble();
        int nextInt(int max) {
            return ((int)(max*nextDouble()))%max;
        }

        void save(VCRand::State &state);
        void save(int *stream);
        void init(int *stream);
    private:
        int idum1;
        int idum2;
        int iy;
        int iv[NTAB];
};

typedef unsigned int BlockID;
typedef unsigned int FaultID;
typedef unsigned int SectionID;

#define UNDEFINED_BLOCK_ID      UINT_MAX
#define UNDEFINED_FAULT_ID      UINT_MAX
#define UNDEFINED_SECTION_ID    UINT_MAX
#define ALL_FAULTS_ID           UINT_MAX-1

// TODO: add VCRand?
struct StateCheckpointData {
    double      slipDeficit;
    double      cff;
    double      stressS;
    double      stressN;
    double      updateField;
};

typedef struct StateCheckpointData StateCheckpointData;

typedef std::map<BlockID, StateCheckpointData> CheckpointSet;

struct BlockVal {
    double      val;
    BlockID     block_id;
};

typedef struct BlockVal BlockVal;

enum BlockValOp {
    BLOCK_VAL_UNDEFINED,
    BLOCK_VAL_MIN,          // Get the minimum value and associated block
    BLOCK_VAL_MAX,          // Get the maximum value and associated block
    BLOCK_VAL_SUM,          // Get the sum over all blocks
};

struct BlockSweepVals {
    double      slip, init_shear, init_normal, final_shear, final_normal;
    BlockID     block_id;
};

typedef struct BlockSweepVals BlockSweepVals;

std::ostream &operator<<(std::ostream &os, const BlockVal &bv);

/*!
 Class to keep dynamic block attributes. Not all of them are needed to save between simulations,
 only slipDeficit, but they are saved anyway to provide debugging and other quickly accessible info.
 */
class State {
    public:
        State(void) : slipDeficit(0), cff(0), cff0(0), stressS0(0), stressN0(0), stressS(NULL), stressN(NULL), updateField(NULL) {};
        //! slip - Vt
        double slipDeficit;
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
        //! ptr to shear stress in simulation
        double *stressS;
        //! ptr to normal stress in simulation
        double *stressN;
        //! ptr to update_field in simulation
        double *updateField;

        //! To have identical mpi and serial results
        VCRand          rand;

        StateCheckpointData readCheckpointData(void) const;
        void storeCheckpointData(const StateCheckpointData &ckpt_data);

        friend std::ostream &operator<<(std::ostream &os, const State &b);
};

/*!
 Class of block attributes. These attributes are static during the simulation except rand.
*/
class Block : public quakelib::SimElement {
    private:
        BlockID         id;                 // block id
        FaultID         fid;                // fault id
        SectionID       sid;                // section id

        // constants during simulation
        double          self_shear;         // Greens function shear by block on itself
        double          self_normal;        // Greens function normal by block on itself
        double          dynamic_val;        // dynamic failure value for this block

        double          stress_drop;        // predetermined stress drop (in pascals)

        double          init_shear_stress;  // initial shear stress
        double          init_normal_stress; // initial normal stress

        double          rhogd;              // normal stress

        double          max_slip;           // Maximum slip of this block based on whatever scaling relationship we're using
        double          friction_val;       // friction value for this block

        double          dynamic_strength;
        double          static_strength;

        std::string     fault_name;
        bool            failed;

    public:

        void get_rake_and_normal_stress_due_to_block(double stresses[2], const Block &source_block) const;

        double getDynamicStrength(void) const {
            return dynamic_strength;
        };
        double getStaticStrength(void) const {
            return static_strength;
        };
        void setStrengths(double new_dynamic, double new_static) {
            dynamic_strength = new_dynamic;
            static_strength = new_static;
        };

        State           state;

        //! Resets the values of this block.
        void clear(void);

        //! Calculates the static friction of this block.
        void calcFriction(void);
        //! Get the static friction of this block.
        double friction(void) const {
            return friction_val;
        };

        //! Whether the block experienced static friction failure.
        //! This occurs if the Coulomb failure function goes over 0.
        bool cffFailure(void) const;
        //! Whether the block experienced dynamic friction failure.
        //! This occurs if the block is on the same fault as the trigger
        //! and it undergoes a significant CFF shift during this event
        //! (where significant is defined in terms of dynamic_val).
        bool dynamicFailure(const FaultID &event_fault) const;
        bool getFailed(void) const {
            return failed;
        };
        void setFailed(bool in_failed) {
            failed = in_failed;
        };

        //! Returns the precalculated CFF of this block.
        double getCFF(void) const {
            return state.cff;
        };
        //! Returns the slip deficit of this block.
        double getSlipDeficit(void) const {
            return state.slipDeficit;
        };
        //! Returns the slip deficit of this block.
        double getShearStress(void) const {
            return state.stressS[0];
        };
        //! Returns the slip deficit of this block.
        double getNormalStress(void) const {
            return state.stressN[0];
        };
        //! Calculates and stores the CFF of this block.
        void calcCFF(void);
        //! Record the current stresses and CFF.
        void saveStresses(void) {
            state.stressS0 = state.stressS[0];
            state.stressN0 = state.stressN[0];
            state.cff0 = state.cff;
            state.slipDeficit0 = state.slipDeficit;
        };
        //! Restore the saved stresses and CFF.
        void restoreStresses(void) {
            *(state.stressS) = state.stressS0;
            *(state.stressN) = state.stressN0;
            state.cff = state.cff0;
            state.slipDeficit = state.slipDeficit0;
        };
        //! Get the recorded shear stress.
        double getStressS0(void) const {
            return state.stressS0;
        };
        //! Get the recorded normal stress.
        double getStressN0(void) const {
            return state.stressN0;
        };
        //! Get the recorded CFF.
        double getCFF0(void) const {
            return state.cff0;
        };
        //! Set pointers into arrays for the stress values of this block
        void setStatePtrs(double *stressS, double *stressN, double *update_field);

        // Functions for manipulation of block parameters
        //! Set the block ID of this block.
        void setBlockID(const BlockID &new_id) {
            id = new_id;
        };
        //! Return the block ID of this block.
        BlockID getBlockID(void) const {
            return id;
        };

        //! Set the fault ID of this block.
        void setFaultID(const FaultID &new_fid) {
            fid = new_fid;
        };
        //! Return the fault ID of this block.
        FaultID getFaultID(void) const {
            return fid;
        };

        //! Set the section ID of this block.
        void setSectionID(const SectionID &new_sid) {
            sid = new_sid;
        };
        //! Return the section ID of this block.
        SectionID getSectionID(void) const {
            return sid;
        };

        double getMaxSlip(void) const {
            return max_slip;
        };
        void setMaxSlip(const double &new_max_slip) {
            max_slip = new_max_slip;
        };
        //! Calculate the expected recurrence of this block in years.
        double getRecurrence(void) const {
            return stress_drop/(_slip_rate*(self_shear-self_normal*friction_val));
        };
        //! Set the Greens function self shear of this block.
        void setSelfStresses(const double &new_self_shear, const double &new_self_normal) {
            self_shear = new_self_shear;
            self_normal = new_self_normal;
        };
        //! Get the Greens function self shear of this block.
        double getSelfStresses(void) const {
            return self_shear - friction()*self_normal;
        };

        //! Get the dynamic triggering value for this block.
        double getDynamicVal(void) const {
            return dynamic_val;
        };
        //! Set the dynamic triggering value for this block.
        void setDynamicVal(const double &new_dyn_val) {
            dynamic_val = new_dyn_val;
        };

        //! Set the initial shear stress of this block in bars.
        void setInitShearStress(const double &new_stress) {
            init_shear_stress = new_stress;
        };
        //! Get the initial shear stress of this block in bars.
        double getInitShearStress(void) const {
            return init_shear_stress;
        };
        //! Set the initial normal stress of this block in bars.
        void setInitNormalStress(const double &new_stress) {
            init_normal_stress = new_stress;
        };
        //! Get the initial normal stress of this block in bars.
        double getInitNormalStress(void) const {
            return init_normal_stress;
        };

        //! Get the stress drop for this block in bars.
        double getStressDrop(void) const {
            return stress_drop;
        };
        //! Set the stress drop for this block in bars.
        void setStressDrop(const double &new_stress_drop) {
            stress_drop = new_stress_drop;
            calcFriction();
        };
        double getRhogd(void) const {
            return rhogd;
        };
        void setRhogd(const double &new_rhogd) {
            rhogd = new_rhogd;
            calcFriction();
        };

        //! Get the fault name for this block.
        std::string getFaultName(void) const {
            return fault_name;
        };
        //! Set the fault name for this block.
        void setFaultName(const std::string &new_fault_name) {
            fault_name = new_fault_name;
        };

        //! Get the string header for state files.
        static std::string getStateHeader(void) {
            return "       cff updateField   stressS   stressN";
        };
};

typedef std::vector<Block> BlockList;
typedef std::vector<BlockID> BlockIDList;
typedef std::map<BlockID, int> BlockIDMap;

std::ostream &operator<<(std::ostream &os, const Block &b);

#endif
