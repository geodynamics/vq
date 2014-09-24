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

#include "VCBlock.h"
#include "QuakeLibIO.h"
#include <map>
#include <set>

#ifndef _VCEVENT_H_
#define _VCEVENT_H_

// Class recording data associated with a block that slipped during an event
// mu is static, but we retain it for use in calculating magnitude.
class EventBlockData {
    public:
        double      slip, area, mu;
        double      shear_init, shear_final;
        double      normal_init, normal_final;

        EventBlockData(void) : slip(0), area(0), mu(0), shear_init(0), shear_final(0), normal_init(0), normal_final(0) {};
};

typedef std::map<BlockID, EventBlockData> EventBlockMap;
typedef std::set<BlockID> BlockIDSet;
typedef std::map<BlockID, int> BlockIDProcMapping;
typedef std::set<FaultID> FaultIDSet;
typedef std::set<SectionID> SectionIDSet;
typedef std::map<FaultID, quakelib::ElementIDSet> FaultBlockMapping;
typedef std::map<SectionID, quakelib::ElementIDSet> SectionBlockMapping;
typedef std::vector<FaultID> FaultIDList;
typedef std::vector<SectionID> SectionIDList;

/*!
 An event consists of one or more sweeps. Each sweep is a list of failed blocks and
 the slip they experienced. This gives the order that the failure propagated
 through the fault(s).
 */
class VCEventSweep {
    private:
        //! Set of block IDs and related data for this sweep (including slip, area, mu for magnitude calculation)
        EventBlockMap       block_vals;

    public:
        typedef EventBlockMap::iterator         iterator;
        typedef EventBlockMap::const_iterator   const_iterator;

        iterator begin(void) {
            return block_vals.begin();
        };
        iterator end(void) {
            return block_vals.end();
        };
        const_iterator begin(void) const {
            return block_vals.begin();
        };
        const_iterator end(void) const {
            return block_vals.end();
        };

        void setSlipAndArea(const BlockID &block_id, const double &slip, const double &area, const double &mu) {
            block_vals[block_id].slip = slip;
            block_vals[block_id].area = area;
            block_vals[block_id].mu = mu;
        };

        void setInitStresses(const BlockID &block_id, const double &shear_init, const double &normal_init) {
            block_vals[block_id].shear_init = shear_init;
            block_vals[block_id].normal_init = normal_init;
        }

        void setFinalStresses(const BlockID &block_id, const double &shear_final, const double &normal_final) {
            block_vals[block_id].shear_final = shear_final;
            block_vals[block_id].normal_final = normal_final;
        }

        double getBlockSlip(const BlockID &block_id) {
            return (block_vals.count(block_id) > 0 ? block_vals[block_id].slip : 0);
        };
        double getBlockArea(const BlockID &block_id) {
            return (block_vals.count(block_id) > 0 ? block_vals[block_id].area : 0);
        };
        double getBlockMu(const BlockID &block_id) {
            return (block_vals.count(block_id) > 0 ? block_vals[block_id].mu : 0);
        };

        void getInvolvedBlocks(BlockIDSet &block_set) const;
        BlockID getIntersectingBlock(const BlockIDSet &block_set) const;

        void clear(void) {
            block_vals.clear();
        };
        unsigned int count(const BlockID &block_id) const {
            return block_vals.count(block_id);
        };
        unsigned int size(void) {
            return block_vals.size();
        };

        friend std::ostream &operator<<(std::ostream &os, const VCEventSweep &es);
};

/*!
 Generalized events are similar to aftershocks but don't have a generation.
 These represent background events in the simulation off of main faults.
 We use floats to save space and improve speed since there can be on the
 order of 1e6 or more general events per VCEvent.
 */
class VCGeneralEvent {
    public:
        //! Magnitude of earthquake
        float mag;
        //! Time of earthquake
        float t;
        //! "x"-position (km) of earthquake
        float x;
        //! "y"-position (km) of earthquake
        float y;
        //! "z"-position (km) of earthquake
        float z;

        VCGeneralEvent(float _m=0.0, float _t=0.0, float _x=0.0, float _y=0.0, float _z=0.0) : mag(_m), t(_t), x(_x), y(_y), z(_z) {};
        void clear(void) {
            mag = t = x = y = z = 0;
        };

        quakelib::Vec<3> loc(void) const {
            return quakelib::Vec<3>(x,y,z);
        };

        // Comparison operator (for sorting)
        bool operator < (const VCGeneralEvent &as) const {
            return (this->t < as.t);
        };
        friend std::ostream &operator<<(std::ostream &os, const VCGeneralEvent &e);
};

typedef std::vector<VCGeneralEvent> VCGeneralEventSet;

/*!
 Aftershocks are events caused by a rupture. Aftershocks have a generation to indicate
 how far removed they are from the original event (aftershocks caused by aftershocks and so on).
 */
class VCEventAftershock : public VCGeneralEvent {
    public:
        unsigned int gen;  // Generation number of earthquake (0=mainshock, N=aftershock generation N)

        // Constructor
        VCEventAftershock(float _m=0.0, float _t=0.0, float _x=0.0, float _y=0.0, unsigned int _g=0)
            : VCGeneralEvent(_m, _t, _x, _y), gen(_g) { };

        void clear(void) {
            VCGeneralEvent::clear();
            gen = 0;
        };
        friend std::ostream &operator<<(std::ostream &os, const VCEventAftershock &e);
};

typedef std::vector<VCEventSweep> EventSweeps;
typedef std::vector<VCEventAftershock> AftershockVector;
typedef std::set<VCEventAftershock> AftershockSet;

/*!
 A VCEvent is a rupture of one or more blocks on one or more faults in a VC model.
 This consists of origin information (trigger block, year) and information about
 how the event propagated through the system (event sweeps).
 */
class VCEvent {
    private:
        //! The event number is used as a simple unique identifier for each event
        unsigned int    event_number;

        //! The year the event occurred
        double          event_year;

        //! Which block triggered the event through a static friction failure
        BlockID         event_trigger;

        //! Whether the event trigger block is on this node or not (used in parallel simulation)
        bool            event_trigger_on_this_node;

        //! The set of event sweeps ordered in time
        EventSweeps     event_sweeps;

        //! Sum of the EventSweeps slips, used to quickly calculate magnitude
        EventBlockMap   total_slip;

        //! Initial and final sum of shear stress on all elements involved in the event
        double shear_stress_init, shear_stress_final;

        //! Initial and final sum of normal stress on all elements involved in the event
        double normal_stress_init, normal_stress_final;

    public:
        typedef EventBlockMap::iterator iterator;
        typedef EventBlockMap::const_iterator const_iterator;

        void setEventStresses(const double &new_init_shear, const double &new_final_shear,
                              const double &new_init_normal, const double &new_final_normal) {
            shear_stress_init = new_init_shear;
            shear_stress_final = new_final_shear;
            normal_stress_init = new_init_normal;
            normal_stress_final = new_final_normal;
        }

        double getShearStressInit(void) {
            return shear_stress_init;
        };
        double getShearStressFinal(void) {
            return shear_stress_final;
        };
        double getNormalStressInit(void) {
            return normal_stress_init;
        };
        double getNormalStressFinal(void) {
            return normal_stress_final;
        };

        //! Set the unique identifier number for this event
        void setEventNumber(const unsigned int &new_event_num) {
            event_number = new_event_num;
        };
        //! Get the unique identifier number for this event
        unsigned int getEventNumber(void) const {
            return event_number;
        };

        //! Set the year this event occurred
        void setEventYear(const double &new_year) {
            event_year = new_year;
        };
        //! Get the year this event occurred
        double getEventYear(void) const {
            return event_year;
        };

        //! Set the block that triggered this event
        void setEventTrigger(const BlockID &new_trigger) {
            event_trigger = new_trigger;
        };
        //! Get the block that triggered this event
        BlockID getEventTrigger(void) const {
            return event_trigger;
        };

        //! Set whether this event occurred on this node
        void setEventTriggerOnThisNode(const bool &new_totn) {
            event_trigger_on_this_node = new_totn;
        };
        //! Get whether this event occurred on this node
        bool getEventTriggerOnThisNode(void) const {
            return event_trigger_on_this_node;
        };

        //! Add a list of event sweeps to this event.
        void addSweeps(EventSweeps &sweep_list);

        //! Get the total amount a given block slipped during this event
        double getEventSlip(const BlockID &block_id) const {
            return (total_slip.count(block_id) > 0 ? total_slip.at(block_id).slip : 0);
        }

        //! Get the area of a given block that slipped during this event
        double getEventArea(const BlockID &block_id) const {
            return (total_slip.count(block_id) > 0 ? total_slip.at(block_id).area : 0);
        }

        //! Get the value of Mu for a given block
        double getBlockMu(const BlockID &block_id) const {
            return (total_slip.count(block_id) > 0 ? total_slip.at(block_id).mu : 0);
        }

        //! Get a set of block IDs for all the blocks that failed in this event.
        void getInvolvedBlocks(BlockIDSet &block_set) const;
        //! Get the magnitude of the earthquake in this event based on the set of specified blocks.
        double getMagnitude(const BlockIDSet &involved_blocks) const;
        //! Get the magnitude of the earthquake in this event.
        double getMagnitude(void) const;

        //! Get the sum of slips for the listed blocks
        void getTotalSlipAndArea(const BlockIDSet &involved_blocks, double &slip_sum, double &area_sum) const;

        //! Order the list of sections by the order in which the rupture spread through them.
        void orderSectionsByRupture(SectionIDList &section_ordering, const SectionBlockMapping &event_sections);

        //! Order the list of faults by the order in which the rupture spread through them.
        void orderFaultsByRupture(FaultIDList &fault_ordering, const FaultBlockMapping &event_faults);
        //! Get the first block among involved_blocks that ruptured in the event.
        BlockID getHypocenter(const BlockIDSet &involved_blocks) const;

        //! Get the total number of blocks that failed in this event.
        unsigned int size(void) const {
            return total_slip.size();
        };
        //! Get the number of sweeps that occurred during this event.
        unsigned int getNumSweeps(void) const {
            return event_sweeps.size();
        };

        void clear(void);

        iterator begin(void) {
            return total_slip.begin();
        };
        iterator end(void) {
            return total_slip.end();
        };

        EventSweeps::iterator sweepBegin(void) {
            return event_sweeps.begin();
        };
        EventSweeps::iterator sweepEnd(void) {
            return event_sweeps.end();
        };

        friend std::ostream &operator<<(std::ostream &os, const VCEvent &e);
};

#endif
