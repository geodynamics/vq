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

#include "VCBlock.h"
#include <map>
#include <set>

#ifndef _VCEVENT_H_
#define _VCEVENT_H_

// Class recording data associated with a block that slipped during an event
// mu is static, but we retain it for use in calculating magnitude.
class EventBlockData {
public:
	double		slip, area, mu;
	double		shear_init, shear_final;
	double		normal_init, normal_final;
	
	EventBlockData(void) : slip(0), area(0), mu(0), shear_init(0), shear_final(0), normal_init(0), normal_final(0) {};
};

typedef std::map<BlockID, EventBlockData> EventBlockMap;
typedef std::set<BlockID> BlockIDSet;
typedef std::set<FaultID> FaultIDSet;
typedef std::set<quakelib::SectionID> SectionIDSet;
typedef std::map<FaultID, BlockIDSet> FaultBlockMapping;
typedef std::map<FaultID, double> FaultFailureAreaMapping;
typedef std::map<quakelib::SectionID, BlockIDSet> SectionBlockMapping;
typedef std::vector<FaultID> FaultIDList;
typedef std::vector<quakelib::SectionID> SectionIDList;

/*!
 An event consists of one or more sweeps. Each sweep is a list of failed blocks and
 the slip they experienced. This gives the order that the failure propagated
 through the fault(s).
 */
class VCEventSweep {
private:
	//! Set of block IDs and related data for this sweep (including slip, area, mu for magnitude calculation)
	EventBlockMap		block_vals;
	
public:
	typedef EventBlockMap::iterator			iterator;
	typedef EventBlockMap::const_iterator	const_iterator;
	
	iterator begin(void) { return block_vals.begin(); };
	iterator end(void) { return block_vals.end(); };
	const_iterator begin(void) const { return block_vals.begin(); };
	const_iterator end(void) const { return block_vals.end(); };
	
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
	
    double getBlockSlip(const BlockID &block_id) { return (block_vals.count(block_id) > 0 ? block_vals[block_id].slip : 0); };
	double getBlockArea(const BlockID &block_id) { return (block_vals.count(block_id) > 0 ? block_vals[block_id].area : 0); };
	double getBlockMu(const BlockID &block_id) { return (block_vals.count(block_id) > 0 ? block_vals[block_id].mu : 0); };

	void getInvolvedBlocks(BlockIDSet &block_set) const;
	BlockID getIntersectingBlock(const BlockIDSet &block_set) const;
	
	void clear(void) { block_vals.clear(); };
	unsigned int count(const BlockID &block_id) const { return block_vals.count(block_id); };
	unsigned int size(void) { return block_vals.size(); };
	
	friend std::ostream& operator<<(std::ostream& os, const VCEventSweep& es);
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
	
	VCGeneralEvent(float _m=0.0, float _t=0.0, float _x=0.0, float _y=0.0) : mag(_m), t(_t), x(_x), y(_y) {};
	void clear(void) { mag = t = x = y = 0; };
	
	// Comparison operator (for sorting)
	bool operator < (const VCGeneralEvent &as) const { return (this->t < as.t); };
	friend std::ostream& operator<<(std::ostream& os, const VCGeneralEvent& e);
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
	
	void clear(void) { VCGeneralEvent::clear(); gen = 0; };
	friend std::ostream& operator<<(std::ostream& os, const VCEventAftershock& e);
};

typedef std::vector<VCEventSweep> EventSweeps;
typedef std::vector<VCEventAftershock> AftershockSet;

/*!
 A VCEvent is a rupture of one or more blocks on one or more faults in a VC model.
 This consists of origin information (trigger block, year) and information about
 how the event propagated through the system (event sweeps). It also contains
 the set of aftershocks associated with this event which may not be on any fault.
 */
class VCEvent {
private:
	//! The event number is used as a simple unique identifier for each event
	unsigned int	event_number;
	
	//! The year the event occurred
	double			event_year;
	
	//! Which block triggered the event through a static friction failure
	BlockID			event_trigger;
	
	//! Whether the event trigger block is on this node or not (used in parallel simulation)
	bool			event_trigger_on_this_node;
	
	//! The set of event sweeps ordered in time
	EventSweeps		event_sweeps;
	
	//! Sum of the EventSweeps slips, used to quickly calculate magnitude
	EventBlockMap	total_slip;
	
	//! Generated event aftershocks
	AftershockSet	aftershocks;
    
    double shear_stress_init, shear_stress_final;
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
    
    double getShearStressInit(void) { return shear_stress_init; };
    double getShearStressFinal(void) { return shear_stress_final; };
    double getNormalStressInit(void) { return normal_stress_init; };
    double getNormalStressFinal(void) { return normal_stress_final; };
    
	//! Set the unique identifier number for this event
	void setEventNumber(const unsigned int &new_event_num) { event_number = new_event_num; };
	//! Get the unique identifier number for this event
	unsigned int getEventNumber(void) const { return event_number; };
	
	//! Set the year this event occurred
	void setEventYear(const double &new_year) { event_year = new_year; };
	//! Get the year this event occurred
	double getEventYear(void) const { return event_year; };
	
	//! Set the block that triggered this event
	void setEventTrigger(const BlockID &new_trigger) { event_trigger = new_trigger; };
	//! Get the block that triggered this event
	BlockID getEventTrigger(void) const { return event_trigger; };
	
	//! Set whether this event occurred on this node
	void setEventTriggerOnThisNode(const bool &new_totn) { event_trigger_on_this_node = new_totn; };
	//! Get whether this event occurred on this node
	bool getEventTriggerOnThisNode(void) const { return event_trigger_on_this_node; };
	
	//! Add a list of event sweeps to this event.
	void addSweeps(EventSweeps &sweep_list);
	//! Add an aftershock to this event.
	void addAftershock(VCEventAftershock &aftershock) { aftershocks.push_back(aftershock); };
	//! Get a pointer to the set of aftershocks associated with this event.
	AftershockSet *getAftershockPtr(void) { return &aftershocks; };
	
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
	unsigned int size(void) const { return total_slip.size(); };
	//! Get the number of sweeps that occurred during this event.
	unsigned int getNumSweeps(void) const { return event_sweeps.size(); };
	//! Get the number of aftershocks resulting from this event.
	unsigned int getNumAftershocks(void) const { return aftershocks.size(); };
	
	void clear(void);
	
	iterator begin(void) { return total_slip.begin(); };
	iterator end(void) { return total_slip.end(); };
	
	EventSweeps::iterator sweepBegin(void) { return event_sweeps.begin(); };
	EventSweeps::iterator sweepEnd(void) { return event_sweeps.end(); };

	AftershockSet::iterator aftershockBegin(void) { return aftershocks.begin(); };
	AftershockSet::iterator aftershockEnd(void) { return aftershocks.end(); };
	
	friend std::ostream& operator<<(std::ostream& os, const VCEvent& e);
};

#endif
