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

#include "VCEvent.h"
#include <limits.h>
#include <iomanip>

void VCEvent::clear(void) {
    event_number = -1;
    event_year = -1;
    event_trigger = UNDEFINED_BLOCK_ID;
    event_sweeps.clear();
    total_slip.clear();
    aftershocks.clear();
}

void VCEvent::addSweeps(EventSweeps &sweep_list) {
    EventSweeps::iterator       lit;
    VCEventSweep::iterator      it;

    for (lit=sweep_list.begin(); lit!=sweep_list.end(); ++lit) {
        for (it=lit->begin(); it!=lit->end(); ++it) {
            total_slip[it->first].slip += it->second.slip;
            total_slip[it->first].area = it->second.area;
            total_slip[it->first].mu = it->second.mu;
        }
    }

    event_sweeps.insert(event_sweeps.end(), sweep_list.begin(), sweep_list.end());
}

BlockID VCEventSweep::getIntersectingBlock(const BlockIDSet &block_set) const {
    BlockIDSet::const_iterator      it;

    for (it=block_set.begin(); it!=block_set.end(); ++it) {
        if (block_vals.count(*it) > 0) return *it;
    }

    return UNDEFINED_BLOCK_ID;
}

void VCEventSweep::getInvolvedBlocks(BlockIDSet &block_id_set) const {
    EventBlockMap::const_iterator   it;

    for (it=block_vals.begin(); it!=block_vals.end(); ++it) {
        block_id_set.insert(it->first);
    }
}

void VCEvent::getInvolvedBlocks(BlockIDSet &block_id_set) const {
    EventSweeps::const_iterator it;

    for (it=event_sweeps.begin(); it!=event_sweeps.end(); ++it) {
        it->getInvolvedBlocks(block_id_set);
    }
}

double VCEvent::getMagnitude(void) const {
    BlockIDSet  block_id_set;
    getInvolvedBlocks(block_id_set);
    return getMagnitude(block_id_set);
}

double VCEvent::getMagnitude(const BlockIDSet &involved_blocks) const {
    BlockIDSet::const_iterator  it;
    double                      moment;

    moment = 0;

    for (it=involved_blocks.begin(); it!=involved_blocks.end(); ++it) {
        moment += getBlockMu(*it)*getEventSlip(*it)*getEventArea(*it);
    }

    return (2.0/3.0)*log10(1e7*moment) - 10.7;
}

void VCEvent::orderSectionsByRupture(SectionIDList &section_ordering, const SectionBlockMapping &event_sections) {
    SectionBlockMapping::const_iterator             it;
    BlockIDSet::const_iterator                      sit;
    EventSweeps::const_iterator                     eit;
    VCEventSweep::const_iterator                    slit;
    std::map<BlockID, SectionID>                    reverse_map;
    SectionIDSet                                    found_sections;
    std::multimap<int, SectionID>                   sweep_sections;
    std::multimap<int, SectionID>::const_iterator sw_it;
    int                                             sweep_num;
    SectionID                                       sid;

    // Build a reverse map from block IDs to fault IDs
    for (it=event_sections.begin(); it!=event_sections.end(); ++it) {
        for (sit=it->second.begin(); sit!=it->second.end(); ++sit) {
            reverse_map[*sit] = it->first;
        }
    }

    // Go through the sweeps
    sweep_num = 0;

    for (eit=event_sweeps.begin(); eit!=event_sweeps.end(); ++eit,++sweep_num) {
        for (slit=eit->begin(); slit!=eit->end(); ++slit) {
            // Get the fault ID of each block
            sid = reverse_map[slit->first];

            // If we haven't seen this fault before, this is the first sweep it has appeared
            if (!found_sections.count(sid)) {
                found_sections.insert(sid);
                sweep_sections.insert(std::make_pair(sweep_num, sid));
            }
        }
    }

    // Sort the faults into section_ordering by their sweep number
    for (sw_it=sweep_sections.begin(); sw_it!=sweep_sections.end(); ++sw_it) {
        section_ordering.push_back(sw_it->second);
    }
}

void VCEvent::orderFaultsByRupture(FaultIDList &fault_ordering, const FaultBlockMapping &event_faults) {
    FaultBlockMapping::const_iterator       it;
    BlockIDSet::const_iterator              sit;
    EventSweeps::const_iterator             eit;
    VCEventSweep::const_iterator            slit;
    std::map<BlockID, FaultID>              reverse_map;
    FaultIDSet                              found_faults;
    std::multimap<int, FaultID>             sweep_faults;
    std::multimap<int, FaultID>::const_iterator sw_it;
    int                                     sweep_num;
    FaultID                                 fid;

    // Build a reverse map from block IDs to fault IDs
    for (it=event_faults.begin(); it!=event_faults.end(); ++it) {
        for (sit=it->second.begin(); sit!=it->second.end(); ++sit) {
            reverse_map[*sit] = it->first;
        }
    }

    // Go through the sweeps
    sweep_num = 0;

    for (eit=event_sweeps.begin(); eit!=event_sweeps.end(); ++eit,++sweep_num) {
        for (slit=eit->begin(); slit!=eit->end(); ++slit) {
            // Get the fault ID of each block
            fid = reverse_map[slit->first];

            // If we haven't seen this fault before, this is the first sweep it has appeared
            if (!found_faults.count(fid)) {
                found_faults.insert(fid);
                sweep_faults.insert(std::make_pair(sweep_num, fid));
            }
        }
    }

    // Sort the faults into fault_ordering by their sweep number
    for (sw_it=sweep_faults.begin(); sw_it!=sweep_faults.end(); ++sw_it) {
        fault_ordering.push_back(sw_it->second);
    }
}

// Finds the hypocenter among the list of blocks, which is defined as the block
// with the earliest sweep
BlockID VCEvent::getHypocenter(const BlockIDSet &involved_blocks) const {
    EventSweeps::const_iterator     it;
    BlockID                         bid;

    for (it=event_sweeps.begin(); it!=event_sweeps.end(); ++it) {
        bid = it->getIntersectingBlock(involved_blocks);

        if (bid != UNDEFINED_BLOCK_ID) return bid;
    }

    return UNDEFINED_BLOCK_ID;
}

void VCEvent::getTotalSlipAndArea(const BlockIDSet &involved_blocks, double &slip_sum, double &area_sum) const {
    BlockIDSet::const_iterator      it;

    slip_sum = area_sum = 0;

    for (it=involved_blocks.begin(); it!=involved_blocks.end(); ++it) {
        slip_sum += total_slip.at(*it).slip;
        area_sum += total_slip.at(*it).area;
    }
}

std::ostream &operator<<(std::ostream &os, const VCEventSweep &es) {
    return os;
}

std::ostream &operator<<(std::ostream &os, const VCEventAftershock &e) {
    os << "M" << std::setprecision(2) << e.mag;
    os << " T" << e.t;
    os << " (" << std::setprecision(1) << e.x << "," << std::setprecision(1) << e.y << ")";
    return os;
}

std::ostream &operator<<(std::ostream &os, const VCEvent &e) {
    os << e.event_number << " " << e.event_year;
    return os;
}
