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

#include "ProgressMonitor.h"
#include <iomanip>
#include <numeric>

void ProgressMonitor::initDesc(const SimFramework *_sim) const {
    const Simulation          *sim = static_cast<const Simulation *>(_sim);

    sim->console() << "# Displaying simulation progress every ";

    if (sim->getProgressPeriod() == 1) sim->console() << "second." << std::endl;
    else sim->console() << sim->getProgressPeriod() << " seconds." << std::endl;
}

/*!
 Print initial information about the simulation environment and model.
 */
void ProgressMonitor::init(SimFramework *_sim) {
    Simulation    *sim = static_cast<Simulation *>(_sim);
    BlockVal        min, avr, max;
    int             width = 30;

    getStats(sim, min, avr, max);

    next_event_count = sim->itersPerSecond();

    if (!sim->isRootNode()) return;

    sim->console() << "#" << std::endl;
    sim->console() << "# **********************************************" << std::endl;
    sim->console() << "# ***" << std::endl;
    sim->console() << std::setw(width) << std::left << "# *** Blocks" << ": " << sim->numGlobalBlocks() << std::endl;
    sim->console() << std::setw(width) << std::left << "# *** Faults" << ": " << sim->numFaults() << std::endl;
    sim->console() << std::setw(width) << std::left << "# ***" << std::endl;
    sim->console() << std::setw(width) << std::left << "# *** Present Time (years)" << ": " << sim->getYear() << std::endl;
    sim->console() << std::setw(width) << std::left << "# *** Min  cff" << ": " << min.val << std::endl;
    sim->console() << std::setw(width) << std::left << "# *** Mean cff" << ": " << avr.val << std::endl;
    sim->console() << std::setw(width) << std::left << "# *** Max  cff" << ": " << max.val << std::endl;
    sim->console() << std::setw(width) << std::left << "# ***" << std::endl;
    sim->console() << std::right;
}

/*!
 Gather statistics over all computing nodes about the CFF values.
 */
void ProgressMonitor::getStats(Simulation *sim, BlockVal &min, BlockVal &avg, BlockVal &max) {
    int         lid;
    BlockID     gid;
    BlockVal    min_val, max_val, sum_val;
    double      cur_cff;

    min_val.val = DBL_MAX;
    max_val.val = -DBL_MAX;
    sum_val.val = 0;
    min_val.block_id = max_val.block_id = sum_val.block_id = UNDEFINED_ELEMENT_ID;

    // Get the minimum, maximum and sum CFF values on this node
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        cur_cff = sim->getCFF(gid);

        if (cur_cff < min_val.val) {
            min_val.val = cur_cff;
            min_val.block_id = gid;
        }

        sum_val.val += cur_cff;

        if (cur_cff > max_val.val) {
            max_val.val = cur_cff;
            max_val.block_id = gid;
        }
    }

    // Reduce to the min/avg/max over all nodes
    sim->allReduceBlockVal(min_val, min, BLOCK_VAL_MIN);
    sim->allReduceBlockVal(max_val, max, BLOCK_VAL_MAX);
    sim->allReduceBlockVal(sum_val, avg, BLOCK_VAL_SUM);
    avg.val /= sim->numGlobalBlocks();
}

/*!
 Periodically (roughly every getProgressPeriod() seconds) print simulation information.
 This includes the number of events and current stress values for the system.
 */
SimRequest ProgressMonitor::run(SimFramework *_sim) {
    Simulation            *sim = static_cast<Simulation *>(_sim);
    BlockVal                min, avg, max;
    int                     max_bid_width;
    std::ios_base::fmtflags fmt_flags;

    // Print the header on the first time
    if (first) {
        sim->console() << "# events        year      minCff[index]      avrCff      maxCff[index]" << std::endl;
        first = false;
    }

    if (sim->getEventCount() >= next_event_count) {
        // Gather simulation-wide statistics
        getStats(sim, min, avg, max);

        // Estimate how many events until next progress report and notify other processes
        next_event_count += (sim->getProgressPeriod()*sim->itersPerSecond())+1;

        // Write the current progress
        max_bid_width = (int)(log10(sim->numGlobalBlocks())+1);

        fmt_flags = sim->console().flags();
        sim->console()  << std::fixed << std::setw(8) << sim->getEventCount()       // number of events
                        << std::setw(12) << std::setprecision(1) << sim->getYear()  // current simulation year
                        << std::scientific << std::setw(17-max_bid_width) << std::setprecision(3) << min.val
                        << "[" << std::setw(max_bid_width) << min.block_id << "]"   // min cff
                        << std::setw(12) << std::setprecision(3) << avg.val         // average cff
                        << std::setw(17-max_bid_width) << std::setprecision(3) << max.val
                        << "[" << std::setw(max_bid_width) << max.block_id << "]"   // max cff
                        << std::endl;
        sim->console().flags(fmt_flags);
    }

    return SIM_STOP_OK;
}
