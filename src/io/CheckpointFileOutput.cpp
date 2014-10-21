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

#include "CheckpointFileOutput.h"
#include <sstream>

void CheckpointFileOutput::initDesc(const SimFramework *_sim) const {
    const Simulation          *sim = static_cast<const Simulation *>(_sim);

    sim->console() << "# Saving checkpoint every " << sim->getCheckpointPeriod() << " events." << std::endl;
}

void CheckpointFileOutput::writeCheckpoint(const std::string &ckpt_file_name, const Simulation *sim) {
    CheckpointSet           checkpoint_set;
    int                     i;
    BlockID                 bid;

    // Get current state for local blocks
    for (i=0; i<sim->numLocalBlocks(); ++i) {
        bid = sim->getGlobalBID(i);
        checkpoint_set[bid] = sim->getBlock(bid).state.readCheckpointData();
    }

#ifdef HDF5_FOUND
    HDF5CheckpointWriter checkpoint_file(ckpt_file_name,
                                         sim->numGlobalBlocks(),
                                         sim->getYear(),
                                         sim->getEventCount(),
                                         checkpoint_set);
#endif
}

SimRequest CheckpointFileOutput::run(SimFramework *_sim) {
    Simulation            *sim = static_cast<Simulation *>(_sim);
    std::stringstream       ss;

    // Periodically checkpoint simulation state to a file
    if (sim->getCheckpointPeriod() > 0 && sim->getEventCount() % sim->getCheckpointPeriod() == 0) {
        ss << sim->getCheckpointPrefix() << checkpoint_num << ".h5";
        writeCheckpoint(ss.str(), sim);
        checkpoint_num++;
    }

    return SIM_STOP_OK;
}

void CheckpointFileOutput::finish(SimFramework *_sim) {
    Simulation            *sim = static_cast<Simulation *>(_sim);
    std::stringstream       ss;

    ss << sim->getCheckpointPrefix() << "final" << ".h5";
    writeCheckpoint(ss.str(), sim);
}
