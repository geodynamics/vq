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

#include "HDF5DataShare.h"
#include "HDF5Data.h"

bool HDF5DataShare::pauseFileExists(void) {
    FILE            *fp;

    if ((fp = fopen(HDF5_PAUSE_FILE_NAME, "r"))) {
        fclose(fp);
        return true;
    }

    return false;
}

void HDF5DataShare::initDesc(const SimFramework *_sim) const {
    const VCSimulation          *sim = static_cast<const VCSimulation *>(_sim);

    sim->console() << "# To access the HDF5 file during the simulation, pause it " << std::endl;
    sim->console() << "# by creating the file " << HDF5_PAUSE_FILE_NAME ".  Delete the file to resume." << std::endl;
}

void HDF5DataShare::init(SimFramework *_sim) {
    VCSimulation                *sim = static_cast<VCSimulation *>(_sim);
    BlockList::const_iterator   it;

    h5_data = NULL;
    next_pause_check = sim->itersPerSecond();

    // Only the root node writes to the output file
    if (sim->isRootNode()) {
        h5_data = new HDF5DataWriter(sim->getHDF5File(), sim->numGlobalBlocks());

        // Set the start and end years of the simulation
        h5_data->setStartEndYears(sim->getYear(), sim->getSimDuration());

        // Record block information
        for (it=sim->begin(); it!=sim->end(); ++it) h5_data->setBlockInfo(*it);
    }
}

/*!
 During the simulation, this framework just writes events to the HDF5 file.
 */
SimRequest HDF5DataShare::run(SimFramework *_sim) {
    VCSimulation        *sim = static_cast<VCSimulation *>(_sim);

    if (h5_data) h5_data->writeEvent(sim->getCurrentEvent(), sim->getBGEvents());

    // Check if the pause file exists each second
    if (sim->getEventCount() >= next_pause_check) {
        next_pause_check += sim->itersPerSecond();

        if (sim->isRootNode() && pauseFileExists()) {
            // Flush out the HDF5 data
            h5_data->flush();
            sim->console() << "# Pausing simulation due to presence of file " << HDF5_PAUSE_FILE_NAME << std::endl;

            while (pauseFileExists()) sleep(1);
        }

        // Other processes wait until the pause file is gone
        sim->barrier();
    }

    return SIM_STOP_OK;
}

/*!
 Finish by unattaching from the shared memory structure.
 */
void HDF5DataShare::finish(SimFramework *_sim) {
    if (h5_data) delete h5_data;
}
