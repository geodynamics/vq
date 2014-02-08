// Copyright (c) 2012-2013 Eric M. Heien, Michael K. Sachs, John B. Rundle
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

#include "EventOutput.h"
#include "HDF5Data.h"

#include <unistd.h>

bool EventOutput::pauseFileExists(void) {
    FILE            *fp;

    if ((fp = fopen(PAUSE_FILE_NAME, "r"))) {
        fclose(fp);
        return true;
    }

    return false;
}

void EventOutput::initDesc(const SimFramework *_sim) const {
    const VCSimulation          *sim = static_cast<const VCSimulation *>(_sim);

    sim->console() << "# To access the event output file during the simulation, pause " << std::endl;
    sim->console() << "# by creating the file " << PAUSE_FILE_NAME ".  Delete the file to resume." << std::endl;
    sim->console() << "# Writing events in format " << sim->getEventOutfileType() << " to file " << sim->getEventOutfile() << std::endl;
}

void EventOutput::init(SimFramework *_sim) {
    VCSimulation                *sim = static_cast<VCSimulation *>(_sim);
    BlockList::const_iterator   it;

    h5_data = NULL;
    next_pause_check = sim->itersPerSecond();

    // Only the root node writes to the output file
    if (sim->isRootNode()) {
        if (sim->getEventOutfileType() == "hdf5") {
#ifdef HDF5_FOUND
            h5_data = new HDF5DataWriter(sim->getEventOutfile());

            // Set the start and end years of the simulation
            h5_data->setStartEndYears(sim->getYear(), sim->getSimDuration());
#else
            std::cerr << "ERROR: HDF5 library not linked, cannot use HDF5 output files." << std::endl;
            exit(-1);
#endif
        } else if (sim->getEventOutfileType() == "text") {
            event_outfile.open(sim->getEventOutfile().c_str());
            sweep_outfile.open(sim->getSweepOutfile().c_str());

            if (!event_outfile.good() || !sweep_outfile.good()) {
                std::cerr << "ERROR: Could not open output file " << sim->getEventOutfile() << std::endl;
                exit(-1);
            }
        } else {
            std::cerr << "ERROR: Unknown output file type " << sim->getEventOutfileType() << std::endl;
            exit(-1);
        }
    }
}

/*!
 During the simulation, this framework writes events to either an HDF5 or text file.
 */
SimRequest EventOutput::run(SimFramework *_sim) {
    VCSimulation        *sim = static_cast<VCSimulation *>(_sim);

    if (sim->getEventOutfileType() == "hdf5") {
#ifdef HDF5_FOUND
        h5_data->writeEvent(sim->getCurrentEvent());
#endif
    } else if (sim->getEventOutfileType() == "text") {
        VCEvent &event = sim->getCurrentEvent();
        BlockIDSet                  involved_blocks;
        EventSweeps::iterator       it;
        VCEventSweep::iterator      eit;
        unsigned int                i, sweep_num, rec_num;

        // Count the number of sweep records
        for (it=event.sweepBegin(),rec_num=0; it!=event.sweepEnd(); ++it) {
            for (eit=it->begin(); eit!=it->end(); ++eit) {
                rec_num++;
            }
        }

        // Write the event details
        event.getInvolvedBlocks(involved_blocks);
        event_outfile << event.getEventNumber() << " " << event.getEventYear() << " ";
        event_outfile << event.getEventTrigger() << " " << event.getMagnitude(involved_blocks) << " ";
        event_outfile << event.getShearStressInit() << " " << event.getNormalStressInit() << " ";
        event_outfile << event.getShearStressFinal() << " " << event.getNormalStressFinal() << "\n";
        event_outfile.flush();

        // Write the sweeps
        for (it=event.sweepBegin(),i=0,sweep_num=0; it!=event.sweepEnd(); ++it,++sweep_num) {
            for (eit=it->begin(); eit!=it->end(); ++eit,++i) {
                sweep_outfile << event.getEventNumber() << " " << sweep_num << " " << eit->first << " " << eit->second.slip << " ";
                sweep_outfile << eit->second.area << " " << eit->second.mu << " ";
                sweep_outfile << eit->second.shear_init << " " << eit->second.shear_final << " ";
                sweep_outfile << eit->second.normal_init << " " << eit->second.normal_final << "\n";
            }
        }

        sweep_outfile.flush();
    }

    // Check if the pause file exists each second
    if (sim->getEventCount() >= next_pause_check) {
        next_pause_check += sim->itersPerSecond();

        if (sim->isRootNode() && pauseFileExists()) {
            // Flush out the HDF5 data
#ifdef HDF5_FOUND
            h5_data->flush();
#endif
            sim->console() << "# Pausing simulation due to presence of file " << PAUSE_FILE_NAME << std::endl;

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
void EventOutput::finish(SimFramework *_sim) {
    VCSimulation        *sim = static_cast<VCSimulation *>(_sim);

#ifdef HDF5_FOUND

    if (h5_data) delete h5_data;

#endif

    if (sim->getEventOutfileType() == "text") {
        event_outfile.flush();
        event_outfile.close();
        sweep_outfile.flush();
        sweep_outfile.close();
    }
}
