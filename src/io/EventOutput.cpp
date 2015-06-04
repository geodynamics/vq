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

#include "EventOutput.h"
#include "HDF5Data.h"

bool EventOutput::pauseFileExists(void) {
    FILE            *fp;

    if ((fp = fopen(PAUSE_FILE_NAME, "r"))) {
        fclose(fp);
        return true;
    }

    return false;
}

void EventOutput::initDesc(const SimFramework *_sim) const {
    const Simulation          *sim = static_cast<const Simulation *>(_sim);

    sim->console() << "# To access the event output file during the simulation, pause " << std::endl;
    sim->console() << "# by creating the file " << PAUSE_FILE_NAME ".  Delete the file to resume." << std::endl;
    sim->console() << "# Writing events in format " << sim->getEventOutfileType() << " to file " << sim->getEventOutfile() << std::endl;
}

/*!
 Initialize the HDF5 writer using the specified model dimensions.
 */
#ifdef HDF5_FOUND
void EventOutput::open_hdf5_file(const std::string &hdf5_file_name, const double &start_year, const double &end_year) {
    hid_t   plist_id;
    herr_t  status;
    double  tmp[2];
    hid_t   sim_years_set;
    hid_t   pair_val_dataspace;
    hsize_t dimsf[2];

    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    if (plist_id < 0) exit(-1);

    // Create the data file, overwriting any old files
    data_file = H5Fcreate(hdf5_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    if (data_file < 0) exit(-1);

    // Create dataspace for pairs of values
    dimsf[0] = 2;
    pair_val_dataspace = H5Screate_simple(1, dimsf, NULL);

    // Create entries for the simulation start/stop years and base longitude/latitude
    sim_years_set = H5Dcreate2(data_file, "sim_years", H5T_NATIVE_DOUBLE, pair_val_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (sim_years_set < 0) exit(-1);

    // Create the event table
    quakelib::ModelEvent::setup_event_hdf5(data_file);

    // Create the event sweeps table
    quakelib::ModelSweeps::setup_sweeps_hdf5(data_file);

    // Record the simulation start/end years
    tmp[0] = start_year;
    tmp[1] = end_year;
    status = H5Dwrite(sim_years_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
    herr_t      res;

    // Close the handles we've used
    res = H5Dclose(sim_years_set);

    if (res < 0) exit(-1);

    H5Pclose(plist_id);
}
#endif

void EventOutput::init(SimFramework *_sim) {
    Simulation                *sim = static_cast<Simulation *>(_sim);
    BlockList::const_iterator   it;

#ifdef HDF5_FOUND
    data_file = 0;
#endif
    sweep_count = 0;
    next_pause_check = sim->itersPerSecond();

    // Only the root node writes to the output file
    if (sim->isRootNode()) {
        if (sim->getEventOutfileType() == "hdf5") {
#ifdef HDF5_FOUND
            open_hdf5_file(sim->getEventOutfile(), sim->getYear(), sim->getSimDuration());
#else
            sim->errConsole() << "ERROR: HDF5 library not linked, cannot use HDF5 output files." << std::endl;
            exit(-1);
#endif
        } else if (sim->getEventOutfileType() == "text") {
            event_outfile.open(sim->getEventOutfile().c_str());
            sweep_outfile.open(sim->getSweepOutfile().c_str());

            if (!event_outfile.good() || !sweep_outfile.good()) {
                sim->errConsole() << "ERROR: Could not open output file " << sim->getEventOutfile() << std::endl;
                exit(-1);
            }

            quakelib::ModelEvent::write_ascii_header(event_outfile);
            quakelib::ModelSweeps::write_ascii_header(sweep_outfile);
        } else {
            sim->errConsole() << "ERROR: Unknown output file type " << sim->getEventOutfileType() << std::endl;
            exit(-1);
        }
    }
}

/*!
 During the simulation, this framework writes events to either an HDF5 or text file.
 */
SimRequest EventOutput::run(SimFramework *_sim) {
    Simulation        *sim = static_cast<Simulation *>(_sim);

    unsigned int num_sweeps = sim->getCurrentEvent().getSweeps().size();
    sim->getCurrentEvent().setStartEndSweep(sweep_count, sweep_count+num_sweeps);
    sweep_count += num_sweeps;

    if (sim->isRootNode()) {
        if (sim->getEventOutfileType() == "hdf5") {
#ifdef HDF5_FOUND
            sim->getCurrentEvent().append_event_hdf5(data_file);
            sim->getCurrentEvent().getSweeps().append_sweeps_hdf5(data_file);
#endif
        } else if (sim->getEventOutfileType() == "text") {
            // Write the event details
            sim->getCurrentEvent().write_ascii(event_outfile);
            event_outfile.flush();

            // Write the sweep details
            sim->getCurrentEvent().getSweeps().write_ascii(sweep_outfile);
            sweep_outfile.flush();
        }
    }

    // Check if the pause file exists each second
    if (sim->getEventCount() >= next_pause_check) {
        next_pause_check += sim->itersPerSecond();

        if (sim->isRootNode() && pauseFileExists()) {
            // Flush out the HDF5 data
#ifdef HDF5_FOUND
            H5Fflush(data_file, H5F_SCOPE_GLOBAL);
#endif
            sim->console() << "# Pausing simulation due to presence of file " << PAUSE_FILE_NAME << std::endl;

            while (pauseFileExists()) {
#ifdef VQ_HAVE_USLEEP_FUNC
                usleep(1000000);
#elif defined VQ_HAVE_SLEEP_FUNC
                sleep(1);
#endif
            }
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
    Simulation        *sim = static_cast<Simulation *>(_sim);
    //
    // debugging and runtime note: this code will execute on all nodes. however, as vq is parameterized now, if(data_file)==true only
    // on the root node, so the h5/txt actions will occur only on the head node. in the event that we allow parallel writing,
    // we may need to be more careful about flushing output from chile nodes and closing only from the root node after the child nodes finish.
#ifdef HDF5_FOUND
	//
    if (data_file) {
    	// should this only happen on the root node? is this only called by the root node?
    	//printf("**Debug(%d/%d): Closing h5 data file.\n", sim->getNodeRank(), getpid());
        herr_t res = H5Fclose(data_file);

        if (res < 0) exit(-1);
    }

#endif
    //
    if (sim->getEventOutfileType() == "text") {
        event_outfile.flush();
        event_outfile.close();
        sweep_outfile.flush();
        sweep_outfile.close();
    }
    //printf("**Debug(%d/%d): Finish EventOuput::finish().\n", sim->getNodeRank(), getpid());
}
