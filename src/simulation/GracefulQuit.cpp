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

#include "GracefulQuit.h"

#ifdef VQ_HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef VQ_HAVE_SIGNAL_H
#include <signal.h>
#endif

GracefulQuit *gq_obj;

void catch_quit_signal(int) {
    gq_obj->mark_signal_quit();
}

bool GracefulQuit::quitFileExists(void) {
    FILE            *fp;

    if ((fp = fopen(GRACEFUL_QUIT_FILE_NAME, "r"))) {
        fclose(fp);
        return true;
    }

    return false;
}

void GracefulQuit::initDesc(const SimFramework *_sim) const {
    const Simulation  *sim = static_cast<const Simulation *>(_sim);
    sim->console() << "# To gracefully quit, create the file "
                   << GRACEFUL_QUIT_FILE_NAME << " in the run directory";
#ifdef VQ_HAVE_SIGNAL_H
    sim->console() << " or use a SIGINT (Control-C)";
#endif
    sim->console() << "." << std::endl;
}

void GracefulQuit::init(SimFramework *_sim) {
    Simulation    *sim = static_cast<Simulation *>(_sim);

    next_check_event = sim->itersPerSecond();
    signal_quit = false;

    if (sim->isRootNode()) {
        if (quitFileExists()) {
            sim->errConsole() << "ERROR: File " << GRACEFUL_QUIT_FILE_NAME <<
                              " already exists. Delete this file and run again." << std::endl;
            exit(-1);
        }

#ifdef VQ_HAVE_SIGNAL_H
        gq_obj = this;
        signal(SIGINT, catch_quit_signal);
#endif
    }
}

/*
 Have the root node check if the quit file exists, and broadcast whether
 to quit or not to other processes. This prevents the possibility of the
 file being created during the check and one process seeing it but another
 not seeing it.
 */
SimRequest GracefulQuit::run(SimFramework *_sim) {
    Simulation    *sim = static_cast<Simulation *>(_sim);
    int             quit_sim = 0, all_quit;
    unsigned int    cur_event;

    cur_event = sim->getEventCount();

    if (cur_event >= next_check_event) {
        next_check_event = cur_event + sim->itersPerSecond();

        if (sim->isRootNode()) {
            quit_sim = (quitFileExists() || signal_quit ? 1 : 0);
        }

        all_quit = sim->broadcastValue(quit_sim);

        if (all_quit) return SIM_STOP_REQUIRED;
    }

    return SIM_STOP_OK;
}
