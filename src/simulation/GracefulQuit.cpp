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
#include <stdlib.h>

bool GracefulQuit::quitFileExists(void) {
	FILE			*fp;
	if ((fp = fopen(GRACEFUL_QUIT_FILE_NAME, "r"))) {
		fclose(fp);
		return true;
	}
	return false;
}

void GracefulQuit::initDesc(const SimFramework *_sim) const {
	const VCSimulation	*sim = static_cast<const VCSimulation*>(_sim);
	sim->console() << "# To gracefully quit, create the file "
	<< GRACEFUL_QUIT_FILE_NAME << " in the run directory." << std::endl;
}

void GracefulQuit::init(SimFramework *_sim) {
	VCSimulation	*sim = static_cast<VCSimulation*>(_sim);
	
	next_check_event = sim->itersPerSecond();
	
	if (sim->isRootNode()) {
		if (quitFileExists()) {
			sim->errConsole() << "ERROR: File " << GRACEFUL_QUIT_FILE_NAME <<
			" already exists. Delete this file and run again." << std::endl;
			exit(-1);
		}
	}
}

/*
 Have the root node check if the quit file exists, and broadcast whether
 to quit or not to other processes. This prevents the possibility of the
 file being created during the check and one process seeing it but another
 not seeing it.
 */
SimRequest GracefulQuit::run(SimFramework *_sim) {
	VCSimulation	*sim = static_cast<VCSimulation*>(_sim);
	int				quit_sim = 0, all_quit;
	unsigned int	cur_event;
	
	cur_event = sim->getEventCount();
	if (cur_event >= next_check_event) {
		next_check_event = cur_event + sim->itersPerSecond();
		
		if (sim->isRootNode()) {
			quit_sim = (quitFileExists() ? 1 : 0);
		}
		
		all_quit = sim->broadcastValue(quit_sim);
		if (all_quit) return SIM_STOP_REQUIRED;
	}
	
	return SIM_STOP_OK;
}
