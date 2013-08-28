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
