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

//#ifdef HAVE_CONFIG_H
#include "config.h"
//#endif

#include "SystemFileOutput.h"
#include "SimError.h"
#include <fstream>
#include <iomanip>

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

void SystemFileOutput::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# System output file: " << sim->getSystemOutfile() << std::endl;
}

void SystemFileOutput::init(SimFramework *_sim) {
	VCSimulation				*sim = static_cast<VCSimulation*>(_sim);
	std::string					file_name = sim->getSystemOutfile();
	std::ofstream				out_file(file_name.c_str());
	BlockList::const_iterator	it;
	Block						tmp_block;
	char						time_buf[256];
	struct tm					*cur_date;
	time_t						cur_time;
	
	// Only output the system info at the root node
	if (!_sim->isRootNode()) return;
	
	// Check that the state file opened correctly
	assertThrow(out_file.good(), "Couldn't open file " + file_name + ".");
	
	// Get the current formatted time
	cur_time = time(NULL);
	cur_date = localtime(&cur_time);
	strftime(time_buf, 255, "%F %T", cur_date);
	time_buf[255] = 0;
	
	// Write the system file header information
	// TODO: Add SVN revision number
	out_file << "#####################################################\n";
#ifdef HAVE_CONFIG_H
	out_file << "# Virtual California " << VERSION << "\n";
#else
	out_file << "# Virtual California\n";
#endif
	out_file << "# Simulation start: " << time_buf << "\n";
	out_file << "#####################################################\n";
	out_file << "#" << tmp_block.header() << "\n";
	out_file << "#" << tmp_block.headerUnits() << "\n";
	
	// Write the block information
	for (it=sim->begin();it!=sim->end();++it) {
		out_file << *it << "\n";
	}
	
	out_file.flush();
	out_file.close();
}
