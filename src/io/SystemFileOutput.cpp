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
