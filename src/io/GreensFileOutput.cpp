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

#include "GreensFileOutput.h"
#include "SimError.h"
#include <fstream>

void GreensFileOutput::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# Greens output file: " << sim->getGreensOutfile() << std::endl;
}

/*!
 Writes Greens function values to the specified file.
 This is written to an HDF5 style file, so multiple processes
 can write in parallel.
 */
void GreensFileOutput::init(SimFramework *_sim) {
	VCSimulation			*sim = static_cast<VCSimulation*>(_sim);
	std::string				file_name = sim->getGreensOutfile();
	unsigned int			green_dim;
	BlockID					row, col, global_row;
	double					*shear_vals, *norm_vals;
	HDF5GreensDataWriter	*h5_greens_data;
	
	// Open the file
	green_dim = sim->numGlobalBlocks();
	h5_greens_data = new HDF5GreensDataWriter(file_name, green_dim);
	
	// Record the Greens function values
	shear_vals = new double[green_dim];
	norm_vals = new double[green_dim];
	for (row=0;row<sim->numLocalBlocks();++row) {
		global_row = sim->getGlobalBID(row);
		for (col=0;col<green_dim;++col) {
			shear_vals[col] = sim->getGreenShear(global_row, col);
			norm_vals[col] = sim->getGreenNormal(global_row, col);
		}
		h5_greens_data->setGreensVals(global_row, shear_vals, norm_vals);
	}
	delete shear_vals;
	delete norm_vals;
	delete h5_greens_data;
}
