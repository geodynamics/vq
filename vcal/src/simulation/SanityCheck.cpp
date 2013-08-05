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

#include "SanityCheck.h"
#include <cmath>
#include <iomanip>

void SanityCheck::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# Performing sanity checks." << std::endl;
}

bool SanityCheck::assertCFFValueCorrectness(VCSimulation *sim) {
	BlockID		gid;
	int			lid;
	bool		failed = false;
	
	for (lid=0;lid<sim->numLocalBlocks();++lid) {
		gid = sim->getGlobalBID(lid);
        double val = sim->getBlock(gid).getCFF();
        if(std::isnan(val) || fabs(val) > 1e20) {
			failed_cffs.push_back(std::make_pair(gid, val));
			failed = true;
        }
    }
	return failed;
}

bool SanityCheck::assertUpdateFieldCorrectness(VCSimulation *sim) {
	BlockID		gid;
	int			lid;
	bool		failed = false;
	
	for (lid=0;lid<sim->numLocalBlocks();++lid) {
		gid = sim->getGlobalBID(lid);
		double val = sim->getUpdateField(gid);
		if(std::isnan(val) || fabs(val) > 1e20) {
			failed_update_fields.push_back(std::make_pair(gid, val));
			failed = true;
		}
    }
	return failed;
}

/*!
 Checks the local CFF and updateField values for the simulation.
 If any values are unusual they are recorded and the simulation is halted.
 */
SimRequest SanityCheck::run(SimFramework *_sim) {
	VCSimulation		*sim = static_cast<VCSimulation*>(_sim);
	if (assertCFFValueCorrectness(sim) || assertUpdateFieldCorrectness(sim)) return SIM_STOP_REQUIRED;
	else return SIM_STOP_OK;
}

/*!
 Print any unreasonable values and their associated blocks
 at the end of the simulation.
 */
void SanityCheck::finish(SimFramework *_sim) {
	std::vector<std::pair<BlockID, double> >::const_iterator	it;
	if (!failed_cffs.empty()) {
		_sim->console() << "CFFs with unreasonable values" << std::endl;
		_sim->console() << "Block ID\tValue" << std::endl;
		for (it=failed_cffs.begin();it!=failed_cffs.end();++it) {
			_sim->console() << it->first << "\t\t\t" << std::scientific << std::setw(6) << it->second << std::endl;
		}
	}
	
	if (!failed_update_fields.empty()) {
		_sim->console() << "Update fields with unreasonable values" << std::endl;
		_sim->console() << "Block ID\tValue" << std::endl;
		for (it=failed_update_fields.begin();it!=failed_update_fields.end();++it) {
			_sim->console() << it->first << "\t\t\t" << std::scientific << std::setw(6) << it->second << std::endl;
		}
	}
}
