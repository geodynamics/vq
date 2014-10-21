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

#include "SanityCheck.h"

void SanityCheck::initDesc(const SimFramework *_sim) const {
    const Simulation          *sim = static_cast<const Simulation *>(_sim);

    sim->console() << "# Performing sanity checks." << std::endl;
}

bool SanityCheck::assertCFFValueCorrectness(Simulation *sim) {
    BlockID     gid;
    int         lid;
    bool        failed = false;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        double val = sim->getBlock(gid).getCFF();

        if (std::isnan(val) || fabs(val) > 1e20) {
            failed_cffs.push_back(std::make_pair(gid, val));
            failed = true;
        }
    }

    return failed;
}

bool SanityCheck::assertUpdateFieldCorrectness(Simulation *sim) {
    BlockID     gid;
    int         lid;
    bool        failed = false;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        double val = sim->getUpdateField(gid);

        if (std::isnan(val) || fabs(val) > 1e20) {
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
    Simulation        *sim = static_cast<Simulation *>(_sim);

    if (assertCFFValueCorrectness(sim) || assertUpdateFieldCorrectness(sim)) return SIM_STOP_REQUIRED;
    else return SIM_STOP_OK;
}

/*!
 Print any unreasonable values and their associated blocks
 at the end of the simulation.
 */
void SanityCheck::finish(SimFramework *_sim) {
    std::vector<std::pair<BlockID, double> >::const_iterator    it;

    if (!failed_cffs.empty()) {
        _sim->console() << "CFFs with unreasonable values" << std::endl;
        _sim->console() << "Block ID\tValue" << std::endl;

        for (it=failed_cffs.begin(); it!=failed_cffs.end(); ++it) {
            _sim->console() << it->first << "\t\t\t" << std::scientific << std::setw(6) << it->second << std::endl;
        }
    }

    if (!failed_update_fields.empty()) {
        _sim->console() << "Update fields with unreasonable values" << std::endl;
        _sim->console() << "Block ID\tValue" << std::endl;

        for (it=failed_update_fields.begin(); it!=failed_update_fields.end(); ++it) {
            _sim->console() << it->first << "\t\t\t" << std::scientific << std::setw(6) << it->second << std::endl;
        }
    }
}
