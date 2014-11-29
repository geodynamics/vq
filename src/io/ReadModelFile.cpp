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

#include "ReadModelFile.h"

#include <iostream>
#include <iomanip>

void ReadModelFile::initDesc(const SimFramework *_sim) const {
}

/*!
 For dry runs the EqSim files are parsed as usual.
 */
void ReadModelFile::dryRun(SimFramework *_sim) {
    init(_sim);
}

/*!
 Parse the geometry file and optionally the condition and friction files.
 */
void ReadModelFile::init(SimFramework *_sim) {
    Simulation                  *sim = static_cast<Simulation *>(_sim);
    std::string                 file_name, file_type;
    quakelib::ModelWorld        world;
    quakelib::siterator         sit;
    quakelib::eiterator         eit;
    int                         err;
    quakelib::Conversion        c;

    file_name = sim->getModelFile();
    file_type = sim->getModelFileType();

    // Read the file into the world
    if (file_type == "text") {
        err = world.read_file_ascii(file_name);
    } else if (file_type == "hdf5") {
        err = world.read_file_hdf5(file_name);
    } else {
        sim->errConsole() << "ERROR: unknown file type " << file_type << std::endl;
        return;
    }

    // If there was an error then exit
    if (err) {
        sim->errConsole() << "ERROR: could not read file " << file_name << std::endl;
        return;
    }

    // Convert input world to simulation elements
    for (eit=world.begin_element(); eit!=world.end_element(); ++eit) {
        Block                   new_block;
        quakelib::SimElement    new_element;

        new_block.clear();
        new_element.clear();
        new_element = world.create_sim_element(eit->id());

        // Set SimElement values
        for (unsigned int i=0; i<3; ++i) new_block.set_vert(i, new_element.vert(i));

        new_block.set_is_quad(new_element.is_quad());
        new_block.set_rake(new_element.rake());
        new_block.set_slip_rate(new_element.slip_rate());
        new_block.set_aseismic(new_element.aseismic());
        new_block.set_lame_mu(new_element.lame_mu());
        new_block.set_lame_lambda(new_element.lame_lambda());
        new_block.set_max_slip(new_element.max_slip());

        // Set VC specific values
        //new_block.setFaultID(eit->fault_id());
        new_block.setSectionID(eit->section_id());    // TODO: add sections?

        BlockID bid = sim->addBlock(new_block);

        sim->setInitShearNormalStress(bid, 0, 0);
    }
}
