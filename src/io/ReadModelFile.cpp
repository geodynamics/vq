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
    std::string                 file_name, file_type, stress_filename, stress_file_type, stress_index_filename;
    quakelib::ModelWorld        world;
    quakelib::ModelStressSet    stress_set;
    quakelib::ModelStress       stress;
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

    stress_file_type = sim->getStressInfileType();
    stress_filename = sim->getStressInfile();
    stress_index_filename = sim->getStressIndexInfile();
    
    if (stress_filename != "" && stress_file_type != "" && stress_index_filename != "") {
        // Read the stress input file for initial stress conditions
        if (stress_file_type == "text") {
            err |= stress_set.read_file_ascii(stress_index_filename, stress_filename);
        } else if (stress_file_type == "hdf5") {
            // TODO: Schultz, add hdf5 reading
            sim->errConsole() << "ERROR: only text files supported right now " << file_type << std::endl;
            return;
        } else {
            sim->errConsole() << "ERROR: unknown file type " << file_type << std::endl;
            return;
        }
    }
    
    
    // If there was an error then exit
    if (err) {
        sim->errConsole() << "ERROR: could not read file " << file_name << std::endl;
        return;
    }
    
    // Schultz: Currently we only support a single stress state. We may want to keep writing stress
    // states every N events, then just load the last event saved in the stress state file.
    if (stress_filename != "" && stress_file_type != "" && stress_index_filename != "" && !err) {
        // Grab the stress values of the first (only) stress state
        stress = stress_set[0].stresses();
        // Also set the sim year to the year the stresses were saved
        sim->setYear(stress_set[0].getYear());
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

        // We also need distance along strike data for computing section length for computing stress drops
        for (unsigned int i=0; i<3; ++i) new_block.set_das(i, new_element.das(i));

        new_block.set_is_quad(new_element.is_quad());
        new_block.set_rake(new_element.rake());
        new_block.set_slip_rate(new_element.slip_rate());
        new_block.set_aseismic(new_element.aseismic());
        new_block.set_lame_mu(new_element.lame_mu());
        new_block.set_lame_lambda(new_element.lame_lambda());
        new_block.set_max_slip(new_element.max_slip());
        new_block.set_stress_drop(new_element.stress_drop());

        // Set VQ specific values
        // Uncommenting the line setting fault_id (Kasey)

        new_block.setFaultID(world.section(eit->section_id()).fault_id());
        new_block.setSectionID(eit->section_id());    // TODO: add sections?

        BlockID bid = sim->addBlock(new_block);

        // If given an initial stress state, set those stresses
        if (stress_filename != "" && stress_file_type != "" && stress_index_filename != "") {
            assert(stress[bid]._element_id == bid);
            sim->console() << bid << "  "  << stress[bid]._shear_stress << "  " << stress[bid]._normal_stress << "  " << std::endl;
            sim->setInitShearNormalStress(bid, stress[bid]._shear_stress, stress[bid]._normal_stress);
        } else {
            sim->setInitShearNormalStress(bid, 0, 0);
        }
    }
}
