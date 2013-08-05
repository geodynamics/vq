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

#include "EqSimFileParse.h"

#include <iostream>
#include <iomanip>

void EqSimFileParse::initDesc(const SimFramework *_sim) const {
	const VCSimulation			*sim = static_cast<const VCSimulation*>(_sim);
	
	sim->console() << "# EqSim condition file: " << (sim->getEqSimConditionFile().empty() ? "none" : sim->getEqSimConditionFile()) << std::endl;
	sim->console() << "# EqSim friction file: " << sim->getEqSimFrictionFile() << std::endl;
	sim->console() << "# EqSim geometry file: " << sim->getEqSimGeometryFile() << std::endl;
	sim->console() << "# Section params file: " << (sim->getSectionParamsFile().empty() ? "none" : sim->getSectionParamsFile()) << std::endl;
}

/*!
 For dry runs the EqSim files are parsed as usual.
 */
void EqSimFileParse::dryRun(SimFramework *_sim) {
	init(_sim);
}

/*!
 Parse the geometry file and optionally the condition and friction files.
 */
void EqSimFileParse::init(SimFramework *_sim) {
	VCSimulation								*sim = static_cast<VCSimulation*>(_sim);
	std::string									condition_file_name = sim->getEqSimConditionFile();
	std::string									friction_file_name = sim->getEqSimFrictionFile();
	std::string									geometry_file_name = sim->getEqSimGeometryFile();
	std::string									section_params_file_name = sim->getSectionParamsFile();
	quakelib::EQSimGeomSectionMap::const_iterator	sit;
	quakelib::EQSimGeomRectangleMap::const_iterator	it;
	quakelib::Conversion						convert(sim->getBaseLatLon());
	bool										read_cond_file, read_section_params_file;
	std::map<FaultID, std::map<double, int> >	layer_map;
    std::map<FaultID, int>						fault_size_map;
	quakelib::EQSimErrors						error_report;
    
	read_cond_file = false;
	// Parse and verify the condition file
	read_cond_file = condition_data.parse_file(condition_file_name);
	if (read_cond_file) {
		condition_data.validate(error_report);
		error_report.write(_sim->console());
		error_report.clear();
	}
	
	// Parse and verify the friction file
	friction_data.parse_file(friction_file_name);
	friction_data.validate(error_report);
	error_report.write(_sim->console());
	error_report.clear();
	
	// Parse and verify the geometry file
	geometry_data.parse_file(geometry_file_name);
	geometry_data.validate(error_report);
	error_report.write(_sim->console());
	error_report.clear();
	
	read_section_params_file = false;
	// Parse the section params file
	// TODO: verify the section params file
	read_section_params_file = section_params_data.parse_file(section_params_file_name);
	/*
	if (read_cond_file) {
		condition_data.validate(error_report);
		error_report.write(_sim->console());
		error_report.clear();
	}
	*/
	
	// Triangle elements are currently not supported
	if (geometry_data.num_triangles() > 0)
		_sim->errConsole() << "ERROR: Currently cannot handle triangle elements. These elements will be ignored." << std::endl;
	
	for (sit=geometry_data.sections.begin();sit!=geometry_data.sections.end();++sit) {
		// Assuming aligned rectangular elements
		for (it=sit->second.rectangles.begin();it!=sit->second.rectangles.end();++it) {
			if (it->second.slip_rate() != 0) {
				Block					new_block;
				quakelib::EQSimGeometryVertex	v[4];
				quakelib::Vec<3>			v_xy;
				double					stress_drop;
				unsigned int			i;
				
				new_block.clear();
				
				new_block.setFaultID(sit->second.fid());
				new_block.setSectionID(sit->second.sid());
				v[0] = sit->second.vertices.find(it->second.vertex(0))->second;
				v[1] = sit->second.vertices.find(it->second.vertex(1))->second;
				v[2] = sit->second.vertices.find(it->second.vertex(2))->second;
				v[3] = sit->second.vertices.find(it->second.vertex(3))->second;
				
				// Take the 4 vertices that define a block
				for (i=0;i<4;++i) {
					// Convert input vertex lat/lon into km basis
					v_xy = convert.convert2xyz(quakelib::LatLonDepth(v[i].loc().lat(), v[i].loc().lon(), v[i].loc().altitude()));
					new_block.set_vert(i, v_xy);
					new_block.set_das(i, v[i].das());
					
					//quakelib::LatLonDepth convert_back = convert.convert2LatLon(v_xy);
					
					//_sim->console() << std::fixed
					//	<< " "  << std::setprecision (10)<< v[i].loc().lat() << " " << std::setprecision (10) << convert_back.lat()
					//	<< " " << std::setprecision (10) << v[i].loc().lon() << " " << std::setprecision (10)<< convert_back.lon() << std::endl;
					
					switch (v[i].trace_flag()) {
						case 0:
							new_block.set_trace_flag(i, quakelib::NOT_ON_TRACE);
							break;
						case 1:
							new_block.set_trace_flag(i, quakelib::MIDDLE_TRACE);
							break;
						case 2:
							new_block.set_trace_flag(i, quakelib::BEGINNING_TRACE);
							break;
						case 3:
							new_block.set_trace_flag(i, quakelib::END_TRACE);
							break;
						default:
							_sim->errConsole() << "WARNING: Unknown vertex trace value: " << v[i].trace_flag() << std::endl;
							break;
					}
				}
				
				new_block.setUnitSlip(1.0);
				
				/*!*/
				
				new_block.set_rake(convert.deg2rad(it->second.rake()));
				
				//new_block.setDip(convert.deg2rad(it->second.dip));		// TODO: Figure out why some blocks have negative dip
				//new_block.setDipDir(convert.deg2rad(90+it->second.strike));		// Convert strike angle to dip direction
				
				new_block.set_slip_rate(it->second.slip_rate());
				
				layer_map[sit->second.sid()][v[2].loc().altitude()] += 1;
				
				fault_size_map[sit->second.fid()] += 1;
				
				// this is a hack to make the recurrence times look good. use either this or the slip rate hack not both
				double stat_str, dyn_str;
				stat_str = friction_data.get_static_strength(it->first);
				dyn_str = friction_data.get_dynamic_strength(it->first);
				new_block.setLambda(friction_data.get_lame_lambda());
				new_block.setMu(friction_data.get_lame_mu());
				//	dyn_str = 0;
				//	stat_str = 2.2e5;
				
				stress_drop = dyn_str - stat_str;
				new_block.setStrengths(dyn_str, stat_str);
				//new_block.setStressDrop(2.09753 * stress_drop);
				new_block.setStressDrop(stress_drop);
				
				if (read_section_params_file) {
					new_block.setDynamicVal(section_params_data.get_dyn_val(sit->second.sid()));
					new_block.setSlipScalingThreshold(section_params_data.get_st_val(sit->second.sid()));
				} else {
					new_block.setDynamicVal(sim->getDynamic());
					new_block.setSlipScalingThreshold(sim->getSlipScalingThreshold());
				}
				
				new_block.setFaultName(sit->second.name());
				new_block.set_aseismic(it->second.aseismic());
				new_block.setValid(true);
				if (read_cond_file) {
					double	shear_stress, normal_stress;
					shear_stress = condition_data.get_shear_stress(it->first);
					normal_stress = condition_data.get_normal_stress(it->first);
					new_block.setRhogd(normal_stress);
					new_block.setInitShearStress(shear_stress);
					new_block.setInitNormalStress(normal_stress);
				} else {
					new_block.setRhogd(155700000);
					new_block.setInitShearStress(0);
					new_block.setInitNormalStress(0);
				}
				
				new_block.setFailed(false);
				new_block.setFailedThisSweep(false);
				
				sim->addBlock(new_block);
			}
		}
	}
    
    BlockList::iterator			bit;
	for (bit=sim->begin();bit!=sim->end();++bit) {
		bit->setSectionLayers(  layer_map[bit->getSectionID()].size() );
        //double new_dynamic = bit->getDynamicVal()/bit->getSectionLayers();
        //bit->setDynamicVal(new_dynamic);
        //sim->console() << fault_size_map[bit->getFaultID()] << std::endl;
	
        bit->setFaultSize(fault_size_map[bit->getFaultID()]);
	}
	
}
