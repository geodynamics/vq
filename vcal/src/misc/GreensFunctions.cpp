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

#include "GreensFunctions.h"
#include "HDF5Data.h"
#include "SimError.h"
#include <iomanip>
#include <set>
#include <cmath>

void GreensFuncCalc::progressBar(VCSimulation *sim, const int &thread_num, const int &num_done_blocks) {
	if (thread_num == 0 && sim->curTime() > last_update+1) {
		outcnt++;
		last_update = sim->curTime();
		if ((outcnt%5)==0) {
			sim->console() << int(100*float(num_done_blocks)/sim->numLocalBlocks()) << "%";
		} else {
			sim->console() << ".";
		}
		sim->console() << std::flush;
	}
}

// **************************************************************************
// *** 2011 Okada code
// **************************************************************************
void GreensFuncCalc2011::CalculateGreens(VCSimulation *sim) {
    std::vector<int>		row_sizes;
    int                     num_blocks, t_num, n;
    
    num_blocks = sim->numGlobalBlocks();
    
    row_sizes.clear();
	for (n=0;n<sim->numGlobalBlocks();++n) row_sizes.push_back(sim->isLocalBlockID(n)?num_blocks:sim->numLocalBlocks());
	GreensValsSparseMatrix ssh = GreensValsSparseMatrix(row_sizes);
	
	row_sizes.clear();
	for (n=0;n<sim->numGlobalBlocks();++n) row_sizes.push_back(sim->isLocalBlockID(n)?num_blocks:0);
	GreensValsSparseMatrix snorm = GreensValsSparseMatrix(row_sizes);
    
    // Get the current thread # for OpenMP to avoid printing multiple progress bars.
#ifdef _OPENMP
	t_num = omp_get_thread_num();
#else
	t_num = 0;
#endif
    
    // Use OpenMP to parallelize the loop, with each thread calculating
	// the Greens function for a block at a time
#pragma omp parallel for schedule(static,1)
	for (n=0;n<sim->numLocalBlocks();++n) {
		progressBar(sim, t_num, n);
        
        InnerCalc2011(sim, sim->getGlobalBID(n), ssh, snorm);
	}
    
	// Symmetrize the shear stress matrix
	symmetrizeMatrix(sim, ssh);
	
	// Set the simulation Greens values to the local ones calculated
	for (int r=0;r<sim->numLocalBlocks();++r) {
		for (int c=0;c<num_blocks;++c) {
			int global_r = sim->getGlobalBID(r);
			sim->setGreens(global_r, c, ssh[global_r][c], snorm[global_r][c]);
		}
	}
}

void GreensFuncCalc2011::InnerCalc2011(VCSimulation *sim,
									   const BlockID &bnum,
									   GreensValsSparseMatrix &ssh,
									   GreensValsSparseMatrix &snorm) {
	Block source_block = sim->getBlock(bnum);
	BlockIDList							target_blocks;
	BlockIDList::const_iterator			bit;
	BlockList::const_iterator			it;
	double								stress_values[2];
	
	// Get a list of block IDs we want to calculate the Greens function for
	for (it=sim->begin();it!=sim->end();++it) target_blocks.push_back(it->getBlockID());
	
	// Convert the strains to Greens values and record them
	for (bit=target_blocks.begin();bit!=target_blocks.end();++bit) {
        Block target_block = sim->getBlock(*bit);
        
        target_block.get_rake_and_normal_stress_due_to_block(stress_values, source_block);
        
        ssh[bnum][*bit] = stress_values[0];
		snorm[bnum][*bit] = stress_values[1];
	}
}

/*!
 Output the contents of the sparse matrix.
 */
std::ostream& operator<<(std::ostream& os, const GreensValsSparseMatrix& m) {
	/*GreensMatrix::const_iterator		it;
	os << "|";
	for (int i=0;i<m.num_rows;++i) os << "-";
	os << "|" << std::endl;
	for (it=m.matrix.begin();it!=m.matrix.end();++it) os << *(it->second);
	os << "|";
	for (int i=0;i<m.num_rows;++i) os << "-";
	os << "|" << std::endl;*/
	return os;
}

// Read the Greens function values in from a specified file
void GreensFuncFileParse::CalculateGreens(VCSimulation* sim) {
	HDF5GreensDataReader	*greens_file_reader;
	BlockID					gid;
	int						i, j, num_global_blocks;
	double					*in_shear_green, *in_normal_green;
	
	// Open the Greens data file and initialize arrays to read in Greens values
	if (sim->getGreensInputfile().empty()) {
		sim->errConsole() << "ERROR: Greens input file undefined. Quitting." << std::endl;
		exit(-1);
	}
	num_global_blocks = sim->numGlobalBlocks();
	greens_file_reader = new HDF5GreensDataReader(sim->getGreensInputfile());
	in_shear_green = new double[num_global_blocks];
	in_normal_green = new double[num_global_blocks];
	
	// Read the Greens function shear and normal values
	for (i=0;i<sim->numLocalBlocks();++i) {
		progressBar(sim, 0, i);
		
		gid = sim->getGlobalBID(i);
		greens_file_reader->getGreensVals(gid, in_shear_green, in_normal_green);
		for (j=0;j<num_global_blocks;++j) {
			sim->setGreens(gid, j, in_shear_green[j], in_normal_green[j]);
		}
	}
	
	delete greens_file_reader;
	delete in_shear_green;
	delete in_normal_green;
}

void GreensFuncCalc::symmetrizeMatrix(VCSimulation *sim, GreensValsSparseMatrix &ssh) {
    double		sxrl, sxru;
    int			ir, ic, n;
	
	// If we're using MPI, exchange Greens values between nodes
	// in order to symmetrize the shear stress matrix
#ifdef HAVE_MPI
	int						i, world_size, local_rank, root_node, num_local_blocks, num_global_blocks;
	int						*local_counts, *local_ids, *global_ids, *displs;
	GREEN_VAL				*send_buf, *recv_buf;
	MPI_Datatype			data_type;
	
	// MPI and decomposition values
	world_size = sim->getWorldSize();
	local_rank = sim->getNodeRank();
	
	// Get the size of the data we will transmit
	if (sizeof(GREEN_VAL)==4) data_type = MPI_FLOAT;
	else if (sizeof(GREEN_VAL)==8) data_type = MPI_DOUBLE;
	
	// Find out how many local objects each node has
	local_counts = new int[world_size];
	num_local_blocks = sim->numLocalBlocks();
	num_global_blocks = sim->numGlobalBlocks();
	MPI_Allgather(&num_local_blocks, 1, MPI_INT, local_counts, 1, MPI_INT, MPI_COMM_WORLD);
	
	// Setup array for displacements
	displs = new int[world_size];
	displs[0] = 0;
	for (n=1;n<world_size;++n) displs[n] = displs[n-1]+local_counts[n-1];
	
	// Find out the global block IDs for each local block on other nodes
	local_ids = new int[num_local_blocks];
	for (n=0;n<num_local_blocks;++n) local_ids[n] = sim->getGlobalBID(n);
	global_ids = new int[num_global_blocks];
	MPI_Allgatherv(local_ids, num_local_blocks, MPI_INT, global_ids, local_counts, displs, MPI_INT, MPI_COMM_WORLD);
	
	// Setup send and receive buffers
	send_buf = new GREEN_VAL[num_global_blocks];
	recv_buf = new GREEN_VAL[num_local_blocks];
	
	// Go through each of the Greens function rows
	// If we computed the row, send pieces of it to other processors
	// If we need the row for symmetrization, receives pieces of it from another processor
	for (i=0;i<num_global_blocks;++i) {
		root_node = sim->getBlockNode(i);
		
		// Fill the send buffer with the appropriate values
		if (root_node == local_rank) {
			for (n=0;n<num_global_blocks;++n) {
				send_buf[n] = ssh[i][global_ids[n]];
			}
		}
		
		// Scatter each of the required matrix sections
		// There is no need to send snorm since it isn't used outside the local node
		MPI_Scatterv(send_buf, local_counts, displs, data_type,
					 recv_buf, num_local_blocks, data_type,
					 root_node, MPI_COMM_WORLD);
		
		// Copy the receive buffer to the matrix for non-local transfers
		for (n=0;n<num_local_blocks;++n) {
			if (root_node != local_rank) {
				ssh[i][n] = recv_buf[n];
			}
		}
	}
	delete displs;
	delete local_counts;
	delete send_buf;
	delete recv_buf;
	delete global_ids;
	delete local_ids;
#endif
	
	// Symmetrize the shear stress matrices
	for (ir=0;ir<sim->numGlobalBlocks();++ir) {
		bool full_row = sim->isLocalBlockID(ir);
		int num_elems = (full_row ? sim->numGlobalBlocks() : sim->numLocalBlocks());
		for (n=0;n<num_elems;++n) {
			ic = (full_row ? n : sim->getGlobalBID(n));
			double ir_area, ic_area;
			ir_area = sim->getBlock(ir).get_area();
			ic_area = sim->getBlock(ic).get_area();
			assertThrow(ir_area > 0 && ic_area > 0, "Blocks cannot have negative area.");
			
			int local_ir = (sim->isLocalBlockID(ic) ? ir : sim->getLocalInd(ir));
			sxru=ssh[ir][n]*ic_area;
			sxrl=ssh[ic][local_ir]*ir_area;
			
			ssh[ir][n]=0.5*(sxrl + sxru)/ic_area;
			ssh[ic][local_ir]=0.5*(sxrl + sxru)/ir_area;
		}
	}
}

void GreensFuncCalcBarnesHut::CalculateGreens(VCSimulation *sim) {
	BlockList::iterator		bit;
	quakelib::Octree<3>		*tree;
	quakelib::RectBound<3>	total_bound;
	unsigned int			n;
	int						t_num;
	
	// Get a list of 3D midpoints of all segments
	for (bit=sim->begin();bit!=sim->end();++bit) {
		total_bound.extend_bound(bit->center());
	}
	
	// Set up an octree of the model space
	tree = new quakelib::Octree<3>(total_bound);
	
	// Fill the octree with center points of the faults
	for (bit=sim->begin();bit!=sim->end();++bit) {
		tree->add_point(bit->center(), bit->getBlockID());
	}
	
	// Get the current thread # for OpenMP to avoid printing multiple progress bars.
#ifdef _OPENMP
	t_num = omp_get_thread_num();
#else
	t_num = 0;
#endif
	
	// TO FIX: Problem with shared global variables (run_bounds)
//#pragma omp parallel for shared(last_update) schedule(static,1)
	for (n=0;n<sim->numLocalBlocks();++n) {
		progressBar(sim, t_num, n);
		
		bhInnerCalc(sim, tree, sim->getGlobalBID(n));
	}
	
	// TODO: symmetrize Greens matrix
}

void GreensFuncCalcBarnesHut::bhInnerCalc(VCSimulation *sim, quakelib::Octree<3> *tree, const BlockID &bid) {
	Block								&source_block = sim->getBlock(bid);
	BlockIDList							target_blocks;
	BlockIDList::const_iterator			bit;
	BlockList::const_iterator			it;
	double								stress_values[2];
	quakelib::Vec<3>					target_point;
	quakelib::Octree<3>					*local_node;
	quakelib::RectBound<3>				local_bound;
	BlockID								start_id, end_id;
	std::set<std::pair<BlockID, BlockID> >::const_iterator	sbit;
	unsigned int						i;
	quakelib::Vec<3>					mean_center, mean_normal, rake_vec, up_vec, v0, v1, v2, v3;
	double								mean_dip, mean_rake, mean_area, side_length;
	double								mean_lambda, mean_mu, mean_unit_slip;
	Block								repr_block;
	
	// Allocate row for this block (if using compressible matrices)
	sim->allocateShearNormalRows(bid);
	
	target_point = sim->getBlock(bid).center();
	local_node = tree->get_leaf_containing_point(target_point);
	local_bound = local_node->bound();
	
	quakelib::BHAnalyzer<3>				bh_analysis(local_bound, sim->getBarnesHutTheta());
	tree->traverse(&bh_analysis);
	
	// Ensure that we got all the fault segments
	assertThrow(bh_analysis.num_checks()==sim->numGlobalBlocks(), "Didn't find all necessary points in the octree.");
	
	sbit = bh_analysis.run_bounds.begin();
	
	for (;sbit!=bh_analysis.run_bounds.end();++sbit) {
		start_id = sbit->first;
		end_id = sbit->second;
		
		// Make a list of target blocks
		target_blocks.clear();
		for (i=start_id;i!=end_id;++i) target_blocks.push_back(i);
		
		// Reset mean values/vectors
		mean_center = rake_vec = mean_normal = quakelib::Vec<3>();
		mean_dip = mean_rake = mean_area = mean_lambda = mean_mu = mean_unit_slip = 0;
		// Convert the strains to Greens values and record them
		for (bit=target_blocks.begin();bit!=target_blocks.end();++bit) {
			Block target_block = sim->getBlock(*bit);
			
			mean_center += target_block.center();
			rake_vec += target_block.rake_vector();
			mean_rake += target_block.rake();
			mean_area += target_block.get_area();
			mean_normal += target_block.normal();
			mean_lambda += target_block.getLambda();
			mean_mu += target_block.getMu();
			mean_unit_slip += target_block.getUnitSlip();
		}
		
		mean_center /= target_blocks.size();
		rake_vec /= target_blocks.size();
		mean_rake /= target_blocks.size();
		mean_area /= target_blocks.size();
		mean_lambda	/= target_blocks.size();
		mean_mu	/= target_blocks.size();
		mean_unit_slip	/= target_blocks.size();
		mean_normal /= target_blocks.size();
		
		side_length = sqrt(mean_area);
		up_vec = -mean_normal.cross(rake_vec);
		
		v0 = mean_center + rake_vec*(side_length/2) + up_vec*(side_length/2);
		v1 = mean_center + rake_vec*(side_length/2) - up_vec*(side_length/2);
		v2 = mean_center - rake_vec*(side_length/2) - up_vec*(side_length/2);
		v3 = mean_center - rake_vec*(side_length/2) + up_vec*(side_length/2);
		
		repr_block.set_rake(mean_rake);
		repr_block.set_vert(0, v0);
		repr_block.set_vert(1, v1);
		repr_block.set_vert(2, v2);
		repr_block.set_vert(3, v3);
		repr_block.setLambda(mean_lambda);
		repr_block.setMu(mean_mu);
		repr_block.setUnitSlip(mean_unit_slip);
		
		repr_block.get_rake_and_normal_stress_due_to_block(stress_values, source_block);
		
		// Set Green's values to averaged values
		for (bit=target_blocks.begin();bit!=target_blocks.end();++bit) {
			sim->setGreens(*bit, bid, stress_values[0], stress_values[1]);
		}
	}	
	
	// Compress the computed matrix row if savings are greater than 30%
	sim->compressShearRow(bid, 0.7);
	sim->compressNormalRow(bid, 0.7);
}
