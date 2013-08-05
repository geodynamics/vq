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

#include "VCComm.h"
#include "SimFramework.h"

#ifdef HAVE_MPI

// Note: We assume these operations are only called for the BlockVal MPI datatype
/*!
 Calculates the minimum value and its associated block among the specified input array.
 */
void BlockValMinimum(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
	BlockVal		*in, *out;
	int				i;
	
	in = (BlockVal*)invec;
	out = (BlockVal*)inoutvec;
	
	for (i=0;i<*len;++i) {
		if (in[i].val < out[i].val) {
			out[i].val = in[i].val;
			out[i].block_id = in[i].block_id;
		}
		// else leave the out array as it is
	}
}

/*!
 Calculates the maximum value and its associated block among the specified input array.
 */
void BlockValMaximum(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
	BlockVal		*in, *out;
	int				i;
	
	in = (BlockVal*)invec;
	out = (BlockVal*)inoutvec;
	
	for (i=0;i<*len;++i) {
		if (in[i].val > out[i].val) {
			out[i].val = in[i].val;
			out[i].block_id = in[i].block_id;
		}
		// else leave the out array as it is
	}
}

/*!
 Calculates the sum of all block values and the number of blocks among the specified input array.
 */
void BlockValSum(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
	BlockVal		*in, *out;
	int				i;
	
	in = (BlockVal*)invec;
	out = (BlockVal*)inoutvec;
	
	for (i=0;i<*len;++i) {
		out[i].val += in[i].val;
		out[i].block_id = *len;
	}
}
#endif

/*!
 Performs an all-reduce of the specified operation on the given input.
 This can be a minimum, maximum or sum over all BlockVals.
 */
void VCComm::allReduceBlockVal(BlockVal &in_val, BlockVal &out_val, const BlockValOp &op) {
#ifdef HAVE_MPI
#ifdef DEBUG
	startTimer(reduce_comm_timer);
#endif
	switch (op) {
		case BLOCK_VAL_MIN:
			MPI_Allreduce(&in_val, &out_val, 1, block_val_type, bv_min_op, MPI_COMM_WORLD);
			break;
		case BLOCK_VAL_MAX:
			MPI_Allreduce(&in_val, &out_val, 1, block_val_type, bv_max_op, MPI_COMM_WORLD);
			break;
		case BLOCK_VAL_SUM:
			MPI_Allreduce(&in_val, &out_val, 1, block_val_type, bv_sum_op, MPI_COMM_WORLD);
			break;
		default:
			std::cerr << "Unknown reduce operation. Quitting." << std::endl;
			break;
	}
#ifdef DEBUG
	stopTimer(reduce_comm_timer);
#endif
#else
	out_val.block_id = in_val.block_id;
	out_val.val = in_val.val;
#endif
}

/*!
 Performs a global synchronization check to see if any nodes have more blocks to fail.
 Returns the number of nodes with blocks to fail.
 */
int VCComm::blocksToFail(const bool &local_fail) {
	int		global_fail, my_fail;
	
	my_fail = local_fail;
#ifdef HAVE_MPI
#ifdef DEBUG
	startTimer(fail_comm_timer);
#endif
	MPI_Allreduce(&my_fail, &global_fail, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#ifdef DEBUG
	stopTimer(fail_comm_timer);
#endif
#else
	global_fail = my_fail;
#endif
	
	return global_fail;
}

/*!
 Broadcasts the specified value from the root node to all other nodes.
 */
int VCComm::broadcastValue(const int &bval) {
#ifdef HAVE_MPI
	int bcast_val = bval;
	MPI_Bcast(&bcast_val, 1, MPI_INT, ROOT_NODE_RANK, MPI_COMM_WORLD);
	return bcast_val;
#else
	return bval;
#endif
}

/*!
 Register the block ID/value MPI datatype and block sweep datatype.
 This must exactly match the contents of BlockVal and BlockSweepVals.
 */
VCComm::VCComm(void) {
#ifdef HAVE_MPI
	int				block_lengths[2];
	MPI_Aint		displacements[2];
	MPI_Datatype	datatypes[2];
	
	updateFieldCounts = updateFieldDisps = NULL;
	updateFieldSendBuf = updateFieldRecvBuf = NULL;
	updateFieldSendIDs = updateFieldRecvIDs = NULL;
	failBlockSendBuf = failBlockRecvBuf = NULL;
	
	// Register BlockVal datatype
	block_lengths[0] = block_lengths[1] = 1;	// 1 member for each block
	displacements[0] = 0;
	displacements[1] = sizeof(double);
	datatypes[0] = MPI_DOUBLE;
	datatypes[1] = MPI_INT;
	
	MPI_Type_struct(2, block_lengths, displacements, datatypes, &block_val_type);
	MPI_Type_commit(&block_val_type);
	
	// Register BlockVal related operations 
	MPI_Op_create(BlockValMinimum, true, &bv_min_op);
	MPI_Op_create(BlockValMaximum, true, &bv_max_op);
	MPI_Op_create(BlockValSum, true, &bv_sum_op);
	
	// Register BlockSweepVals datatype
	block_lengths[0] = 5;
	block_lengths[1] = 1;
	displacements[0] = 0;
	displacements[1] = 5*sizeof(double);
	datatypes[0] = MPI_DOUBLE;
	datatypes[1] = MPI_INT;
	
	MPI_Type_struct(2, block_lengths, displacements, datatypes, &block_sweep_type);
	MPI_Type_commit(&block_sweep_type);
#endif
}

VCComm::~VCComm(void) {
#ifdef HAVE_MPI
	if (updateFieldCounts) delete updateFieldCounts;
	if (updateFieldDisps) delete updateFieldDisps;
	if (updateFieldSendBuf) delete updateFieldSendBuf;
	if (updateFieldRecvBuf) delete updateFieldRecvBuf;
	if (updateFieldSendIDs) delete updateFieldSendIDs;
	if (updateFieldRecvIDs) delete updateFieldRecvIDs;
	if (failBlockSendBuf) delete failBlockSendBuf;
	if (failBlockRecvBuf) delete failBlockRecvBuf;
	if (failBlockCounts) delete failBlockCounts;
	if (failBlockDisps) delete failBlockDisps;
	
	MPI_Type_free(&block_val_type);
	MPI_Op_free(&bv_min_op);
	MPI_Op_free(&bv_max_op);
	MPI_Op_free(&bv_sum_op);
	MPI_Type_free(&block_sweep_type);
#endif
}
