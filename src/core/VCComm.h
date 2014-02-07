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

#include "SimTimer.h"
#include "VCBlock.h"

//#ifdef HAVE_CONFIG_H
#include "config.h"
//#endif

#ifdef MPI_C_FOUND
#include "mpi.h"
#endif

#ifndef _VC_COMM_H_
#define _VC_COMM_H_

class VCComm : virtual public SimTimer {
protected:
	//! Debugging timers for communication
	int							reduce_comm_timer, fail_comm_timer, dist_comm_timer, sweep_comm_timer;
	
#ifdef MPI_C_FOUND
	//! Information regarding the update field used for parallel data synchronization
	int							*updateFieldCounts, *updateFieldDisps;
	
	//! Buffers to send/receive update field values
	double						*updateFieldSendBuf, *updateFieldRecvBuf;
	BlockID						*updateFieldSendIDs, *updateFieldRecvIDs;
	
	int							*failBlockSendBuf, *failBlockRecvBuf;
	int							*failBlockCounts, *failBlockDisps;
	
	//! Registered MPI datatype for the block-value combination structure
	MPI_Datatype				block_val_type;
	
	//! Registered MPI operations to get minimum, maximum and sum of block-value structures
	MPI_Op						bv_min_op, bv_max_op, bv_sum_op;
	
	//! Registered MPI datatype for the block-sweep data structure
	MPI_Datatype				block_sweep_type;
#endif
	
public:
	VCComm(void);
	~VCComm(void);
	
	void allReduceBlockVal(BlockVal &in_val, BlockVal &out_val, const BlockValOp &op);
	int blocksToFail(const bool &local_fail);
	int broadcastValue(const int &bval);
};

#endif
