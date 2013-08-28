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
