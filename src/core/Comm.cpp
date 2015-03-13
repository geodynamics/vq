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

#include "Comm.h"
#include "SimFramework.h"

#ifdef MPI_C_FOUND

// Note: We assume these operations are only called for the BlockVal MPI datatype
/*!
 Calculates the minimum value and its associated block among the specified input array.
 If two blocks have the same value, order by block ID.
 */
void BlockValMinimum(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
    BlockVal        *in, *out;
    int             i;

    in = (BlockVal *)invec;
    out = (BlockVal *)inoutvec;

    for (i=0; i<*len; ++i) {
        if (in[i].val < out[i].val) {
            out[i].val = in[i].val;
            out[i].block_id = in[i].block_id;
        } else if (in[i].val == out[i].val) {
            out[i].block_id = (in[i].block_id < out[i].block_id ? in[i].block_id : out[i].block_id);
        }

        // else leave the out array as it is
    }
}

/*!
 Calculates the maximum value and its associated block among the specified input array.
 If two blocks have the same value, order by block ID.
 */
void BlockValMaximum(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
    BlockVal        *in, *out;
    int             i;

    in = (BlockVal *)invec;
    out = (BlockVal *)inoutvec;

    for (i=0; i<*len; ++i) {
        if (in[i].val > out[i].val) {
            out[i].val = in[i].val;
            out[i].block_id = in[i].block_id;
        } else if (in[i].val == out[i].val) {
            out[i].block_id = (in[i].block_id < out[i].block_id ? in[i].block_id : out[i].block_id);
        }

        // else leave the out array as it is
    }
}

/*!
 Calculates the sum of all block values and the number of blocks among the specified input array.
 */
void BlockValSum(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
    BlockVal        *in, *out;
    int             i;

    in = (BlockVal *)invec;
    out = (BlockVal *)inoutvec;

    for (i=0; i<*len; ++i) {
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
#ifdef MPI_C_FOUND
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
    int     global_fail, my_fail;

    my_fail = local_fail;
#ifdef MPI_C_FOUND
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
#ifdef MPI_C_FOUND
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
#ifdef MPI_C_FOUND
    int             block_lengths[3];
    MPI_Aint        displacements[3];
    MPI_Datatype    datatypes[3];

    updateFieldCounts = updateFieldDisps = NULL;
    updateFieldSendBuf = updateFieldRecvBuf = NULL;
    updateFieldSendIDs = updateFieldRecvIDs = NULL;
    failBlockSendBuf = failBlockRecvBuf = NULL;

    // Register BlockVal datatype
    block_lengths[0] = block_lengths[1] = 1;    // 1 member for each block
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
    block_lengths[2] = 1;
    displacements[0] = 0;
    displacements[1] = 5*sizeof(double);
    displacements[2] = 5*sizeof(double)+sizeof(unsigned int);
    datatypes[0] = MPI_DOUBLE;
    datatypes[1] = MPI_UNSIGNED;
    datatypes[2] = MPI_INT;

    MPI_Type_struct(3, block_lengths, displacements, datatypes, &element_sweep_type);
    MPI_Type_commit(&element_sweep_type);
#endif
}

VCComm::~VCComm(void) {
#ifdef MPI_C_FOUND

    // these are declared with "new type[]", so use delete []:
    //
    if (updateFieldCounts) delete [] updateFieldCounts;

    if (updateFieldDisps) delete [] updateFieldDisps;

    if (updateFieldSendBuf) delete [] updateFieldSendBuf;

    if (updateFieldRecvBuf) delete [] updateFieldRecvBuf;

    if (updateFieldSendIDs) delete [] updateFieldSendIDs;

    if (updateFieldRecvIDs) delete [] updateFieldRecvIDs;

    if (failBlockSendBuf) delete [] failBlockSendBuf;

    if (failBlockRecvBuf) delete [] failBlockRecvBuf;

    if (failBlockCounts) delete [] failBlockCounts;

    if (failBlockDisps) delete [] failBlockDisps;

    MPI_Type_free(&block_val_type);
    MPI_Op_free(&bv_min_op);
    MPI_Op_free(&bv_max_op);
    MPI_Op_free(&bv_sum_op);
    MPI_Type_free(&element_sweep_type);
#endif
}
