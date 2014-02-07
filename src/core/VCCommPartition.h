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

#include "VCBlock.h"

#ifndef _VC_COMM_PARTITION_H_
#define _VC_COMM_PARTITION_H_

class VCCommPartition {
    protected:
        //! Map of local array indices to global block IDs for blocks assigned to this node
        BlockIDList                 local_block_ids;

        //! Map of global block IDs to array indices for blocks assigned to this node
        BlockIDMap                  global_block_ids;

        //! Map of which blocks belong to which nodes
        std::map<BlockID, int>      block_node_map;

        //! Map of which nodes manage which blocks (reverse of block_node_map)
        std::multimap<int, BlockID> node_block_map;

    public:
        unsigned int numLocalBlocks(void) const {
            return local_block_ids.size();
        };

        BlockID getGlobalBID(const int &local_id) const {
            return local_block_ids.at(local_id);
        };
        int getLocalInd(const BlockID &global_id) const {
            return global_block_ids.at(global_id);
        };
        int getBlockNode(const BlockID &global_id) const {
            return block_node_map.at(global_id);
        };
        bool isLocalToNode(const BlockID &global_id) const {
            return (global_block_ids.count(global_id) > 0);
        };
};

#endif
