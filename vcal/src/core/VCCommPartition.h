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

#include "VCBlock.h"

#ifndef _VC_COMM_PARTITION_H_
#define _VC_COMM_PARTITION_H_

class VCCommPartition {
protected:
	//! Map of local array indices to global block IDs for blocks assigned to this node
	BlockIDList					local_block_ids;
	
	//! Map of global block IDs to array indices for blocks assigned to this node
	BlockIDMap					global_block_ids;
	
	//! Map of which blocks belong to which nodes
	std::map<BlockID, int>		block_node_map;
	
	//! Map of which nodes manage which blocks (reverse of block_node_map)
	std::multimap<int, BlockID>	node_block_map;
	
public:
	unsigned int numLocalBlocks(void) const { return local_block_ids.size(); };
	
	BlockID getGlobalBID(const int &local_id) const { return local_block_ids.at(local_id); };
	int getLocalInd(const BlockID &global_id) const { return global_block_ids.at(global_id); };
	int getBlockNode(const BlockID &global_id) const { return block_node_map.at(global_id); };
	bool isLocalToNode(const BlockID &global_id) const { return (global_block_ids.count(global_id) > 0); };
};

#endif
