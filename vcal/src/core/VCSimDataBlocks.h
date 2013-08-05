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

#ifndef _VCSIM_DATA_BLOCKS_H_
#define _VCSIM_DATA_BLOCKS_H_

class VCSimDataBlocks {
private:
	//! Set of simulation model blocks
	BlockList					blocks;

	//! Temporarily store recurrence intervals during Greens function calculation
	std::map<BlockID, double>	recurrences;
	
public:
	BlockList::iterator begin(void) { return blocks.begin(); };
	BlockList::iterator end(void) { return blocks.end(); };
	BlockList::const_iterator begin(void) const { return blocks.begin(); };
	BlockList::const_iterator end(void) const { return blocks.end(); };
	
	BlockID addBlock(const Block &new_block);
	Block &getBlock(const BlockID &block_num) { return blocks[block_num]; };
	const Block &getBlock(const BlockID &block_num) const { return blocks[block_num]; };
	BlockInfo getBlockInfo(const Block &block);
	
	void setRecurrence(const BlockID &bid, const double &rec) { recurrences[bid] = rec; };
	double getRecurrence(const BlockID &bid) { return recurrences[bid]; };
	
	unsigned int numGlobalBlocks(void) const { return blocks.size(); };
};

#endif
