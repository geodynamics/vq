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

#include "VCSimDataBlocks.h"
#include "SimError.h"

/*!
 Add a block to the block list and perform some simple correctness checks on it.
 Returns the simulation assigned ID of the added block.
 */
BlockID VCSimDataBlocks::addBlock(const Block &new_block) {
	BlockID		new_block_id;
	
	// Perform block sanity checks
	// Blocks may not have a negative slip rate (instead, use rake)
	//assertThrow(new_block.dip()>=0&&new_block.dip()<=M_PI/2.0, "Block dip must be between 0 and 90 degrees.");
	assertThrow(new_block.getFaultID()>=0, "Block fault ID must be non-negative.");
	assertThrow(new_block.min_depth()>new_block.max_depth(), "Block top must be higher than block bottom.");
	assertThrow(new_block.slip_rate()>=0, "Blocks may not have a negative slip rate (use rake to specify direction instead).");
	assertThrow(new_block.getRhogd()>0, "Blocks must have a positive rhogd.");
	
	new_block_id = blocks.size();
	blocks.push_back(new_block);
	blocks.back().setBlockID(new_block_id);
	
	return new_block_id;
}
