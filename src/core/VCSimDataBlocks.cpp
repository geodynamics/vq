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

#include "VCSimDataBlocks.h"
#include "SimError.h"

/*!
 Add a block to the block list and perform some simple correctness checks on it.
 Returns the simulation assigned ID of the added block.
 */
BlockID VCSimDataBlocks::addBlock(const Block &new_block) {
    BlockID     new_block_id;

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
