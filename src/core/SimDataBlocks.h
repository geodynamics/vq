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

#include "Block.h"

#ifndef _SIM_DATA_BLOCKS_H_
#define _SIM_DATA_BLOCKS_H_

class VCSimDataBlocks {
    private:
        //! Set of simulation model blocks
        BlockList                   blocks;

    public:
        BlockList::iterator begin(void) {
            return blocks.begin();
        };
        BlockList::iterator end(void) {
            return blocks.end();
        };
        BlockList::const_iterator begin(void) const {
            return blocks.begin();
        };
        BlockList::const_iterator end(void) const {
            return blocks.end();
        };

        BlockID addBlock(const Block &new_block);
        Block &getBlock(const BlockID &block_num) {
            assertThrow(block_num<blocks.size(), std::domain_error("Invalid block number"));
            return blocks[block_num];
        };
        const Block &getBlock(const BlockID &block_num) const {
            return blocks[block_num];
        };

        unsigned int numGlobalBlocks(void) const {
            return blocks.size();
        };
};

#endif
