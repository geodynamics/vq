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

#include <string>

#include "config.h"

#ifdef MPI_C_FOUND
#include <mpi.h>
#endif

#ifndef _MPI_DEBUG_OUTPUT_STREAM_H_
#define _MPI_DEBUG_OUTPUT_STREAM_H_

/*
 MPI Debug Output Class
 Problem: Using normal cerr or cout on MPI causes output to mix, you don't know which machine wrote the output, it's hard to merge output from different systems
 Solution: Provide a unified interface to output on an MPI system
 What we need:
 Should act like a standard stream, during initialization user specifies a standard stream it will output to
 Should provide options to clarify output (node rank, node name?, time stamps, etc)
 Should provide debug levels to ignore output for non-debug runs
 Should provide gathering function to collect values over multiple nodes and output in a specified order
*/

class MPIDebugOutputStream {
private:
#ifdef MPI_C_FOUND
	MPI_Comm		comm_world;
#endif
	bool			print_node_rank, print_timestamp;
	std::string		output_buf;
	int				output_node_rank, my_node_rank, world_size;
	
public:
	MPIDebugOutputStream(const int &output_node);
	void setOutputFormat(const bool &print_rank, const bool &print_time);
	void write(const std::string &out_str);
	void flush(void);
};

#endif
