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

#include <string>

//#ifdef HAVE_CONFIG_H
#include "config.h"
//#endif

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
