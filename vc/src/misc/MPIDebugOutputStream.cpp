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

#include "MPIDebugOutputStream.h"
#include <sstream>

MPIDebugOutputStream::MPIDebugOutputStream(const int &output_node) : print_node_rank(true), print_timestamp(false), output_node_rank(output_node) {
#ifdef MPI_C_FOUND
	comm_world = MPI_COMM_WORLD;
	MPI_Comm_rank(comm_world, &my_node_rank);
	MPI_Comm_size(comm_world, &world_size);
#endif
}

void MPIDebugOutputStream::setOutputFormat(const bool &print_rank, const bool &print_time) {
	print_node_rank = print_rank;
	print_timestamp = print_time;
}

void MPIDebugOutputStream::write(const std::string &out_str) {
	std::stringstream		ss;
	
	if (print_node_rank) ss << "[N" << my_node_rank << "] ";
	if (print_timestamp) ss << "(" << "time" << ") ";
	ss << out_str << "\n";
	output_buf += ss.str();
}

// "Flush" the output stream by collecting all messages on the output node
// and writing them to the specified output stream
void MPIDebugOutputStream::flush(void) {
#ifdef MPI_C_FOUND
	int		*msg_lens, *msg_disps, my_msg_len, i, total_len;
	char	*my_msg_storage, *all_msg_storage;
	
	my_msg_len = output_buf.size();
	if (my_node_rank == output_node_rank) {
		msg_lens = new int[world_size];
		msg_disps = new int[world_size];
	}
	
	// Get the message sizes on the output node
	MPI_Gather(&my_msg_len, 1, MPI_INT,
			   msg_lens, 1, MPI_INT,
			   output_node_rank, comm_world);
	
	if (my_node_rank == output_node_rank) {
		for (i=0;i<world_size;++i) {
			msg_disps[i] = total_len;
			total_len += msg_lens[i];
		}
		all_msg_storage = new char[total_len];
	}
	
	// Gather the actual messages on the output node
	MPI_Gatherv(my_msg_storage, my_msg_len, MPI_CHAR,
				all_msg_storage, msg_lens, msg_disps, MPI_CHAR,
				output_node_rank, comm_world);
	
	// Write the messages to the output stream
	for (i=0;i<world_size;++i) {
		std::cerr << std::string(&all_msg_storage[msg_disps[i]], msg_lens[i]);
	}
	// Flush everything before continuing
	std::cerr << std::flush;
	
	// Cleanup the arrays
	delete my_msg_storage;
	if (my_node_rank == output_node_rank) {
		delete msg_lens;
		delete msg_disps;
		delete all_msg_storage;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}
