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
    std::stringstream       ss;

    if (print_node_rank) ss << "[N" << my_node_rank << "] ";

    if (print_timestamp) ss << "(" << "time" << ") ";

    ss << out_str << "\n";
    output_buf += ss.str();
}

// "Flush" the output stream by collecting all messages on the output node
// and writing them to the specified output stream
void MPIDebugOutputStream::flush(void) {
#ifdef MPI_C_FOUND
    int     *msg_lens, *msg_disps, my_msg_len, i, total_len;
    char    *my_msg_storage, *all_msg_storage;

    total_len = 0;
    my_msg_storage = all_msg_storage = NULL;
    msg_lens = NULL;
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
        for (i=0; i<world_size; ++i) {
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
    for (i=0; i<world_size; ++i) {
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
