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

#include <vector>
#include <string>
#include <stack>
#include <iostream>

#include "config.h"

#ifdef VQ_HAVE_MATH_H
#include <math.h>
#endif

#ifdef MPI_C_FOUND
#include <mpi.h>
#endif

#ifndef _SIM_TIMER_H_
#define _SIM_TIMER_H_

class SimTimer {
    private:
        std::vector<std::string>    timer_names;
        std::vector<double>         accumulated_times;
        std::vector<double>         start_times;
        std::vector<bool>           use_barrier;
        std::vector<bool>           can_preempt;
        std::stack<int>             paused_timers;
        std::vector<int>            num_timings;

    public:
        double curTime(const bool &use_barrier=false);
        int initTimer(const std::string &timer_name, const bool &barrier, const bool &preemptable);
        void startTimer(const int &timer_id);
        void stopTimer(const int &timer_id);
        void printAllTimers(std::ostream &out_stream, int world_size, int node_rank, int root_node_rank);
};

#endif
