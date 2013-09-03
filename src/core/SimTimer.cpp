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

#include "SimTimer.h"
#include "SimError.h"
#include "QuakeLibUtil.h"

#include <iomanip>
#include <math.h>
#include <stdexcept>
#include <valarray>

#ifdef VC_HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

/*!
 Initializes a timer within this simulation framework. The timer consists of a name,
 whether it should have a global barrier before checking time, and whether
 it is preemptable by other times (e.g. communication timers).
 */
int SimTimer::initTimer(const std::string &timer_name, const bool &barrier, const bool &preemptable) {
	int		next_num;
	
	next_num = accumulated_times.size();
	timer_names.push_back(timer_name);
	accumulated_times.push_back(0.0);
	start_times.push_back(0.0);
	use_barrier.push_back(barrier);
	can_preempt.push_back(preemptable);
	num_timings.push_back(0);
	
	return next_num;
}

/*!
 Start a timer given the timer ID.
 */
void SimTimer::startTimer(const int &timer_id) {
	unsigned int		i;
	
	// Make sure the timer isn't already running
	assertThrow(!start_times[timer_id], "Timer already started.");
	
	// If any preemptable timers are running, suspend them
	for (i=0;i<accumulated_times.size();++i) {
		if (can_preempt[i] && start_times[i]) {
			paused_timers.push(i);
		}
	}
	// Record the start time of this timer
	start_times[timer_id] = curTime(use_barrier[timer_id]);
}

/*!
 Stop a timer given the timer ID.
 */
void SimTimer::stopTimer(const int &timer_id) {
	double accum_time;
	int paused_id;
	
	// Make sure the timer is running before we stop it
	assertThrow(start_times[timer_id], "Timer stopped without first being started.");
	
	accum_time = curTime(use_barrier[timer_id]) - start_times[timer_id];
	// For any paused timers, subtract this time interval from their running total
	while (!paused_timers.empty()) {
		paused_id = paused_timers.top();
		paused_timers.pop();
		start_times[paused_id] += accum_time;
	}
	// Add this time interval to the accumulated time
	accumulated_times[timer_id] += accum_time;
	start_times[timer_id] = 0;
	num_timings[timer_id]++;
}

/*!
 Check the current time using whatever function is
 appropriate to the environment.
 */
double SimTimer::curTime(const bool &use_barrier) {
#ifdef MPI_C_FOUND
	if (use_barrier) MPI_Barrier(MPI_COMM_WORLD);
	return MPI_Wtime();
#else
    struct timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_sec + (tv.tv_usec/1.e6);
#endif
}

/*!
 Gathers global timing statistics and prints the formatted values.
 */
void SimTimer::printAllTimers(std::ostream &out_stream, int world_size, int node_rank, int root_node_rank) {
	int		num_timers, i, max_name_len, ind_col_width, num_col_width, num_prec, num_times, *count_recv_buf, avg_counts;
	double	my_time, *time_recv_buf, max_time, min_time, sum_time, mean_time, stddev_time;
	
	num_timers = accumulated_times.size();
	
	// Calculate column widths for formatting
	for (max_name_len=0,i=0;i<num_timers;i++)
		max_name_len = (int)fmax((long unsigned int)max_name_len, (long unsigned int)timer_names[i].size()+5);
	ind_col_width = (int)(log(num_timers) + 1);
	num_prec = 3;
	num_col_width = (int)fmax(num_prec+2, 11);
	
	// Print the header line
	if (node_rank == root_node_rank) {
		out_stream.setf(std::ios::right);
		out_stream << std::setw(ind_col_width) << "#";
		out_stream.setf(std::ios::right);
		out_stream << std::setw(max_name_len+1) << "Timer Name ";
		out_stream.setf(std::ios::right);
		out_stream << std::setw(num_col_width) << "AvgCounts ";
		out_stream.setf(std::ios::right);
		out_stream << std::setw(num_col_width) << "Min ";
		out_stream.setf(std::ios::right);
		out_stream << std::setw(num_col_width) << "Max ";
		out_stream.setf(std::ios::right);
		out_stream << std::setw(num_col_width) << "Mean ";
		out_stream.setf(std::ios::right);
		out_stream << std::setw(num_col_width) << "StdDev ";
		out_stream << std::endl;
		time_recv_buf = new double[world_size];
		count_recv_buf = new int[world_size];
	} else {
		time_recv_buf = NULL;
		count_recv_buf = NULL;
	}
	
	// For each timer, get statistics over all nodes in comm and print them to out_stream on node proc_num
	for (i=0;i<num_timers;++i) {
		my_time = accumulated_times[i];
		num_times = num_timings[i];
#ifdef MPI_C_FOUND
		MPI_Gather(&my_time, 1, MPI_DOUBLE, time_recv_buf, 1, MPI_DOUBLE, root_node_rank, MPI_COMM_WORLD);
		MPI_Gather(&num_times, 1, MPI_INT, count_recv_buf, 1, MPI_INT, root_node_rank, MPI_COMM_WORLD);
#else
		time_recv_buf[0] = my_time;
		count_recv_buf[0] = num_times;
#endif
		if (node_rank == root_node_rank) {
			std::valarray<double> all_times(time_recv_buf, world_size);
			std::valarray<int> all_counts(count_recv_buf, world_size);
			max_time = all_times.max();
			min_time = all_times.min();
			sum_time = all_times.sum();
			avg_counts = all_counts.sum()/world_size;
			mean_time = sum_time/world_size;
			all_times -= mean_time;
			all_times *= all_times;
			stddev_time = sqrt(all_times.sum()/world_size);
			
			out_stream.setf(std::ios::right);
			out_stream << std::setw(ind_col_width) << i;
			out_stream.setf(std::ios::right);
			out_stream << std::setw(max_name_len) << timer_names[i];
			out_stream.setf(std::ios::fixed);
			out_stream << std::setprecision(num_prec);
			out_stream << std::setw(num_col_width) << avg_counts;
			out_stream << std::setw(num_col_width) << min_time;
			out_stream << std::setw(num_col_width) << max_time;
			out_stream << std::setw(num_col_width) << mean_time;
			out_stream << std::setw(num_col_width) << stddev_time;
			out_stream << std::endl;
			out_stream.unsetf(std::ios::right);
			out_stream.unsetf(std::ios::fixed);
		}
	}
	if (time_recv_buf) delete time_recv_buf;
	if (count_recv_buf) delete count_recv_buf;
}
