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

#include "SimFramework.h"
#include <iomanip>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

/*
 * Simple Parallel Debugging Tool
 * Debugging parallel programs is HARD!
 * Rather than forking over hundreds of dollars for TotalView, you can use this flag to debug with GDB
 * Uncomment the #defines below and run the program. You will get a message like:
 *      PID 55848 on Frinkiac-9.local ready for attach
 * Open GDB and "attach 55848"
 * The program is now paused, waiting for your signal
 * To let it know to continue, set the i variable to something non-zero and start the program,
 * like "up 3" then "set var i = 5" then "cont"
 * You can now debug this process like a normal serial program
 * To change the target MPI process rank, use the MPI_GDB_DEBUG_RANK #define
 */
//#define MPI_GDB_DEBUG
#define MPI_GDB_DEBUG_RANK              0

/*!
 Initialize SimFramework given the command line arguments.
 If any argument is "-R" or "--dry_run" then turn dry run mode on.
 */
SimFramework::SimFramework(int argc, char **argv) {
	int		i;
	
	total_timer = initTimer("Total Time", true, false);
	init_finish_timer = initTimer("Init and Cleanup", false, true);
	barrier_timer = initTimer("Comm Barrier", false, false);
	
#ifdef MPI_C_FOUND
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &node_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
#ifdef MPI_GDB_DEBUG
	if (getNodeRank() == MPI_GDB_DEBUG_RANK) {
		volatile int i = 0;
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		printf("PID %d on %s ready for attach\n", getpid(), hostname);
		fflush(stdout);
		while (0 == i)
			sleep(5);
	}
	barrier();
#endif
#else
	node_rank = ROOT_NODE_RANK;
	world_size = 1;
#endif
	
	// Initialize the PAPI library
#ifdef HAVE_PAPI_H
	int		retval;
	retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT && retval > 0) {
		std::cerr << "PAPI library version mismatch" << std::endl;
		exit(1);
	} else if (retval < 0) {
		std::cerr << "PAPI error: " << PAPI_strerror(retval) << std::endl;
		exit(1);
	}
#endif
	
	dry_run = false;
	for (i=0;i<argc;++i) {
		if (!strcmp(argv[i], "--dry_run") || !strcmp(argv[i], "-R")) {
			dry_run = true;
		}
	}
}

/*!
 Cleanup the SimFramework by releasing MPI.
 */
SimFramework::~SimFramework(void) {
	// Finish using the PAPI library and free associated resources
#ifdef HAVE_PAPI_H
	PAPI_shutdown();
#endif
	
#ifdef MPI_C_FOUND
	MPI_Finalize();
#endif
}

/*!
 Arranges the plugins in order based on their dependencies.
 This uses a simple topological sort algorithm to do so.
 If there are mutually exclusive dependencies it errors out
 if both plugins are enabled.  If there are required dependencies
 it errors out if they are not correctly satisfied.
 */
void SimFramework::orderPlugins(void) {
	std::map<PluginID, SimPlugin*>::iterator		mit;
	DependenceMap::iterator				it, kit, jit;
	std::map<PluginID, int>				num_deps;
	std::map<PluginID, int>::iterator	eit;
	int									i;
	PluginID							removed_id;
	bool								had_error;
	std::stringstream					ss;
	
	had_error = false;
	ordered_plugins.clear();
	
	for (mit=plugin_objs.begin();mit!=plugin_objs.end();++mit) {
		num_deps[mit->first] = 0;
	}
	
	// Create a list of how many unresolved dependencies each plugin has
	for (it=dep_types.begin();it!=dep_types.end();++it) {
		if (it->second != DEP_EXCLUSIVE) {
			num_deps[it->first.first]++;
		}
	}
	
	// For each of the plugins
	for (i=0;i<plugin_objs.size();++i) {
		// Find a plugin with no unresolved dependencies
		for (eit=num_deps.begin();eit!=num_deps.end();++eit) {
			if (eit->second == 0) {
				// removed_id will be the next plugin to be executed
				removed_id = eit->first;
				
				// Decrement the dependency count for plugins attached to the outgoing edges
				for (it=dep_types.begin();it!=dep_types.end();++it) {
					if (it->first.second == removed_id) {
						if (it->second != DEP_EXCLUSIVE) num_deps[it->first.first]--;
						// If this dependency is required, make sure the dependent plugin is also active
						if (it->second == DEP_REQUIRE) {
							// If the target is active, the source must be active as well
							if (plugin_active[it->first.first] && !plugin_active[it->first.second]) {
								ss << "Dependency error: \"" << plugin_objs[it->first.first]->name()
									<< "\" is active so \"" << plugin_objs[it->first.second]->name()
									<< "\" must also be active." << std::endl;
								throw std::logic_error(ss.str());
							}
						}
						// For exclusive dependencies if the target is active then
						// the source must not be active (and vice versa)
						if (it->second == DEP_EXCLUSIVE) {
							if (plugin_active[it->first.first]) {
								if (plugin_active[it->first.second]) {
									ss << "Dependency error: \"" << plugin_objs[it->first.first]->name()
										<< "\" is active so \"" << plugin_objs[it->first.second]->name()
										<< "\" must not be active." << std::endl;
									throw std::logic_error(ss.str());
								}
							} else if (plugin_active[it->first.second]) {
								if (plugin_active[it->first.first]) {
									ss << "Dependency error: \"" << plugin_objs[it->first.second]->name()
										<< "\" is active so \"" << plugin_objs[it->first.first]->name()
										<< "\" must not be active." << std::endl;
									throw std::logic_error(ss.str());
								}
							}
						}
					}
				}
				
				// Push this plugin on the list
				if (plugin_active[removed_id]) ordered_plugins.push_back(removed_id);
				break;
			}
		}
		// If we couldn't find any 0 dependency objects, we have a cycle in the graph
		if (eit == num_deps.end()) {
			ss << "Found a cycle in the plugin dependencies." << std::endl;
			throw std::logic_error(ss.str());
		}
		num_deps.erase(eit);
	}
	
	is_ordered = true;
}

/*!
 Broadcasts a value from the root node to all other nodes.
 Used for loop iteration timing.
 */
int SimFramework::internalBroadcast(const unsigned int &val) {
#ifdef MPI_C_FOUND
	unsigned int bcast_val = val;
	MPI_Bcast(&bcast_val, 1, MPI_UNSIGNED, ROOT_NODE_RANK, MPI_COMM_WORLD);
	return bcast_val;
#else
	return val;
#endif
}

/*!
 Registers a new plugin and whether it is active during this simulation.
 Returns a PluginID to refer to the plugin when assigning dependencies.
 */
PluginID SimFramework::registerPlugin(SimPlugin *new_plugin, const bool &is_active) {
	PluginID		next_id;
	
	next_id = plugin_active.size();
	plugin_active[next_id] = is_active;
	plugin_objs[next_id] = new_plugin;
	return next_id;
}

/*!
 Specifies dependencies between plugins. There are three dependency types - optional,
 required and exclusive. Because this new dependency may change the execution order,
 we must recalculate plugin ordering after adding new dependencies.
 */
void SimFramework::registerDependence(const PluginID &func_a, const PluginID &func_b, const DependenceType &dep_type) {
	PluginPair		dep_pair(func_a, func_b);
	
	is_ordered = false;
	dep_types.insert(std::make_pair(dep_pair, dep_type));
}

/*!
 Writes a DOT style graph file showing the plugin dependencies.
 This can be useful to understand execution flow.
 */
void SimFramework::writeDOT(std::ostream &os) const {
	DependenceMap::const_iterator		it;
	
	os << "digraph Dependencies {" << std::endl;
	for (it=dep_types.begin();it!=dep_types.end();++it) {
		os << "\"" << plugin_objs.at(it->first.second)->name() << "\" -> \"";
		os << plugin_objs.at(it->first.first)->name() << "\" [";
		switch (it->second) {
			case DEP_REQUIRE:
				os << "style=solid";
				break;
			case DEP_OPTIONAL:
				os << "style=dashed";
				break;
			case DEP_EXCLUSIVE:
				os << "dir=none style=bold";
				break;
		}
		os << "];" << std::endl;
	}
	os << "};" << std::endl;
}

/*!
 Initialize the simulation by setting up timers, changing settings for a
 parallel environment, and running the init() function for active plugins.
 */
void SimFramework::init(void) {
	std::vector<PluginID>::iterator	it;
	SimPlugin		*cur_plugin;
	int				width;
	
	// Start the total and initialization timers
	startTimer(total_timer);
	startTimer(init_finish_timer);
	
	next_bcast_iter = iters_per_sec = 10;
	cur_iter = last_bcast_iter = 0;
	
	// Only allow output at the root node
	if (!isRootNode()) disableConsole();
	
	// If the plugins have not been ordered, do so now
	if (!is_ordered) orderPlugins();
	
	// Output multiprocessor information
	width = 30;
#ifdef MPI_C_FOUND
	console() << std::setw(width) << std::left << "# *** MPI CPU count" << ": " << getWorldSize() << std::endl;
#endif
#ifdef _OPENMP
	console() << std::setw(width) << std::left << "# *** OpenMP Threads" << ": " << omp_get_max_threads() << std::endl;
#endif
	
	// Do the dry run or normal initialization
	for(it=ordered_plugins.begin();it!=ordered_plugins.end();++it) {
		cur_plugin = plugin_objs[*it];
		if (dry_run) {
			cur_plugin->dryRun(this);
		} else {
			cur_plugin->initDesc(this);
			cur_plugin->init(this);
			plugin_timers[*it] = (cur_plugin->needsTimer() ? initTimer(cur_plugin->name(), false, true) : -1);
		}
	}
	
	stopTimer(init_finish_timer);
}

/*!
 Run the simulation in its entirety. This involves calling the run() function
 of each plugin until either a plugin actively requests a stop (SIM_STOP_REQUIRED)
 or none of the plugins require the simulation to continue (SIM_STOP_OK).
 For dry runs this currently does nothing.
 */
void SimFramework::run(void) {
	std::vector<PluginID>::iterator	it;
	SimRequest						sim_status, end_request;
	double							cur_time, time_diff;
	
	if (dry_run) return;
	if (!is_ordered) orderPlugins();
	
	// Initialize loop timing mechanism. Start by timing the first 10 iterations
	if (isRootNode()) last_time = curTime();
	
	sim_status = SIM_CONTINUE;
	// While there are no active or passive stop requests
    while((sim_status & SIM_CONTINUE) && !(sim_status & SIM_STOP_REQUIRED) && !(sim_status & SIM_STOP_ERR)) {
		sim_status = SIM_STOP_OK;
		// Execute each plugin in order
		for(it=ordered_plugins.begin();it!=ordered_plugins.end();++it) {
			// If there is a timer associated with the plugin, run it during the plugin execution
			if (plugin_timers[*it] >= 0) startTimer(plugin_timers[*it]);
			try {
				end_request = plugin_objs[*it]->run(this);
			} catch (std::exception &e) {
				end_request = SIM_STOP_ERR;
				console() << "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
				console() << "# Plugin (" << plugin_objs[*it]->name() << ") had error: " << e.what() << std::endl;
				console() << "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			}
			if (plugin_timers[*it] >= 0) stopTimer(plugin_timers[*it]);
			
			// If a plugin actively requested a stop, notify the user which plugin it was
			if (end_request == SIM_STOP_REQUIRED) {
				console() << "# Plugin requested simulation end: " << plugin_objs[*it]->name() << std::endl;
			}
			sim_status = (SimRequest) (sim_status | end_request);
		}
		
		// Check if it's time to update the loop iteration timer
		cur_iter++;
		if (cur_iter == next_bcast_iter) {
			// Check how long the past X iterations took
			if (isRootNode()) {
				cur_time = curTime();
				time_diff = cur_time - last_time;
				last_time = cur_time;
				iters_per_sec = fmax(1, ((double)(cur_iter - last_bcast_iter))/time_diff);
			}
			// Share the iteration rate with other nodes
			iters_per_sec = internalBroadcast((unsigned int)iters_per_sec);
			last_bcast_iter = cur_iter;
			next_bcast_iter = cur_iter + iters_per_sec;
		}
	}
	
	// If the simulation ended due to a passive stop request, notify the user
	if (!(end_request | SIM_STOP_REQUIRED)) {
		console() << "# No plugins request continue. Simulation ending." << std::endl;
	}
}

/*!
 Finish the simulation by calling the finish() function for each
 plugin and stopping necessary timers.
 */
void SimFramework::finish(void) {
	std::vector<PluginID>::iterator	it;
	
	if (dry_run) return;
	if (!is_ordered) orderPlugins();
	
	// Start the init/finish timer and run finish() for all plugins()
	startTimer(init_finish_timer);
	for(it=ordered_plugins.begin();it!=ordered_plugins.end();++it) {
		plugin_objs[*it]->finish(this);
	}
	
	// Stop the timers
	stopTimer(init_finish_timer);
	stopTimer(total_timer);
}

/*!
 Returns a random double value in [0, 1).
 */
double SimFramework::randDouble(void) const {
	return ((double)random())*pow(2,-31);
}

/*!
 Returns a random float value in [0, 1).
 */
float SimFramework::randFloat(void) const {
	return ((float)random())*powf(2,-31);
}

/*!
 Returns a random integer in [0, max_int).
 */
int SimFramework::randInt(const int &max_int) const {
	return random()%max_int;
}
