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

#include <map>

//#ifdef HAVE_CONFIG_H
#include "config.h"
//#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_PAPI_H
#include "papi.h"
#include "papiStdEventDefs.h"
#endif

#include "SimTimer.h"

#ifndef _SIM_FRAMEWORK_H_
#define _SIM_FRAMEWORK_H_

class SimFramework;

#define ROOT_NODE_RANK		0

enum SimRequest {
	SIM_STOP_OK			= 0x00,		// The simulation may stop if no other plugins require a continue
	SIM_STOP_REQUIRED	= 0x01,		// The simulation must stop at this time step
	SIM_CONTINUE		= 0x02,		// The simulation should continue if no other plugins require a stop
	SIM_STOP_ERR		= 0x04		// The simulation had an error and must stop
};

/*!
 Represents an action that is taken during the simulation.
 Simulations can be written by subclassing this and rewriting key functions.
 Plugins are added and ordered in SimFramework so that their
 execution will fulfill required dependency ordering.
 By default, each function does nothing and requests a passive simulation stop.
 */
class SimPlugin {
public:
    virtual ~SimPlugin(void) {}
	
	//! A short name of the plugin
    virtual std::string name() const = 0;
	
	//! A longer description of the plugin shown at startup
	virtual void initDesc(const SimFramework *_sim) const = 0;
	
	//! Whether the plugin should be timed while it is running
	virtual bool needsTimer(void) const { return false; }
	
	//! Estimate how this plugin will behave in a real simulation
	//! (memory used, calculation speed, etc) given the parameters
	virtual void dryRun(SimFramework *_sim) {};
	
	//! Initialize the plugin given a simulation with full blocks
    virtual void init(SimFramework *_sim) {};
	
	//! Run the plugin during the course of the simulation
    virtual SimRequest run(SimFramework *_sim) { return SIM_STOP_OK; };
	
	//! Finish and cleanup the plugin and any related memory
    virtual void finish(SimFramework *_sim) {};
};

enum DependenceType {
	DEP_REQUIRE,		// For A -> B, if B is active then A must be active and A must be executed before B.
	DEP_OPTIONAL,		// For A -> B, if A and B are active then A must be executed before B.
	DEP_EXCLUSIVE		// For A -> B, if B is active then A must not be active
};

typedef int PluginID;
typedef std::pair<PluginID, PluginID> PluginPair;
typedef std::map<PluginPair, DependenceType> DependenceMap;

/*!
 A framework specifying which components exist in a simulation and how they interact.
 This is the main engine to drive simulations.  It provides functions to register
 plugins and their dependencies, simple functions for parallel computing environments,
 timing functions for plugins, random number generation and output functionality.
 */
class SimFramework : public virtual SimTimer {
private:
	//! A map of plugin-plugin dependencies and their types
	DependenceMap					dep_types;
	
	//! A map of which plugins are active
	std::map<PluginID, bool>		plugin_active;
	
	//! Map of plugin IDs to actual objects
	std::map<PluginID, SimPlugin*>	plugin_objs;
	
	//! Order of plugin execution determined by dependence and activation
	std::vector<PluginID>			ordered_plugins;
	
	//! Whether ordered_plugins is correctly sorted currently
	bool							is_ordered;
	
	//! Timer IDs for each plugin (if needed)
	std::map<PluginID, int>			plugin_timers;
	
	// Determines the correct order for plugins and writes to ordered_plugins
	void orderPlugins(void);
	
	//! Broadcast function used internally for framework
	int internalBroadcast(const unsigned int &val);
	
protected:
	//! Whether this simulation is a dry run or not
	bool							dry_run;
	
	int								node_rank, world_size;
	int								total_timer, init_finish_timer, barrier_timer;
	
	//! Recent average number of loop iterations per second and related values
	//! This is used for plugins to ensure periodic events are synchronized between
	//! parallel nodes by loop iteration count rather than time,
	//! since time can vary depending on the node
	double							last_time;
	unsigned int					cur_iter, last_bcast_iter, next_bcast_iter;
	double							iters_per_sec;
	
public:
	SimFramework(int argc, char **argv);
    ~SimFramework(void);
	
	// Dependency related functions
	PluginID registerPlugin(SimPlugin *new_plugin, const bool &is_active);
	void registerDependence(const PluginID &func_a, const PluginID &func_b, const DependenceType &dep_type);
	void writeDOT(std::ostream &os) const;

	// Simulation execution related functions
	void init(void);
	void run(void);
	void finish(void);
	
	// Multiprocessor related functions
	int getNodeRank(void) const { return node_rank; };
	int getWorldSize(void) const { return world_size; };
	bool isRootNode(void) const { return (node_rank == ROOT_NODE_RANK); };
#ifdef HAVE_MPI
	void barrier(void) { startTimer(barrier_timer); MPI_Barrier(MPI_COMM_WORLD); stopTimer(barrier_timer); };
#else
	void barrier(void) {};
#endif
	
	// Random number generation functions
	double randDouble(void) const;
	float randFloat(void) const;
	int randInt(const int &max_int) const;
	
	// Output related functions
	std::ostream &console(void) const { return std::cout; };
	std::ostream &errConsole(void) const { return std::cerr; };
	void disableConsole(void) { std::cout.setstate(std::ios_base::badbit); };
	
	// Loop timing functions
	unsigned int itersPerSecond(void) const { return (unsigned int)iters_per_sec; };
};

#endif
