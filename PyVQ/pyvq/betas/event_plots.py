#!/usr/bin/env python

# ---------------------- IMPORTS -------------------
from __future__ import print_function
import math
import sys
import argparse
import quakelib
import gc
import operator
scipy_available = True
try:
    import scipy.stats
except ImportError:
    scipy_available = False
matplotlib_available = True
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib_available = False
numpy_available = True
try:
    import numpy as np
except ImportError:
    numpy_available = False
h5py_available = True
try:
    import h5py
except ImportError:
    h5py_available = False
# ----------------------         -------------------


# ======= SINGLE EVENT I/O ============================================
def read_events(sim_file, event_numbers=None):
    # TODO: Add event filters
    with h5py.File(sim_file) as vq_data:
        events = vq_data['events'][()]
    # If event_numbers specified, only return those events
    if event_numbers is not None:
        # Handle single events separately
        if isinstance(event_numbers, int): 
            events = np.core.records.fromarrays(zip(*filter(lambda x: x['event_number'] == event_numbers, events)), dtype=events.dtype)
        else:
            d_type = events.dtype
            events = np.core.records.fromarrays(zip(*filter(lambda x: x['event_number'] in event_numbers, events)), dtype=d_type)	
	return events

def read_sweeps(sim_file, event_number=0, block_ids=None):
	# Read sweeps sequence for multiple blocks (unless block_id specified) in a single event.
	with h5py.File(sim_file) as vq_data:
		sweep_range = [vq_data['events'][event_number]['start_sweep_rec'],
                       vq_data['events'][event_number]['end_sweep_rec']]
		sweeps = vq_data['sweeps'][sweep_range[0]:sweep_range[1]][()]
	# If block_id specified, only return those sweeps for that block
	if block_ids is not None:
		d_type = sweeps.dtype
		sweeps = np.core.records.fromarrays(zip(*filter(lambda x: x['block_id'] in block_ids, sweeps)), dtype=d_type)	
	return sweeps

def parse_sweeps(sim_file=None, block_id=None, event_number=0, do_print=True, sweeps=None):
    # Read sweep data if not provided
	if sweeps is None: sweeps = read_sweeps(sim_file, block_id=block_id, event_number=event_number)
	# Grab data
	data = [[rw['sweep_number'], rw['block_id'], rw['block_slip'], rw['shear_init'],
             rw['shear_final'], rw['normal_init'],rw['normal_final'], 
             (rw['shear_final']-rw['shear_init'])/rw['shear_init'], 
             (rw['normal_final']-rw['normal_init'])/rw['normal_init']] for rw in sweeps]
	if do_print:
		for rw in data: print(rw)
	cols = ['sweep_number', 'block_id', 'block_slip', 'shear_init', 
            'shear_final', 'normal_init', 'normal_final', 'shear_change', 'normal_change']
	return np.core.records.fromarrays(zip(*data), names=cols, formats = [type(x).__name__ for x in data[0]])


# ======= SIM DATA CLASSES ===========================================
class Sweeps(object):
    def __init__(self, sim_file, event_number=0, block_ids=None):
        self.sweeps = read_sweeps(sim_file, event_number=event_number, block_ids=block_ids)
        self.sweep_data = parse_sweeps(sweeps=self.sweeps, do_print=False)
        self.block_ids = self.sweep_data['block_id'].tolist()
        self.mag = read_events(sim_file,event_numbers=event_number)['event_magnitude'][0]
        self.event_number = event_number
        print("Read event {} sweeps from {}".format(event_number,sim_file))
        # we could also, at this point, parse out the individual block sequences, maybe make a class Block().
    #
    def plot_event_block_slips(self, block_ids=None, fignum=0):
        block_ids = self.check_block_ids_list(block_ids)
        plt.figure(fignum)
        plt.clf()
        for block_id in block_ids:
            rws = np.core.records.fromarrays(zip(*filter(lambda x: x['block_id']==block_id, self.sweep_data)), dtype=self.sweep_data.dtype)
            plt.semilogy(rws['sweep_number'], rws['block_slip'], '.-', label=block_id)
        if len(block_ids) <= 10:
            plt.legend(loc='best', numpoints=1,fontsize=8,ncol=3,handlelength=2,handletextpad=1)
        plt.title('Event {} (M={:.2f}) slips for {} blocks'.format(self.event_number,self.mag,len(block_ids)))
        plt.xlabel('sweep number')
        plt.ylabel('slip [m]')
        min_sweep = 0
        max_sweep = int(max(self.sweep_data['sweep_number']))
        if max(self.sweep_data['sweep_number']) < 3:
            max_sweep += 1
        ticks = range(max_sweep+1)
        plt.xticks(ticks,[str(tick) for tick in ticks])
        plt.xlim(min_sweep, max_sweep)
    #
    def plot_stress_drop(self, block_ids=None, fignum=0, shear=True):
        block_ids = self.check_block_ids_list(block_ids)
        #
        plt.figure(fignum)
        plt.clf()
        #
        for block_id in block_ids:
            rws = np.core.records.fromarrays(zip(*filter(lambda x: x['block_id']==block_id, self.sweep_data)), dtype=self.sweep_data.dtype)
            if shear: 
                plt.plot(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
            else: 
                plt.plot(rws['sweep_number'], rws['normal_change'], '.-', label='block_id: %d' % block_id)
		plt.plot([min(self.sweep_data['sweep_number']), max(self.sweep_data['sweep_number'])], [0., 0.], 'k-')
		plt.legend(loc=0, numpoints=1)
        if shear: 
            plt.title('Block shear_stress drop sequences')
        else: 
            plt.title('Block normal_stress drop sequences')
        plt.xlabel('sweep number')
        plt.ylabel('fractional stress change')
    #    
    def check_block_ids_list(self, block_ids):
        # Make sure the block_ids are a list
        if block_ids is None: block_ids=self.block_ids
        if isinstance(block_ids, float): block_ids=[int(block_ids)]
        if isinstance(block_ids, int): block_ids = [block_ids]
        return block_ids
    #
class Events(object):
    def __init__(self, sim_file):
        # TODO: Add event filters
        self.events = read_events(sim_file)
        self.block_ids = self.sweep_data['block_id'].tolist()
        print("Reading event {} sweeps from {}".format(event_number,sim_file))
        # we could also, at this point, parse out the individual block sequences, maybe make a class Block().


# ============================ TEMP. RUNNER ===================

# TODO: Change this to a small sim and add it to github
SIM_FILE = "../../../../Desktop/RUNNING/UCERF2/events_ALLCAL2_VQmeshed_3km_EQSim_StressDrops_1600yr_22June2015.h5"
EVENT_NUM = 504
BLOCK_IDS = None
sim_sweeps = Sweeps(SIM_FILE, event_number=EVENT_NUM, block_ids=BLOCK_IDS)
sim_sweeps.plot_event_block_slips()
savename = "../../../../VQScripts/event_{}_slips{}.png".format(sim_sweeps.event_number,SIM_FILE.split("/")[-1].split(".")[0].split("events")[-1])
plt.savefig(savename,dpi=100)













