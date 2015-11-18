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
    import matplotlib.colors as mcolor
    import matplotlib.animation as manimation
    import matplotlib.colorbar as mcolorbar
    from mpl_toolkits.axes_grid1 import make_axes_locatable
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

# ======= h5py I/O ============================================
def read_events_h5(sim_file, event_numbers=None):
    # TODO: Add event filters
    with h5py.File(sim_file) as vq_data:
        events = vq_data['events'][()]
    # If event_numbers specified, only return those events
    if event_numbers is not None:
        # Handle single events separately
        if isinstance(event_numbers, int): 
            events = np.core.records.fromarrays(zip(*filter(lambda x: x['event_number'] == event_numbers, events)), dtype=events.dtype)
        else:
            events = np.core.records.fromarrays(zip(*filter(lambda x: x['event_number'] in event_numbers, events)), dtype=events.dtype)	
	return events

def read_sweeps_h5(sim_file, event_number=0, block_ids=None):
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

def parse_sweeps_h5(sim_file=None, block_id=None, event_number=0, do_print=True, sweeps=None):
    # Read sweep data if not provided
	if sweeps is None: sweeps = read_sweeps_h5(sim_file, block_id=block_id, event_number=event_number)
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
class Events:
    def __init__(self, sim_file):
        # TODO: Add event filters
        self.events = read_events_h5(sim_file)
        print("Reading events from {}".format(sim_file))


class Sweeps:
    def __init__(self, sim_file, event_number=0, block_ids=None):
        self.sweeps = read_sweeps_h5(sim_file, event_number=event_number, block_ids=block_ids)
        self.sweep_data = parse_sweeps_h5(sweeps=self.sweeps, do_print=False)
        self.block_ids = self.sweep_data['block_id'].tolist()
        self.mag = read_events_h5(sim_file,event_numbers=event_number)['event_magnitude'][0]
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
    def plot_stress_changes(self, block_ids=None, fignum=0, shear=True,log=False,max_val=None):
        block_ids = self.check_block_ids_list(block_ids)
        #
        plt.figure(fignum)
        plt.clf()
        #
        for block_id in block_ids:
            rws = np.core.records.fromarrays(zip(*filter(lambda x: x['block_id']==block_id, self.sweep_data)), dtype=self.sweep_data.dtype)
            if shear: 
                if not log:
                    plt.plot(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
                else:
                    plt.semilogy(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
            else: 
                if not log:
                    plt.plot(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
                else:
                    plt.semilogy(rws['sweep_number'], rws['shear_change'], '.-', label=block_id)
		plt.plot([min(self.sweep_data['sweep_number']), max(self.sweep_data['sweep_number'])], [0., 0.], 'k-')
        if len(block_ids) <= 10:
            plt.legend(loc='best', numpoints=1,fontsize=8,ncol=3,handlelength=2,handletextpad=1)
        if shear: 
            plt.title('Event {} (M={:.2f}) shear stress changes for {} blocks'.format(self.event_number,self.mag,len(block_ids)))
        else: 
            plt.title('Event {} (M={:.2f}) normal stress changes for {} blocks'.format(self.event_number,self.mag,len(block_ids)))
        plt.xlabel('sweep number')
        plt.ylabel('fractional stress change')
        min_sweep = 0
        max_sweep = int(max(self.sweep_data['sweep_number']))
        if max(self.sweep_data['sweep_number']) < 3:
            max_sweep += 1
        ticks = range(max_sweep+1)
        plt.xticks(ticks,[str(tick) for tick in ticks])
        plt.xlim(min_sweep, max_sweep)
        if max_val is not None: plt.ylim(-max_val,max_val)
    #    
    def check_block_ids_list(self, block_ids):
        # Make sure the block_ids are a list
        if block_ids is None: block_ids=self.block_ids
        if isinstance(block_ids, float): block_ids=[int(block_ids)]
        if isinstance(block_ids, int): block_ids = [block_ids]
        return block_ids



def event_stress_movie(model, event_output_text_file, savefile, plotting, FPS=3, DPI=100):

    event_data = np.genfromtxt(event_output_text_file, dtype=[('sweep_num','int'),('block_id','int'), ('shear_stress','f8'), ('normal_stress','f8'), ('cff','f8'), ('stress_drop','f8')])
    split_data = np.split(event_data, np.unique(event_data['sweep_num']).shape[0])

    # Currently only works for perfectly rectangular faults
    # Currently only plotting the elements on the triggering section
    triggerID = split_data[0][np.where(split_data[0]['sweep_num']==0)][0]['block_id']
    num_sweeps = len(split_data)
    print("Read "+str(num_sweeps)+" sweeps.")
    
    sectionID = model.element(triggerID).section_id()
    ele_length = np.sqrt(model.create_sim_element(triggerID).area())
    triggerSecElements = [id for id in range(model.num_elements()) if model.element(id).section_id() == sectionID]
    sec_name = model.section(sectionID).name()
    min_id    = triggerSecElements[0]
      
    max_val = max(event_data[plotting]/np.abs(np.mean(event_data['stress_drop'])))
    min_val = min(event_data[plotting]/np.abs(np.mean(event_data['stress_drop'])))
    shear_bound = max(np.abs(max_val), np.abs(min_val))
            
    section_length = model.section_length(sectionID)
    section_depth = abs(model.section_max_depth(sectionID))
    num_elements_down = int(round(section_depth/ele_length))
    num_elements_across = int(round(section_length/ele_length))
    assert(len(triggerSecElements) == num_elements_across*num_elements_down)
    element_grid = np.zeros((num_elements_down,num_elements_across))
    fig = plt.figure()
    ax = plt.gca()
    if min_val > 0:
        cmap = plt.get_cmap('Reds')
        norm = mcolor.Normalize(vmin=0, vmax=shear_bound)
    else: 
        cmap = plt.get_cmap('seismic')
        norm = mcolor.Normalize(vmin=-shear_bound, vmax=shear_bound)
    
    # Initialize movie writing stuff
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='VQ event', artist='Matplotlib',comment='Testing.')
    writer = FFMpegWriter(fps=FPS, metadata=metadata)
    
    plt.xlabel("along strike")
    plt.ylabel("down dip")
    plt.title("Virtual Quake Event",fontsize=11)
    plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    plt.figtext(0.97, 0.6, '{} [units are mean stress drop]'.format(plotting), rotation='vertical')
    
    # Colorbar
    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.1)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
    
    # Draw the arrow in the rake direction
    mean_rake = 0
    for id in triggerSecElements: mean_rake += model.element(id).rake()/len(triggerSecElements) 
    arrow_tail = np.array([0.13, 0.1])
    arrow_length = 0.08
    arrow_head = np.array([arrow_length*np.cos(mean_rake), arrow_length*np.sin(mean_rake)])
    arrow_head += arrow_tail  #vector addition
    plt.annotate("", xy=arrow_head, xytext=arrow_tail, arrowprops=dict(arrowstyle="->", lw=2), xycoords="figure fraction")
    plt.figtext(0.03, 0.05, 'Rake Direction\n\n\n', bbox={'facecolor':'cyan', 'pad':8, 'alpha':0.3})

    with writer.saving(fig, savefile, DPI):
        # Create the first frame of zero slip
        this_plot = ax.imshow(element_grid, cmap=cmap,origin='upper',interpolation='none',norm=norm)
        writer.grab_frame()
        for sweep_num in range(num_sweeps):
            if sweep_num == 0: sys.stdout.write("Generating frames...")
            if num_sweeps>10: 
                if sweep_num%int(num_sweeps/10.0)==0: sys.stdout.write("...{:.1f}%".format(100*sweep_num/float(num_sweeps-1)))
            sys.stdout.flush()
            # Using here the fact that VQ element numbering goes from top (near surface) to bottom,
            #   then makes a step down the strike (length) once you reach the bottom.
            this_sweep = event_data[ np.where(event_data['sweep_num']==sweep_num) ]
            for row in this_sweep:
                ele_id = int(row['block_id'])
                # Only plotting the elements on the triggering fault
                if model.element(ele_id).section_id() == sectionID:
                    grid_row = int((ele_id-min_id)%num_elements_down)
                    grid_col = int((ele_id-min_id)/num_elements_down)
                    element_grid[grid_row,grid_col] = row[plotting]/np.abs(row['stress_drop'])
                else:
                    sys.stdout.write("\nElement {} involved but not on triggering fault.".format(ele_id))
            # Update the colors
            this_plot.set_data(element_grid)
            # Time stamp
            plt.figtext(0.03, 0.9, 'Sweep: {:03d}'.format(sweep_num), bbox={'facecolor':'yellow', 'pad':8})
            writer.grab_frame()
    sys.stdout.write("\n>> Movie saved to {}\n".format(savefile))


# ============================ TEMP. RUNNER ===================
# Example usage below
"""
SIM_FILE = "../../../../Desktop/RUNNING/UCERF2/events_ALLCAL2_VQmeshed_3km_EQSim_StressDrops_1600yr_22June2015.h5"
EVENT_NUM = 1541 #13, 948, 504, 1541
BLOCK_IDS = None
sim_sweeps = Sweeps(SIM_FILE, event_number=EVENT_NUM, block_ids=BLOCK_IDS)
# ---- plot slips ---------

sim_sweeps.plot_event_block_slips()
savename = "../../../../VQScripts/event_{}_slips{}.png".format(sim_sweeps.event_number,SIM_FILE.split("/")[-1].split(".")[0].split("events")[-1])
plt.savefig(savename,dpi=100)

# ---- plot stresses ---------
sim_sweeps.plot_stress_changes()
savename = "../../../../VQScripts/event_{}_shear_changes{}.png".format(sim_sweeps.event_number,SIM_FILE.split("/")[-1].split(".")[0].split("events")[-1])
plt.savefig(savename,dpi=100)


# ---- plot CFF, stresses, from simulation internal data ---------
# ---- To use this, one must hack VQ to output the following during every sweep for every element
# ----  sweep_num  element_num  shear_stress  normal_stress  CFF   stress_drop

model_file = "/Users/kasey/Desktop/RUNNING/single_small_fault_3000m.txt"
plotting = 'cff'  # 'cff' or 'shear_stress' or 'normal_stress'
event_output_text_file = "/Users/kasey/Desktop/RUNNING/small_8x4_fault_eq_output.txt"
savefile = "/Users/kasey/VQScripts/small_8x4_"+plotting+"_movie.mp4"
model = quakelib.ModelWorld()
model.read_file_ascii(model_file)

event_stress_movie(model, event_output_text_file, savefile, plotting, FPS=1, DPI=100)
"""






