#!/usr/bin/env python

from __future__ import print_function

import math
import sys
import argparse
import quakelib
import gc
import operator
import os

scipy_available = True
try:
    import scipy.stats
except ImportError:
    scipy_available = False

matplotlib_available = True
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import matplotlib.font_manager as mfont
    import matplotlib.colors as mcolor
    import matplotlib.colorbar as mcolorbar
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches
    from PIL import Image 
    #TODO: Move this guy

    # we only want to execute this in the __main__ part of the script, so we can also run plotting scripts interactively.
    #plt.switch_backend('agg') #Required for map plots

    #plt.switch_backend('agg') #Required for map plots
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
    
# ----------------- Global constants -------------------------------------------
# Kasey: These are only relevent for few-element field plots
#LAT_LON_DIFF_FACTOR = 1.333 
#MIN_LON_DIFF = 0.01   # 1 corresponds to ~ 100km at lat,lon = (40.35, -124.85)
#MIN_LAT_DIFF = MIN_LON_DIFF/LAT_LON_DIFF_FACTOR   # 0.8 corresponds to ~ 100km at lat,lon = (40.35, -124.85)
#MIN_FIT_MAG  = 5.0     # lower end of magnitude for fitting freq_mag plot with b=1 curve

#-------------------------------------------------------------------------------
# Given a set of maxes and mins return a linear value betweem them.
# Used to compute cutoff for field value evaluation, cutoff scales with
# number of involved elements for FieldPlotter instances.
#-------------------------------------------------------------------------------
def linear_interp(x, x_min, x_max, y_min, y_max):
    return ((y_max - y_min)/(x_max - x_min) * (x - x_min)) + y_min
    
def calculate_averages(x,y,log_bin=False,num_bins=None):
    if num_bins is None:
        num_bins = math.floor(len(x)/100)
        if num_bins < 20:
            num_bins = 20
        elif num_bins > 100:
            num_bins = 100
    x = np.array(x)
    y = np.array(y)
    #if np.min(x) == 0:
    #    bin_min = 1
    #else:
    if log_bin: bin_min = math.floor(math.log(np.min(x),10))
    else: bin_min = math.floor(np.min(x))
    if log_bin: bin_max = math.ceil(math.log(np.max(x),10))
    else: bin_max = math.ceil(np.max(x))
    if log_bin: bins = np.logspace(bin_min,bin_max,num=num_bins)
    else: bins = np.linspace(bin_min,bin_max,num=num_bins)
    inds = np.digitize(x, bins)
    binned_data = {}
    for n, i in enumerate(inds):
        try:
            binned_data[i].append(y[n])
        except KeyError:
            binned_data[i] = [y[n]]
    x_ave = []
    y_ave = []
    for k in sorted(binned_data.keys()):
        if k != 0:
            x_ave.append(0.5*(bins[k-1]+bins[k]))
            y_ave.append(sum(binned_data[k])/float(len(binned_data[k])))
    return x_ave, y_ave

class SaveFile:
    def event_plot(self, event_file, plot_type, min_mag):
        min_mag = str(min_mag)
        # Remove any folders in front of model_file name
        if len(event_file.split("/")) > 1:
            event_file = event_file.split("/")[-1]
        # Add tags to convey the subsets/cuts being made
        add=""
        if min_year is not None: add+="_yearMin"+str(int(min_year))    
        if max_year is not None: add+="_yearMax"+str(int(max_year))
        if args.use_sections is not None:
            for sec in args.use_sections:
                add+="_"+geometry.model.section(sec).name()
        if min_mag is not None: 
            # e.g. min_mag = 7.5, filename has '7-5'
            if len(min_mag.split(".")) > 1:
                add += "_minMag_"+min_mag.split(".")[0]+"-"+min_mag.split(".")[1]
            else:
                add += "_minMag_"+min_mag

        return plot_type+add+"_"+event_file.split(".")[0]+".png"
        
    def field_plot(self, model_file, field_type, uniform_slip, event_id):
        # Remove any folders in front of model_file name
        if len(model_file.split("/")) > 1:
            model_file = model_file.split("/")[-1]
        if uniform_slip is None and event_id is not None:
            return model_file.split(".")[0]+"_"+field_type+"_event"+str(event_id)+".png"
        elif uniform_slip is not None and event_id is None:
            return model_file.split(".")[0]+"_"+field_type+"_uniform_slip"+str(int(uniform_slip))+"m.png"
        else:
            raise "Must specify either uniform_slip or event_id"
            
    def greens_plot(self, name, field_type, slip):
        return "greens_"+field_type+"_"+name+"_slip"+str(int(slip))+"m.png"
            
    def trace_plot(self, model_file):
        # Remove any folders in front of model_file name
        if len(model_file.split("/")) > 1:
            model_file = model_file.split("/")[-1]
        return "traces_"+model_file.split(".")[0]+".png"
        
    def diagnostic_plot(self, event_file, plot_type, min_year=None, max_year=None):
        # Remove any folders in front of model_file name
        if len(event_file.split("/")) > 1:
            event_file = event_file.split("/")[-1]
        add=""
        if min_year is not None: add+="_yearMin"+str(int(min_year))    
        if max_year is not None: add+="_yearMax"+str(int(max_year))
        if args.use_sections is not None:
            for sec in args.use_sections:
                add+="_"+geometry.model.section(sec).name()
            
        return plot_type+"_diagnostic"+add+"_"+event_file.split(".")[0]+".png"
    
    def event_kml_plot(self, event_file, event_id):
        if len(event_file.split("/")) > 1:
            event_file = event_file.split("/")[-1]
        event_file = event_file.split("events")[-1]
        return "event_"+str(event_id)+event_file.split(".")[0]+".kml"
    

class MagFilter:
    def __init__(self, min_mag=None, max_mag=None):
        self._min_mag = min_mag if min_mag is not None else -float("inf")
        self._max_mag = max_mag if max_mag is not None else float("inf")

    def test_event(self, event):
        return (event.getMagnitude() >= self._min_mag and event.getMagnitude() <= self._max_mag)

    def plot_str(self):
        label_str = ""
# TODO: change to <= character
        if self._min_mag != -float("inf"): label_str += str(self._min_mag)+"<"
        if self._max_mag != float("inf"): label_str += "M<"+str(self._max_mag)
        return label_str

class YearFilter:
    def __init__(self, min_year=None, max_year=None):
        self._min_year = min_year if min_year is not None else -float("inf")
        self._max_year = max_year if max_year is not None else float("inf")

    def test_event(self, event):
        return (event.getEventYear() >= self._min_year and event.getEventYear() <= self._max_year)

    def plot_str(self):
        label_str = ""
# TODO: change to <= character
        if self._min_year != -float("inf"): label_str += str(self._min_year)+"<"
        if self._max_year != float("inf"): label_str += "year<"+str(self._max_year)
        return label_str

class EventNumFilter:
    def __init__(self, min_event_num=None, max_event_num=None):
        self._min_event_num = min_event_num if min_event_num is not None else -sys.maxint
        self._max_event_num = max_event_num if max_event_num is not None else sys.maxint

    def test_event(self, event):
        return (event.getEventNumber() >= self._min_event_num and event.getEventNumber() <= self._max_event_num)

    def plot_str(self):
        label_str = ""
# TODO: change to <= character
        if self._min_event_num != -sys.maxint: label_str += str(self._min_event_num)+"<"
        if self._max_event_num != sys.maxint: label_str += "event num<"+str(self._max_event_num)
        return label_str

class SectionFilter:
    def __init__(self, geometry, section_list):
        self._section_list = section_list
        self._elem_to_section_map = {elem_num: geometry.model.element(elem_num).section_id() for elem_num in range(geometry.model.num_elements())}

    def test_event(self, event):
        event_elements = event.getInvolvedElements()
        for elem_num in event_elements:
            elem_section = self._elem_to_section_map[elem_num]
            if elem_section in self._section_list: return True

        return False

    def plot_str(self):
        return ""

class TriggerSectionFilter:
    def __init__(self, geometry, section_list):
        self._section_list = section_list
        self._elem_to_section_map = {elem_num: geometry.model.element(elem_num).section_id() for elem_num in range(geometry.model.num_elements())}

    def test_event(self, event):
        triggerID = event.getEventTrigger()
        elem_section = self._elem_to_section_map[triggerID]
        if elem_section in self._section_list: return True

        return False

    def plot_str(self):
        return ""
        
class SlipFilter:
    def __init__(self, min_slip=None, max_slip=None):
        self._min_slip = min_slip if min_slip is not None else -float("inf")
        self._max_slip = max_slip if max_slip is not None else float("inf")

    def test_event(self, event):
        return (event.calcMeanSlip() >= self._min_slip and event.calcMeanSlip() <= self._max_slip)

    def plot_str(self):
        label_str = ""
# TODO: change to <= character
        if self._min_slip != -float("inf"): label_str += str(self._min_slip)+"<"
        if self._max_slip != float("inf"): label_str += "year<"+str(self._max_slip)
        return label_str
        
class Geometry:
    def __init__(self, model_file=None, model_file_type=None):
        if model_file is not None:
            self.model = quakelib.ModelWorld()
            if model_file_type =='text' or model_file.split(".")[-1] == 'txt':
                self.model.read_file_ascii(model_file)
            elif model_file_type == 'hdf5' or model_file.split(".")[-1] == 'h5' or model_file.split(".")[-1] == 'hdf5':
                self.model.read_file_hdf5(model_file)
            else:
                raise "Must specify --model_file_type, either hdf5 or text"
            self._elem_to_section_map = {elem_num: self.model.element(elem_num).section_id() for elem_num in self.model.getElementIDs()}
        else:
            if args.use_sections:
                raise "Model file required if specifying fault sections."
            return None

    def get_fault_traces(self):
        traces_lat_lon = {}
        ele_ids = self.model.getElementIDs()
        for eid in ele_ids:
            sid      = self._elem_to_section_map[eid]
            vids     = [self.model.element(eid).vertex(i) for i in range(3)]
            vertices = [self.model.vertex(vid) for vid in vids]
            for vert in vertices:
                if vert.is_trace():
                    lat = vert.lld().lat()
                    lon = vert.lld().lon()
                    try:
                        traces_lat_lon[sid].append((lat,lon))
                    except KeyError:
                        traces_lat_lon[sid] = [(lat,lon)]
                    #break
        return traces_lat_lon
        
    def get_slip_rates(self, elements):
        # Convert slip rates from meters/second to meters/(decimal year)
        CONVERSION = 3.15576*pow(10,7) 
        return {id:self.model.element(id).slip_rate()*CONVERSION for id in elements}
        
    def get_slip_time_series(self, events, elements=None, min_year=None, max_year=None, DT=None):
        # slip_time_series    = dictionary indexed by block_id with entries being arrays of absolute slip at each time step
        # Get slip rates for the elements
        slip_rates = self.get_slip_rates(elements)
        #Initialize blocks with 0.0 slip at time t=0.0
        slip_time_series  = {id:[0.0] for id in elements}
        # Grab the events data
        event_years = events.event_years()
        event_numbers = events.event_numbers()
        #Initialize time steps to evaluate slip    
        time_values = np.arange(min_year+DT, max_year+DT, DT)
        for k in range(len(time_values)):
            if k>0:
                # current time in simulation
                right_now = time_values[k]
                # back slip all elements by subtracting the slip_rate*dt
                for block_id in slip_time_series.keys():
                    last_slip = slip_time_series[block_id][k-1]
                    this_slip = slip_rates[block_id]*DT
                    slip_time_series[block_id].append(last_slip-this_slip)
                # check if any elements slip as part of simulated event in the window of simulation time
                # between (current time - DT, current time), add event slips to the slip at current time 
                # for elements involved
                for j in range(len(event_numbers)):
                    evid    = event_numbers[j]
                    ev_year = event_years[j]
                    if right_now-DT < ev_year <= right_now:
                        event_element_slips = events.get_event_element_slips(evid)
                        for block_id in event_element_slips.keys():
                            try:
                                slip_time_series[block_id][k] += event_element_slips[block_id]
                                #sys.stdout.write("element {} slips {} in event {}\n".format(block_id,event_element_slips[block_id],evid))
                                #sys.stdout.flush()
                            except KeyError:
                                pass # Ignore event elements that we are not asked for (in elements)                            
        return slip_time_series


# ======= h5py I/O ============================================
def read_events_h5(sim_file, event_numbers=None):
    # TODO: Add event filters
    with h5py.File(sim_file) as vq_data:
        events = vq_data['events'][()]
    # If event_numbers specified, only return those events
    if event_numbers is not None:
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
    
    
class Events:
    def __init__(self, event_file, sweep_file = None):
        filetype = event_file.split('.')[-1].lower()
        event_file_type = "text" # default
        if filetype == 'h5' or filetype == 'hdf5': event_file_type = "hdf5"
        if event_file_type == "hdf5":
            # Reading in via QuakeLib
            #if not h5py_available:
            self._events = quakelib.ModelEventSet()
            self._events.read_file_hdf5(event_file)
            print("Read in events via QuakeLib from {}".format(event_file))
            # Reading via h5py
            #else:
            #    self._events = read_events_h5(event_file)
            #    print("Read in events via h5py from {}".format(event_file))
        elif event_file_type == "text" and sweep_file != None:
            self._events = quakelib.ModelEventSet()
            self._events.read_file_ascii(event_file, sweep_file)
        else:
            raise "event_file_type must be hdf5 or text. If text, a sweep_file is required."
            
        self._filtered_events = range(len(self._events))
        self._plot_str = ""

    def plot_str(self):
        return self._plot_str

    def set_filters(self, filter_list):
        self._filtered_events = [evnum for evnum in range(len(self._events))]
        self._plot_str = ""
        for cur_filter in filter_list:
            new_filtered_events = [evnum for evnum in self._filtered_events if cur_filter.test_event(self._events[evnum])]
            self._filtered_events = new_filtered_events
            self._plot_str += cur_filter.plot_str()
        if len(self._filtered_events) == 0:
            raise "No events matching filters found!"

    def interevent_times(self):
        event_times = [self._events[evnum].getEventYear() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]
        return [event_times[i+1]-event_times[i] for i in xrange(len(event_times)-1)]

    def event_years(self):
        return [self._events[evnum].getEventYear() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]

    def event_rupture_areas(self):
        return [self._events[evnum].calcEventRuptureArea() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]

    def event_magnitudes(self):
        return [self._events[evnum].getMagnitude() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]
    # TODO: Handle  NaN magnitudes on the C++ side

    def event_numbers(self):
        return [evnum for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]
        
    def event_mean_slip(self):
        return [self._events[evnum].calcMeanSlip() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]
        
    def get_event_element_slips(self, evnum):
        element_ids = self._events[evnum].getInvolvedElements()
        return {ele_id:self._events[evnum].getEventSlip(ele_id) for ele_id in element_ids}
        
    def get_event_sections(self, evnum, geometry):
        sec_ids = [geometry.model.element(eid).section_id() for eid in self._events[evnum].getInvolvedElements()]
        # Get unique section ids by converting to a set, then back to a list for ease of use
        return list(set(sec_ids))
        
    def get_ids_largest_events(self, num_events):
        mags = {evnum:self._events[evnum].getMagnitude() for evnum in self._filtered_events if self._events[evnum].getMagnitude() != float("-inf")}
        # Sort by decreasing magnitude
        mags_sorted = list(reversed(sorted(mags.items(), key=operator.itemgetter(1))))
        ev_ids = [mags_sorted[i][0] for i in range(len(mags_sorted))]
        return ev_ids[:num_events]
        
    def event_summary(self, evnums):
        mags = [self._events[evnum].getMagnitude() for evnum in evnums if self._events[evnum].getMagnitude() != float("-inf")]
        areas = [self._events[evnum].calcEventRuptureArea() for evnum in evnums]
        times = [self._events[evnum].getEventYear() for evnum in evnums]
        slips = [self._events[evnum].calcMeanSlip() for evnum in evnums]
        print("=======================================================")
        print("evid\tyear\t\tmag\tarea[km^2]\tslip[m]")
        print("-------------------------------------------------------")
        for k in range(len(evnums)):
            print("{}\t{:.1f}\t\t{:.3f}\t{:.4f}\t{:.4f}".format(evnums[k],times[k],mags[k],areas[k]*pow(10,-6),slips[k]))
        print("-------------------------------------------------------\n")
            
    def largest_event_summary(self, num_events):
        evnums = self.get_ids_largest_events(num_events)
        self.event_summary(evnums)
    
    def event_initial_shear_stresses(self):
        return [self._events[evnum].getShearStressInit() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]

    def event_final_shear_stresses(self):
        return [self._events[evnum].getShearStressFinal() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]        
                        
    def event_initial_normal_stresses(self):
        return [self._events[evnum].getNormalStressInit() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]

    def event_final_normal_stresses(self):
        return [self._events[evnum].getNormalStressFinal() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]  
        
    def number_of_sweeps(self):
        return [self._events[evnum].getNumRecordedSweeps() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())] 

class Sweeps:
    # A class for reading/analyzing data from the event sweeps
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
        
class GreensPlotter:
    # Plot Okubo Greens functions for a single fault element
    def __init__(self, field_type, cbar_max=None, levels=None, Nx=690, Ny=422, Xmin=-5000, Xmax=15000, Ymin=-10000, Ymax=10000, L=10000, W=10000, DTTF=1000, slip=5, dip=90, _lambda=3.2e10, _mu=3.0e10, rake=0, g0=None):
        if g0 is None: self.g0 = 9.80665
        else: self.g0 = g0
        self.field_type = field_type.lower()
        self.block = quakelib.Okada()
        self.slip = slip
        self.dip = dip*np.pi/180.0
        self.C = DTTF + W*np.sin(self.dip)
        self.cbar_max = cbar_max
        self.levels = levels
        self.L = L
        self.W = W
        self.rake = rake*np.pi/180.0
        self.US = slip*np.cos(self.rake)
        self.UD = slip*np.sin(self.rake)
        self.UT = 0.0
        self.X = np.linspace(Xmin,Xmax,num=Nx)
        self.Y = np.linspace(Ymin,Ymax,num=Ny)
        self.XX, self.YY = np.meshgrid(self.X, self.Y)
        self.field = np.zeros(self.XX.shape)
        self._lambda = _lambda
        self._mu = _mu
        self.cmap = plt.get_cmap('seismic')
        
    def compute_field(self):
        if self.field_type == 'gravity':
            for i in range(self.XX.shape[0]):
                for j in range(self.XX.shape[1]):
                    loc       = quakelib.Vec2(self.XX[i][j], self.YY[i][j])
                    self.field[i][j] = self.block.calc_dg(loc, self.C, self.dip, self.L, self.W, self.US, self.UD, self.UT, self._lambda, self._mu)*pow(10,8)
        elif self.field_type == 'dilat_gravity':
            for i in range(self.XX.shape[0]):
                for j in range(self.XX.shape[1]):
                    loc       = quakelib.Vec2(self.XX[i][j], self.YY[i][j])
                    self.field[i][j] = self.block.calc_dg_dilat(loc, self.C, self.dip, self.L, self.W, self.US, self.UD, self.UT, self._lambda, self._mu)*pow(10,8)
        elif self.field_type == 'potential':
            for i in range(self.XX.shape[0]):
                for j in range(self.XX.shape[1]):
                    loc       = quakelib.Vec3(self.XX[i][j], self.YY[i][j], 0.0)
                    self.field[i][j] = self.block.calc_dV(loc, self.C, self.dip, self.L, self.W, self.US, self.UD, self.UT, self._lambda, self._mu)
        elif self.field_type == 'geoid':
            for i in range(self.XX.shape[0]):
                for j in range(self.XX.shape[1]):
                    loc       = quakelib.Vec3(self.XX[i][j], self.YY[i][j], 0.0)
                    self.field[i][j] = -self.block.calc_dV(loc, self.C, self.dip, self.L, self.W, self.US, self.UD, self.UT, self._lambda, self._mu)/self.g0
        elif self.field_type == 'displacement':
            for i in range(self.XX.shape[0]):
                for j in range(self.XX.shape[1]):
                    loc       = quakelib.Vec3(self.XX[i][j], self.YY[i][j], 0.0)
                    self.field[i][j] = self.block.calc_displacement_vector(loc, self.C, self.dip, self.L, self.W, self.US, self.UD, self.UT, self._lambda, self._mu)[2]
    
    def plot_field(self, output_file, no_labels=False, cbar_loc='top', tick_font=18, frame_font=18, x_ticks=True):
        ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=tick_font)
        framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=frame_font)
        # PAD 40 for 12pt font, 52 for 18pt
        PAD = 52
        
        if self.field_type == 'gravity': cbar_lab = r'total gravity changes $[\mu gal]$'
        elif self.field_type == 'dilat_gravity': cbar_lab = r'dilat. gravity changes $[\mu gal]$'
        elif self.field_type == 'displacement': cbar_lab = r'$\Delta h \ [m]$'
        elif self.field_type == 'potential': cbar_lab = 'Grav. potential changes'
        elif self.field_type == 'geoid': cbar_lab = r'Geoid height changes $[m]$'
        else: sys.exit("Field type not supported.")
        
        if self.cbar_max is not None:
            self.norm = mcolor.Normalize(vmin=-self.cbar_max, vmax=self.cbar_max)
        else:
            self.norm = None
        fig = plt.figure()
        fig_axes = plt.subplot(111)
        if self.levels is not None:
            img = plt.contourf(self.field, self.levels, cmap=self.cmap, norm=self.norm, extend='both', extent=[self.X.min()/1000.0,self.X.max()/1000.0,self.Y.min()/1000.0,
                                 self.Y.max()/1000.0])
        else:
            img = plt.imshow(self.field, origin = 'lower',interpolation='nearest',
                         extent=[self.X.min()/1000.0,self.X.max()/1000.0,self.Y.min()/1000.0,
                                 self.Y.max()/1000.0], cmap=self.cmap, norm=self.norm)
        img_ax = fig.gca()
        if not no_labels:
            img_ax.set_xlabel(r'along fault [$km$]',labelpad=-1, fontproperties=framelabelfont)
            img_ax.set_ylabel(r'[$km$]',labelpad=-5, fontproperties=framelabelfont)
        divider = make_axes_locatable(fig_axes)
        if cbar_loc=='top':
            cbar_ax      = divider.append_axes("top", size="5%",pad=0.02)                       
        else:
            cbar_ax      = divider.append_axes("bottom", size="5%",pad=0.02)
        cb = mcolorbar.ColorbarBase(cbar_ax, cmap=self.cmap, norm=self.norm, orientation='horizontal')
        if not no_labels:
            cbar_ax.set_xlabel(cbar_lab,labelpad=-PAD, fontproperties=framelabelfont)
        if cbar_loc=='bottom':
            PAD = 2.5 
            TOP = False
            BOTTOM = True 
        else:
            PAD = -.5        
            TOP = True
            BOTTOM = False
        cbar_ax.tick_params(axis='x',labelbottom=BOTTOM,labeltop=TOP,
                        bottom='off',top='off',right='off',left='off',pad=PAD)
        if self.field_type == "gravity" or self.field_type == "dilat_gravity":
            forced_ticks  = [int(num) for num in np.linspace(-self.cbar_max, self.cbar_max, len(cbar_ax.xaxis.get_ticklabels()))]
        else:
            forced_ticks  = [round(num, 3) for num in np.linspace(-self.cbar_max, self.cbar_max, len(cbar_ax.xaxis.get_ticklabels()))]
        cb_tick_labs    = [str(num) for num in forced_ticks]
        cb_tick_labs[0] = '<'+cb_tick_labs[0]
        cb_tick_labs[-1]= '>'+cb_tick_labs[-1]
        cbar_ax.set_xticklabels(cb_tick_labs)
        for label in img_ax.xaxis.get_ticklabels()+img_ax.yaxis.get_ticklabels():
            label.set_fontproperties(framelabelfont)  
        for label in cbar_ax.xaxis.get_ticklabels()+cbar_ax.yaxis.get_ticklabels():
            label.set_fontproperties(ticklabelfont)
        if not x_ticks:
            plt.setp(img_ax.xaxis.get_ticklabels(),visible=False)
        W_proj      = self.W*np.cos(self.dip)  #projected width of fault due to dip angle
        fault_proj  = mpl.patches.Rectangle((0.0,0.0),self.L/1000.0,W_proj/1000.0,
                                        ec='k',fc='none',fill=False,
                                        ls='solid',lw=4.0)
        fig_axes.add_patch(fault_proj)
        plt.savefig(output_file, dpi=100)
        print("----Greens function plot saved: "+output_file)
        plt.clf()

class TracePlotter:
    # Plot fault traces on a map
    def __init__(self, geometry, output_file, use_sections=None, small_model=False):
        self.small_model = small_model
        plot_height = 768.0
        max_map_width = 690.0
        max_map_height = 658.0
        map_res  = 'i'
        padding  = 0.08
        map_proj = 'cyl'
        # Read elements and slips into the SlippedElementList
        involved_sections = geometry.model.getSectionIDs()
        self.elements = quakelib.SlippedElementList()
        element_ids   = geometry.model.getElementIDs()
        for ele_id in element_ids:
            new_ele = geometry.model.create_slipped_element(ele_id)
            new_ele.set_slip(0.0)
            self.elements.append(new_ele)
        # Grab base Lat/Lon from fault model, used for lat/lon <-> xyz conversion
        base = geometry.model.get_base()
        self.base_lat = self.min_lat = base[0]
        self.base_lon = self.min_lon = base[1]
        self.min_lat, self.max_lat, self.min_lon, self.max_lon = geometry.model.get_latlon_bounds()
        # Expand lat/lon range in the case of plotting a few elements
        if self.small_model:
            self.min_lat = self.min_lat - MIN_LAT_DIFF*10
            self.max_lat = self.max_lat + MIN_LAT_DIFF*10
            self.min_lon = self.min_lon - MIN_LON_DIFF*10
            self.max_lon = self.max_lon + MIN_LON_DIFF*10
        # Adjust bounds for good framing on plot
        lon_range = self.max_lon - self.min_lon
        lat_range = self.max_lat - self.min_lat
        max_range = max((lon_range, lat_range))
        self.min_lon = self.min_lon - lon_range*padding
        self.min_lat = self.min_lat - lat_range*padding
        self.max_lon = self.max_lon + lon_range*padding
        self.max_lat = self.max_lat + lat_range*padding
        self.lat0, self.lon0  = (self.max_lat+self.min_lat)/2.0, (self.max_lon+self.min_lon)/2.0
        self.llcrnrlat = self.min_lat
        self.llcrnrlon = self.min_lon
        self.urcrnrlat = self.max_lat
        self.urcrnrlon = self.max_lon
        # We need a map instance to calculate the aspect ratio
        map = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        # Using the aspect ratio (h/w) to find the actual map width and height in pixels
        if map.aspect > max_map_height/max_map_width:
            map_height = max_map_height
            map_width = max_map_height/map.aspect
        else:
            map_width = max_map_width
            map_height = max_map_width*map.aspect
        # A conversion instance for doing the lat-lon to x-y conversions
        base_lld = quakelib.LatLonDepth(self.base_lat, self.base_lon, 0.0)
        self.convert = quakelib.Conversion(base_lld)
        self.lons_1d = np.linspace(self.min_lon, self.max_lon, num=int(map_width))
        self.lats_1d = np.linspace(self.min_lat, self.max_lat, num=int(map_height))
        _lons_1d = quakelib.FloatList()
        _lats_1d = quakelib.FloatList()
        for lon in self.lons_1d:
            _lons_1d.append(lon)
        for lat in self.lats_1d:
            _lats_1d.append(lat)
        # Set up the points for field evaluation, convert to xyz basis
        self.grid_1d = self.convert.convertArray2xyz(_lats_1d,_lons_1d)
        self.fault_traces_latlon = geometry.get_fault_traces()
        # Font/color presets
        font = mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='normal')
        font_bold = mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='bold')
        cmap = plt.get_cmap('seismic')
        water_color = '#4eacf4'
        boundary_color = '#000000'
        boundary_width = 1.0
        coastline_color = '#000000'
        coastline_width = 1.0
        country_color = '#000000'
        country_width = 1.0
        land_color = cmap(0.5)
        state_color = '#000000'
        state_width = 1.0
        river_width = 0.25
        fault_color = '#000000'
        event_fault_color = '#ff0000'
        fault_width = 0.5
        fault_width_bold = 2.5
        grid_color = '#000000'
        grid_width = 0.5
        num_grid_lines = 5
        map_resolution = map_res
        map_projection = map_proj
        plot_resolution = 72.0
        map_tick_color = '#000000'
        map_frame_color = '#000000'
        map_frame_width = 1
        map_fontsize = 26.0  # default 12   THIS IS BROKEN
        #---------------------------------------------------------------------------
        # m1, fig1 is all of the boundary data plus fault traces.
        #---------------------------------------------------------------------------
        self.m1 = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, 
            lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        mw = self.lons_1d.size
        mh = self.lats_1d.size

        if mh > mw:
            ph = 768.0
            pw = mw + 70.0 + 40.0
        else:
            pw = 790.0
            ph = mh + 70.0 + 40.0

        width_frac = mw/pw
        height_frac = mh/ph
        left_frac = 70.0/pw
        bottom_frac = 70.0/ph

        pwi = pw/plot_resolution
        phi = ph/plot_resolution
        
        #-----------------------------------------------------------------------
        # Set the map dimensions
        #-----------------------------------------------------------------------
        mw = self.lons_1d.size
        mh = self.lats_1d.size
        mwi = mw/plot_resolution
        mhi = mh/plot_resolution

        #-----------------------------------------------------------------------
        # Fig1 is the background land, ocean, and fault traces.
        #-----------------------------------------------------------------------
        fig1 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        #self.m1.ax = fig1.add_axes((0,0,1,1))
        self.m1.ax = fig1.add_axes((left_frac,bottom_frac,width_frac,height_frac))

        self.m1.drawmapboundary(
            color=boundary_color,
            linewidth=boundary_width,
            fill_color=water_color
        )
        self.m1.fillcontinents(
            color=land_color,
            lake_color=water_color
        )
        
        # draw coastlines, edge of map.
        self.m1.drawcoastlines(color=coastline_color, linewidth=coastline_width)

        # draw countries
        self.m1.drawcountries(linewidth=country_width, color=country_color)

        # draw states
        self.m1.drawstates(linewidth=state_width, color=state_color)

        # draw parallels.
        parallels = np.linspace(self.lats_1d.min(), self.lats_1d.max(), num_grid_lines+1)
        m1_parallels = self.m1.drawparallels(parallels, fontsize=map_fontsize, labels=[1,0,0,0], color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])

        # draw meridians
        meridians = np.linspace(self.lons_1d.min(), self.lons_1d.max(), num_grid_lines+1)
        m1_meridians = self.m1.drawmeridians(meridians, fontsize=map_fontsize, labels=[0,0,1,0], color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
        
        # Plot faults on lon-lat plot
        for sid, sec_trace in self.fault_traces_latlon.iteritems():
            sec_trace_lons = [lat_lon[1] for lat_lon in sec_trace]
            sec_trace_lats = [lat_lon[0] for lat_lon in sec_trace]
            
            trace_Xs, trace_Ys = self.m1(sec_trace_lons, sec_trace_lats)
            
            if use_sections is not None:
                if sid in use_sections:
                    linewidth = fault_width_bold
                else:
                    linewidth = fault_width
            else:  
                linewidth = fault_width_bold

            self.m1.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')

        fig1.savefig(output_file, format='png', dpi=plot_resolution)
        sys.stdout.write('Plot saved: {}\n'.format(output_file))
        sys.stdout.flush()
        
class FieldPlotter:
    def __init__(self, geometry, field_type, element_slips=None, event_id=None, event=None,
                cbar_max=None, levels=None, small_model=False, g0=None):
        if g0 is None: 
            self.g0 = 9.80665
        else:
            self.g0 = g0
        self.field_type = field_type.lower()
        self.small_model = small_model
        self.levels = levels
        self.event_id = event_id
        plot_height = 768.0
        max_map_width = 690.0
        max_map_height = 658.0
        map_res  = 'i'
        padding  = 0.08
        map_proj = 'cyl'
        self.norm = None
        # Define how the cutoff value scales if it is not explitly set.
        # Cutoff is the max distance away from elements to compute
        # the field given in units of element length.
        if self.field_type == 'gravity' or self.field_type == 'dilat_gravity' or self.field_type == 'potential' or self.field_type == 'geoid':
            self.cutoff_min_size = 20.0
            self.cutoff_min = 20.0
            self.cutoff_p2_size = 65.0
            self.cutoff_p2 = 90.0
        elif self.field_type == 'displacement' or self.field_type == 'insar':
            self.look_azimuth = None
            self.look_elevation = None
            self.cutoff_min_size = 20.0
            self.cutoff_min = 46.5
            self.cutoff_p2_size = 65.0
            self.cutoff_p2 = 90.0    
            self.dX = None
            self.dY = None
            self.dZ = None
            self.dX_min = sys.float_info.max
            self.dY_min = sys.float_info.max
            self.dZ_min = sys.float_info.max
            self.wavelength = 0.03
        # Read elements and slips into the SlippedElementList
        self.elements = quakelib.SlippedElementList()
        if event_id is None and event is None and element_slips is None:
            raise "Must specify event_id for event fields or element_slips (dictionary of slip indexed by element_id) for custom field."
        else:
            self.element_ids = element_slips.keys()
            self.element_slips = element_slips
        if len(self.element_slips) != len(self.element_ids):
            raise "Must specify slip for all elements."
        for ele_id in self.element_ids:
            new_ele = geometry.model.create_slipped_element(ele_id)
            new_ele.set_slip(self.element_slips[ele_id])
            self.elements.append(new_ele)
        self.slip_map = quakelib.SlipMap()
        self.slip_map.add_elements(self.elements)
        # Grab base Lat/Lon from fault model, used for lat/lon <-> xyz conversion
        base = geometry.model.get_base()
        self.base_lat = self.min_lat = base[0]
        self.base_lon = self.min_lon = base[1]
        self.min_lat, self.max_lat, self.min_lon, self.max_lon = geometry.model.get_latlon_bounds()
        # Expand lat/lon range in the case of plotting a few elements
        if self.small_model:
            self.min_lat = self.min_lat - MIN_LAT_DIFF*10
            self.max_lat = self.max_lat + MIN_LAT_DIFF*10
            self.min_lon = self.min_lon - MIN_LON_DIFF*10
            self.max_lon = self.max_lon + MIN_LON_DIFF*10
        # Adjust bounds for good framing on plot
        lon_range = self.max_lon - self.min_lon
        lat_range = self.max_lat - self.min_lat
        max_range = max((lon_range, lat_range))
        self.min_lon = self.min_lon - lon_range*padding
        self.min_lat = self.min_lat - lat_range*padding
        self.max_lon = self.max_lon + lon_range*padding
        self.max_lat = self.max_lat + lat_range*padding
        self.lat0, self.lon0  = (self.max_lat+self.min_lat)/2.0, (self.max_lon+self.min_lon)/2.0
        self.llcrnrlat = self.min_lat
        self.llcrnrlon = self.min_lon
        self.urcrnrlat = self.max_lat
        self.urcrnrlon = self.max_lon
        # We need a map instance to calculate the aspect ratio
        map = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        # Using the aspect ratio (h/w) to find the actual map width and height in pixels
        if map.aspect > max_map_height/max_map_width:
            map_height = max_map_height
            map_width = max_map_height/map.aspect
        else:
            map_width = max_map_width
            map_height = max_map_width*map.aspect
        # A conversion instance for doing the lat-lon to x-y conversions
        base_lld = quakelib.LatLonDepth(self.base_lat, self.base_lon, 0.0)
        self.convert = quakelib.Conversion(base_lld)
        self.lons_1d = np.linspace(self.min_lon, self.max_lon, num=int(map_width))
        self.lats_1d = np.linspace(self.min_lat, self.max_lat, num=int(map_height))
        _lons_1d = quakelib.FloatList()
        _lats_1d = quakelib.FloatList()
        for lon in self.lons_1d:
            _lons_1d.append(lon)
        for lat in self.lats_1d:
            _lats_1d.append(lat)
        # Set up the points for field evaluation, convert to xyz basis
        self.grid_1d = self.convert.convertArray2xyz(_lats_1d,_lons_1d)
        self.fault_traces_latlon = geometry.get_fault_traces()
        self._plot_str = ""
        #-----------------------------------------------------------------------
        # Gravity map configuration  #TODO: Put in switches for field_type
        #-----------------------------------------------------------------------
        self.dmc = {
            'font':               mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='normal'),
            'font_bold':          mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='bold'),
        #water
            'water_color':          '#4eacf4',
            'water_color_f':        '#4eacf4',
        #map boundaries
            'boundary_color':       '#000000',
            'boundary_width':       1.0,
            'coastline_color':      '#000000',
            'coastline_width':      1.0,
            'country_color':        '#000000',
            'country_width':        1.0,
            'state_color':          '#000000',
            'state_width':          1.0,
        #rivers
            'river_width':          0.25,
        #faults
            'fault_color':          '#000000',
            'event_fault_color':    '#ff0000',
            'fault_width':          0.5,
        #lat lon grid
            'grid_color':           '#000000',
            'grid_width':           0.0,
            'num_grid_lines':       5,
        #map props
            'map_resolution':       map_res,
            'map_projection':       map_proj,
            'plot_resolution':      72.0,
            'map_tick_color':       '#000000',
            'map_frame_color':      '#000000',
            'map_frame_width':      1,
        #map_fontsize = 12
            'map_fontsize':         26.0,   # 12   THIS IS BROKEN
        #cb_fontsize = 12
            'cb_fontcolor':         '#000000',
            'cb_height':            20.0,
            'cb_margin_t':          2.0, # 10
        }
        # Set field-specific plotting arguments
        if cbar_max is None:
            if self.field_type == 'gravity' or self.field_type == 'dilat_gravity':
                cbar_max = 20
            elif self.field_type == 'potential':
                cbar_max = 0.002
            elif self.field_type == 'geoid':
                cbar_max = 0.00015
                
        if self.field_type == 'gravity' or self.field_type == 'dilat_gravity' or self.field_type == 'potential' or self.field_type == 'geoid':
            self.dmc['cmap'] = plt.get_cmap('seismic')
            self.dmc['cbar_min'] = -cbar_max
            self.dmc['cbar_max'] = cbar_max
            if self.field_type == 'gravity' or self.field_type == 'dilat_gravity':
                self.dmc['cb_fontsize'] = 20.0
            elif self.field_type == 'potential':
                self.dmc['cb_fontsize'] = 16.0
            elif self.field_type == 'geoid':
                self.dmc['cb_fontsize'] = 16.0
                
        if self.field_type == 'displacement' or self.field_type == 'insar':
            self.dmc['boundary_color_f'] = '#ffffff'
            self.dmc['coastline_color_f'] = '#ffffff'
            if self.levels is None:
                self.dmc['cmap'] = plt.get_cmap('YlOrRd')
            else:
                self.dmc['cmap'] = plt.get_cmap('seismic')
            self.dmc['cmap_f'] = plt.get_cmap('jet')
            self.dmc['country_color_f'] = '#ffffff'
            self.dmc['state_color_f'] = '#ffffff'
            self.dmc['fault_color_f'] = '#ffffff'
            self.dmc['event_fault_color_f'] = '#ffffff'
            self.dmc['grid_color_f'] = '#ffffff'
            self.dmc['map_tick_color_f'] = '#000000'
            self.dmc['map_frame_color_f'] = '#000000'
            self.dmc['arrow_inset'] = 18.0
            self.dmc['arrow_fontsize'] = 14.0
            self.dmc['cb_fontcolor_f'] = '#000000'
            self.dmc['cb_margin_t'] = 4.0
            self.dmc['cb_fontsize'] = 20.0
            if self.levels:
                self.dmc['cbar_min'] = -cbar_max
                self.dmc['cbar_max'] = cbar_max
        #-----------------------------------------------------------------------
        # m1, fig1 is the oceans and the continents. This will lie behind the
        # masked data image.
        #-----------------------------------------------------------------------
        self.m1 = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, 
            lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        #-----------------------------------------------------------------------
        # m2, fig2 is the plotted deformation data.
        #-----------------------------------------------------------------------
        self.m2 = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, 
            lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        #-----------------------------------------------------------------------
        # m3, fig3 is the ocean land mask.
        #-----------------------------------------------------------------------
        self.m3 = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, 
            lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        
    def compute_field(self, cutoff=None):
        self.lame_lambda = 3.2e10
        self.lame_mu     = 3.0e10
        #-----------------------------------------------------------------------
        # If the cutoff is none (ie not explicitly set) calculate the cutoff for
        # this event. 
        #-----------------------------------------------------------------------
        num_involved_elements = float(len(self.element_slips.keys()))
        if cutoff is None:
            if  num_involved_elements >= self.cutoff_min_size:
                cutoff = linear_interp(
                    num_involved_elements,
                    self.cutoff_min_size,
                    self.cutoff_p2_size,
                    self.cutoff_min,
                    self.cutoff_p2
                    )
            else:
                cutoff = self.cutoff_min
        sys.stdout.write('{:0.2f} cutoff [units of element length] : '.format(cutoff))
        self.fringes = False
        if self.field_type == "gravity":
            sys.stdout.write(" Computing gravity field :")
            self.field_1d = self.slip_map.gravity_changes(self.grid_1d, self.lame_lambda, self.lame_mu, cutoff)
            # Reshape field
            self.field = np.array(self.field_1d).reshape((self.lats_1d.size,self.lons_1d.size))
        if self.field_type == "dilat_gravity":
            sys.stdout.write(" Computing dilatational gravity field :")
            self.field_1d = self.slip_map.dilat_gravity_changes(self.grid_1d, self.lame_lambda, self.lame_mu, cutoff)
            self.field = np.array(self.field_1d).reshape((self.lats_1d.size,self.lons_1d.size))
        if self.field_type == "potential":
            sys.stdout.write(" Computing gravitational potential field :")
            self.field_1d = self.slip_map.potential_changes(self.grid_1d, self.lame_lambda, self.lame_mu, cutoff)
            self.field = np.array(self.field_1d).reshape((self.lats_1d.size,self.lons_1d.size))
        elif self.field_type == "geoid":
            sys.stdout.write(" Computing geoid height change field :")
            self.field_1d = self.slip_map.potential_changes(self.grid_1d, self.lame_lambda, self.lame_mu, cutoff)
            self.field = np.array(self.field_1d).reshape((self.lats_1d.size,self.lons_1d.size))
            # To convert from potential to geoid height, divide by mean surface gravity
            self.field /= -1*self.g0
            sys.stdout.write(" g0 {} :".format(self.g0))
        elif self.field_type == "displacement" or self.field_type == "insar":
            if self.field_type == "displacement": 
                sys.stdout.write(" Computing displacement field :")
            else:
                sys.stdout.write(" Computing InSAR field :")
                self.fringes = True
            self.field_1d = self.slip_map.displacements(self.grid_1d, self.lame_lambda, self.lame_mu, cutoff)
            disp = np.array(self.field_1d).reshape((self.lats_1d.size,self.lons_1d.size,3))
            # Parse returned VectorList into separate dX,dY,dZ 2D arrays
            self.dX = np.empty((self.lats_1d.size, self.lons_1d.size))
            self.dY = np.empty((self.lats_1d.size, self.lons_1d.size))
            self.dZ = np.empty((self.lats_1d.size, self.lons_1d.size))
            it = np.nditer(self.dX, flags=['multi_index'])
            while not it.finished:
                self.dX[it.multi_index] = disp[it.multi_index][0]
                self.dY[it.multi_index] = disp[it.multi_index][1]
                self.dZ[it.multi_index] = disp[it.multi_index][2]
                it.iternext()
                
        sys.stdout.flush()
        
    def plot_str(self):
        return self._plot_str
        
    def create_field_image(self, angles=None):
        #-----------------------------------------------------------------------
        # Set all of the plotting properties
        #-----------------------------------------------------------------------
        if self.field_type == 'displacement' or self.field_type == 'insar':
            if self.field_type == 'insar':
                cmap            = self.dmc['cmap_f']
                water_color     = self.dmc['water_color_f']
                boundary_color  = self.dmc['boundary_color_f']
            else:
                cmap            = self.dmc['cmap']
                water_color     = self.dmc['water_color']
                boundary_color  = self.dmc['boundary_color']
            land_color      = cmap(0)
            if angles is not None:
                self.look_azimuth = angles[0]
                self.look_elevation = angles[1]
            else:
                if self.field_type == 'insar':
                    # Typical angles for InSAR are approx 30 deg and 40 deg respectively
                    self.look_azimuth = 30.0*np.pi/180.0
                    self.look_elevation = 40.0*np.pi/180.0
                else:
                    self.look_azimuth = 0.0
                    self.look_elevation = 0.0
            sys.stdout.write("Displacements projected along azimuth={:.1f}deg and elevation={:.1f}deg : ".format(self.look_azimuth*180.0/np.pi, self.look_elevation*180.0/np.pi))
            
        if self.field_type == 'gravity' or self.field_type == 'dilat_gravity' or self.field_type == 'potential' or self.field_type == 'geoid':
            cmap            = self.dmc['cmap']
            water_color     = self.dmc['water_color']
            boundary_color  = self.dmc['boundary_color']
            land_color      = cmap(0.5)
        sys.stdout.flush()
        plot_resolution = self.dmc['plot_resolution']
        #-----------------------------------------------------------------------
        # Set the map dimensions
        #-----------------------------------------------------------------------
        mw = self.lons_1d.size
        mh = self.lats_1d.size
        mwi = mw/plot_resolution
        mhi = mh/plot_resolution
        #-----------------------------------------------------------------------
        # Fig1 is the background land and ocean.
        #-----------------------------------------------------------------------
        fig1 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        self.m1.ax = fig1.add_axes((0,0,1,1))
        self.m1.drawmapboundary(
            color=boundary_color,
            linewidth=0,
            fill_color=water_color
        )
        self.m1.fillcontinents(
            color=land_color,
            lake_color=water_color
        )
        #-----------------------------------------------------------------------
        # Fig2 is the deformations.
        #-----------------------------------------------------------------------
        fig2 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        self.m2.ax = fig2.add_axes((0,0,1,1))
        
        if self.field_type == 'displacement' or self.field_type == 'insar':
            # Use observing angles to compute projection (field_proj) along the observing direction
            self.field_proj = -self.dX * math.sin(self.look_azimuth) * math.cos(self.look_elevation) - self.dY * math.cos(self.look_azimuth) * math.cos(self.look_elevation) + self.dZ * math.sin(self.look_elevation)
            
            # Make sure field values are at correct map location
            self.field_transformed = self.m2.transform_scalar(self.field_proj, self.lons_1d, self.lats_1d, self.lons_1d.size, self.lats_1d.size)
                        
            if self.fringes:
                # prepare the colors for the InSAR plot and do the plot
                self.insar = np.empty((self.field_transformed.shape[0],self.field_transformed.shape[1],4))
                r,g,b,a = cmap(0)
                self.insar[:,:,0].fill(r)
                self.insar[:,:,1].fill(g)
                self.insar[:,:,2].fill(b)
                self.insar[:,:,3].fill(a)
                non_zeros = self.field_transformed.nonzero()
                for n,i in enumerate(non_zeros[0]):
                    j = non_zeros[1][n]
                    r,g,b,a = cmap(math.modf(abs(self.field_transformed[i,j])/self.wavelength)[0])
                    self.insar[i, j, 0] = r
                    self.insar[i, j, 1] = g
                    self.insar[i, j, 2] = b
                    self.insar[i, j, 3] = a
                if self.norm is None:
                    self.norm = mcolor.Normalize(vmin=0, vmax=self.wavelength)
                self.m2.imshow(self.insar, interpolation='spline36')
            else:
                # Prepare the displacement plot
                if self.levels is None:
                    self.insar = np.empty(self.field_transformed.shape)
                    non_zeros = self.field_transformed.nonzero()
                    self.insar.fill(5e-4)
                    self.insar[non_zeros] = np.fabs(self.field_transformed[non_zeros])
                    vmax = np.amax(self.insar)
                    if vmax <= 1:
                        mod_vmax = 1
                    elif vmax > 1 and vmax <= 10:
                        mod_vmax = 10
                    elif vmax > 10 and vmax <= 100:
                        mod_vmax = 100
                    elif vmax > 100 and vmax <= 1000:
                        mod_vmax = 1000
                    elif vmax > 1000:
                        mod_vmax = 1000
                    if self.norm is None:
                        self.norm = mcolor.LogNorm(vmin=5e-4, vmax=mod_vmax, clip=True)
                    self.m2.imshow(self.insar, cmap=cmap, norm=self.norm)
                else:
                    map_x, map_y = self.m2(self.lons_1d, self.lats_1d)
                    XX,YY = np.meshgrid(map_x, map_y)
                    self.norm = mcolor.Normalize(vmin=self.dmc['cbar_min'], vmax=self.dmc['cbar_max'])
                    self.m2.contourf(XX, YY, self.field_transformed, self.levels, cmap=cmap, norm=self.norm, extend='both')
        
        else:
            # make sure the values are located at the correct location on the map
            self.field_transformed = self.m2.transform_scalar(self.field, self.lons_1d, self.lats_1d, self.lons_1d.size, self.lats_1d.size)
            
            if self.norm is None:
                self.norm = mcolor.Normalize(vmin=self.dmc['cbar_min'], vmax=self.dmc['cbar_max'])
        
            # Changed units to microgals (multiply MKS unit by 10^8)
            if self.field_type == 'gravity': self.field_transformed *= float(pow(10,8))
            if self.field_type == 'dilat_gravity': self.field_transformed *= float(pow(10,8))
            if self.field_type == 'geoid': self.field_transformed *= float(pow(10,2))
            
            # Plot the field on the map
            if self.levels is None:
                self.m2.imshow(self.field_transformed, cmap=cmap, norm=self.norm)
            else:
                map_x, map_y = self.m2(self.lons_1d, self.lats_1d)
                XX,YY = np.meshgrid(map_x, map_y)
                self.m2.contourf(XX, YY, self.field_transformed, self.levels, cmap=cmap, norm=self.norm, extend='both')
        #-----------------------------------------------------------------------
        # Composite fig 1 - 2 together
        #-----------------------------------------------------------------------
        # FIGURE 1 draw the renderer
        fig1.canvas.draw()
        
        # FIGURE 1 Get the RGBA buffer from the figure
        w,h = fig1.canvas.get_width_height()
        buf = np.fromstring ( fig1.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 1 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im1 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        
        # FIGURE 2 draw the renderer
        fig2.canvas.draw()
        
        # FIGURE 2 Get the RGBA buffer from the figure
        w,h = fig2.canvas.get_width_height()
        buf = np.fromstring ( fig2.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 2 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im2 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        # Clear all three figures
        fig1.clf()
        fig2.clf()
        plt.close('all')
        gc.collect()
        return im2

    def plot_field(self, output_file=None, angles=None):
        map_image = self.create_field_image(angles=angles)
        
        sys.stdout.write('map overlay : ')
        sys.stdout.flush()
        #---------------------------------------------------------------------------
        # Plot all of the geographic info on top of the displacement map image.
        #---------------------------------------------------------------------------
        # Grab all of the plot properties that we will need.
        # properties that are fringes dependent
        if self.field_type == 'insar':
            cmap            = self.dmc['cmap_f']
            coastline_color = self.dmc['coastline_color_f']
            country_color   = self.dmc['country_color_f']
            state_color     = self.dmc['state_color_f']
            fault_color     = self.dmc['fault_color_f']
            map_tick_color  = self.dmc['map_tick_color_f']
            map_frame_color = self.dmc['map_frame_color_f']
            grid_color      = self.dmc['grid_color_f']
            cb_fontcolor    = self.dmc['cb_fontcolor_f']
            arrow_inset     = self.dmc['arrow_inset']
            arrow_fontsize  = self.dmc['arrow_fontsize']
        else:
            cmap            = self.dmc['cmap']
            coastline_color = self.dmc['coastline_color']
            country_color   = self.dmc['country_color']
            state_color     = self.dmc['state_color']
            fault_color     = self.dmc['fault_color']
            map_tick_color  = self.dmc['map_tick_color']
            map_frame_color = self.dmc['map_frame_color']
            grid_color      = self.dmc['grid_color']
            cb_fontcolor    = self.dmc['cb_fontcolor']
        if self.field_type == 'displacement':             
            arrow_inset     = self.dmc['arrow_inset']
            arrow_fontsize  = self.dmc['arrow_fontsize']
        
        boundary_width  = self.dmc['boundary_width']
        coastline_width = self.dmc['coastline_width']
        country_width   = self.dmc['country_width']
        state_width     = self.dmc['state_width']
        river_width     = self.dmc['river_width']
        fault_width     = self.dmc['fault_width']
        map_frame_width = self.dmc['map_frame_width']
        map_fontsize    = self.dmc['map_fontsize']
        cb_fontsize     = self.dmc['cb_fontsize']
        cb_height       = self.dmc['cb_height']
        cb_margin_t     = self.dmc['cb_margin_t']
        grid_width      = self.dmc['grid_width']
        num_grid_lines  = self.dmc['num_grid_lines']
        font            = self.dmc['font']
        font_bold       = self.dmc['font_bold']
        map_resolution  = self.dmc['map_resolution']
        map_projection  = self.dmc['map_projection']
        plot_resolution = self.dmc['plot_resolution']

        # The sizing for the image is tricky. The aspect ratio of the plot is fixed,
        # so we cant set all of margins to whatever we want. We will set the anchor
        # to the top, left margin position. Then scale the image based on the
        # bottom/right margin, whichever is bigger.

        mw = self.lons_1d.size
        mh = self.lats_1d.size

        if mh > mw:
            ph = 768.0
            pw = mw + 70.0 + 40.0
        else:
            pw = 790.0
            ph = mh + 70.0 + 40.0

        width_frac = mw/pw
        height_frac = mh/ph
        left_frac = 70.0/pw
        bottom_frac = 70.0/ph

        pwi = pw/plot_resolution
        phi = ph/plot_resolution

        fig_res = plot_resolution

        fig4 = plt.figure(figsize=(pwi, phi), dpi=fig_res)

        #---------------------------------------------------------------------------
        # m4, fig4 is all of the boundary data.
        #---------------------------------------------------------------------------
        m4 = Basemap(
            llcrnrlon=self.min_lon,
            llcrnrlat=self.min_lat,
            urcrnrlon=self.max_lon,
            urcrnrlat=self.max_lat,
            lat_0=(self.max_lat+self.min_lat)/2.0,
            lon_0=(self.max_lon+self.min_lon)/2.0,
            resolution=map_resolution,
            projection=map_projection,
            suppress_ticks=True
        )
        m4.ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))

        # draw coastlines, edge of map.
        m4.drawcoastlines(color=coastline_color, linewidth=coastline_width)

        # draw countries
        m4.drawcountries(linewidth=country_width, color=country_color)

        # draw states
        m4.drawstates(linewidth=state_width, color=state_color)

        # draw parallels.
        parallels = np.linspace(self.lats_1d.min(), self.lats_1d.max(), num_grid_lines+1)
        m4_parallels = m4.drawparallels(parallels, fontsize=map_fontsize, labels=[1,0,0,0], color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])

        # draw meridians
        meridians = np.linspace(self.lons_1d.min(), self.lons_1d.max(), num_grid_lines+1)
        m4_meridians = m4.drawmeridians(meridians, fontsize=map_fontsize, labels=[0,0,1,0], color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])

        if self.field_type == 'displacement' or self.field_type == 'insar':
            box_size = 70.0
            # draw the azimuth look arrow
            az_width_frac    = box_size/pw
            az_height_frac   = box_size/ph
            az_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
            az_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac)/ph
            az_ax = fig4.add_axes((az_left_frac,az_bottom_frac,az_width_frac,az_height_frac))

            az_ax.set_xlim((0,1.0))
            az_ax.set_ylim((0,1.0))
            for item in az_ax.yaxis.get_ticklabels() + az_ax.xaxis.get_ticklabels() + az_ax.yaxis.get_ticklines() + az_ax.xaxis.get_ticklines():
                item.set_alpha(0)

            az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(self.look_azimuth)
            az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(self.look_azimuth)
            az_arrow_dx      = 0.8*math.sin(self.look_azimuth)
            az_arrow_dy      = 0.8*math.cos(self.look_azimuth)

            az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
            az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
            az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - self.convert.rad2deg(self.look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
            az_ax.text(1.0, 1.0, 'az = {:0.1f}{}'.format(self.convert.rad2deg(self.look_azimuth),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')

            # draw the altitude look arrow
            al_width_frac    = box_size/pw
            al_height_frac   = box_size/ph
            al_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
            al_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac - ph*al_height_frac)/ph
            al_ax = fig4.add_axes((al_left_frac,al_bottom_frac,al_width_frac,al_height_frac))

            al_ax.set_xlim((0,1.0))
            al_ax.set_ylim((0,1.0))
            for item in al_ax.yaxis.get_ticklabels() + al_ax.xaxis.get_ticklabels() + al_ax.yaxis.get_ticklines() + al_ax.xaxis.get_ticklines():
                item.set_alpha(0)

            al_arrow_start_x = 0.1 + 0.8*math.cos(self.look_elevation)
            al_arrow_start_y = 0.1 + 0.8*math.sin(self.look_elevation)
            al_arrow_dx      = -0.8*math.cos(self.look_elevation)
            al_arrow_dy      = -0.8*math.sin(self.look_elevation)

            al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
            al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
            al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=self.convert.rad2deg(self.look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
            al_ax.text(1.0, 1.0, 'al = {:0.1f}{}'.format(self.convert.rad2deg(self.look_elevation),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')
            
            # draw the box with the magnitude
            mag_width_frac    = box_size/pw
            if self.fringes:
                mag_height_frac   = 25.0/ph    # originally 10.0/ph
            else:
                mag_height_frac   = 15.0/ph 
            mag_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
            mag_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac  - ph*az_height_frac - ph*mag_height_frac)/ph
            mag_ax = fig4.add_axes((mag_left_frac,mag_bottom_frac,mag_width_frac,mag_height_frac))

            mag_ax.set_xlim((0,1.0))
            mag_ax.set_ylim((0,1.0))
            for item in mag_ax.yaxis.get_ticklabels() + mag_ax.xaxis.get_ticklabels() + mag_ax.yaxis.get_ticklines() + mag_ax.xaxis.get_ticklines():
                item.set_alpha(0)
            
            if self.event_id is not None:
                mag_ax.text(0.5, 0.5, 'm = {:0.3f}'.format(float(events._events[self.event_id].getMagnitude())), fontproperties=font_bold, size=arrow_fontsize, ha='center', va='center')
            else:
                avg_slip = np.average([x[1] for x in self.element_slips.items()])
                mag_ax.text(0.5, 0.5, 'mean slip \n{:0.3f}m'.format(avg_slip), fontproperties=font_bold, size=arrow_fontsize-1, ha='center', va='center')

        # add the map image to the plot
        m4.imshow(map_image, origin='upper')
        
        # If plotting event field, get involved sections
        if self.event_id is not None:
            involved_sections = events.get_event_sections(self.event_id, geometry)
            sys.stdout.write(" Event slips on {} sections out of {} : ".format(len(involved_sections), len(geometry.model.getSectionIDs()) ))
        else:
            involved_sections = geometry.model.getSectionIDs()

        # print faults on lon-lat plot
        for sid, sec_trace in self.fault_traces_latlon.iteritems():
            sec_trace_lons = [lat_lon[1] for lat_lon in sec_trace]
            sec_trace_lats = [lat_lon[0] for lat_lon in sec_trace]
            
            trace_Xs, trace_Ys = m4(sec_trace_lons, sec_trace_lats)
            
            if sid in involved_sections:
                linewidth = fault_width + 2.5
            else:
                linewidth = fault_width

            m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')

        #plot the cb
        left_frac = 70.0/pw
        bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
        width_frac = mw/pw
        height_frac = cb_height/ph

        cb_ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
        norm = self.norm
        cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
               norm=norm,
               orientation='horizontal')
        if self.field_type == 'displacement' or self.field_type == 'insar':
            if self.fringes:
                cb_title = 'Displacement [m]'
            else:
                cb_title = 'Total displacement [m]'
                if self.levels:
                    # Make first and last ticks on colorbar be <MIN and >MAX.
                    # Values of colorbar min/max are set in FieldPlotter init.
                    cb_tick_labs    = [item.get_text() for item in cb_ax.get_xticklabels()]
                    cb_tick_labs[0] = '<'+cb_tick_labs[0]
                    cb_tick_labs[-1]= '>'+cb_tick_labs[-1]
                    cb_ax.set_xticklabels(cb_tick_labs)
    
        else:
            if self.field_type == 'gravity' or self.field_type == 'dilat_gravity':
                cb_title = r'Gravity changes [$\mu gal$]'
            elif self.field_type == 'potential':
                cb_title = r'Gravitational potential changes [$m^2$/$s^2$]'
            elif self.field_type == 'geoid':
                cb_title = 'Geoid height change [cm]'
            # Make first and last ticks on colorbar be <MIN and >MAX.
            # Values of colorbar min/max are set in FieldPlotter init.
            cb_tick_labs    = [item.get_text() for item in cb_ax.get_xticklabels()]
            cb_tick_labs[0] = '<'+cb_tick_labs[0]
            cb_tick_labs[-1]= '>'+cb_tick_labs[-1]
            cb_ax.set_xticklabels(cb_tick_labs)

        cb_ax.set_title(cb_title, fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )

        for label in cb_ax.xaxis.get_ticklabels():
            label.set_fontproperties(font)
            label.set_fontsize(cb_fontsize)
            label.set_color(cb_fontcolor)
        for line in cb_ax.xaxis.get_ticklines():
            line.set_alpha(0)

        fig4.savefig(output_file, format='png', dpi=fig_res)
        sys.stdout.write('\nPlot saved: {}'.format(output_file))
        sys.stdout.write('\ndone\n')
        sys.stdout.flush()
        
        
        


# Evaluate an event field at specified lat/lon coords
# Currently only for displacement field
class FieldEvaluator:
    def __init__(self, geometry, event_id, event, element_slips, LLD_file):
        # LLD file contains columns of lat/lon/depth for the points we wish to evaluate
        self.LLDdata = np.genfromtxt(LLD_file, dtype=[('lat','f8'),('lon','f8'), ('z','f8')],skip_header=4)
        # Set field and event data
        self.event_id = event_id
        self.LLD_file = LLD_file
        self.elements = quakelib.SlippedElementList()
        self.element_ids = element_slips.keys()    
        self.slip_map = quakelib.SlipMap()
        # Assign the slips from element_slips
        for ele_id in self.element_ids:
            new_ele = geometry.model.create_slipped_element(ele_id)
            new_ele.set_slip(element_slips[ele_id])
            self.elements.append(new_ele)
        # Add the elements to the slip map
        self.slip_map.add_elements(self.elements)
        # A conversion instance for doing the lat-lon to x-y conversions
        base = geometry.model.get_base()
        base_lld = quakelib.LatLonDepth(base[0], base[1], 0.0)
        self.convert = quakelib.Conversion(base_lld)
        _lons_1d = quakelib.FloatList()
        _lats_1d = quakelib.FloatList()
        self.lons_1d = self.LLDdata['lon']
        self.lats_1d = self.LLDdata['lat']
        # Set up the points for field evaluation, convert to xyz basis
        self.grid_1d = quakelib.VectorList()
        for i in range(len(self.lons_1d)):
            self.grid_1d.append(self.convert.convert2xyz(quakelib.LatLonDepth(self.lats_1d[i],self.lons_1d[i])))
    #            
    def compute_field(self):        
        self.lame_lambda = 3.2e10
        self.lame_mu     = 3.0e10
        self.field_1d    = self.slip_map.displacements(self.grid_1d, self.lame_lambda, self.lame_mu, 1e9)
        outname = self.LLD_file.split(".tx")[0]+"_dispField_event"+str(self.event_id)+".txt"
        outfile = open(outname,'w')
        # Write the header with the number of points
        outfile.write("#### number of points ####\n")
        outfile.write("{}\n".format(len(self.field_1d)))
        outfile.write("##########################\n")
        for i in range(len(self.field_1d)):
            outfile.write("{}\t{}\t{}\n".format(self.lats_1d[i], self.lons_1d[i], self.field_1d[i][2]))
        outfile.close()
        sys.stdout.write("\n---> Event displacements written to "+outname)
        sys.stdout.write("\n")



class BasePlotter:
    def create_plot(self, plot_type, log_y, x_data, y_data, plot_title, x_label, y_label, filename):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        if plot_type == "scatter":
            ax.scatter(x_data, y_data, color='g')
        elif plot_type == "line":
            ax.plot(x_data, y_data, color='g')
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

    def multi_line_plot(self, x_data, y_data, labels, linewidths, plot_title, x_label, y_label, legend_str, filename, colors=None, linestyles=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if linestyles is None: linestyles = ["-" for each in x_data]
        fig.suptitle(plot_title, fontsize=10)
        if colors is not None:
            if not (len(x_data) == len(y_data) and len(x_data) == len(colors) and len(colors) == len(labels) and len(linewidths) == len(colors)):
                raise "These lists must be the same length: x_data, y_data, colors, labels, linewidths."
            for i in range(len(x_data)):
                ax.plot(x_data[i], y_data[i], color=colors[i], label=labels[i], linewidth=linewidths[i], ls=linestyles[i])
        else:
            if not (len(x_data) == len(y_data) and len(x_data) == len(labels) and len(linewidths) == len(y_data)):
                raise "These lists must be the same length: x_data, y_data, labels, linewidths."
            for i in range(len(x_data)):
                ax.plot(x_data[i], y_data[i], label=labels[i], linewidth=linewidths[i], ls=linestyles[i])
        ax.legend(title=legend_str, loc='best')
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

    def t0_vs_dt_plot(self, t0_dt_plot, wait_75, filename):
# TODO: Set fonts explicitly
        t0_dt_main_line_color   = '#000000'
        t0_dt_sub_line_color    = '#737373'
        t0_dt_main_line_width   = 2.0
        t0_dt_sub_line_width    = 1.0
        t0_dt_range_color       = plt.get_cmap('autumn')(0.99)
        years_since_line_color  = 'blue'
        legend_loc              = 'best'
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(r't$_0$ [years]')
        ax.set_ylabel(r'$\Delta$t [years]')
        percents = t0_dt_plot.keys()
        ax.fill_between(t0_dt_plot[min(percents)]['x'], t0_dt_plot[min(percents)]['y'], y2=t0_dt_plot[max(percents)]['y'], linewidth=0, facecolor=t0_dt_range_color)
        for percent in t0_dt_plot.iterkeys():
            if percent == min(percents):
                linewidth = t0_dt_sub_line_width
                color = t0_dt_sub_line_color
                linestyle = '--'
            elif percent == max(percents):
                linewidth = t0_dt_sub_line_width
                color = t0_dt_sub_line_color
                linestyle = ':'
            else:
                linewidth = t0_dt_main_line_width
                color = t0_dt_main_line_color
                linestyle = '-'
            ax.plot(t0_dt_plot[percent]['x'], t0_dt_plot[percent]['y'], color=color, linewidth=linewidth, linestyle=linestyle, label='{}%'.format(percent))
        if wait_75 is not None:
            # Draw vertical dotted line where "today" is denoted by years_since
            ax.axvline(x=years_since,ymin=0,ymax=wait_75,color=years_since_line_color,linewidth=t0_dt_main_line_width,linestyle='--')
        ax.legend(title='event prob.', loc=legend_loc, handlelength=5)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

    def scatter_and_errorbar(self, log_y, x_data, y_data, err_x, err_y, y_error, err_label, plot_title, x_label, y_label, filename, add_x = None, add_y = None, add_label = None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        ax.scatter(x_data, y_data)
        ax.errorbar(err_x, err_y, yerr = y_error, label=err_label, ecolor='r')
        if add_x is not None:
            if log_y: ax.semilogy(add_x, add_y, label = add_label, c = 'k')
            if not log_y: ax.plot(add_x, add_y, label = add_label, c = 'k')
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        ax.legend(loc = "best")
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

    def scatter_and_line(self, log_y, x_data, y_data, line_x, line_y, line_label, plot_title, x_label, y_label, filename):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        ax.scatter(x_data, y_data, color='g')
        if line_x is not None and line_y is not None:
            ax.plot(line_x, line_y, label = line_label, ls='-', c = 'k', lw=2)
            ax.legend(loc = "best")
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        
        if args.zoom: plt.ylim(-5,5)
        
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))
        
    def scatter_and_multiline(self, log_y, x_data, y_data, lines_x, lines_y, line_labels, line_widths, line_styles, colors, plot_title, x_label, y_label, filename):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y: ax.set_yscale('log')
        ax.scatter(x_data, y_data, color='g')
        for i in range(len(lines_x)):
            ax.plot(lines_x[i], lines_y[i], label = line_labels[i], ls=line_styles[i], lw=line_widths[i], c = colors[i])
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        ax.legend(loc = "lower right")
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

class MagnitudeRuptureAreaPlot(BasePlotter):
    def plot(self, events, filename, WC94=False):
        ra_list = events.event_rupture_areas()
        mag_list = events.event_magnitudes()
        ra_renorm_list = [quakelib.Conversion().sqm2sqkm(ra) for ra in ra_list]
        min_mag, max_mag = min(mag_list), max(mag_list)
        if WC94:
            scale_x, scale_y = Distributions().wells_coppersmith('area')
            scale_label = "Wells & Coppersmith 1994"
            full_x, full_y = Distributions().wells_coppersmith('area', min_mag=min_mag, max_mag=max_mag)
            lines_x = [scale_x, full_x]
            lines_y = [scale_y, full_y]
            line_labels = [scale_label, None]
            line_widths = [2.0, 1.0]
            line_styles = ['-', '--']
            colors = ['k', 'k']
            self.scatter_and_multiline(True, mag_list, ra_renorm_list, lines_x, lines_y, line_labels, line_widths, line_styles, colors, "",   "Magnitude", "Rupture Area (square km)", filename)
        else:
            self.create_plot("scatter", True, mag_list, ra_renorm_list, events.plot_str(), "Magnitude", "Rupture Area (square km)", filename)

class MagnitudeMeanSlipPlot(BasePlotter):
    def plot(self, events, filename, WC94=False):
        slip_list = events.event_mean_slip()
        mag_list = events.event_magnitudes()
        min_mag, max_mag = min(mag_list), max(mag_list)
        if WC94:
            scale_x, scale_y = Distributions().wells_coppersmith('slip')
            scale_label = "Wells & Coppersmith 1994"
            full_x, full_y = Distributions().wells_coppersmith('slip', min_mag=min_mag, max_mag=max_mag)
            lines_x = [scale_x, full_x]
            lines_y = [scale_y, full_y]
            line_labels = [scale_label, None]
            line_widths = [2.0, 1.0]
            line_styles = ['-', '--']
            colors = ['k', 'k']
            self.scatter_and_multiline(True, mag_list, slip_list, lines_x, lines_y, line_labels, line_widths, line_styles, colors, "",   "Magnitude", "Mean Slip (meters)", filename)
        else:
            self.create_plot("scatter", True, mag_list, slip_list, events.plot_str(), "Magnitude", "Mean Slip (meters)", filename)

class FrequencyMagnitudePlot(BasePlotter):
    def plot(self, events, filename, UCERF2 = False, b1 = False):
        # California observed seismicity rates and errorbars (UCERF2)
        x_UCERF = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
        y_UCERF = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]
        y_error_UCERF = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
        add_x, add_y, add_label = None, None, None
        mag_list = events.event_magnitudes()
        cum_freq = {}
        freq_x, freq_y = [], []
        num_events = len(mag_list)
        years = events.event_years()
        if min(years) < max(years)*.01:
            # In most sims, it takes 30-100 years for first sim to occur, but the sim started at year 0. So if the first event occurs at a year that's before the first 1% of sim time, consider the filtered events to represent year=0 and onward. Needed for accurate # of events/yr
            year_range = max(years)
        else:
            year_range = max(years) - min(years)
        for num, mag in enumerate(sorted(mag_list)):
            cum_freq[mag] = num_events - (num + 1)
        for mag in sorted(cum_freq.iterkeys()):
            freq_x.append(mag)
            freq_y.append(float(cum_freq[mag])/year_range)
        if b1:
            add_x = np.linspace(min(freq_x),max(freq_x),10)
            fit_point = freq_x[(np.abs(np.array(freq_x)-MIN_FIT_MAG)).argmin()]
            add_y = 10**(math.log(fit_point,10)+freq_x[0]-add_x)
            add_label = "b==1"
        if UCERF2:
            self.scatter_and_errorbar(True, freq_x, freq_y, x_UCERF, y_UCERF, y_error_UCERF, "UCERF2", events.plot_str(), "Magnitude (M)", "# events/year with mag > M", filename, add_x=add_x, add_y=add_y, add_label=add_label)
        if b1 and not UCERF2:
            self.scatter_and_line(True, freq_x, freq_y, add_x, add_y, add_label, events.plot_str(), "Magnitude (M)", "# events/year with mag > M", filename)
        if not UCERF2 and not b1:
            self.create_plot("scatter", True, freq_x, freq_y, events.plot_str(), "Magnitude (M)", "# events/year with mag > M", filename)

class StressHistoryPlot(BasePlotter):
    def plot(self, stress_set, elements):
        stress_histories = {}
        for element in elements:
            stress_histories[element] = []
        for stress_state in stress_set:
            if stress_state.getSweepNum() == 0:
                stresses = stress_state.stresses()
                for stress in stresses:
                    if stress._element_id in elements:
                        stress_histories[element].append(stress._shear_stress)
        for element in elements:
            print(stress_histories[element])
        #self.create_plot("scatter", True, mag_vals, mag_norm, events.plot_str(), "Shear Stress", "Year")
        
class DiagnosticPlot(BasePlotter):
    def plot_shear_stress_changes(self, events, filename):
        shear_init = np.array(events.event_initial_shear_stresses())
        shear_final = np.array(events.event_final_shear_stresses())
        years = events.event_years()
        stress_changes = (shear_final-shear_init)/shear_init
        # Generate the binned averages too
        x_ave, y_ave = calculate_averages(years,stress_changes,log_bin=False,num_bins=20)
        self.scatter_and_line(False, years, stress_changes, x_ave, y_ave, "binned average", "Event shear stress changes", "simulation time [years]", "fractional change", filename)
        
    def plot_normal_stress_changes(self, events, filename):
        normal_init = np.array(events.event_initial_normal_stresses())
        normal_final = np.array(events.event_final_normal_stresses())
        years = events.event_years()
        stress_changes = (normal_final-normal_init)/normal_init
        # Generate the binned averages too
        x_ave, y_ave = calculate_averages(years,stress_changes,log_bin=False,num_bins=20)
        self.scatter_and_line(False, years, stress_changes, x_ave, y_ave, "binned average", "Event normal stress changes", "simulation time [years]", "fractional change", filename)
        
    def plot_number_of_sweeps(self, events, filename):
        num_sweeps = np.array(events.number_of_sweeps())
        years = events.event_years()
        # Generate the binned averages too
        x_ave, y_ave = calculate_averages(years,num_sweeps,log_bin=False,num_bins=20)
        self.scatter_and_line(True, years, num_sweeps, x_ave, y_ave, "binned average", " ", "simulation time [years]", "number of event sweeps", filename)
        
    def plot_mean_slip(self, events, filename):
        slips = np.array(events.event_mean_slip())
        years = events.event_years()
        # Generate the binned averages too
        x_ave, y_ave = calculate_averages(years,slips,log_bin=False,num_bins=20)
        self.scatter_and_line(True, years, slips, x_ave, y_ave, "binned average", " ", "simulation time [years]", "event mean slip [m]", filename)

class ProbabilityPlot(BasePlotter):
    def plot_p_of_t(self, events, filename):
        # Cumulative probability P(t) as a function of interevent time t
        intervals = np.array(events.interevent_times())
        prob = {}
        prob['x'] = np.sort(intervals)
        prob['y'] = np.arange(float(intervals.size))/float(intervals.size)
        self.create_plot("line", False, prob['x'], prob['y'], events.plot_str(),"t [years]", "P(t)", filename)

    def plot_conditional_fixed_dt(self, events, filename, fixed_dt=30.0):
        # P(t0 + dt, t0) vs. t0 for fixed dt
        intervals = np.array(events.interevent_times())
        prob_dt = {'x':[],'y':[]}
        t0_to_eval = np.arange(0.0,int(intervals.max())+.01,1.0)
        for t0 in t0_to_eval:
            int_t0_dt = intervals[np.where( intervals > t0+fixed_dt)]
            int_t0 = intervals[np.where( intervals > t0)]
            if int_t0.size != 0:
                prob_dt['x'].append(t0)
                prob_dt['y'].append(1.0 - float(int_t0_dt.size)/float(int_t0.size))
        self.create_plot("line", False, prob_dt['x'], prob_dt['y'], events.plot_str(),"t0 [years]", "P(t0 + dt, t0)", filename)

    def plot_p_of_t_multi(self, events, filename, beta=None, tau=None, num_t0=4, numPoints=200):
        # Cumulative conditional probability P(t,t0) as a function of
        # interevent time t, computed for multiple t0. Beta/Tau are Weibull parameters
        line_colormap = plt.get_cmap('autumn')
        intervals = np.array(events.interevent_times())
        conditional = {}
        weibull = {}
        max_t0 = int(intervals.max())
        t0_to_eval = list(np.linspace(0, max_t0, num=numPoints))
        t0_to_plot = [int(t) for t in np.linspace(0, int(max_t0/2.0), num=num_t0)]
        # To get the lines of P(t,t0) evaluated at integer values of t0
        t0_to_eval = np.sort(t0_to_eval+t0_to_plot)
        t0_to_plot = np.array(t0_to_plot)
        for t0 in t0_to_eval:
            int_t0 = intervals[np.where( intervals > t0)]
            if int_t0.size != 0:
                conditional[t0] = {'x':[],'y':[]}
                weibull[t0]     = {'x':[],'y':[]}
                for dt in range(max_t0-int(t0)):
                    int_t0_dt  = intervals[np.where( intervals > t0+dt)]
                    prob_t0_dt = 1.0 - float(int_t0_dt.size)/float(int_t0.size)
                    conditional[t0]['x'].append(t0+dt)
                    conditional[t0]['y'].append(prob_t0_dt)
                    if beta is not None and tau is not None:
                        weibull[t0]['x'].append(t0+dt)
                        weibull_t0_dt = Distributions().cond_weibull(weibull[t0]['x'][-1],t0,beta,tau)
                        weibull[t0]['y'].append(weibull_t0_dt)
            else:
                conditional[t0] = None
                weibull[t0] = None
        x_data_prob = [conditional[t0]['x'] for t0 in t0_to_plot]
        y_data_prob = [conditional[t0]['y'] for t0 in t0_to_plot]
        t0_colors   = [line_colormap(float(t0*.8)/t0_to_plot.max()) for t0 in t0_to_plot]
        prob_lw     = [2 for t0 in t0_to_plot]
        if beta is not None and tau is not None:
            x_data_weib = [weibull[t0]['x'] for t0 in t0_to_plot]
            y_data_weib = [weibull[t0]['y'] for t0 in t0_to_plot]
            weib_colors = ['k' for t0 in t0_to_plot]
            weib_labels = [None for t0 in t0_to_plot]
            weib_lw     = [1 for t0 in t0_to_plot]
            # List concatenation, not addition
            colors = t0_colors + weib_colors
            x_data = x_data_prob + x_data_weib
            y_data = y_data_prob + y_data_weib
            labels = [t0 for t0 in t0_to_plot] + weib_labels
            linewidths = prob_lw + weib_lw
        else:
            colors = t0_colors
            x_data = x_data_prob
            y_data = y_data_prob
            labels = [t0 for t0 in t0_to_plot]
            linewidths = prob_lw
        legend_string = r't$_0$='
        y_lab         = r'P(t, t$_0$)'
        x_lab         = r't = t$_0$ + $\Delta$t [years]'
        plot_title    = ""
        self.multi_line_plot(x_data, y_data, labels, linewidths, plot_title, x_lab, y_lab, legend_string, filename, colors=colors)

    def plot_dt_vs_t0(self, events, filename, years_since=None):
        # Plot the waiting times corresponding to 25/50/75% conditional probabilities
        # as a function of t0 (time since last earthquake on the selected faults).
        # years_since is the number of years since the last observed (real) earthquake
        # on the selected faults.
        intervals = np.array(events.interevent_times())
        conditional = {}
        wait_75 = None
        max_t0 = int(intervals.max())
        # t0_to_eval used to evaluate waiting times with 25/50/75% probability given t0=years_since
        t0_to_eval = np.arange(0, max_t0, 1.0)
        # t0_to_plot is "smoothed" so that the plots aren't as jagged
        t0_to_plot = np.linspace(0, int(max_t0), num=10)
        t0_to_plot = [int(t0) for t0 in t0_to_plot]
        t0_dt      = {}
        t0_dt_plot = {}
        # First generate the conditional distributions P(t,t0) for each t0
        for t0 in t0_to_eval:
            int_t0 = intervals[np.where( intervals > t0)]
            if int_t0.size != 0:
                conditional[t0] = {'x':[],'y':[]}
                for dt in range(max_t0-int(t0)):
                    int_t0_dt = intervals[np.where( intervals > t0+dt)]
                    conditional[t0]['x'].append(t0+dt)
                    prob_t0_dt    = 1.0 - float(int_t0_dt.size)/float(int_t0.size)
                    conditional[t0]['y'].append(prob_t0_dt)
        # Loop over the probabilities whose waiting times we want to plot, invert P(t,t0)
        for percent in [0.25, 0.5, 0.75]:
            t0_dt[int(percent*100)]      = {'x':[],'y':[]}
            t0_dt_plot[int(percent*100)] = {'x':[],'y':[]}
            for t0 in t0_to_eval:
                if conditional[t0] is not None:
                    # Invert the conditional probabilities, find the recurrence time closest to
                    # the current percent
                    index   = (np.abs(np.array(conditional[t0]['y'])-percent)).argmin()
                    dt      = conditional[t0]['x'][index]-t0
                    t0_dt[int(percent*100)]['x'].append(t0)
                    t0_dt[int(percent*100)]['y'].append(dt)
                    if t0 in t0_to_plot:
                        t0_dt_plot[int(percent*100)]['x'].append(t0)
                        t0_dt_plot[int(percent*100)]['y'].append(dt)
        if years_since is not None:
            # Print out the "Forecast", the 25/50/75% probability given t0=years_since
            ind_25 = (np.abs(np.array(t0_dt[25]['x'])-years_since)).argmin()
            ind_50 = (np.abs(np.array(t0_dt[50]['x'])-years_since)).argmin()
            ind_75 = (np.abs(np.array(t0_dt[75]['x'])-years_since)).argmin()
            wait_25 = t0_dt[25]['y'][ind_25]
            wait_50 = t0_dt[50]['y'][ind_50]
            wait_75 = t0_dt[75]['y'][ind_75]
            sys.stdout.write('For t0 = {:.2f} years'.format(year_eval))
            sys.stdout.write('\n25% waiting time: {:.2f} years'.format(wait_25))
            sys.stdout.write('\n50% waiting time: {:.2f} years'.format(wait_50))
            sys.stdout.write('\n75% waiting time: {:.2f} years'.format(wait_75))
            sys.stdout.write('\n=======================================\n\n')
        self.t0_vs_dt_plot(t0_dt_plot, wait_75, filename)

class Distributions:
    def weibull(self, X, beta, tau):
        # Return the Weibull distribution at a point
        return 1-np.exp( -(X/float(tau))**beta)

    def cond_weibull(self, X, t0, beta, tau):
        # Return the conditional Weibull distribution at a single point
        return 1-np.exp( (t0/float(tau))**beta - (X/float(tau))**beta)

    def wells_coppersmith(self, type, min_mag=None, max_mag=None, num=5):
        # Return empirical scaling relations from Wells & Coppersmith 1994
        log_10 = np.log(10)
        if type.lower() == 'area':
            if min_mag is None: min_mag, max_mag = 4.8, 7.9
            a = -3.49 
            b =  0.91
        elif type.lower() == 'slip':
            if min_mag is None: min_mag, max_mag = 5.6, 8.1
            a = -4.80
            b =  0.69
        else:
            raise "Must specify rupture area or mean slip"
        x_data = np.linspace(min_mag, max_mag, num=num)
        y_data = np.array([pow(10,a+b*m) for m in x_data])
        #y_err  = np.array([log_10*y_data[i]*np.sqrt(sig_a**2 + sig_b**2 * x_data[i]**2) for i in range(len(x_data))])
        return x_data, y_data

if __name__ == "__main__":
    # yoder:
    # when run as a command-line, switch pyplot to "background" mode, but permit interactive mode for... well, interactive
    # mode. needs to be double-checked on all counts.
    plt.switch_backend('agg') #Required for map plots
    #
    #
    # Specify arguments
    parser = argparse.ArgumentParser(description="PyVQ.")

    # Event/model file arguments
    parser.add_argument('--event_file', required=False,
            help="Name of event file to analyze.")
    parser.add_argument('--sweep_file', required=False,
            help="Name of sweep file to analyze.")
    parser.add_argument('--model_file', required=False,
            help="Name of model (geometry) file to use in analysis.")
    parser.add_argument('--model_file_type', required=False,
            help="Model file type, either hdf5 or text.")
    parser.add_argument('--stress_index_file', required=False,
            help="Name of stress index file to use in analysis.")
    parser.add_argument('--stress_file', required=False,
            help="Name of stress file to use in analysis.")
    parser.add_argument('--summary', type=int, required=False,
            help="Specify the number of largest magnitude EQs to summarize.")

    # Event filtering arguments
    parser.add_argument('--min_magnitude', type=float, required=False,
            help="Minimum magnitude of events to process.")
    parser.add_argument('--max_magnitude', type=float, required=False,
            help="Maximum magnitude of events to process.")
    parser.add_argument('--min_year', type=float, required=False,
            help="Minimum year of events to process.")
    parser.add_argument('--max_year', type=float, required=False,
            help="Maximum year of events to process.")
    parser.add_argument('--min_slip', type=float, required=False,
            help="Minimum mean slip of events to process.")
    parser.add_argument('--max_slip', type=float, required=False,
            help="Maximum mean slip of events to process.")
    parser.add_argument('--min_event_num', type=float, required=False,
            help="Minimum event number of events to process.")
    parser.add_argument('--max_event_num', type=float, required=False,
            help="Maximum event number of events to process.")
    parser.add_argument('--use_sections', type=int, nargs='+', required=False,
            help="List of model sections to use (all sections used if unspecified).")
    parser.add_argument('--use_trigger_sections', type=int, nargs='+', required=False,
            help="List of model triggering sections to use for subsetting events.")

    # Statisical plotting arguments
    parser.add_argument('--plot_freq_mag', required=False, type=int,
            help="Generate frequency magnitude plot. 1: Only event data, 2: Plot b=1 Gutenberg-Richter relation, 3: Plot UCERF2 observed seismicity rates, 4: Plot UCERF2 and the b=1 line.")
    parser.add_argument('--plot_mag_rupt_area', required=False, action='store_true',
            help="Generate magnitude vs rupture area plot.")
    parser.add_argument('--plot_mag_mean_slip', required=False, action='store_true',
            help="Generate magnitude vs mean slip plot.")
    parser.add_argument('--all_stat_plots', required=False, action='store_true',
            help="Generate frequency-magnitude, magnitude vs rupture area, and magnitude vs mean slip plots.")  
    parser.add_argument('--wc94', required=False, action='store_true',
            help="Plot Wells and Coppersmith 1994 scaling relations.")             

    # Probability plotting arguments
    parser.add_argument('--plot_prob_vs_t', required=False, action='store_true',
            help="Generate earthquake recurrence probability at time t plot.")
    parser.add_argument('--plot_prob_vs_t_fixed_dt', required=False, action='store_true',
            help="Generate earthquake recurrence probability at time t + dt vs t plot.")
    parser.add_argument('--plot_cond_prob_vs_t', required=False, action='store_true',
            help="Generate earthquake recurrence conditional probabilities at time t = t0 + dt for multiple t0.")
    parser.add_argument('--plot_waiting_times', required=False, action='store_true',
            help="Generate waiting times until the next earthquake as function of time since last earthquake.")
    parser.add_argument('--beta', required=False, type=float,
            help="Beta parameter for the Weibull distribution, must also specify Tau")
    parser.add_argument('--tau', required=False, type=float,
            help="Tau parameter for the Weibull distribution, must also specify Beta")
            
    # Field plotting arguments
    parser.add_argument('--field_plot', required=False, action='store_true',
            help="Plot surface field for a specified event, e.g. gravity changes or displacements.")
    parser.add_argument('--field_type', required=False, help="Field type: gravity, dilat_gravity, displacement, insar, potential, geoid")
    parser.add_argument('--colorbar_max', required=False, type=float, help="Max unit for colorbar")
    parser.add_argument('--event_id', required=False, type=int, help="Event number for plotting event fields")
    parser.add_argument('--uniform_slip', required=False, type=float, help="Amount of slip for each element in the model_file, in meters.")
    parser.add_argument('--angles', type=float, nargs='+', required=False,
            help="Observing angles (azimuth, elevation) for InSAR or displacement plots, in degrees.")
    parser.add_argument('--levels', type=float, nargs='+', required=False,
            help="Levels for contour plot.")
    parser.add_argument('--small_model', required=False, action='store_true', help="Small fault model, used to specify map extent.")    
    parser.add_argument('--traces', required=False, action='store_true', help="Plot the fault traces from a fault model on a map.") 
    parser.add_argument('--field_eval', required=False, action='store_true', help="Evaluate an event field at specified lat/lon. Must provide the file, --lld_file")
    parser.add_argument('--lld_file', required=False, help="File containing lat/lon columns to evaluate an event field.")
    
    # Greens function plotting arguments
    parser.add_argument('--greens', required=False, action='store_true', help="Plot single element Okubo Green's functions. Field type also required.")            
    parser.add_argument('--plot_name', required=False, help="Name for saving the plot to file.")
    parser.add_argument('--Nx', required=False, type=int, help="Number of points along x axis to evaluate function (default 690).")
    parser.add_argument('--Ny', required=False, type=int, help="Number of points along y axis to evaluate function. (default 422)")
    parser.add_argument('--Xmin', required=False, type=float, help="Minimum value of x in meters (along strike direction) for plotting. (default -5km)")
    parser.add_argument('--Xmax', required=False, type=float, help="Maximum value of x in meters (along strike direction) for plotting. (default 15km)")
    parser.add_argument('--Ymin', required=False, type=float, help="Minimum value of y in meters (distance from fault direction) for plotting. (default -10km)")
    parser.add_argument('--Ymax', required=False, type=float, help="Maximum value of y in meters (distance from fault direction) for plotting. (default 10km)")
    parser.add_argument('--L', required=False, type=float, help="Length of the fault in meters. (default 10km)")
    parser.add_argument('--W', required=False, type=float, help="Down-dip Width of the fault in meters. (default 10km)")
    parser.add_argument('--DTTF', required=False, type=float, help="Distance to the top of the fault in meters (i.e. distance below ground that the fault is buried). (default 1km)")
    parser.add_argument('--dip', required=False, type=float, help="Dip angle of fault in degrees.")
    parser.add_argument('--rake', required=False, type=float, help="Rake angle of fault in degrees.")
    parser.add_argument('--g', required=False, type=float, help="Mean surface gravity in meters/s^2, default is 9.81.")
    parser.add_argument('--_lambda', required=False, type=float, help="Lame's first parameter, default 3.2e10.")
    parser.add_argument('--mu', required=False, type=float, help="Shear modulus, default 3.0e10.")
    
    # Stress plotting arguments
    parser.add_argument('--stress_elements', type=int, nargs='+', required=False,
            help="List of elements to plot stress history for.")
            
    # Diagnostic plots
    parser.add_argument('--diagnostics', required=False, action='store_true',
            help="Plot all diagnostic plots")
    parser.add_argument('--num_sweeps', required=False, action='store_true',
            help="Plot the number of sweeps for events")
    parser.add_argument('--event_shear_stress', required=False, action='store_true',
            help="Plot shear stress changes for events")
    parser.add_argument('--event_normal_stress', required=False, action='store_true',
            help="Plot normal stress changes for events")
    parser.add_argument('--event_mean_slip', required=False, action='store_true',
            help="Plot the mean slip for events")
    parser.add_argument('--zoom', required=False, action='store_true',
            help="Force zoomed bounds on scatter and line plots")
            
    # Geometry
    parser.add_argument('--slip_rates', required=False, action='store_true',
            help="Print element id and slip rate for all elements.")
    parser.add_argument('--elements', type=int, nargs='+', required=False,
            help="List of elements for filtering.")
    parser.add_argument('--slip_time_series', required=False, action='store_true',
            help="Return the slip time series for all specified --elements.")
    parser.add_argument('--dt', required=False, type=float, help="Time step for slip rate plots, unit is decimal years.")
    parser.add_argument('--event_kml', required=False, action='store_true',
            help="Save a KML (Google Earth) file of the event elements, colored by event slip.")

    # Validation/testing arguments
    parser.add_argument('--validate_slip_sum', required=False,
            help="Ensure the sum of mean slip for all events is within 1 percent of the specified value.")
    parser.add_argument('--validate_mean_interevent', required=False,
            help="Ensure the mean interevent time for all events is within 2 percent of the specified value.")

    args = parser.parse_args()
    
    # ------------------------------------------------------------------------
    # Catch these errors before reading events to save unneeded computation
    if args.uniform_slip:
        if float(args.uniform_slip) < 0: raise "Slip must be positive"
    
    if args.field_plot:
        if args.model_file is None:
            raise "Must specify --model_file for field plots"
        elif args.field_type is None:
            raise "Must specify --field_type for field plots"
            
    if args.traces:
        if args.model_file is None:
            raise "Must specify --model_file for fault trace plots"
            
    # Check that if either beta or tau is given then the other is also given
    if (args.beta and not args.tau) or (args.tau and not args.beta):
        raise "Must specify both beta and tau."
        
    # Check that field_type is one of the supported types
    if args.field_type:
        type = args.field_type.lower()
        if type != "gravity" and type != "dilat_gravity" and type != "displacement" and type != "insar" and type!= "potential" and type != "geoid":
            raise "Field type is one of gravity, dilat_gravity, displacement, insar, potential, geoid"
    # ------------------------------------------------------------------------

    # Read the event and sweeps files
    if args.event_file:
        if not os.path.isfile(args.event_file):
            raise "Event file does not exist: "+args.event_file
        else:
            events = Events(args.event_file, args.sweep_file)

    # Read the geometry model if specified
    if args.model_file:
        if args.model_file_type:
            geometry = Geometry(model_file=args.model_file, model_file_type=args.model_file_type)
        else:
            geometry = Geometry(model_file=args.model_file)

    # Read the stress files if specified
    if args.stress_index_file and args.stress_file:
        stress_set = quakelib.ModelStressSet()
        stress_set.read_file_ascii(args.stress_index_file, args.stress_file)
    else:
        stress_set = None
        
    if args.all_stat_plots:
        args.plot_freq_mag = 3
        args.plot_mag_rupt_area = True
        args.plot_mag_mean_slip = True
        args.wc94 = True

    # Set up filters
    event_filters = []
    if args.min_magnitude or args.max_magnitude:
        event_filters.append(MagFilter(min_mag=args.min_magnitude, max_mag=args.max_magnitude))

    if args.min_year or args.max_year:
        event_filters.append(YearFilter(min_year=args.min_year, max_year=args.max_year))

    if args.min_slip or args.max_slip:
        # Setting default lower limit on mean slip of 1cm
        #if args.min_slip is None: args.min_slip = 0.01
        event_filters.append(SlipFilter(min_slip=args.min_slip, max_slip=args.max_slip))

    if args.min_event_num or args.max_event_num:
        event_filters.append(EventNumFilter(min_mag=args.min_event_num, max_mag=args.max_event_num))

    if args.use_sections:
        if not args.model_file:
            raise "Must specify --model_file for --use_sections to work."
        event_filters.append(SectionFilter(geometry, args.use_sections))
        # Also grab all the elements from this section in case this is being used to grab element ids
        if args.elements is None:
            args.elements = [elem_num for elem_num in range(geometry.model.num_elements()) if geometry.model.element(elem_num).section_id() in args.use_sections]
        
        
    if args.use_trigger_sections:
        if not args.model_file:
            raise "Must specify --model_file for --use_trigger_sections to work."
        event_filters.append(TriggerSectionFilter(geometry, args.use_trigger_sections))

    if args.event_file:
        events.set_filters(event_filters)
        
    # Print out event summary data if requested
    if args.summary:
        print("\n"+args.event_file)
        events.largest_event_summary(args.summary)

    # Generate plots
    if args.diagnostics:
        args.num_sweeps = True
        args.event_shear_stress = True
        args.event_normal_stress = True
        args.event_mean_slip = True
        args.plot_freq_mag = 3
        args.plot_mag_mean_slip = True
        args.plot_mag_rupt_area = True
        args.wc94 = True
    if args.plot_freq_mag:
        filename = SaveFile().event_plot(args.event_file, "freq_mag", args.min_magnitude)
        if args.plot_freq_mag == 1: UCERF2,b1 = False, False
        if args.plot_freq_mag == 2: UCERF2,b1 = False, True
        if args.plot_freq_mag == 3: UCERF2,b1 = True, False
        if args.plot_freq_mag == 4: UCERF2,b1 = True, True
        FrequencyMagnitudePlot().plot(events, filename, UCERF2=UCERF2, b1=b1)
    if args.plot_mag_rupt_area:
        filename = SaveFile().event_plot(args.event_file, "mag_rupt_area", args.min_magnitude)
        MagnitudeRuptureAreaPlot().plot(events, filename, WC94=args.wc94)
    if args.plot_mag_mean_slip:
        filename = SaveFile().event_plot(args.event_file, "mag_mean_slip", args.min_magnitude)
        MagnitudeMeanSlipPlot().plot(events, filename, WC94=args.wc94)
    if args.plot_prob_vs_t:
        filename = SaveFile().event_plot(args.event_file, "prob_vs_time", args.min_magnitude)
        ProbabilityPlot().plot_p_of_t(events, filename)
    if args.plot_prob_vs_t_fixed_dt:
        filename = SaveFile().event_plot(args.event_file, "p_vs_t_fixed_dt", args.min_magnitude)
        ProbabilityPlot().plot_conditional_fixed_dt(events, filename)
    if args.plot_cond_prob_vs_t:
        filename = SaveFile().event_plot(args.event_file, "cond_prob_vs_t", args.min_magnitude)
        if args.beta:
            ProbabilityPlot().plot_p_of_t_multi(events, filename, beta=args.beta, tau=args.tau)
        else:
            ProbabilityPlot().plot_p_of_t_multi(events, filename)
    if args.plot_waiting_times:
        filename = SaveFile().event_plot(args.event_file, "waiting_times", args.min_magnitude)
        ProbabilityPlot().plot_dt_vs_t0(events, filename)
    if args.field_plot:
        type = args.field_type.lower()
        if args.colorbar_max: cbar_max = args.colorbar_max
        else: cbar_max = None
        if args.levels: levels = args.levels
        else: levels = None
        filename = SaveFile().field_plot(args.model_file, type, args.uniform_slip, args.event_id)
        if args.angles: 
            if len(args.angles) != 2:
                raise "Must specify 2 angles"
            else:
                angles = np.array(args.angles)*np.pi/180.0
        else: angles = None
        if args.event_id is None:
            element_ids = geometry.model.getElementIDs()
            ele_slips = {}
            if args.uniform_slip is None: uniform_slip = 5.0
            else: uniform_slip = args.uniform_slip
            sys.stdout.write(" Computing field for uniform slip {}m :".format(int(uniform_slip)))
            for ele_id in element_ids:
                ele_slips[ele_id] = uniform_slip
            event = None
        else:
            sys.stdout.write(" Processing event {}, M={:.2f} : ".format(args.event_id, events._events[args.event_id].getMagnitude()))
            ele_slips = events.get_event_element_slips(args.event_id)
            event = events._events[args.event_id]
        
        if len(ele_slips.keys()) == 0:
            raise "Error in processing slips."
        else:
            sys.stdout.write(" Loaded slips for {} elements :".format(len(ele_slips.keys()))) 
        sys.stdout.flush()
        
        FP = FieldPlotter(geometry, args.field_type, element_slips=ele_slips, event=event, event_id=args.event_id, cbar_max=cbar_max, levels=levels, small_model=args.small_model, g0=args.g)
        FP.compute_field(cutoff=1000)
        FP.plot_field(output_file=filename, angles=angles)
    if args.field_eval:
        filename = SaveFile().field_plot(args.model_file, "displacement", args.uniform_slip, args.event_id)
        sys.stdout.write(" Processing event {}, M={:.2f} : ".format(args.event_id, events._events[args.event_id].getMagnitude()))
        ele_slips = events.get_event_element_slips(args.event_id)
        event = events._events[args.event_id]
        if len(ele_slips.keys()) == 0:
            raise "Error in processing slips."
        else:
            sys.stdout.write(" Loaded slips for {} elements :".format(len(ele_slips.keys()))) 
        sys.stdout.flush()
        FE = FieldEvaluator(geometry, args.event_id, event, ele_slips, args.lld_file)
        FE.compute_field()
    if args.greens:
        # Set default values
        if args.dip is None: sys.exit("Must specify --dip")
        if args.rake is None: sys.exit("Must specify --rake")
        if args.Nx is None: args.Nx = 690
        if args.Ny is None: args.Ny = 422
        if args.Xmin is None: args.Xmin = -5000
        if args.Xmax is None: args.Xmax = 15000
        if args.Ymin is None: args.Ymin = -10000
        if args.Ymax is None: args.Ymax = 10000
        if args.L is None: args.L = 10000
        if args.W is None: args.W = 10000
        if args.DTTF is None: args.DTTF = 1000
        if args.g is None: args.g = 9.81
        if args._lambda is None: args._lambda = 3.2e10
        if args.mu is None: args.mu = 3.2e10
        filename = SaveFile().greens_plot(args.plot_name, args.field_type, args.uniform_slip)
        GP = GreensPlotter(args.field_type, cbar_max=args.colorbar_max, levels=args.levels, Nx=args.Nx, Ny=args.Ny, Xmin=args.Xmin, Xmax=args.Xmax, Ymin=args.Ymin, Ymax=args.Ymax, L=args.L, W=args.W, DTTF=args.DTTF, slip=args.uniform_slip, dip=args.dip, _lambda=args._lambda, _mu=args.mu, rake=args.rake, g0=args.g)
        GP.compute_field()
        GP.plot_field(filename)
    if args.traces:
        filename = SaveFile().trace_plot(args.model_file)
        if args.small_model is None: args.small_model = False
        TP = TracePlotter(geometry, filename, use_sections=args.use_sections, small_model=args.small_model)

    if args.slip_rates:
        if args.elements is None: args.elements = geometry.model.getElementIDs()
        slip_rates = geometry.get_slip_rates(args.elements)
        for id in slip_rates.keys():
            sys.stdout.write("{}  {}\n".format(id,slip_rates[id]))
            
    if args.slip_time_series:
        if args.elements is None: raise "Must specify element ids, e.g. --elements 0 1 2"
        if args.min_year is None: args.min_year = 0.0
        if args.max_year is None: args.max_year = 20.0
        if args.dt is None: args.dt = 0.5  # Unit is decimal years
        if args.use_sections is not None:
            if len(args.use_sections) > 1:
                section_name = ""
                for sec in args.use_sections:
                    section_name += geometry.model.section(sec).name()+", "
            else:
                section_name = geometry.model.section(args.use_sections[0]).name()+", "
        time_series = geometry.get_slip_time_series(events, elements=args.elements, min_year=args.min_year, max_year=args.max_year, DT=args.dt)
        if len(time_series.keys()) < 10: 
            labels = time_series.keys()+[""]
        else:
            labels = [None for each in range(len(time_series.keys())+1)]
        x_data = [list(np.arange(args.min_year+args.dt, args.max_year+args.dt, args.dt)) for key in time_series.keys()]+[[args.min_year,args.max_year]]
        linewidths = [0.8 for key in time_series.keys()]+[1]
        styles = ["-" for key in time_series.keys()]+["--"]
        y_data = time_series.values()+[[0,0]]
        if args.use_sections is not None:
            plot_title = "Slip time series for {}from years {} to {} with step {}\n{}".format(section_name, args.min_year,args.max_year,args.dt,args.event_file.split("/")[-1])
        else:
            plot_title = "Slip time series for {} elements, from years {} to {} with step {}\n{}".format(len(args.elements), args.min_year,args.max_year,args.dt,args.event_file.split("/")[-1])
        filename = SaveFile().diagnostic_plot(args.event_file, "slip_time_series", min_year=args.min_year, max_year=args.max_year)
        BasePlotter().multi_line_plot(x_data, y_data, labels, linewidths, plot_title, "sim time [years]", "cumulative slip [m]", "", filename, linestyles=styles)

    if args.event_kml:
        if args.event_id is None or args.event_file is None or args.model_file is None:
            raise "Must specify an event to plot with --event_id and provide an --event_file and a --model_file."
        else:
            event = events._events[args.event_id]
            filename = SaveFile().event_kml_plot(args.event_file, args.event_id)
            geometry.model.write_event_kml(filename, event)

    # Generate stress plots
    if args.stress_elements:
# TODO: check that stress_set is valid
        StressHistoryPlot().plot(stress_set, args.stress_elements)
        
    if args.num_sweeps:
        filename = SaveFile().diagnostic_plot(args.event_file, "num_sweeps", min_year=args.min_year, max_year=args.max_year)
        DiagnosticPlot().plot_number_of_sweeps(events, filename)
    if args.event_shear_stress:
        filename = SaveFile().diagnostic_plot(args.event_file, "shear_stress", min_year=args.min_year, max_year=args.max_year)
        DiagnosticPlot().plot_shear_stress_changes(events, filename)
    if args.event_normal_stress:
        filename = SaveFile().diagnostic_plot(args.event_file, "normal_stress", min_year=args.min_year, max_year=args.max_year)
        DiagnosticPlot().plot_normal_stress_changes(events, filename)
    if args.event_mean_slip:
        filename = SaveFile().diagnostic_plot(args.event_file, "mean_slip", min_year=args.min_year, max_year=args.max_year)
        DiagnosticPlot().plot_mean_slip(events, filename)

    # Validate data if requested
    err = False
    if args.validate_slip_sum:
        mean_slip = sum(events.event_mean_slip())
        if abs(mean_slip-args.validate_slip_sum)/args.validate_slip_sum > 0.01: err = True
        print("Calculated mean slip:", mean_slip, "vs. expected:", args.validate_slip_sum)

    if args.validate_mean_interevent:
        ie_times = events.interevent_times()
        mean_ie = sum(ie_times)/len(ie_times)
        if abs(mean_ie-args.mean_interevent)/args.mean_interevent > 0.02: err = True
        print("Calculated mean interevent:", mean_interevent, "vs. expected:", args.mean_interevent)

    if err: exit(1)

