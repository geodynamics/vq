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
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import matplotlib.font_manager as mfont
    import matplotlib.colors as mcolor
    import matplotlib.colorbar as mcolorbar
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches
    from PIL import Image 
    import matplotlib.animation as manimation
    # we only want to execute this in the __main__ part of the script, so we can also run plotting scripts interactively.
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
    
cPickle_available = True
try:
    import cPickle as pickle
except ImportError:
    cPickle_available = False
    
    
# ----------------- Global constants -------------------------------------------
# Kasey: These are only relevent for few-element field plots
LAT_LON_DIFF_FACTOR = 1.333 
MIN_LON_DIFF = 0.4   # 1 corresponds to ~ 100km at lat,lon = (40.35, -124.85)
MIN_LAT_DIFF = MIN_LON_DIFF/LAT_LON_DIFF_FACTOR   # 0.8 corresponds to ~ 100km at lat,lon = (40.35, -124.85)
#MIN_FIT_MAG  = 5.0     # lower end o08 f magnitude for fitting freq_mag plot with b=1 curve

STAT_COLOR_CYCLE = ['k','b','cyan','purple','g']
SCATTER_ALPHA = 0.5
SCATTER_SIZE = 10

#-------------------------------------------------------------------------------
# Given a set of maxes and mins return a linear value betweem them.
# Used to compute cutoff for field value evaluation, cutoff scales with
# number of involved elements for FieldPlotter instances.
#-------------------------------------------------------------------------------
def linear_interp(x, x_min, x_max, y_min, y_max):
    return ((y_max - y_min)/(x_max - x_min) * (x - x_min)) + y_min
    
def fit_to_weibull(prob_x, prob_y, beta_guess, tau_guess):
    # Least squares fitting to the Weibull distribution
    from scipy.optimize import leastsq
    errorFunction = lambda beta_tau,x,y: Distributions().weibull(x, beta_tau[0], beta_tau[1])-y
    beta_tau_guess = (beta_guess,tau_guess)
    fit_params, success = leastsq(errorFunction, beta_tau_guess[:], args=(prob_x,prob_y))
    return fit_params
    
def calculate_averages(x,y,log_bin=False,num_bins=None):
    if num_bins is None:
        num_bins = math.floor(len(x)/100)
        if num_bins < 20:
            num_bins = 20
        elif num_bins > 100:
            num_bins = 100
    x = np.array(x)
    y = np.array(y)
    if log_bin: bin_min = math.floor(math.log(np.min(x),10))
    else: bin_min = math.floor(np.min(x))
    if log_bin: bin_max = math.ceil(math.log(np.max(x),10))
    else: bin_max = math.ceil(np.max(x))
    if log_bin: bins = np.logspace(bin_min,bin_max,num=num_bins)
    else: bins = np.linspace(bin_min,bin_max,num=num_bins)
    inds = np.digitize(x, bins)
    binned_data = {}
    for n, i in enumerate(inds):
        # i-> i-1, bins are 1-indexed, want arrays storing 0-indexed stuff
        try:
            binned_data[i-1].append(y[n])
        except KeyError:
            binned_data[i-1] = [y[n]]
    x_ave = []
    y_ave = []
    for k in sorted(binned_data.keys()):
        if k != 0:
            x_ave.append(0.5*(bins[k-1]+bins[k]))
            y_ave.append(sum(binned_data[k])/float(len(binned_data[k])))
    return x_ave, y_ave


def standardize_time_series(time_series):
    # This function subtracts the mean from the time series 
    #    and scales it by the standard deviation.
    time_series = np.array(time_series)
    series_mean = np.mean(time_series)
    series_std  = np.std(time_series)
    return (time_series - series_mean)/series_std


def get_fault_group_average_time_series_from_pickle(file_list):
    if not cPickle_available:
        raise BaseException("\nRequested time series cannot be fetched, cPickle module is not available.")
    else:
        num_faults = len(file_list)
        average_time_series = [] ## create the object used for averaging
        num_points = 0
        for i, file_name in enumerate(file_list):
            this_file = open(file_name, 'r')
            time_data, time_series = pickle.load(this_file)
            this_file.close()
            sys.stdout.write("Read time series from {}\n".format(file_name))
            if i==0:
                ## Being slick here and creating the average time series as soon as we know the number of points per time series.
                ## This should make it faster for very large numbers time series points.
                num_points = len(time_series)
                average_time_series = np.zeros(num_points)
            else:
                ## If we have read a time series already, make sure we read the same number of data points as the previous series.
                assert(len(time_series)==num_points)
            average_time_series += time_series/float(num_faults)
        return [time_data, average_time_series]

class SaveFile:
    def __init__(self):
        if args.pdf: self.file_type = '.pdf'
        elif args.eps: self.file_type = '.eps'
        else: self.file_type = '.png'

    def event_plot(self, event_file, plot_type, min_mag, min_year, max_year, combine):
        # Add tags to convey the subsets/cuts being made
        add=""
        if len(event_file) > 1: 
            add += "_MULTI_EVENT_FILE"
        event_file = event_file[0]
        min_mag = str(min_mag)
        # Remove any folders in front of model_file name
        if len(event_file.split("/")) > 1:
            event_file = event_file.split("/")[-1]
        if min_year is not None: add+="_yearMin"+str(int(min_year))    
        if max_year is not None: add+="_yearMax"+str(int(max_year))
        if args.use_sections is not None:
            add+="_triggeredBySection_"
            for sec in args.use_sections:
                add+="-"+geometry.model.section(sec).name()
        if args.use_faults is not None:
            add+="_triggeredByFault_"
            for fault in args.use_faults:
                add+="-"+str(fault)
        if min_mag is not None: 
            # e.g. min_mag = 7.5, filename has '7-5'
            if len(min_mag.split(".")) > 1:
                add += "_minMag_"+min_mag.split(".")[0]+"-"+min_mag.split(".")[1]
            else:
                add += "_minMag_"+min_mag
        if combine is not None:
            add+="_combined"

        return plot_type+add+"_"+event_file.split(".")[0]+self.file_type
        
    def field_plot(self, model_file, field_type, uniform_slip, event_id, wavelength):
        # Remove any folders in front of model_file name
        if len(model_file.split("/")) > 1:
            model_file = model_file.split("/")[-1]
        if wavelength is None:
            wave = ""
        else:
            wave = "_"+str(int(round(wavelength*100,0)))+"cm"
        if uniform_slip is None and event_id is not None:
            return model_file.split(".")[0]+"_"+field_type+"_event"+str(event_id)+wave+self.file_type
        elif uniform_slip is not None and event_id is None:
            return model_file.split(".")[0]+"_"+field_type+"_uniform_slip"+str(int(uniform_slip))+"m"+wave+self.file_type
        else:
            raise BaseException("\nMust specify either uniform_slip or event_id")
            
    def greens_plot(self, name, field_type, slip):
        return "greens_"+field_type+"_"+name+"_slip"+str(int(slip))+"m"+self.file_type
            
    def trace_plot(self, model_file):
        # Remove any folders in front of model_file name
        if len(model_file.split("/")) > 1:
            model_file = model_file.split("/")[-1]
        return "traces_"+model_file.split(".")[0]+self.file_type

    def distribution_plot(self, model_file, type):
        # Remove any folders in front of model_file name
        if len(model_file.split("/")) > 1:
            model_file = model_file.split("/")[-1]
        return type+"_"+model_file.split(".")[0]+self.file_type
        
    def diagnostic_plot(self, event_file, plot_type, min_year=None, max_year=None, min_mag=None, combine=None):
        # Add tags to convey the subsets/cuts being made
        add=""
        if len(event_file) > 1: 
            add += "_MULTI_EVENT_FILE"
        event_file = event_file[0]
        # Remove any folders in front of model_file name
        if len(event_file.split("/")) > 1:
            event_file = event_file.split("/")[-1]
        if min_year is not None: add+="_yearMin"+str(int(min_year))    
        if max_year is not None: add+="_yearMax"+str(int(max_year))
        if args.use_sections is not None:
            for sec in args.use_sections:
                add+="_"+geometry.model.section(sec).name()
        if min_mag is not None:
            min_mag = str(min_mag)
            # e.g. min_mag = 7.5, filename has '7-5'
            if len(min_mag.split(".")) > 1:
                add += "_minMag_"+min_mag.split(".")[0]+"-"+min_mag.split(".")[1]
            else:
                add += "_minMag_"+min_mag
        if combine is not None:
            add += "_combined"
                        
        return plot_type+"_diagnostic"+add+"_"+event_file.split(".")[0]+self.file_type

    def time_series_plot(self, event_file, plot_type, min_year=None, max_year=None, min_mag=None, combine=None, dt=None):
        # Add tags to convey the subsets/cuts being made
        add=""
        if len(event_file) > 1: 
            add += "_MULTI_EVENT_FILE"
        event_file = event_file[0]
        # Remove any folders in front of model_file name
        if len(event_file.split("/")) > 1:
            event_file = event_file.split("/")[-1]
        if min_year is not None: add+="_yearMin"+str(int(min_year))    
        if max_year is not None: add+="_yearMax"+str(int(max_year))
        if min_mag is not None:
            min_mag = str(min_mag)
            # e.g. min_mag = 7.5, filename has '7-5'
            if len(min_mag.split(".")) > 1:
                add += "_minMag_"+min_mag.split(".")[0]+"-"+min_mag.split(".")[1]
            else:
                add += "_minMag_"+min_mag
        if combine is not None:
            add += "_combined"
        if dt is not None:
            add += "_dt{}yrs".format(dt)
        return plot_type+"_timeseries"+add+"_"+event_file.split(".")[0]+self.file_type
        
    def fault_time_series_pickle(self, event_file, fault_id, min_year=None, max_year=None, min_mag=None, combine=None, dt=None, standardized=False):
        # Add tags to convey the subsets/cuts being made
        add=""
        event_file = event_file[0]
        # Remove any folders in front of model_file name
        if len(event_file.split("/")) > 1:
            event_file = event_file.split("/")[-1]
        if standardized:
            add += "_standardized"
        if min_year is not None: add+="_yearMin"+str(int(min_year))    
        if max_year is not None: add+="_yearMax"+str(int(max_year))
        if min_mag is not None:
            min_mag = str(min_mag)
            # e.g. min_mag = 7.5, filename has '7-5'
            if len(min_mag.split(".")) > 1:
                add += "_minMag_"+min_mag.split(".")[0]+"-"+min_mag.split(".")[1]
            else:
                add += "_minMag_"+min_mag
        if combine is not None:
            add += "_combined"
        if dt is not None:
            add += "_dt{}yrs".format(dt)
        return "fault_"+str(fault_id)+"_timeseries"+add+"_"+event_file.split(".")[0]+".pickle"

    def event_movie(self, event_file, event_id):
        # Remove any folders in front of model_file name
        if len(event_file.split("/")) > 1:
            event_file = event_file.split("/")[-1]
        return "movie_event_{}_{}.mp4".format(event_id, event_file.split(".")[0])
    
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
        label_str = "  "
# TODO: change to <= character
        if self._min_mag != -float("inf"): label_str += str(self._min_mag)+"<"
        label_str += "M"
        if self._max_mag != float("inf"): label_str += "<"+str(self._max_mag)
        return label_str

class YearFilter:
    def __init__(self, min_year=None, max_year=None):
        self._min_year = min_year if min_year is not None else -float("inf")
        self._max_year = max_year if max_year is not None else float("inf")

    def test_event(self, event):
        return (event.getEventYear() >= self._min_year and event.getEventYear() <= self._max_year)

    def plot_str(self):
        label_str = "  "
# TODO: change to <= character
        if self._min_year != -float("inf"): label_str += str(self._min_year)+"<"
        label_str += "year"
        if self._max_year != float("inf"): label_str += "<"+str(self._max_year)
        return label_str

class EventNumFilter:
    def __init__(self, min_event_num=None, max_event_num=None):
        self._min_event_num = min_event_num if min_event_num is not None else -sys.maxint
        self._max_event_num = max_event_num if max_event_num is not None else sys.maxint

    def test_event(self, event):
        return (event.getEventNumber() >= self._min_event_num and event.getEventNumber() <= self._max_event_num)

    def plot_str(self):
        label_str = "  "
# TODO: change to <= character
        if self._min_event_num != -sys.maxint: label_str += str(self._min_event_num)+"<"
        label_str += "event num"
        if self._max_event_num != sys.maxint: label_str += "<"+str(self._max_event_num)
        return label_str

class NumElementsFilter:
    def __init__(self, min_num_elements=None, max_num_elements=None):
        self._min_num_elements = min_num_elements if min_num_elements is not None else 0
        self._max_num_elements = max_num_elements if max_num_elements is not None else sys.maxint
    
    def test_event(self, event):
        return (len(event.getInvolvedElements()) >= self._min_num_elements and len(event.getInvolvedElements()) <= self._max_num_elements)
    
    def plot_str(self):
        label_str = "  "
        # TODO: change to <= character
        if self._min_num_elements != 0: label_str += str(self._min_num_elements)+"<"
        label_str += "num elements"
        if self._max_num_elements != sys.maxint: label_str += "<"+str(self._max_num_elements)
        return label_str

        
class SlipFilter:
    def __init__(self, min_slip=None, max_slip=None):
        self._min_slip = min_slip if min_slip is not None else -float("inf")
        self._max_slip = max_slip if max_slip is not None else float("inf")

    def test_event(self, event):
        return (event.calcMeanSlip() >= self._min_slip and event.calcMeanSlip() <= self._max_slip)

    def plot_str(self):
        label_str = "   "
# TODO: change to <= character
        if self._min_slip != -float("inf"): label_str += str(self._min_slip)+"<"
        label_str += "slip"
        if self._max_slip != float("inf"): label_str += "<"+str(self._max_slip)
        return label_str
        
class AreaFilter:
    def __init__(self, min_area=None, max_area=None):
        # Convert from the input km^2 to the unit of calcEventRuptureArea() which is m^2
        self._min_area = quakelib.Conversion().sqkm2sqm(min_area) if min_area is not None else -float("inf")
        self._max_area = quakelib.Conversion().sqkm2sqm(max_area) if max_area is not None else float("inf")

    def test_event(self, event):
        return (event.calcEventRuptureArea() >= self._min_area and event.calcEventRuptureArea() <= self._max_area)

    def plot_str(self):
        label_str = "  "
# TODO: change to <= character
        if self._min_area != -float("inf"): label_str += str(self._min_area)+"<"
        label_str+="area"
        if self._max_area != float("inf"): label_str += "<"+str(self._max_area)
        return label_str
        
class Geometry:
    def __init__(self, model_file=None, model_file_type=None):
        if model_file is not None:
            self.model = quakelib.ModelWorld()
            if model_file_type =='text' or model_file.split(".")[-1] == 'txt':
                self.model.read_file_ascii(model_file)
                sys.stdout.write("Read fault model from {}\n".format(model_file))
            elif model_file_type == 'hdf5' or model_file.split(".")[-1] == 'h5' or model_file.split(".")[-1] == 'hdf5':
                self.model.read_file_hdf5(model_file)
                sys.stdout.write("Read fault model from {}\n".format(model_file))
            else:
                raise BaseException("\nMust specify --model_file_type, either hdf5 or text")
            self._elem_to_section_map = {elem_num: self.model.element(elem_num).section_id() for elem_num in self.model.getElementIDs()}
            self._elem_to_fault_map = {elem_num: self.model.section(self.model.element(elem_num).section_id()).fault_id() for elem_num in self.model.getElementIDs()}
            ###!!!!!
            #sys.stdout.write(list(self.model.getElementIDs()))
            #sys.stdout.write("Elem to Section : {}".format(self._elem_to_section_map))
            #sys.stdout.write("Elem to Fault: {}".format(self._elem_to_fault_map))
            ###!!!!!
            
        else:
            if args.use_sections:
                raise BaseException("\nModel file required if specifying fault sections.")
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
        
    def get_aseismics(self, elements):
        return {id:self.model.element(block_id).aseismic() for id in elements}
        
    def get_slip_time_series(self, events, elements=None, max_year=None, DT=None):
        if max_year is None:
            raise BaseException("Must specify --max_year.")
        # slip_time_series    = dictionary indexed by block_id with entries being arrays of absolute slip at each time step
        # Get slip rates for the elements
        slip_rates = self.get_slip_rates(elements)
        aseismic_fracs = self.get_aseismics(elements)
        #Initialize blocks with 0.0 slip at time t=0.0
        slip_time_series  = {id:[0.0] for id in elements}
        # Grab the events data
        event_years = events.event_years()
        event_numbers = events.event_numbers()
        #Initialize time steps to evaluate slip    
        time_values = np.arange(DT, max_year+DT, DT)
        for k in range(len(time_values)):
            if k>0:
                # current time in simulation
                right_now = time_values[k]
                # back slip all elements by subtracting the slip_rate*dt
                for block_id in slip_time_series.keys():
                    last_slip = slip_time_series[block_id][k-1]
                    this_slip = (1-geometry.model.element(block_id).aseismic())*slip_rates[block_id]*DT
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
        
        
    def get_fault_averaged_slip_time_series(self, events, fault_id=None, max_year=None, DT=1.0, standardized=False):
        if max_year is None: raise BaseException("Must specify --max_year.")
        sys.stdout.write("{}...".format(fault_id)) # Update the console as to which fault we are processing now.
        sys.stdout.flush()  ## Make sure the status prints right now and doesn't lag due to hi CPU demand
        # slip_time_series    = dictionary indexed by block_id with entries being arrays of slip deficit at each time step
        # Get slip rates for the elements
        all_elements = geometry.model.getElementIDs()
        elements     = [elem_num for elem_num in all_elements if geometry.model.section(geometry.model.element(elem_num).section_id()).fault_id() == fault_id]
        num_elements = float(len(elements))
        slip_rates   = self.get_slip_rates(elements)
        # Grab the events data
        event_years = list(events.event_years())
        # Restrict the time series values to end after the last EQ
        event_numbers = list(events.event_numbers())
        #Initialize time steps to evaluate slip    
        time_values = np.arange(0.0, max_year+DT, DT)
        #Pre-allocate the full time series array for speed 
        slip_time_series  = {id:np.zeros(len(time_values)) for id in elements}
        for k in range(len(time_values)):
            # Skip the k=0 case, we already initialized our slip
            if k>0:
                # current time in simulation (think of this as the right edge of a bin with width DT)
                right_now = time_values[k]
                # back slip all elements by subtracting the slip_rate*dt
                for block_id in slip_time_series.keys():
                    last_slip = slip_time_series[block_id][k-1]
                    # Must reduce the slip rate by the aseismic fraction.
                    this_slip = (1-geometry.model.element(block_id).aseismic())*slip_rates[block_id]*DT
                    slip_time_series[block_id][k] = last_slip-this_slip
                # check if any elements slip as part of simulated event in the window of simulation time
                # between (current time - DT, current time), add event slips to the slip at current time 
                # for elements involved
                for j, evid in enumerate(event_numbers):
                    ev_year = event_years[j]
                    if right_now-DT < ev_year <= right_now:
                        event_element_slips = events.get_event_element_slips(evid)
                        for block_id in event_element_slips.keys():
                            try:
                                slip_time_series[block_id][k] += event_element_slips[block_id]
                                #sys.stdout.write("element {} slips {} in event {}\n".format(block_id,event_element_slips[block_id],evid))
                                #sys.stdout.flush()
                            except KeyError:
                                pass # Ignore event elements that we are not asking for
                                
        ### Average the time series over all elements in the fault
        fault_time_series = np.zeros(len(time_values))
        for ele_id in slip_time_series.keys():
            fault_time_series += slip_time_series[ele_id]/num_elements
        
        ## Standardize the averaged time series
        if standardized:
            return standardize_time_series(fault_time_series)
        else:
            return fault_time_series
            
        


    def get_stress_drops(self):
        return [self.model.element(ele).stress_drop() for ele in self.model.getElementIDs()]
    
    def get_stress_drop_factor(self):
        return self.model.stressDropFactor()



class InvolvedSectionFilter:
    def __init__(self, geometry, event_file, section_list):
        self._section_list = section_list
        self._elem_to_section_map = {elem_num: geometry.model.element(elem_num).section_id() for elem_num in geometry.model.getElementIDs()}
        self._event_file = event_file
        self._all_sweeps = AllSweeps(self._event_file)
        #self._elem_to_section_map = geometry._elem_to_section_map

    def test_event(self, event):
        #event_sweeps = Sweeps(self._event_file, event.getEventNumber())
        event_id = event.getEventNumber()
        event_secs = [self._elem_to_section_map[elID] for elID in self._all_sweeps.event_elements[event_id]]
        for sec in self._section_list:
            if sec in event_secs: return True
        return False

    def plot_str(self):
        label_str = "  partSections"
        for sec in self._section_list:
            label_str += "-"+geometry.model.section(sec).name()
        return label_str


class TriggerSectionFilter:
    def __init__(self, geometry, section_list):
        self._section_list = section_list
        self._elem_to_section_map = {elem_num: geometry.model.element(elem_num).section_id() for elem_num in geometry.model.getElementIDs()}
        #self._elem_to_section_map = geometry._elem_to_section_map

    def test_event(self, event):
        triggerID = event.getEventTrigger()
        elem_section = self._elem_to_section_map[triggerID]
        if elem_section in self._section_list: return True

        return False

    def plot_str(self):
        label_str = "  triggerSections"
        for sec in self._section_list:
            label_str += "-"+geometry.model.section(sec).name()
        return label_str


class TriggerFaultFilter:
    def __init__(self, geometry, fault_list):
        self._fault_list = fault_list
        #self._elem_to_fault_map = geometry._elem_to_fault_map
        self._elem_to_fault_map = {elem_num: geometry.model.section(geometry.model.element(elem_num).section_id()).fault_id() for elem_num in geometry.model.getElementIDs()}

    def test_event(self, event):
        triggerID = event.getEventTrigger()
        elem_fault = self._elem_to_fault_map[triggerID]
        if elem_fault in self._fault_list: 
            #sys.stdout.write("Event {} was triggered on fault {}\n".format(event.getEventNumber(),elem_fault))
            return True
        return False

    def plot_str(self):
        label_str = "  triggerFaults"
        for fault in self._fault_list:
            # TODO: Add fault names once ModelFault::name() exists
            #label_str += "-"+geometry.model.fault(fault).name()
            label_str += "-"+str(fault)
        return label_str


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

def read_all_sweeps_h5(sim_file, block_ids=None):
    # Read sweeps sequence for multiple blocks (unless block_id specified)
    with h5py.File(sim_file) as vq_data:
        sweeps = vq_data['sweeps'][()]
    # If block_id specified, only return those sweeps for that block
    if block_ids is not None:
        d_type = sweeps.dtype
        sweeps = np.core.records.fromarrays(zip(*filter(lambda x: x['block_id'] in block_ids, sweeps)), dtype=d_type)	
    return sweeps

def parse_all_sweeps_h5(sim_file=None, block_id=None, do_print=True, sweeps=None):
    # Read sweep data if not provided
    if sweeps is None: sweeps = read_all_sweeps_h5(sim_file, block_id=block_id)
    
    data = [[rw['event_number'], rw['sweep_number'], rw['block_id'], rw['block_slip'], rw['shear_init'],
             rw['shear_final'], rw['normal_init'],rw['normal_final'], 
             (rw['shear_final']-rw['shear_init'])/rw['shear_init'], 
             (rw['normal_final']-rw['normal_init'])/rw['normal_init']] for rw in sweeps]
    if do_print:
        for rw in data: sys.stdout.write(rw)
    cols = ['event_number', 'sweep_number', 'block_id', 'block_slip', 'shear_init', 
            'shear_final', 'normal_init', 'normal_final', 'shear_change', 'normal_change']
    return np.core.records.fromarrays(zip(*data), names=cols, formats = [type(x).__name__ for x in data[0]])

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
    
    data = [[rw['sweep_number'], rw['block_id'], rw['block_slip'], rw['shear_init'],
             rw['shear_final'], rw['normal_init'],rw['normal_final'], 
             (rw['shear_final']-rw['shear_init'])/rw['shear_init'], 
             (rw['normal_final']-rw['normal_init'])/rw['normal_init']] for rw in sweeps]
    if do_print:
        for rw in data: sys.stdout.write(rw)
    cols = ['sweep_number', 'block_id', 'block_slip', 'shear_init', 
            'shear_final', 'normal_init', 'normal_final', 'shear_change', 'normal_change']
    return np.core.records.fromarrays(zip(*data), names=cols, formats = [type(x).__name__ for x in data[0]])
    
    
class Events:
    def __init__(self, event_file, sweep_file = None, combine_file=None, stress_file=None, stress_index_file=None):
        filetype = event_file.split('.')[-1].lower()
        event_file_type = "text" # default
        if filetype == 'h5' or filetype == 'hdf5': event_file_type = "hdf5"
        if event_file_type == "hdf5":
            # Reading in via QuakeLib
            #if not h5py_available:
            self._events = quakelib.ModelEventSet()
            sys.stdout.write("Reading the HDF5 event file.....".format(event_file))
            self._events.read_file_hdf5(event_file)
            sys.stdout.write("Read in events via QuakeLib from {}\n".format(event_file))
            # Reading via h5py
            #else:
            #    self._events = read_events_h5(event_file)
            #    sys.stdout.write("Read in events via h5py from {}".format(event_file))
        elif event_file_type == "text" and sweep_file != None:
            self._events = quakelib.ModelEventSet()
            sys.stdout.write("Reading the TEXT event file.....".format(event_file))
            self._events.read_file_ascii(event_file, sweep_file)
            sys.stdout.write("Read in events via QuakeLib from {}\n".format(event_file))
        else:
            raise BaseException("\nevent_file_type must be hdf5 or text. If text, a sweep_file is required.")
            
        if combine_file is not None and event_file_type == 'hdf5' and stress_file is not None and not stress_file.split(".")[-1]=="txt":
            if not os.path.isfile(stress_file) or not os.path.isfile(combine_file):
                raise BaseException("\nOne or more files does not exist!")
            # If stress state was saved as hdf5
            with h5py.File(stress_file) as state_data:
                stress_state = state_data['stress_state'][()]
            stress_state = np.core.records.fromarrays(stress_state[-1], dtype=stress_state.dtype)
            add_year = float(stress_state['year'])
            add_evnum = int(stress_state['event_num'])
            self._events.append_from_hdf5(combine_file, add_year, add_evnum)
            sys.stdout.write("## Combined with: "+combine_file+"\n")
        elif combine_file is not None and event_file_type == 'hdf5' and stress_file is not None and stress_index_file is not None:
            if not os.path.isfile(stress_file) or not os.path.isfile(combine_file) or not os.path.isfile(stress_index_file):
                raise BaseException("\nOne or more files does not exist!")
            # If stress state was saved as text
            add_year, add_evnum, start_rec, end_rec = np.genfromtxt(stress_index_file)
            self._events.append_from_hdf5(combine_file, add_year, int(add_evnum))
            sys.stdout.write("## Combined with: "+combine_file+"\n")
            
        self._filtered_events = range(len(self._events))
        self._plot_str = ""

    def plot_str(self):
        return self._plot_str

    def set_filters(self, filter_list):
        #self._filtered_events = [evnum for evnum in range(len(self._events))]
        #self._plot_str = ""
        for cur_filter in filter_list:
            new_filtered_events = [evnum for evnum in self._filtered_events if cur_filter.test_event(self._events[evnum])]
            self._filtered_events = new_filtered_events
            self._plot_str += cur_filter.plot_str()
        if len(self._filtered_events) == 0:
            raise BaseException("\nNo events matching filters found!")

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
        
    def event_summary(self, evnums, geometry):
        mags = [self._events[evnum].getMagnitude() for evnum in evnums if self._events[evnum].getMagnitude() != float("-inf")]
        areas = [self._events[evnum].calcEventRuptureArea() for evnum in evnums]
        times = [self._events[evnum].getEventYear() for evnum in evnums]
        slips = [self._events[evnum].calcMeanSlip() for evnum in evnums]
        triggers = [self._events[evnum].getEventTrigger() for evnum in evnums]
        trigger_fault_names = [geometry.model.section( geometry.model.element(triggerID).section_id() ).name() for triggerID in triggers]
        if min(slips) > 1e-4:
            sys.stdout.write("==============================================================================\n")
            sys.stdout.write("evid\tyear\t\tmag\tarea[km^2]\tslip[m]\ttrigger\ttrigger fault\n")
            sys.stdout.write("------------------------------------------------------------------------------\n")
            for k in range(len(evnums)):
                sys.stdout.write("{}\t{:>.1f}\t\t{:>.3f}\t{:>.4f}\t{:>.4f}\t{}\t{}\n".format(evnums[k],times[k],mags[k],areas[k]*pow(10,-6),slips[k],triggers[k], trigger_fault_names[k]))
            sys.stdout.write("------------------------------------------------------------------------------\n")
        else:
            sys.stdout.write("==============================================================================\n")
            sys.stdout.write("evid\tyear\t\tmag\tarea[km^2]\tslip[m]\t\ttrigger\ttrigger fault\n")
            sys.stdout.write("------------------------------------------------------------------------------\n")
            for k in range(len(evnums)):
                sys.stdout.write("{}\t{:>.1f}\t\t{:>.3f}\t{:>.4f}\t{:>.4e}\t{}\t{}\n".format(evnums[k],times[k],mags[k],areas[k]*pow(10,-6),slips[k],triggers[k], trigger_fault_names[k]))
            sys.stdout.write("------------------------------------------------------------------------------\n")
            
    def largest_event_summary(self, num_events, geometry):
        evnums = self.get_ids_largest_events(num_events)
        self.event_summary(evnums, geometry)
    
    def event_initial_shear_stresses(self):
        return [self._events[evnum].getShearStressInit() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]

    def event_final_shear_stresses(self):
        return [self._events[evnum].getShearStressFinal() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]        
                        
    def event_initial_normal_stresses(self):
        return [self._events[evnum].getNormalStressInit() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]

    def event_final_normal_stresses(self):
        return [self._events[evnum].getNormalStressFinal() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())]  
        
    def number_of_sweep_records(self):
        return [self._events[evnum].getNumRecordedSweeps() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())] 

    def get_num_sweep_records(self, evnum):
        return self._events[evnum].getNumRecordedSweeps()
        
    def number_of_sweeps(self):
        return [self._events[evnum].getMaxSweepNum() for evnum in self._filtered_events if not np.isnan(self._events[evnum].getMagnitude())] 

    def get_num_sweeps(self, evnum):
        return self._events[evnum].getMaxSweepNum()

class AllSweeps:
    # A class for reading/analyzing data from all sweeps
    def __init__(self, sim_file, block_ids=None):
        self.sweeps = read_all_sweeps_h5(sim_file, block_ids=block_ids)
        self.sweep_data = parse_all_sweeps_h5(sweeps=self.sweeps, do_print=False)
        self.block_ids = self.sweep_data['block_id'].tolist()
        sys.stdout.write("\nRead all sweeps from {}".format(sim_file))
        self.event_elements = {0:[]}
        cur_event = 0
        for rw in self.sweep_data:
            if rw['event_number'] != cur_event:
                cur_event = rw['event_number']
                self.event_elements[rw['event_number']]=[rw['block_id']]
            else:
                self.event_elements[rw['event_number']].append(rw['block_id'])
        

class Sweeps:
    # A class for reading/analyzing data from the event sweeps
    def __init__(self, sim_file, event_number=0, block_ids=None):
        self.sweeps = read_sweeps_h5(sim_file, event_number=event_number, block_ids=block_ids)
        self.sweep_data = parse_sweeps_h5(sweeps=self.sweeps, do_print=False, event_number=event_number)
        self.block_ids = self.sweep_data['block_id'].tolist()
        self.mag = read_events_h5(sim_file,event_numbers=event_number)['event_magnitude'][0]
        self.event_number = event_number
        sys.stdout.write("\nRead event {} sweeps from {}".format(event_number,sim_file))
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
        
    def event_movie(self, geometry, events, savefile, FPS=3, DPI=100):
        # Currently only works for perfectly rectangular faults
        # Currently only plotting the elements on the triggering section
        # TODO: Change the plot to the triggering fault.
        # TODO: Better element ID and DAS handling so we can plot entire faults with sections that have different depths.
        #          e.g. A section that is 5 elements deep next to a section that's 4 elements deep.
        triggerID = int(self.sweep_data[ np.where(self.sweep_data['sweep_number']==0) ]['block_id'][0])
        num_sweeps = max([sweep_num for sweep_num in self.sweep_data['sweep_number'] ])+1
        sectionID = geometry.model.element(triggerID).section_id()
        ele_length = np.sqrt(geometry.model.create_sim_element(triggerID).area())
        triggerSecElements = [id for id in geometry.model.getElementIDs() if geometry.model.element(id).section_id() == sectionID]       
        sec_name = geometry.model.section(sectionID).name()
        min_id    = triggerSecElements[0]
        magnitude = events._events[self.event_number].getMagnitude()
        mean_slip = events._events[self.event_number].calcMeanSlip()
        ele_slips = events.get_event_element_slips(self.event_number)
        max_slip = max(ele_slips.values())
        min_slip = min(ele_slips.values())
        section_length = geometry.model.section_length(sectionID)
        section_depth = abs(geometry.model.section_max_depth(sectionID))
        num_elements_down = int(round(section_depth/ele_length))
        num_elements_across = int(round(section_length/ele_length))
        assert(len(triggerSecElements) == num_elements_across*num_elements_down)
        element_grid = np.zeros((num_elements_down,num_elements_across))
        fig = plt.figure()
        ax = plt.gca()
        if min_slip > 0:
            cmap = plt.get_cmap('Reds')
            norm = mcolor.Normalize(vmin=0, vmax=max_slip)
        else: 
            cmap = plt.get_cmap('seismic')
            norm = mcolor.Normalize(vmin=-max_slip, vmax=max_slip)
        
        # Initialize movie writing stuff
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='VQ event {}'.format(self.event_number), artist='Matplotlib',comment='Testing.')
        writer = FFMpegWriter(fps=FPS, metadata=metadata)
        
        plt.xlabel("along strike")
        plt.ylabel("down dip")
        plt.title("Virtual Quake Event {}, M={:.2f}, Fault: {}\n mean slip = {:.2f}m, max slip = {:.2f}m".format(self.event_number,magnitude,sec_name,mean_slip,max_slip),fontsize=11)
        plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
        plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
        plt.figtext(0.96, 0.6, r'cumulative slip $[m]$', rotation='vertical')
        
        # Draw the arrow in the rake direction
        mean_rake = 0
        for id in triggerSecElements: mean_rake += geometry.model.element(id).rake()/len(triggerSecElements) 
        arrow_tail = np.array([0.13, 0.1])
        arrow_length = 0.08
        arrow_head = np.array([arrow_length*np.cos(mean_rake), arrow_length*np.sin(mean_rake)])
        arrow_head += arrow_tail  #vector addition
        plt.annotate("", xy=arrow_head, xytext=arrow_tail, arrowprops=dict(arrowstyle="->", lw=2.), xycoords="figure fraction")
        plt.figtext(0.03, 0.05, 'Rake Direction\n\n\n', bbox={'facecolor':'cyan', 'pad':8, 'alpha':0.3})
        
        # Colorbar
        divider = make_axes_locatable(ax)
        cbar_ax = divider.append_axes("right", size="5%",pad=0.1)
        cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)

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
                this_sweep = self.sweep_data[ np.where(self.sweep_data['sweep_number']==sweep_num) ]
                for row in this_sweep:
                    ele_id = int(row['block_id'])
                    # Only plotting the elements on the triggering fault
                    if geometry.model.element(ele_id).section_id() == sectionID:
                        grid_row = int((ele_id-min_id)%num_elements_down)
                        grid_col = int((ele_id-min_id)/num_elements_down)
                        element_grid[grid_row,grid_col] += row['block_slip']
                    else:
                        sys.stdout.write("\nElement {} involved but not on triggering fault.".format(ele_id))
                # Update the colors
                this_plot.set_data(element_grid)
                # Time stamp
                plt.figtext(0.03, 0.9, 'Sweep: {:03d}'.format(sweep_num), bbox={'facecolor':'yellow', 'pad':8})
                writer.grab_frame()
        sys.stdout.write("\n>> Movie saved to {}\n".format(savefile))

    


class SpaceTimePlot:
    def __init__(self, geometry=None, min_year=None, max_year=None, event_file=None, trigger_fault=None, title=None):
        if min_year is None or max_year is None: 
            raise BaseException("\nMust specify --min_year and --max_year")
        elif geometry is None:
            raise BaseException("\nMust specify --model_file")
        elif event_file is None:
            raise BaseException("\nMust specify --event_file")
        elif trigger_fault is None :
            raise BaseException("\nMust specify a single fault id (example fault #N) with --use_faults N")
        else:
            self.trigger_fault = trigger_fault
            self.num_events = len(events[0]._filtered_events)
            self.event_ids = events[0]._filtered_events;
            self.geometry = geometry
            self.title = title
            self.min_year = min_year
            self.max_year = max_year
            if args.min_magnitude is None: self.min_mag = min(events[0].event_magnitudes())
            else: self.min_mag = args.min_magnitude
            self.lower_min_mag = min(events[0].event_magnitudes())*.85
            if args.max_magnitude is None: self.max_mag = max(events[0].event_magnitudes())
            else: self.max_mag = args.max_magnitude
            self.mag_slope = 1.0/(self.max_mag - self.lower_min_mag)
            self.cmap = plt.get_cmap('Reds')
            
    def get_color(self, magnitude):
         return self.cmap(self.mag_slope * (magnitude - self.lower_min_mag))
            
            
    def plot(self, fig): 
        ax = plt.gca()
        ax.set_xlabel("Distance along strike [km]")
        ax.set_ylabel("Simulation time [years]")
        ax.set_title(self.title)
        global_max_das, global_min_das = 0.0, 0.0
        ax.plot([],[])
        
        sys.stdout.write("----Plotting {:d} events from simulation years {:.1f} to {:.1f} on Fault {}\n".format(self.num_events,self.min_year,self.max_year,self.trigger_fault))
        
        for id in self.event_ids:
            event_elements = events[0]._events[id].getInvolvedElements()
            event_time = events[0]._events[id].getEventYear()
            mag = events[0]._events[id].getMagnitude()
            trigger_elem = events[0]._events[id].getEventTrigger()
            sys.stdout.write("Time {:.1f}\tEvent {:d}\tM = {:.3f}\tNum. Elements {:d}\n".format(event_time,id,mag,len(event_elements)))
            for element in event_elements:
                # Only plot slips on the specified fault
                if self.geometry._elem_to_fault_map[element] == self.trigger_fault:
                    das_min = self.geometry.model.element_min_das(element)/1000.0
                    das_max = self.geometry.model.element_max_das(element)/1000.0
                    # Keep track of the total max DAS
                    global_max_das = max(global_max_das,das_max)
                    global_min_das = min(global_min_das,das_min)
                    (r,g,b,a) = self.get_color(mag)
                    ax.plot([das_min, das_max],[event_time, event_time], color=(r,g,b), zorder=1)
                    # ===== Add a black star denoting the initial failure for the triggering element ======
                    if element == trigger_elem:
                        mean_das = 0.5*(das_max+das_min)
                        ax.scatter([mean_das],[event_time], s=60, color='k', marker='*', zorder=10)
                
        fault_length = self.geometry.model.fault(self.trigger_fault).length()/1000.0
        sys.stdout.write("---Fault {:d} has length {:.2f}, and the min/max event element DAS is {}, {}\n".format(self.trigger_fault,fault_length,global_min_das,global_max_das))
        plt.xlim(0.0, global_max_das)
        plt.ylim(args.min_year-3, args.max_year+3)
        #plt.gca().invert_yaxis()
        norm = mcolor.Normalize(vmin=self.lower_min_mag, vmax=self.max_mag)
        divider = make_axes_locatable(ax)
        cbar_ax = divider.append_axes("right", size="3%",pad=0.2)
        cb = mcolorbar.ColorbarBase(cbar_ax, cmap=self.cmap, norm=norm, orientation='vertical')
        cb.set_label('Magnitude $M_W$')
        
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
    
    def save_greens(self, output_file):
        output_file = output_file.split('.png')[0]+'.txt'
        outfile = open(output_file,'w')
        # self.field[i][j] is the field value matrix
        # self.XX[i][j] is the x coordinate
        # self.YY[i][j] is the y coordinate
        sys.stdout.write("Writing greens function values to {} in the format of\nX[m]\tY[m]\tfield value\n".format(output_file))
        for i in range(self.XX.shape[0]):
            for j in range(self.XX.shape[1]):
                outfile.write("{:.6f}\t{:.6f}\t{:.6f}\n".format(self.XX[i][j],self.YY[i][j],self.field[i][j]))
        
    
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
        plt.savefig(output_file, dpi=args.dpi)
        sys.stdout.write("----Greens function plot saved: "+output_file)
        plt.clf()

class TracePlotter:
    # Plot fault traces on a map
    def __init__(self, geometry, output_file, use_sections=None):
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
        if ((self.max_lat - self.min_lat) < MIN_LAT_DIFF) or ((self.max_lon - self.min_lon) < MIN_LON_DIFF):
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
                cbar_max=None, levels=None, g0=None, wavelength=0.3):
        if g0 is None: 
            self.g0 = 9.80665
        else:
            self.g0 = g0
        self.field_type = field_type.lower()
        self.levels = levels
        self.event_id = event_id
        plot_height = 768.0
        max_map_width = 690.0
        max_map_height = 658.0
        map_res  = 'i'
        padding  = 0.5
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
            self.wavelength = wavelength
        # Read elements and slips into the SlippedElementList
        self.elements = quakelib.SlippedElementList()
        if event_id is None and event is None and element_slips is None:
            raise BaseException("\nMust specify event_id for event fields or element_slips (dictionary of slip indexed by element_id) for custom field.")
        else:
            self.element_ids = element_slips.keys()
            self.element_slips = element_slips
        if len(self.element_slips) != len(self.element_ids):
            raise BaseException("\nMust specify slip for all elements.")
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
        if ((self.max_lat - self.min_lat) < MIN_LAT_DIFF) or ((self.max_lon - self.min_lon) < MIN_LON_DIFF):
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
                    
            #-----------------------------------------------------------------------
            # Fig3 is the land/sea mask.
            #-----------------------------------------------------------------------
            fig3 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
            self.m3.ax = fig3.add_axes((0,0,1,1))
            self.m3.fillcontinents(color='#000000', lake_color='#ffffff')
        
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
        
        if self.field_type == 'displacement' or self.field_type == 'insar':
            # FIGURE 3 draw the renderer for the sea mask
            fig3.canvas.draw()
            # FIGURE 3 Get the RGBA buffer from the figure
            w,h = fig3.canvas.get_width_height()
            buf = np.fromstring ( fig3.canvas.tostring_argb(), dtype=np.uint8 )
            buf.shape = ( w, h,4 )
     
            # FIGURE 3 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
            buf = np.roll ( buf, 3, axis = 2 )
            im3 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
            mask = im3.convert('L')
        
        # Clear all three figures
        fig1.clf()
        fig2.clf()
        plt.close('all')
        gc.collect()
        if self.field_type == 'displacement' or self.field_type == 'insar': 
            fig3.clf()
            return Image.composite(im1, im2, mask)
        else:
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
                mag_ax.text(0.5, 0.5, 'm = {:0.3f}'.format(float(events[0]._events[self.event_id].getMagnitude())), fontproperties=font_bold, size=arrow_fontsize, ha='center', va='center')
            else:
                avg_slip = np.average([x[1] for x in self.element_slips.items()])
                mag_ax.text(0.5, 0.5, 'mean slip \n{:0.3f}m'.format(avg_slip), fontproperties=font_bold, size=arrow_fontsize-1, ha='center', va='center')

        # add the map image to the plot
        m4.imshow(map_image, origin='upper')
        
        # If plotting event field, get involved sections
        if self.event_id is not None:
            involved_sections = events[0].get_event_sections(self.event_id, geometry)
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
    def create_plot(self, fig, color_index, plot_type, log_y, x_data, y_data, plot_title, x_label, y_label, label):
        #fig = plt.figure()
        ax = plt.gca()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        if plot_type == "scatter":
            ax.scatter(x_data, y_data, color = STAT_COLOR_CYCLE[color_index%len(STAT_COLOR_CYCLE)], label=label, alpha=SCATTER_ALPHA, s=SCATTER_SIZE)
        elif plot_type == "line":
            ax.plot(x_data, y_data, color = STAT_COLOR_CYCLE[color_index%len(STAT_COLOR_CYCLE)])
        elif plot_type == "loglog":
            ax.loglog(x_data, y_data, color = STAT_COLOR_CYCLE[color_index%len(STAT_COLOR_CYCLE)], lw=2.0)
        elif plot_type == "hist":
            if len(x_data) > 200: BINS=100
            elif len(x_data) < 60: BINS=20
            else: BINS=100
            ax.hist(x_data, bins=BINS, color = STAT_COLOR_CYCLE[color_index%len(STAT_COLOR_CYCLE)], histtype='stepfilled', log=log_y, normed=True)
        if plot_type != "loglog": plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        #plt.savefig(filename,dpi=100)
        #sys.stdout.write("Plot saved: {}\n".format(filename))

    def multi_line_plot(self, fig, x_data, y_data, labels, linewidths, plot_title, x_label, y_label, legend_str, filename, colors=None, linestyles=None):
        #fig = plt.figure()
        ax = plt.gca()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if linestyles is None: linestyles = ["-" for each in x_data]
        fig.suptitle(plot_title, fontsize=10)
        if colors is not None:
            if not (len(x_data) == len(y_data) and len(x_data) == len(colors) and len(colors) == len(labels) and len(linewidths) == len(colors)):
                raise BaseException("\nThese lists must be the same length: x_data, y_data, colors, labels, linewidths.")
            for i in range(len(x_data)):
                ax.plot(x_data[i], y_data[i], color=colors[i], label=labels[i], linewidth=linewidths[i], ls=linestyles[i])
        else:
            if not (len(x_data) == len(y_data) and len(x_data) == len(labels) and len(linewidths) == len(y_data)):
                raise BaseException("\nThese lists must be the same length: x_data, y_data, labels, linewidths.")
            for i in range(len(x_data)):
                ax.plot(x_data[i], y_data[i], label=labels[i], linewidth=linewidths[i], ls=linestyles[i])
        #ax.legend(title=legend_str, loc='best')
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        #plt.savefig(filename,dpi=100)
        #sys.stdout.write("Plot saved: {}\n".format(filename))

    def t0_vs_dt_plot(self, fig, t0_dt_plot, wait_75, filename):
# TODO: Set fonts explicitly
        t0_dt_main_line_color   = '#000000'
        t0_dt_sub_line_color    = '#737373'
        t0_dt_main_line_width   = 2.0
        t0_dt_sub_line_width    = 1.0
        t0_dt_range_color       = plt.get_cmap('autumn')(0.99)
        years_since_line_color  = 'blue'
        legend_loc              = 'best'
        #fig = plt.figure()
        ax = plt.gca()
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
        #ax.legend(title='event prob.', loc=legend_loc, handlelength=5)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        #plt.savefig(filename,dpi=100)
        #sys.stdout.write("Plot saved: {}\n".format(filename))

    def scatter_and_errorbar(self, fig, log_y, x_data, y_data, err_x, err_y, y_error, err_label, plot_title, x_label, y_label, label, add_x = None, add_y = None, add_label = None):
        #fig = plt.figure()
        ax = plt.gca()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        ax.scatter(x_data, y_data, label=label, alpha=SCATTER_ALPHA, color=STAT_COLOR_CYCLE[0], s=SCATTER_SIZE)
        ax.errorbar(err_x, err_y, yerr = y_error, label=err_label, ecolor='r', color='r')
        if add_x is not None:
            if log_y: ax.semilogy(add_x, add_y, label = add_label, c = 'r')
            if not log_y: ax.plot(add_x, add_y, label = add_label, c = 'r')
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        #ax.legend(loc = "best")
        #plt.savefig(filename,dpi=100)
        #sys.stdout.write("Plot saved: {}\n".format(filename))
        
    def scatter_and_error_polygon(self, fig, log_y, x_data, y_data, err_x, err_y, y_error, err_label, plot_title, x_label, y_label, label, add_x = None, add_y = None, add_label = None):
        #fig = plt.figure()
        ax = plt.gca()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        y_low = list(np.array(err_y)-np.array(y_error[0]))
        y_hi = list(np.array(err_y)+np.array(y_error[1]))
        # ===== UCERF3 rates with polygon, vertices are rate bin centers =======
        """
        x_vals = list(err_x)+list(reversed(err_x))
        y_vals = y_hi+list(reversed(y_low))
        assert(len(y_vals)==len(x_vals))
        x_vals = x_vals+[x_vals[0]]
        y_vals = y_vals+[y_vals[0]]
        xy_vals = np.array(zip(x_vals, y_vals))
        ax.add_patch(mpatches.Polygon(xy_vals, label='{} 95% confidence bounds'.format(err_label), alpha=0.55, facecolor='r', zorder=1))
        """
        
        # ===== UCERF3 rates with boxes =======
        x_bounds = [[x0-.25,x0+.25] for x0 in err_x]
        y_bounds = zip(y_low,y_hi)
        x0_y0 = [[x_bounds[i][0], y_bounds[i][0]] for i in range(len(x_bounds))]
        widths = [x_bounds[i][1] - x_bounds[i][0]for i in range(len(x_bounds))]
        heigths = [y_bounds[i][1] - y_bounds[i][0]for i in range(len(y_bounds))]
        for i in range(len(x_bounds)):
            if i==0: ax.add_patch(mpatches.Rectangle(x0_y0[i], widths[i], heigths[i], alpha=0.4, facecolor='r', label='{} 95% confidence'.format(err_label)))
            else: ax.add_patch(mpatches.Rectangle(x0_y0[i], widths[i], heigths[i], alpha=0.4, facecolor='r'))
        
        ax.scatter(x_data, y_data, label=label, alpha=SCATTER_ALPHA, color=STAT_COLOR_CYCLE[0], s=SCATTER_SIZE, zorder=10)
        if add_x is not None:
            if log_y: ax.semilogy(add_x, add_y, label = add_label, c = 'r', zorder=5)
            if not log_y: ax.plot(add_x, add_y, label = add_label, c = 'r', zorder=5)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        #ax.legend(loc = "best")
        #plt.savefig(filename,dpi=100)
        #sys.stdout.write("Plot saved: {}\n".format(filename))

    def scatter_and_line(self, fig, color_index, log_y, x_data, y_data, line_x, line_y, line_label, plot_title, x_label, y_label, label, legend_loc ='upper left'):
        #fig = plt.figure()
        ax = plt.gca()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        ax.scatter(x_data, y_data, label=label, color = STAT_COLOR_CYCLE[color_index%len(STAT_COLOR_CYCLE)], alpha=SCATTER_ALPHA, s=SCATTER_SIZE)
        if line_x is not None and line_y is not None:
            ax.plot(line_x, line_y, label = line_label, ls='-', color = 'r', lw=3.0)
            #ax.legend(loc = legend_loc)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        
        if args.zoom: plt.ylim(-5,5)
        
        #plt.savefig(filename,dpi=100)
        #sys.stdout.write("Plot saved: {}\n".format(filename))
        
    def scatter_and_multiline(self, fig, log_y, x_data, y_data, lines_x, lines_y, line_labels, line_widths, line_styles, colors, plot_title, x_label, y_label, label, legend_loc='upper left'):
        #fig = plt.figure()
        ax = plt.gca()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y: ax.set_yscale('log')
        ax.scatter(x_data, y_data, label=label, s=SCATTER_SIZE)
        for i in range(len(lines_x)):
            ax.plot(lines_x[i], lines_y[i], label = line_labels[i], ls=line_styles[i], lw=line_widths[i], c = colors[i])
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        #ax.legend(loc = legend_loc)

        y_label_words = [s.lower() for s in y_label.split(" ")]
        if "slip" in y_label_words and min(y_data) > 0.95e-2 and max(y_data) < 1.05e1: plt.ylim(1e-2,1e1)
        if "area" in y_label_words and max(y_data) < 2e4 and max(y_data) < 1.05e4: plt.ylim(1,1e4)
        
        #plt.savefig(filename,dpi=100)
        #sys.stdout.write("Plot saved: {}\n".format(filename))

class MagnitudeRuptureAreaPlot(BasePlotter):
    def plot(self, fig, color_index, events, filename, WC94=False, leonard=False, label=None):
        ra_list = events.event_rupture_areas()
        mag_list = events.event_magnitudes()
        ra_renorm_list = [quakelib.Conversion().sqm2sqkm(ra) for ra in ra_list]
        min_mag, max_mag = min(mag_list), max(mag_list)
        PLOT_TITLE = events.plot_str()
        if args.no_titles: PLOT_TITLE = " "
        if label is None: label = filename
        if WC94 and not leonard and color_index == 0:
            scale_x, scale_y = Distributions().wells_coppersmith('area')
            scale_label = "Wells & Coppersmith 1994"
            full_x, full_y = Distributions().wells_coppersmith('area', min_mag=min_mag, max_mag=max_mag)
            lines_x = [scale_x, full_x]
            lines_y = [scale_y, full_y]
            line_labels = [scale_label, None]
            line_widths = [2.0, 1.0]
            line_styles = ['-', '--']
            colors = ['k', 'k']
            self.scatter_and_multiline(fig, True, mag_list, ra_renorm_list, lines_x, lines_y, line_labels, line_widths, line_styles, colors, PLOT_TITLE,   "Magnitude", "Rupture Area (square km)", label)
        elif leonard and not WC94 and color_index == 0:
            scale_label = "Leonard 2010"
            full_x, full_y = Distributions().leonard_2010('area', min_mag=min_mag, max_mag=max_mag)
            lines_x = full_x
            lines_y = full_y
            line_labels = scale_label
            self.scatter_and_line(fig, color_index, True, mag_list, ra_renorm_list, lines_x, lines_y, line_labels, PLOT_TITLE,   "Magnitude", "Rupture Area (square km)", label)
        elif leonard and WC94 and color_index == 0:
            wc_x, wc_y = Distributions().wells_coppersmith('area', min_mag=min_mag, max_mag=max_mag)
            wc_label = "Wells & Coppersmith 1994"
            leo_x, leo_y = Distributions().leonard_2010('area', min_mag=min_mag, max_mag=max_mag)
            leo_label = "Leonard 2010"
            lines_x = [wc_x, leo_x]
            lines_y = [wc_y, leo_y]
            line_labels = [wc_label, leo_label]
            line_widths = [1.0, 1.0]
            line_styles = ['-', '-']
            colors = ['k', 'r']
            self.scatter_and_multiline(fig, True, mag_list, ra_renorm_list, lines_x, lines_y, line_labels, line_widths, line_styles, colors, PLOT_TITLE,   "Magnitude", "Rupture Area (square km)", label)
        else:
            self.create_plot(fig, color_index, "scatter", True, mag_list, ra_renorm_list, PLOT_TITLE, "Magnitude", "Rupture Area (square km)", label)

class MagnitudeMeanSlipPlot(BasePlotter):
    def plot(self, fig, color_index, events, filename, WC94=False, leonard=False, label=None):
    # Color index is an index for the event_file number, 0 is the first file, 1 is the second file
        slip_list = events.event_mean_slip()
        mag_list = events.event_magnitudes()
        min_mag, max_mag = min(mag_list), max(mag_list)
        PLOT_TITLE = events.plot_str()
        if args.no_titles: PLOT_TITLE = " "
        if label is None: label = filename
        if WC94 and not leonard and color_index == 0:
            scale_x, scale_y = Distributions().wells_coppersmith('slip')
            scale_label = "Wells & Coppersmith 1994"
            full_x, full_y = Distributions().wells_coppersmith('slip', min_mag=min_mag, max_mag=max_mag)
            lines_x = [scale_x, full_x]
            lines_y = [scale_y, full_y]
            line_labels = [scale_label, None]
            line_widths = [2.0, 1.0]
            line_styles = ['-', '--']
            colors = ['k', 'k']
            self.scatter_and_multiline(fig, True, mag_list, slip_list, lines_x, lines_y, line_labels, line_widths, line_styles, colors,PLOT_TITLE,   "Magnitude", "Mean Slip (meters)", label)
        elif leonard and not WC94 and color_index == 0:
            scale_label = "Leonard 2010"
            full_x, full_y = Distributions().leonard_2010('slip', min_mag=min_mag, max_mag=max_mag)
            lines_x = full_x
            lines_y = full_y
            line_labels = scale_label
            self.scatter_and_line(fig, color_index, True, mag_list, slip_list, lines_x, lines_y, line_labels, PLOT_TITLE,   "Magnitude", "Mean Slip (meters)", label)
        elif leonard and WC94 and color_index == 0:
            wc_x, wc_y = Distributions().wells_coppersmith('slip', min_mag=min_mag, max_mag=max_mag)
            wc_label = "Wells & Coppersmith 1994"
            leo_x, leo_y = Distributions().leonard_2010('slip', min_mag=min_mag, max_mag=max_mag)
            leo_label = "Leonard 2010"
            lines_x = [wc_x, leo_x]
            lines_y = [wc_y, leo_y]
            line_labels = [wc_label, leo_label]
            line_widths = [1.0, 1.0]
            line_styles = ['-', '-']
            colors = ['k', 'r']
            self.scatter_and_multiline(fig, True, mag_list, slip_list, lines_x, lines_y, line_labels, line_widths, line_styles, colors, PLOT_TITLE,   "Magnitude", "Mean Slip (meters)", label)
        else:
            self.create_plot(fig, color_index, "scatter", True, mag_list, slip_list, PLOT_TITLE, "Magnitude", "Mean Slip (meters)", label)

class FrequencyMagnitudePlot(BasePlotter):
    def plot(self, fig, color_index, events, filename, UCERF2 = False, UCERF3 = False, label=None):
        PLOT_TITLE = events.plot_str()
        if args.no_titles: PLOT_TITLE = " "
        # California observed seismicity rates and errorbars (UCERF2)
        x_UCERF2 = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
        y_UCERF2 = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]
        y_error_UCERF2 = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
        # California observed seismicity rates and errorbars (UCERF3, uses years 1932-2011)
        # From table L12 in Appendix L of UCERF3 Time-Independent, K.R. Felzer
        x_UCERF3 = [5.25, 5.75, 6.25, 6.75, 7.25, 7.75]
        y_UCERF3 = [4.0, 1.4, 0.45, 0.2, 0.0625, 0.0125]
        y_error_UCERF3 = [[0.4, 0.3, 0.09, 0.08, .0375, .0112],[1.5, 0.3, .14, .12, .0855, .0563]]
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
        #if b1 and color_index == 0:
        #    add_x = np.linspace(min(freq_x),max(freq_x),10)
        #    fit_point = freq_x[(np.abs(np.array(freq_x)-MIN_FIT_MAG)).argmin()]
        #    add_y = 10**(math.log(fit_point,10)+freq_x[0]-add_x)
        #    add_label = "b==1"
        if label is None: label = filename
        if UCERF2 and color_index == 0:
            self.scatter_and_error_polygon(fig, True, freq_x, freq_y, x_UCERF2, y_UCERF2, y_error_UCERF2, "UCERF2", PLOT_TITLE, "Magnitude (M)", "cumulative number of events per year with mag > M", label, add_x=add_x, add_y=add_y, add_label=add_label)
        elif UCERF3 and color_index == 0:
            self.scatter_and_error_polygon(fig, True, freq_x, freq_y, x_UCERF3, y_UCERF3, y_error_UCERF3, "UCERF3", PLOT_TITLE, "Magnitude (M)", "cumulative number of events per year with mag > M", label, add_x=add_x, add_y=add_y, add_label=add_label)
        elif not UCERF2 and not UCERF3 and color_index == 0:
            self.scatter_and_line(fig, color_index, True, freq_x, freq_y, add_x, add_y, add_label, PLOT_TITLE, "Magnitude (M)", "cumulative number of events per year with mag > M", label)
        else:
            self.create_plot(fig, color_index, "scatter", True, freq_x, freq_y, PLOT_TITLE, "Magnitude (M)", "cumulative number of events per year with mag > M", label)

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
            sys.stdout.write(stress_histories[element])
        #self.create_plot("scatter", True, mag_vals, mag_norm, events.plot_str(), "Shear Stress", "Year")
        
class DiagnosticPlot(BasePlotter):
    def plot_shear_stress_changes(self, fig, color_index, events, filename):
        shear_init = np.array(events.event_initial_shear_stresses())
        shear_final = np.array(events.event_final_shear_stresses())
        years = events.event_years()
        stress_changes = (shear_final-shear_init)/shear_init
        # Generate the binned averages too
        # Don't bin it if we're looking at a single fault model with a single repeating earthquake
        if len(np.unique(np.array(stress_changes))) > 1:
            x_ave, y_ave = calculate_averages(years,stress_changes,log_bin=False,num_bins=20)
            ave_label = "binned average"
        else:
            x_ave, y_ave = None, None
            ave_label = ""
        self.scatter_and_line(fig, color_index, False, years, stress_changes, x_ave, y_ave, ave_label, "Event shear stress changes", "simulation time [years]", "fractional change", filename)
        
    def plot_normal_stress_changes(self, fig, color_index, events, filename):
        normal_init = np.array(events.event_initial_normal_stresses())
        normal_final = np.array(events.event_final_normal_stresses())
        years = events.event_years()
        stress_changes = (normal_final-normal_init)/normal_init
        # Generate the binned averages too
        # Don't bin it if we're looking at a single fault model with a single repeating earthquake
        if len(np.unique(np.array(stress_changes))) > 1:
            x_ave, y_ave = calculate_averages(years,stress_changes,log_bin=False,num_bins=20)
            ave_label = "binned average"
        else:
            x_ave, y_ave = None, None
            ave_label = ""
        self.scatter_and_line(fig, color_index, False, years, stress_changes, x_ave, y_ave, ave_label, "Event normal stress changes", "simulation time [years]", "fractional change", filename)
        
    def plot_shear_stress_changes_vs_magnitude(self, fig, color_index, events, filename):
        shear_init = np.array(events.event_initial_shear_stresses())
        shear_final = np.array(events.event_final_shear_stresses())
        mags = events.event_magnitudes()
        stress_changes = (shear_final-shear_init)/shear_init
        # Generate the binned averages too
        # Don't bin it if we're looking at a single fault model with a single repeating earthquake
        if len(np.unique(np.array(stress_changes))) > 1:
            x_ave, y_ave = calculate_averages(mags,stress_changes,log_bin=False,num_bins=20)
            ave_label = "binned average"
        else:
            x_ave, y_ave = None, None
            ave_label = ""
        self.scatter_and_line(fig, color_index, False, mags, stress_changes, x_ave, y_ave, ave_label, "Event shear stress changes", "Magnitude", "fractional change", filename)
        
    def plot_number_of_sweeps(self, fig, color_index, events, filename):
        num_sweeps = np.array(events.number_of_sweeps())
        years = events.event_years()
        # Generate the binned averages too, if we have more than 1 distinct value for the number of sweeps.
        # i.e. don't bin it if we're looking at a single fault model with a single repeating earthquake
        if len(np.unique(np.array(num_sweeps))) > 1:
            x_ave, y_ave = calculate_averages(years,num_sweeps,log_bin=False,num_bins=20)
            ave_label = "binned average"
        else:
            x_ave, y_ave = None, None
            ave_label = ""
        self.scatter_and_line(fig, color_index, True, years, num_sweeps, x_ave, y_ave, ave_label, " ", "simulation time [years]", "number of event sweeps", filename)
        
    def plot_mean_slip(self, fig, color_index, events, filename):
        slips = np.array(events.event_mean_slip())
        years = events.event_years()
        # Generate the binned averages too
        # Don't bin it if we're looking at a single fault model with a single repeating earthquake
        if len(np.unique(np.array(slips))) > 1:
            x_ave, y_ave = calculate_averages(years,slips,log_bin=False,num_bins=20)
            ave_label = "binned average"
        else:
            x_ave, y_ave = None, None
            ave_label = ""
        self.scatter_and_line(fig, color_index, True, years, slips, x_ave, y_ave, ave_label, " ", "simulation time [years]", "event mean slip [m]", filename)

class ProbabilityPlot(BasePlotter):
    def plot_p_of_t(self, fig, events, filename, fit_weibull):
        PLOT_TITLE = events.plot_str()
        if args.no_titles: PLOT_TITLE = " "
        # Cumulative probability P(t) as a function of interevent time t
        intervals = np.array(events.interevent_times())
        prob = {}
        weibull = {}
        prob['x'] = np.sort(intervals)
        prob['y'] = np.arange(float(intervals.size))/float(intervals.size)
        
        if fit_weibull:
            mean_t0 = np.mean(intervals)
            t0_to_eval = list(np.linspace(0, int(intervals.max()), num=len(intervals)))
            fit_beta,fit_tau = fit_to_weibull(prob['x'], prob['y'], 1.0, mean_t0)
            sys.stdout.write("FITTED BETA = {:.3f}\nFITTED TAU = {:.3f}\n".format(fit_beta,fit_tau))
            weibull     = {'x':[],'y':[]}
            for dt in range(int(intervals.max())):
                weibull['x'].append(dt)
                weibull_t0_dt = Distributions().weibull(weibull['x'][-1],fit_beta,fit_tau)
                weibull['y'].append(weibull_t0_dt)          
        
        self.create_plot(fig, 0, "line", False, prob['x'], prob['y'], PLOT_TITLE,"t [years]", "P(t)", filename)
        if fit_weibull:
            fig.gca().plot(weibull['x'],weibull['y'],c='r',ls='--',label=r'$\beta = {:.3f}, \ \ \tau = {:.1f}$'.format(fit_beta,fit_tau))

    def plot_conditional_fixed_dt(self, fig, events, filename, fixed_dt=5.0):
        PLOT_TITLE = events.plot_str()
        if args.no_titles: PLOT_TITLE = " "
        # P(t0 + dt, t0) vs. t0 for fixed dt
        intervals = np.array(events.interevent_times())
        prob_dt = {'x':[],'y':[]}
        prob = {}
        prob['x'] = np.sort(intervals)
        prob['y'] = np.arange(float(intervals.size))/float(intervals.size)
        eq_guaranteed = prob['x'][np.where(prob['y']>0.995)[0][0]]
        sys.stdout.write("EQ probability reaches 99.5% after t={} years".format(eq_guaranteed))
        
        t0_to_eval = np.arange(0.0,int(eq_guaranteed)+1,0.2)
        for t0 in t0_to_eval:
            int_t0_dt = intervals[np.where( intervals > t0+fixed_dt)]
            int_t0 = intervals[np.where( intervals > t0)]
            if int_t0.size != 0:
                prob_dt['x'].append(t0)
                prob_dt['y'].append(1.0 - float(int_t0_dt.size)/float(int_t0.size))
        self.create_plot(fig, 0, "line", False, prob_dt['x'], prob_dt['y'], PLOT_TITLE, "t0 [years]", "P(t0 + {:d}, t0)".format(int(fixed_dt)), filename)

    def plot_p_of_t_multi(self, fig, events, filename, beta=None, tau=None, num_t0=4, numPoints=200, fitWeibull=False):
        PLOT_TITLE = events.plot_str()
        if args.no_titles: PLOT_TITLE = " "
        # Cumulative conditional probability P(t,t0) as a function of
        # interevent time t, computed for multiple t0. Beta/Tau are Weibull parameters
        line_colormap = plt.get_cmap('autumn')
        intervals = np.array(events.interevent_times())
        conditional = {}
        weibull = {}
        max_t0 = int(intervals.max())
        mean_t0 = np.mean(intervals)
        median_t0 = np.median(intervals)
        sys.stdout.write("Mean recurrence interval: {:.2f}\n".format(mean_t0))
        sys.stdout.write("Median recurrence interval: {:.2f}\n".format(median_t0))
        t0_to_eval = list(np.linspace(0, max_t0, num=len(intervals)))
        t0_to_plot = [int(0), int(round(0.6*mean_t0,-1)), int(round(1.25*mean_t0,-1))]
        if t0_to_plot[1]==t0_to_plot[2]: t0_to_plot[2]+=10
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
                    if beta is not None and tau is not None and not fitWeibull:
                        weibull[t0]['x'].append(t0+dt)
                        weibull_t0_dt = Distributions().cond_weibull(weibull[t0]['x'][-1],t0,beta,tau)
                        weibull[t0]['y'].append(weibull_t0_dt)
            else:
                conditional[t0] = None
                weibull[t0] = None
                
        if fitWeibull:
            prob = {}
            prob['x'] = np.sort(intervals)
            prob['y'] = np.arange(float(intervals.size))/float(intervals.size)
            fit_beta,fit_tau = fit_to_weibull(prob['x'], prob['y'], 1.0, mean_t0)
            sys.stdout.write("FITTED BETA = {:.3f}\nFITTED TAU = {:.3f}\n".format(fit_beta,fit_tau))
            for t0 in t0_to_eval:
                int_t0 = intervals[np.where( intervals > t0)]
                if int_t0.size != 0:
                    weibull[t0]     = {'x':[],'y':[]}
                    for dt in range(max_t0-int(t0)):
                        weibull[t0]['x'].append(t0+dt)
                        weibull_t0_dt = Distributions().cond_weibull(weibull[t0]['x'][-1],t0,fit_beta,fit_tau)
                        weibull[t0]['y'].append(weibull_t0_dt)
                else:
                    weibull[t0] = None
                
        x_data_prob = [conditional[t0]['x'] for t0 in t0_to_plot]
        y_data_prob = [conditional[t0]['y'] for t0 in t0_to_plot]
        t0_colors   = [line_colormap(float(t0*.8)/t0_to_plot.max()) for t0 in t0_to_plot]
        prob_lw     = [2 for t0 in t0_to_plot]
        if fitWeibull or (beta is not None and tau is not None):
            x_data_weib = [weibull[t0]['x'] for t0 in t0_to_plot]
            y_data_weib = [weibull[t0]['y'] for t0 in t0_to_plot]
            weib_colors = ['k' for t0 in t0_to_plot]
            if fitWeibull: weib_labels = [r'$\beta = {:.3f}, \ \ \tau = {:.1f}$'.format(fit_beta,fit_tau), None, None]
            else: weib_labels = [r'$\beta = {:.3f}, \tau = {:.1f}$'.format(beta,tau), None, None]
            weib_lw     = [1 for t0 in t0_to_plot]
            # List concatenation, not addition
            colors = t0_colors + weib_colors
            x_data = x_data_prob + x_data_weib
            y_data = y_data_prob + y_data_weib
            labels = ["$t_0 =$ {:d}".format(t0) for t0 in t0_to_plot] + weib_labels
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
        self.multi_line_plot(fig, x_data, y_data, labels, linewidths, PLOT_TITLE, x_lab, y_lab, legend_string, filename, colors=colors)

    def plot_dt_vs_t0(self, fig, events, filename, years_since=None):
        # Plot the waiting times corresponding to 25/50/75% conditional probabilities
        # as a function of t0 (time since last earthquake on the selected faults).
        # years_since is the number of years since the last observed (real) earthquake
        # on the selected faults.
        intervals = np.array(events.interevent_times())
        sys.stdout.write("Found {:d} recurrence times...\n".format(len(intervals)))
        conditional = {}
        wait_75 = None
        max_t0 = int(intervals.max())
        # t0_to_eval used to evaluate waiting times with 25/50/75% probability given t0=years_since
        t0_to_eval = np.arange(0, max_t0, 1.0)
        t0_to_plot = np.linspace(0, int(max_t0), num=min(len(intervals),len(t0_to_eval)-5))
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
        self.t0_vs_dt_plot(fig, t0_dt_plot, wait_75, filename)
        
    def print_prob_table(self, t0_list, dt_list, mag_list, events):
        assert(len(t0_list)==len(mag_list))
        dt_vals       = dt_list # units = years
        Mag_vals      = mag_list
        YEAR_INTERVAL = 0.1
        probabilities = {} # Probabilities for the final table
        num_events    = []
        
        for k,MAG in enumerate(Mag_vals):
            # --------- Set the first event filter for the lower magnitude bound ------------
            new_event_filter = [MagFilter(min_mag=MAG, max_mag=None)]
            sys.stdout.write("Events before this MAG > {} filter {}\n".format(MAG,len(events._filtered_events)))    
            events.set_filters(new_event_filter)
            sys.stdout.write("Events after this MAG > {} filter {}\n".format(MAG,len(events._filtered_events)))  
            num_events.append(len(events._filtered_events))
            
            # == Compute conditional probs for M > MAG earthquakes in the time ranges "dt_vals" ==
            intervals       = np.array(events.interevent_times())
            num_intervals   = len(intervals)
            conditional     = {}
            max_t0          = intervals.max() 
            #t0_vals         = list(t0_list)+list([max_t0])
            #print("max_t0 {}, t0_vals {}".format(max_t0, t0_vals))
            #if (max_t0 != max(t0_vals)):
            #    raise BaseException("\nA specified t0 value exceeds all recurrence intervals for M>{:.1f}, consider changing Mag_vals to a smaller maximum value.\n".format(MAG))
            #else:
            t0_to_eval  = list(np.linspace(0, max_t0+YEAR_INTERVAL, num=num_intervals))
            t0_to_eval  += list(t0_list)
            for t0 in sorted(t0_to_eval):
                t0 = round(t0,1)
                int_t0 = intervals[np.where( intervals > t0 )]
                if int_t0.size != 0:
                    conditional[t0] = {'x':[],'y':[]}
                    for dt in np.arange(0.0,max_t0-int(t0), YEAR_INTERVAL):
                        int_t0_dt = intervals[np.where( intervals > t0+dt )]
                        conditional[t0]['x'].append(t0+dt)
                        prob_t0_dt    = 1.0 - float(int_t0_dt.size)/float(int_t0.size)
                        conditional[t0]['y'].append(prob_t0_dt)
            
            prob_in_t0_dt = []
            
            for i in range(len(dt_vals)):
                t0_rounded = round(t0_list[k],1)
                prob_in_t0_dt.append(conditional[ t0_rounded ]['y'][ int(dt_vals[i]/YEAR_INTERVAL) ] - conditional[ t0_rounded ]['y'][0])
            
            probabilities[MAG] = prob_in_t0_dt
        
        column_headers = "\n\t\t\t\t"
        for dt_value in dt_vals:
            column_headers += "\t{:.1f}yr ".format(dt_value)
        column_headers += "\n"
        
        sys.stdout.write(column_headers)
        for i in range(len(Mag_vals)):
            MAG = Mag_vals[i]
            table_line = "(N={:6d})\tM > {}\t  t0={:.1f}".format(num_events[i],MAG,t0_list[i])
            for j in range(len(dt_vals)):
                table_line += "\t{:.2f}".format(probabilities[MAG][j])
            table_line += "\n"
            sys.stdout.write(table_line)
        sys.stdout.write("\n")
        

class Distributions:
    def weibull(self, X, beta, tau):
        # Return the Weibull distribution at a point
        return 1-np.exp( -(np.array(X)/float(tau))**beta)

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
            raise BaseException("\nMust specify rupture area or mean slip")
        x_data = np.linspace(min_mag, max_mag, num=num)
        y_data = np.array([pow(10,a+b*m) for m in x_data])
        return x_data, y_data
        
    def leonard_2010(self, type, min_mag, max_mag, num=5):
        # Return empirical scaling relations from Mark Leonard 2010 BSSA
        if type.lower() == 'area':
            a = -4.0 
            b =  1.0
        elif type.lower() == 'slip':
            a = -3.417
            b =  0.499
        else:
            raise BaseException("\nMust specify rupture area or mean slip")
        if min_mag < max_mag:
            x_data = np.linspace(min_mag, max_mag, num=num)
        else:
            x_data = np.linspace(4, 8, num=num)
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
    parser.add_argument('--event_file', required=False, type=str, nargs='+',
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
    parser.add_argument('--combine_file', required=False,
            help="Name of events hdf5 file to combine with event_file.")
    parser.add_argument('--label', required=False, type=str, nargs='+',
            help="Custom label to use for plot legends, specify one per event file.")

    # Event filtering arguments
    parser.add_argument('--min_magnitude', type=float, required=False,
            help="Minimum magnitude of events to process.")
    parser.add_argument('--max_magnitude', type=float, required=False,
            help="Maximum magnitude of events to process.")
    parser.add_argument('--min_year', type=float, required=False,
            help="Minimum year of events to process.")
    parser.add_argument('--max_year', type=float, required=False,
            help="Maximum year of events to process.")
    parser.add_argument('--full_time_range', required=False, action='store_true',
            help="A speedup for time series plots. If you want a time series for the full range of the simulation, specify this parameter.")
    parser.add_argument('--min_slip', type=float, required=False,
            help="Minimum mean slip of events to process.")
    parser.add_argument('--max_slip', type=float, required=False,
            help="Maximum mean slip of events to process.")
    parser.add_argument('--min_area', type=float, required=False,
            help="Minimum rupture area of events to process (in km^2).")
    parser.add_argument('--max_area', type=float, required=False,
            help="Maximum rupture area of events to process (in km^2).")
    parser.add_argument('--min_event_num', type=float, required=False,
            help="Minimum event number of events to process.")
    parser.add_argument('--max_event_num', type=float, required=False,
            help="Maximum event number of events to process.")
    parser.add_argument('--min_num_elements', type=float, required=False,
                        help="Minimum number of elements involved in an event")
    parser.add_argument('--max_num_elements', type=float, required=False,
                        help="Maximum number of elements involved in an event")
    parser.add_argument('--involved_sections', type=int, nargs='+', required=False,
                        help="List of model sections to use (all sections used if unspecified). All earthquakes involving specified sections included.")
    parser.add_argument('--use_sections', type=int, nargs='+', required=False,
            help="List of model sections to use (all sections used if unspecified). Earthquakes will have initiated on the specfified sections.")
    parser.add_argument('--use_faults', type=int, nargs='+', required=False,
            help="List of model faults to use (all sections used if unspecified). Earthquakes will have initiated on the specfified faults.")
    parser.add_argument('--group1_ids', type=int, nargs='+', required=False,
            help="List of model faults to use. Earthquakes will have initiated on the specified faults. Must also specify --group2_ids. These subsets are used for computing time series correlations.")
    parser.add_argument('--group2_ids', type=int, nargs='+', required=False,
            help="List of model faults to use. Earthquakes will have initiated on the specified faults. Must also specify --group1_ids. These subsets are used for computing time series correlations.")
    parser.add_argument('--group1_files', type=str, nargs='+', required=False,
            help="List of files containing the fault slip time series. Must also specify --group2_files. These subsets are used for computing time series correlations.")
    parser.add_argument('--group2_files', type=str, nargs='+', required=False,
            help="List of files containing the fault slip time series. Must also specify --group1_files. These subsets are used for computing time series correlations.")

    # Statisical plotting arguments
    parser.add_argument('--plot_freq_mag', required=False, action='store_true',
            help="Generate frequency magnitude plot.")
    parser.add_argument('--UCERF2', required=False, action='store_true',
            help="Add to frequency-magnitude plot the observed rates in California from UCERF2 [Field et al. 2009].")
    parser.add_argument('--UCERF3', required=False, action='store_true',
            help="Add to frequency-magnitude plot the observed rates in California from UCERF3 [Field et al. 2014].")
    parser.add_argument('--plot_mag_rupt_area', required=False, action='store_true',
            help="Generate magnitude vs rupture area plot.")
    parser.add_argument('--plot_mag_mean_slip', required=False, action='store_true',
            help="Generate magnitude vs mean slip plot.")
    parser.add_argument('--all_stat_plots', required=False, action='store_true',
            help="Generate frequency-magnitude, magnitude vs rupture area, and magnitude vs mean slip plots.")  
    parser.add_argument('--wc94', required=False, action='store_true',
            help="Plot Wells and Coppersmith 1994 scaling relations.")
    parser.add_argument('--leonard', required=False, action='store_true',
            help="Plot Leonard 2010 scaling relations.")
    parser.add_argument('--plot_recurrence', required=False, action='store_true',
            help="Plot distribution of recurrence intervals.")

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
    parser.add_argument('--probability_table', required=False, action='store_true',
            help="Generate a table of conditional probabilities for the next large earthquakes, must specify the time since the last M>5.0, M>6.0, and M>7.0 earthquakes on the faults being simulated (use --t0 to specify these times in units of decimal years).")
    parser.add_argument('--t0', type=float, nargs='+', required=False,
            help="List of times [rounded to the nearest 0.1 years] since the last observed earthquakes of a given magnitude on the modeled faults. Must also specify the same number of --magnitudes arguments. Example --t0 4.1 12.5 35.9 and --magnitudes 5 6 6.5 for the probabilities of M>5, M>6 and M>6.5. Can request specific number of years in the future to evaluate the probabilities with --t.")
    parser.add_argument('--t', type=float, nargs='+', required=False,
            help="Number of years into the future to evaluate conditional probabilities. E.g. --t 1 5 10 will evaluate the conditional probabilities for earthquakes with magnitudes given by --magnitudes, with the last-observed earthquake of each magnitude occurring at each value of --t0. Example --t0 4.1 12.5 35.9 --magnitudes 5 6 6.5 --t 1 5 10 for the conditional probabilities of M>5, M>6 and M>6.5 at 1, 5 and 10 years from present if the last earthquakes of each magnitude were --t0 years ago.")
    parser.add_argument('--magnitudes', type=float, nargs='+', required=False,
            help="List of magnitudes to compute conditional probability table for. Must also specify the same number of --t0 arguments, where the i-th --t0 value is the time elapsed since the last observed earthquake with magnitude M = i-th --magnitudes value. Example --t0 4.1 12.5 35.9 and --magnitudes 5 6 6.5 for the probabilities of M>5, M>6 and M>6.5. Can request specific number of years in the future to evaluate the probabilities with --t.")
    parser.add_argument('--fit_weibull', action='store_true', required=False,
            help="Fit a Weibull distribution to the simulated earthquake distribution.")
            
            
    # Field plotting arguments
    parser.add_argument('--field_plot', required=False, action='store_true',
            help="Plot surface field for a specified event, e.g. gravity changes or displacements.")
    parser.add_argument('--field_type', required=False, help="Field type: gravity, dilat_gravity, displacement, insar, potential, geoid")
    parser.add_argument('--colorbar_max', required=False, type=float, help="Max unit for colorbar")
    parser.add_argument('--event_id', required=False, type=int, help="Event number for plotting event fields")
    parser.add_argument('--uniform_slip', required=False, type=float, help="Amount of slip for each element in the model_file, in meters.")
    parser.add_argument('--angles', type=float, nargs='+', required=False,
            help="Observing angles (azimuth, elevation) for InSAR or displacement plots, in degrees.")
    parser.add_argument('--wavelength', type=float, required=False,
            help="Observing wavelength for InSAR or displacement plots, in meters. Default is 0.3m")
    parser.add_argument('--levels', type=float, nargs='+', required=False,
            help="Levels for contour plot.")
    parser.add_argument('--small_model', required=False, action='store_true', help="Small fault model, used to specify map extent.")    
    parser.add_argument('--field_eval', required=False, action='store_true', help="Evaluate an event field at specified lat/lon. Must provide the file, --lld_file")
    parser.add_argument('--lld_file', required=False, help="File containing lat/lon columns to evaluate an event field.")
    
    # Greens function plotting arguments
    parser.add_argument('--greens', required=False, action='store_true', help="Plot single element Okubo Green's functions. Field type also required.") 
    parser.add_argument('--save_greens', required=False, action='store_true', help="Save single element Okubo Green's function values.")            
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
            help="Plot all diagnostic plotsall")
    parser.add_argument('--event_elements', required=False, action='store_true',
            help="Print the involved elements, must specify event id.")
    parser.add_argument('--num_sweeps', required=False, action='store_true',
            help="Plot the number of sweeps for events")
    parser.add_argument('--event_shear_stress', required=False, action='store_true',
            help="Plot shear stress changes for events")
    parser.add_argument('--event_normal_stress', required=False, action='store_true',
            help="Plot normal stress changes for events")
    parser.add_argument('--event_shear_stress_vs_magnitude', required=False, action='store_true',
            help="Plot shear stress changes for events vs magnitude")    
    parser.add_argument('--event_mean_slip', required=False, action='store_true',
            help="Plot the mean slip for events")
    parser.add_argument('--zoom', required=False, action='store_true',
            help="Force zoomed bounds on scatter and line plots")

    # Customization
    parser.add_argument('--dpi', required=False, type=float,
            help="Specify the DPI for plots that are saved.")
    parser.add_argument('--no_titles', required=False, action='store_true',
            help="Specify no titles on plots.")
    parser.add_argument('--pdf', required=False, action='store_true',
            help="Save plots as PDF instead of PNG.")
    parser.add_argument('--eps', required=False, action='store_true',
            help="Save plots as EPS instead of PNG.")
            
    # Geometry
    parser.add_argument('--slip_rates', required=False, action='store_true',
            help="Print element id and slip rate for all elements.")
    parser.add_argument('--elements', type=int, nargs='+', required=False,
            help="List of elements for filtering.")
    parser.add_argument('--slip_time_series', required=False, action='store_true',
            help="Return the slip time series for all specified --elements.")
    parser.add_argument('--standardized', required=False, action='store_true',
            help="Standardize the time series by subtracting the mean and dividing by the standard deviation.")            
    parser.add_argument('--fault_time_series', required=False, action='store_true',
            help="Return the average slip time series for all elements on a fault. To specify fault #1, use --use_faults 1.")
    parser.add_argument('--fault_group_time_series_plot', required=False, action='store_true',
            help="After using --fault_time_series to save fault time series data for each fault, make a plot of time series averaged over groups of faults by specifying the fault time series files for each group of faults with --group1_files and --group2_files.")
    parser.add_argument('--fault_group_time_series_correlate', required=False, action='store_true',
            help="After using --fault_time_series to save fault time series data for each fault, correlate the time series and make a plot by specifying the fault time series files for each group of faults with --group1_files and --group2_files.")
    parser.add_argument('--dt', required=False, type=float, help="Time step for slip rate plots, unit is decimal years.")
    parser.add_argument('--event_kml', required=False, action='store_true',
            help="Save a KML (Google Earth) file of the event elements, colored by event slip.")
    parser.add_argument('--block_area_hist', required=False, action='store_true',
            help="Save a histogram of element areas.")
    parser.add_argument('--block_length_hist', required=False, action='store_true',
            help="Save a histogram of element lengths [sqrt(area)].")
    parser.add_argument('--block_aseismic_hist', required=False, action='store_true',
            help="Save a histogram of element aseismic fraction.")
    parser.add_argument('--block_stress_drop_hist', required=False, action='store_true',
            help="Save a histogram of element stress drops.")
    parser.add_argument('--fault_length_hist', required=False, action='store_true',
            help="Save a histogram of fault lengths in the model.")  
    parser.add_argument('--fault_length_distribution', required=False, action='store_true',
            help="Save the cumulative distribution of fault lengths in the model.") 
    parser.add_argument('--reference', required=False, type=float,
            help="Reference value for numbers relative to some value.")
    parser.add_argument('--traces', required=False, action='store_true', help="Plot the fault traces from a fault model on a map.") 
    parser.add_argument('--fault_group_traces', required=False, action='store_true', help="Plot the fault traces on a map for two groups of faults, specified by two lists of fault ids using --group1_ids and --group2_ids.") 
            
    # --------- Spacetime plots -----------
    parser.add_argument('--spacetime', required=False, action='store_true',
            help="Plot a spacetime plot of ruptures for a single fault. Must specify the fault to use with --use_faults.")
            
    # Event movies
    parser.add_argument('--event_movie', required=False, action='store_true',
            help="Make a movie of a specified event, must use --event_id.")

    # Validation/testing arguments
    parser.add_argument('--validate_slip_sum', required=False,
            help="Ensure the sum of mean slip for all events is within 1 percent of the specified value.")
    parser.add_argument('--validate_mean_interevent', required=False,
            help="Ensure the mean interevent time for all events is within 2 percent of the specified value.")

    args = parser.parse_args()
    
    # ------------------------------------------------------------------------
    # Catch these errors before reading events to save unneeded computation
    if args.uniform_slip:
        if float(args.uniform_slip) < 0: raise BaseException("\nSlip must be positive")
    
    if args.field_plot:
        if args.model_file is None:
            raise BaseException("\nMust specify --model_file for field plots")
        elif args.field_type is None:
            raise BaseException("\nMust specify --field_type for field plots")
            
    if args.traces:
        if args.model_file is None:
            raise BaseException("\nMust specify --model_file for fault trace plots")

    if args.fault_group_traces:
        if args.model_file is None or args.group1_ids is None or args.group2_ids is None:
            raise BaseException("\nMust specify --model_file, --group1, and --group2 for fault group trace plots")
            
    # Check that if either beta or tau is given then the other is also given
    if (args.beta and not args.tau) or (args.tau and not args.beta):
        raise BaseException("\nMust specify both beta and tau.")
        
    # Check that field_type is one of the supported types
    if args.field_type:
        type = args.field_type.lower()
        if type != "gravity" and type != "dilat_gravity" and type != "displacement" and type != "insar" and type!= "potential" and type != "geoid":
            raise BaseException("\nField type is one of gravity, dilat_gravity, displacement, insar, potential, geoid")
            
    # Pre-check for probability table evaluation. We need the same number of --t0 arguments as --magnitudes
    if args.t0 or args.magnitudes:
        if args.magnitudes and args.t0:
            if len(args.t0) != len(args.magnitudes):
                raise BaseException("\nFor probability table, must specify the same number of --t0 and --magnitudes arguments, and --magnitudes should be listed in increasing order.")
        else:
            raise BaseException("\nFor probability table, must specify both --t0 and --magnitudes, and --magnitudes should be listed in increasing order.")
    
    # Apply default values for t if we are making probability table and they are not specified
    if args.probability_table and args.t is None:
        args.t = [1.0,5.0,10.0]
    if args.probability_table and args.t:
        args.t = [float(t) for t in args.t]
    # ------------------------------------------------------------------------
    
    if args.dpi is None:
        args.dpi = 100

     # Read the geometry model if specified
    if args.model_file:
        if args.model_file_type:
            if not os.path.isfile(args.model_file):
                raise BaseException("\nModel file does not exist: "+args.model_file)
            else:
                geometry = Geometry(model_file=args.model_file, model_file_type=args.model_file_type)
        else:
            if not os.path.isfile(args.model_file):
                raise BaseException("\nModel file does not exist: "+args.model_file)
            else:
                geometry = Geometry(model_file=args.model_file)
                             
    # Read the event and sweeps files
    if args.event_file and args.sweep_file is None and args.combine_file is None:
        # If given multiple event files
        # Currently only works for hdf5 files, time consuming to add text file support for every new feature
        events = []
        for file in args.event_file:
            # Check that all files exist
            if not os.path.isfile(file):
                raise BaseException("\nEvent file does not exist: "+file)
            else:
                events.append( Events(file, None) )
    elif args.event_file and len(args.event_file)==1 and ( args.sweep_file or args.combine_file or args.stress_file):
        if not os.path.isfile(args.event_file[0]):
            raise BaseException("\nEvent file does not exist: "+args.event_file[0])
        else:
            events = [Events(args.event_file[0], args.sweep_file, stress_file=args.stress_file, combine_file=args.combine_file, stress_index_file=args.stress_index_file)]


    # Read the stress files if specified
    if args.stress_index_file and args.stress_file:
        stress_set = quakelib.ModelStressSet()
        stress_set.read_file_ascii(args.stress_index_file, args.stress_file)
    else:
        stress_set = None
        
    if args.all_stat_plots:
        args.plot_freq_mag = True
        args.plot_mag_rupt_area = True
        args.plot_mag_mean_slip = True
        args.leonard = True
    
    # Set up filters
    event_filters = []
    if args.min_magnitude or args.max_magnitude:
        event_filters.append(MagFilter(min_mag=args.min_magnitude, max_mag=args.max_magnitude))
        sys.stdout.write("Applying magnitude filter...")

    if args.min_num_elements or args.max_num_elements:
        event_filters.append(NumElementsFilter(min_num_elements=args.min_num_elements, max_num_elements=args.max_num_elements))
        sys.stdout.write("Applying number of elements per event filter...")

    if args.min_year or args.max_year:
        ## To prevent testing all events when plotting a time series and using --min_year,--max_year as the min/max event times,
        ##### don't create a filter in this case.
        if not args.full_time_range:
            event_filters.append(YearFilter(min_year=args.min_year, max_year=args.max_year))
            sys.stdout.write("Applying event year filter...")

    # Detectability threshold, min slip 1cm
    #if args.event_file and args.min_slip is None: 
    #    args.min_slip = 0.01
    #    sys.stdout.write(" >>> Applying detectibility cut, minimum mean event slip 1cm <<< \n")
    if args.event_file and args.min_slip is not None and args.min_slip < 0:
        args.min_slip = None

    if args.min_slip or args.max_slip:
        event_filters.append(SlipFilter(min_slip=args.min_slip, max_slip=args.max_slip))
        sys.stdout.write("Applying mean slip filter...")
        
    if args.min_area or args.max_area:
        event_filters.append(AreaFilter(min_area=args.min_area, max_area=args.max_area))
        sys.stdout.write("Applying rupture area filter...")

    if args.min_event_num or args.max_event_num:
        event_filters.append(EventNumFilter(min_event_num=args.min_event_num, max_event_num=args.max_event_num))
        sys.stdout.write("Applying event number filter...")
    
    if args.involved_sections:
        if not args.model_file: raise BaseException("\nMust specify --model_file for --involved_sections to work.")
        for sec_id in args.involved_sections:
            if sec_id not in geometry._elem_to_section_map.values():
                #sys.stdout.write(geometry._elem_to_section_map)
                raise BaseException("\nSection id {} does not exist.".format(sec_id))
        event_filters.append(InvolvedSectionFilter(geometry, args.event_file[0], args.involved_sections))
        sys.stdout.write("Applying fault section filter involving only fault sections {}...".format(args.involved_sections))
    
    if args.use_sections:
        if not args.model_file: raise BaseException("\nMust specify --model_file for --use_sections to work.")
        for sec_id in args.use_sections:
            if sec_id not in geometry._elem_to_section_map.values():
                #sys.stdout.write(geometry._elem_to_section_map)
                raise BaseException("\nSection id {} does not exist.".format(sec_id))
        event_filters.append(TriggerSectionFilter(geometry, args.use_sections))
        sys.stdout.write("Applying fault section filter using only fault sections {}...".format(args.use_sections))

    if args.use_faults:
        if not args.model_file: raise BaseException("\nMust specify --model_file for --use_faults to work.")
        for fault_id in args.use_faults:
            if fault_id not in geometry._elem_to_fault_map.values():
                #sys.stdout.write(geometry._elem_to_fault_map)
                raise BaseException("\nFault id {} does not exist.".format(fault_id))
        event_filters.append(TriggerFaultFilter(geometry, args.use_faults))
        sys.stdout.write("Applying fault filter using only faults {}...".format(args.use_faults))

    if args.event_file:
        if isinstance(args.event_file, list):
            for event_set in events:
                event_set.set_filters(event_filters)
                sys.stdout.write("..Event filtering complete.\n")
        else:
            events.set_filters(event_filters)
            sys.stdout.write("..Event filtering complete.\n")
            
    # Make sure that events is a list
    if args.event_file: assert(isinstance(events, list))
    
    # Make sure that if labels are specified, the number matches the number of event files
    #if args.label: assert(len(args.label)==len(events))
    
    # Print out event summary data if requested
    if args.summary:
        if args.model_file is None: raise BaseException("\nMust specify --model_file for summary.")
        for i, event in enumerate(events):        
            sys.stdout.write("\n Event summary for: "+ args.event_file[i])
            event.largest_event_summary(args.summary, geometry)

    if args.event_elements:
        if args.event_id is None: raise BaseException("\nMust specify --event_id")
        sys.stdout.write("\nEvent {}\n".format(args.event_id))
        sys.stdout.write([each for each in events._events[args.event_id].getInvolvedElements()])

    # Generate plots
    if args.diagnostics:
        args.num_sweeps = True
        args.event_shear_stress = True
        args.event_shear_stress_vs_magnitude = True
        args.event_normal_stress = True
        args.event_mean_slip = True
        args.plot_freq_mag = True
        args.plot_mag_mean_slip = True
        args.plot_mag_rupt_area = True
        args.leonard = True
        
    if args.plot_freq_mag:
        if args.label is None: LABEL_SIZE = 8
        else: LABEL_SIZE = 10
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().event_plot(args.event_file, "freq_mag", args.min_magnitude, args.min_year, args.max_year, args.combine_file)
        for i, event_set in enumerate(events):
            num_unique_mags = len(np.unique(np.array(event_set.event_magnitudes())))
            if num_unique_mags > 5:
                if args.label: LABEL = args.label[i]
                else: LABEL = None
                FrequencyMagnitudePlot().plot(fig, i, event_set, args.event_file[i].split("events_")[-1].split("/")[-1], UCERF2=args.UCERF2, UCERF3=args.UCERF3, label=LABEL)
            else:
                sys.stdout.write("Not plotting frequency/magnitude for {}, found < 5 unique events.".format(args.event_file[i].split("events_")[-1].split("/")[-1]))
        plt.legend(loc='lower left', fontsize=LABEL_SIZE)
        if args.min_magnitude is not None and args.max_magnitude is not None:
            plt.xlim(args.min_magnitude, args.max_magnitude)
        elif args.min_magnitude is not None:
            plt.xlim(args.min_magnitude, plt.xlim()[1])
        elif args.max_magnitude is not None:
            plt.xlim(plt.xlim()[0], args.max_magnitude)
        plt.ylim(0.99e-5,4)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.plot_mag_rupt_area:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().event_plot(args.event_file, "mag_rupt_area", args.min_magnitude, args.min_year, args.max_year, args.combine_file)
        for i, event_set in enumerate(events):
            if args.label: LABEL = args.label[i]
            else: LABEL = None
            MagnitudeRuptureAreaPlot().plot(fig, i, event_set, args.event_file[i].split("events_")[-1].split("/")[-1], WC94=args.wc94, leonard=args.leonard, label=LABEL)
        if args.min_magnitude is not None and args.max_magnitude is not None:
            plt.xlim(args.min_magnitude, args.max_magnitude)
        elif args.min_magnitude is not None:
            plt.xlim(args.min_magnitude, plt.xlim()[1])
        elif args.max_magnitude is not None:
            plt.xlim(plt.xlim()[0], args.max_magnitude)
        plt.legend(loc='upper left', fontsize=10)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.plot_mag_mean_slip:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().event_plot(args.event_file, "mag_mean_slip", args.min_magnitude, args.min_year, args.max_year, args.combine_file)
        for i, event_set in enumerate(events):
            if args.label: LABEL = args.label[i]
            else: LABEL = None
            MagnitudeMeanSlipPlot().plot(fig, i, event_set, args.event_file[i].split("events_")[-1].split("/")[-1], WC94=args.wc94, leonard=args.leonard, label=LABEL)
        if args.min_magnitude is not None and args.max_magnitude is not None:
            plt.xlim(args.min_magnitude, args.max_magnitude)
        elif args.min_magnitude is not None:
            plt.xlim(args.min_magnitude, plt.xlim()[1])
        elif args.max_magnitude is not None:
            plt.xlim(plt.xlim()[0], args.max_magnitude)
        plt.legend(loc='upper left', fontsize=10)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.plot_prob_vs_t:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().event_plot(args.event_file, "prob_vs_time", args.min_magnitude, args.min_year, args.max_year, args.combine_file)
        for event_set in events:
            ProbabilityPlot().plot_p_of_t(fig, event_set, filename, args.fit_weibull)
        ax.legend(loc='best',fontsize=10)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.plot_prob_vs_t_fixed_dt:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().event_plot(args.event_file, "p_vs_t_fixed_dt", args.min_magnitude, args.min_year, args.max_year, args.combine_file)
        for event_set in events:
            ProbabilityPlot().plot_conditional_fixed_dt(fig, event_set, filename)
        #ax.legend(loc='best')
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.plot_cond_prob_vs_t:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().event_plot(args.event_file, "cond_prob_vs_t", args.min_magnitude, args.min_year, args.max_year, args.combine_file)
        if args.beta:
            for event_set in events:
                ProbabilityPlot().plot_p_of_t_multi(fig, event_set, filename, beta=args.beta, tau=args.tau, fitWeibull=args.fit_weibull)
        else:
            for event_set in events:
                ProbabilityPlot().plot_p_of_t_multi(fig, event_set, filename, fitWeibull=args.fit_weibull)
        ax.legend(loc='best', fontsize=10)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.plot_waiting_times:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().event_plot(args.event_file, "waiting_times", args.min_magnitude, args.min_year, args.max_year, args.combine_file)
        for event_set in events:
            ProbabilityPlot().plot_dt_vs_t0(fig, event_set, filename)
        ax.legend(loc='best')
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.plot_recurrence:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        times = [event_set.interevent_times() for event_set in events]
        filename = SaveFile().event_plot(args.event_file, "recurrence", args.min_magnitude, args.min_year, args.max_year, args.combine_file)
        for time in times:
            BasePlotter().create_plot(fig, 0, "hist", False, time, None, events[0].plot_str(), "interevent time [years]", "", filename)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.probability_table:
        if args.t0 is None: raise BaseException("\nMust specify time since last earthquakes with magnitude values given by --magnitudes. Use --t0 and --magnitudes, they must be the same number of arguments, and --magnitudes should be listed in increasing order.")
        else:
            for event_set in events:
                ProbabilityPlot().print_prob_table(args.t0, args.t, args.magnitudes, event_set)
    if args.field_plot:
        type = args.field_type.lower()
        if args.colorbar_max: cbar_max = args.colorbar_max
        else: cbar_max = None
        if args.levels: levels = args.levels
        else: levels = None
        filename = SaveFile().field_plot(args.model_file, type, args.uniform_slip, args.event_id, args.wavelength)
        if args.angles: 
            if len(args.angles) != 2:
                raise BaseException("\nMust specify 2 angles")
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
            sys.stdout.write(" Processing event {}, M={:.2f} : ".format(args.event_id, events[0]._events[args.event_id].getMagnitude()))
            ele_slips = events[0].get_event_element_slips(args.event_id)
            event = events[0]._events[args.event_id]
        
        if len(ele_slips.keys()) == 0:
            raise BaseException("\nError in processing slips.")
        else:
            sys.stdout.write(" Loaded slips for {} elements :".format(len(ele_slips.keys()))) 
        sys.stdout.flush()
        
        if args.wavelength is None:
            args.wavelength = 0.21
            ## L-band 21cm is default
        
        FP = FieldPlotter(geometry, args.field_type, element_slips=ele_slips, event=event, event_id=args.event_id, cbar_max=cbar_max, levels=levels, g0=args.g, wavelength=args.wavelength)
        FP.compute_field(cutoff=1000)
        FP.plot_field(output_file=filename, angles=angles)
    if args.field_eval:
        filename = SaveFile().field_plot(args.model_file, "displacement", args.uniform_slip, args.event_id)
        sys.stdout.write(" Processing event {}, M={:.2f} : ".format(args.event_id, events[0]._events[args.event_id].getMagnitude()))
        sys.stdout.flush()
        ele_slips = events[0].get_event_element_slips(args.event_id)
        event = events[0]._events[args.event_id]
        if len(ele_slips.keys()) == 0:
            raise BaseException("\nError in processing slips.")
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
        if args.save_greens:
            GP.save_greens(filename)
    if args.traces:
        filename = SaveFile().trace_plot(args.model_file)
        if args.small_model is None: args.small_model = False
        TP = TracePlotter(geometry, filename, use_sections=args.use_sections)

    if args.slip_rates:
        if args.elements is None: args.elements = geometry.model.getElementIDs()
        slip_rates = geometry.get_slip_rates(args.elements)
        for id in slip_rates.keys():
            sys.stdout.write("{}  {}\n".format(id,slip_rates[id]))
            
    if args.slip_time_series:
        # TODO: Add multi-event file compatibility to compare between different sims
        if args.elements is None: raise BaseException("\nMust specify element ids, e.g. --elements 0 1 2")
        if args.dt is None: args.dt = 5.0  # Unit is decimal years
        if args.use_sections is not None:
            if len(args.use_sections) > 1:
                section_name = ""
                for sec in args.use_sections:
                    section_name += geometry.model.section(sec).name()+", "
            else:
                section_name = geometry.model.section(args.use_sections[0]).name()+", "
        time_series = geometry.get_slip_time_series(events[0], elements=args.elements, max_year=args.max_year, DT=args.dt)
        if len(time_series.keys()) < 10: 
            labels = time_series.keys()+[""]
        else:
            labels = [None for each in range(len(time_series.keys())+1)]
        x_data = [list(np.arange(args.dt, args.max_year+args.dt, args.dt)) for key in time_series.keys()]+[[0.0,args.max_year]]
        linewidths = [0.8 for key in time_series.keys()]+[1]
        styles = ["-" for key in time_series.keys()]+["--"]
        y_data = time_series.values()+[[0,0]]
        if args.use_sections is not None:
            plot_title = "Slip time series for {}from years 0 to {} with step {}\n{}".format(section_name,args.max_year,args.dt,args.event_file[0].split("/")[-1])
        else:
            plot_title = "Slip time series for {} elements, from years 0 to {} with step {}\n{}".format(len(args.elements),args.max_year,args.dt,args.event_file[0].split("/")[-1])
        filename = SaveFile().time_series_plot(args.event_file, "slip_", min_year=0, max_year=args.max_year, min_mag=args.min_magnitude,dt=args.dt)
        fig = plt.figure()
        BasePlotter().multi_line_plot(fig, x_data, y_data, labels, linewidths, plot_title, "simulation time [years]", "cumulative slip [m]", "", filename, linestyles=styles)
        plt.legend(loc='best')
        plt.savefig(filename, dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
        
    if args.fault_time_series:
        if args.use_faults is None: raise BaseException("\nMust specify fault id, e.g. --use_faults 33")
        if args.dt is None: args.dt = 5  # Unit is decimal years
        standardized = args.standardized
        if args.standardized is None: standardized=False
        x_data_list = [list(np.arange(0.0, args.max_year+args.dt, args.dt)) for fault_id in args.use_faults]
        x_data = x_data_list+[[0.0,args.max_year]]
        sys.stdout.write("Building fault slip time series for fault ")
        fault_time_series_data = [geometry.get_fault_averaged_slip_time_series(events[0], fault_id=fid, max_year=args.max_year, DT=args.dt, standardized=standardized) for fid in args.use_faults]
        sys.stdout.write("done.\n") ## Write that the fault time series computation is finished.
        sys.stdout.flush()  ## Always flush after you're finished
        ## If the cPickle module is available, pickle each fault time series for easier analysis later.
        ##### This saves the time series to a file that should be read like    time_data, slip_data = pickle.load(file)
        if cPickle_available:
            for i,time_series in enumerate(fault_time_series_data):
                fault_ID = args.use_faults[i]
                time_values = x_data_list[i]
                pickle_file_name = SaveFile().fault_time_series_pickle(args.event_file, fault_ID, min_year=0, max_year=args.max_year, min_mag=args.min_magnitude, combine=False, dt=args.dt, standardized=standardized)
                pickle_file = open(pickle_file_name, 'wb')
                pickle.dump([time_values,time_series], pickle_file)
                pickle_file.close()
                sys.stdout.write("Wrote time series for fault {} to {}.\n".format(fault_ID,pickle_file_name))
        else:
            sys.stdout.write("\n===Tried to save fault time series via Pickling, but cPickle module not available.===\n")
        fault_time_series = [series for series in fault_time_series_data]+[[0,0]]
        styles = ["-" for fid in args.use_faults]+["--"]
        linewidths = [1.0 for fid in args.use_faults]+[1.0]
        fault_label = ""
        for fid in args.use_faults:
            fault_label+="-{:d}".format(fid)
        if args.label:
            labels = ["{}".format(lab) for lab in args.label]+[""]
        else:
            labels = ["Fault {}".format(fid) for fid in args.use_faults]+[""]
        file_label = "Faults{}".format(fault_label)
        if standardized:
            file_label = "standardized_"+file_label
            y_label = "standardized cumulative slip"
        else:
            y_label = "cumulative slip [m]"
        filename = SaveFile().time_series_plot(args.event_file, file_label, min_year=0, max_year=args.max_year, min_mag=args.min_magnitude, dt=args.dt)
        fig = plt.figure()
        BasePlotter().multi_line_plot(fig, x_data, fault_time_series, labels, linewidths, " ", "simulation time [years]", y_label, "", filename, linestyles=styles)
        plt.legend(loc='best', fontsize=10)
        plt.savefig(filename, dpi=args.dpi)
        sys.stdout.write("\nPlot saved: {}\n\n".format(filename))


    if args.fault_group_time_series_plot:
        if args.group1_files is None or args.group2_files is None:
            raise BaseException("\nMust specify the fault time series files for each group of faults with --group1 and --group2.")
        time_data_1, average_time_series_1 = get_fault_group_average_time_series_from_pickle(args.group1_files)
        time_data_2, average_time_series_2 = get_fault_group_average_time_series_from_pickle(args.group2_files)
        x_data = [time_data_1, time_data_2]
        y_data = [average_time_series_1, average_time_series_2]
        labels = ["fault group 1", "fault group 2"]
        if args.label:
            labels = ["{}".format(lab) for lab in args.label]
        styles = ["-", "-"]
        linewidths = [1.0, 1.0]
        colors = ['r','b']
        filename = "fault_group_"+"_".join(args.group1_files[0].split('.pickle')[0].split('fault')[1].split("_")[2:])+".png"
        fig = plt.figure()
        BasePlotter().multi_line_plot(fig, x_data, y_data, labels, linewidths, " ", "simulation time [years]", "standardized slip", "", filename, linestyles=styles, colors=colors)
        plt.legend(loc='best', fontsize=13)
        plt.savefig(filename, dpi=args.dpi)
        sys.stdout.write("\nPlot saved: {}\n\n".format(filename))
        
        
    ### TO-DO: Finish this correlation function.
    ## T is total time, s1(t) and s2(t) are the slip time series for different faults at time t, 
    ##     tau is the time difference or lag, dt is time step. 
    ###   C(tau) = (1/T)*Sum_from_0_to_T( s1(t)*s2(t+tau)*dt )
    if args.fault_group_time_series_correlate:
        if args.group1_files is None or args.group2_files is None:
            raise BaseException("\nMust specify the fault time series files for each group of faults with --group1 and --group2.")
        time_data_1, average_time_series_1 = get_fault_group_average_time_series_from_pickle(args.group1_files)
        time_data_2, average_time_series_2 = get_fault_group_average_time_series_from_pickle(args.group2_files)
        total_time = max(time_data_1)-min(time_data_1)
        


    if args.event_kml:
        '''Currently this only works for a the first event file if a list of event files is given
        '''
        if args.event_id is None or args.event_file is None or args.model_file is None:
            raise BaseException("\nMust specify an event to plot with --event_id and provide an --event_file and a --model_file.")
        else:
            event = events[0]._events[args.event_id]
            filename = SaveFile().event_kml_plot(args.event_file[0], args.event_id)
            geometry.model.write_event_kml(filename, event)
            
    if args.block_area_hist:
        if args.model_file is None:
            raise BaseException("\nMust specify a fault model with --model_file.")
        else:
            units = "km^2"
            fig = plt.figure()
            model_file = args.model_file
            areas = [geometry.model.create_sim_element(elem_num).area()/1e6 for elem_num in geometry.model.getElementIDs()]
            if args.reference: 
                areas = [area/args.reference for area in areas]
                units = "{:.5f}".format(args.reference)+units
            filename = SaveFile().distribution_plot(model_file, "area_hist")
            if len(model_file.split("/")) > 1:
                model_file = model_file.split("/")[-1]
            BasePlotter().create_plot(fig, 0, "hist", False, areas, None, model_file, "element area ["+units+"]", "", filename)
            plt.savefig(filename,dpi=args.dpi)
            sys.stdout.write("Plot saved: {}\n".format(filename))

    if args.fault_length_hist:
        if args.model_file is None:
            raise BaseException("\nMust specify a fault model with --model_file.")
        else:
            units = "km"
            fig = plt.figure()
            model_file = args.model_file
            lengths = [geometry.model.fault(f_id).length()/1e3 for f_id in geometry.model.getFaultIDs()]
            if args.reference: 
                lengths = [length/args.reference for length in lengths]
                units = "{:.5f}".format(args.reference)+units
            filename = SaveFile().distribution_plot(model_file, "fault_length_hist")
            if len(model_file.split("/")) > 1:
                model_file = model_file.split("/")[-1]
            BasePlotter().create_plot(fig, 0, "hist", False, lengths, None, model_file, "fault length ["+units+"]", "", filename)
            plt.savefig(filename,dpi=args.dpi)
            sys.stdout.write("Plot saved: {}\n".format(filename))

    if args.fault_length_distribution:
        if args.model_file is None:
            raise BaseException("\nMust specify a fault model with --model_file.")
        else:
            cum_len = {}
            lens_x, lens_y = [], []
            units = "km"
            fig = plt.figure()
            model_file = args.model_file
            lengths = [geometry.model.fault(f_id).length()/1e3 for f_id in geometry.model.getFaultIDs()]
            num_faults = len(lengths)
            for num, size in enumerate(sorted(lengths)):
                cum_len[size] = num_faults - (num + 1)
            for size in sorted(cum_len.iterkeys()):
                lens_x.append(size)
                lens_y.append(cum_len[size])
            filename = SaveFile().distribution_plot(model_file, "fault_length_distrib")
            if len(model_file.split("/")) > 1:
                model_file = model_file.split("/")[-1]
            BasePlotter().create_plot(fig, 0, "loglog", True, lens_x, lens_y, model_file, "fault length L ["+units+"]", "Cumulative faults with length L or larger", filename)
            plt.scatter([max(lengths)],[1], s=50, color='r', marker='*', label="Max Length = {:.2f} [km]".format(max(lengths)))
            plt.legend(loc='best', scatterpoints=1)
            plt.ylim(0.8, plt.ylim()[1])
            plt.xlim(2,3e3)
            plt.savefig(filename,dpi=args.dpi)
            sys.stdout.write("Plot saved: {}\n".format(filename))        

    if args.block_aseismic_hist:
        if args.model_file is None:
            raise BaseException("\nMust specify a fault model with --model_file.")
        else:
            units = "aseismic fraction"
            fig = plt.figure()
            model_file = args.model_file
            fractions = [geometry.model.element(elem_num).aseismic() for elem_num in geometry.model.getElementIDs()]
            filename = SaveFile().distribution_plot(model_file, "aseismic_hist")
            if len(model_file.split("/")) > 1:
                model_file = model_file.split("/")[-1]
            BasePlotter().create_plot(fig, 0, "hist", False, fractions, None, model_file, units, "", filename)
            plt.savefig(filename,dpi=args.dpi)
            sys.stdout.write("Plot saved: {}\n".format(filename))
            
    if args.block_length_hist:
        if args.model_file is None:
            raise BaseException("\nMust specify a fault model with --model_file.")
        else:
            units = "km"
            fig = plt.figure()
            model_file = args.model_file
            lengths = [np.sqrt(geometry.model.create_sim_element(elem_num).area()/1e6) for elem_num in geometry.model.getElementIDs()]
            if args.reference: 
                lengths = [length/args.reference for length in lengths]
                units = "{:.5f}".format(args.reference)+units
            filename = SaveFile().distribution_plot(model_file, "length_hist")
            if len(model_file.split("/")) > 1:
                model_file = model_file.split("/")[-1]
            BasePlotter().create_plot(fig, 0, "hist", False, lengths, None, model_file, "element length ["+units+"]", "", filename)
            plt.savefig(filename,dpi=args.dpi)
            sys.stdout.write("Plot saved: {}\n".format(filename))

    if args.block_stress_drop_hist:
        if args.model_file is None:
            raise BaseException("\nMust specify a fault model with --model_file.")
        else:
            units = "Pa"
            if args.reference is not None: 
                if float(args.reference) == 1e6: units = "MPa"
            fig = plt.figure()
            model_file = args.model_file
            drops = np.fabs(geometry.get_stress_drops())
            factor = geometry.get_stress_drop_factor()
            if args.reference: 
                drops = [drop/args.reference for drop in drops]
            filename = SaveFile().distribution_plot(model_file, "stress_drop_hist")
            if len(model_file.split("/")) > 1:
                model_file = model_file.split("/")[-1]
            BasePlotter().create_plot(fig, 0, "hist", False, drops, None, "", "stress drop ["+units+"]", "", filename)
            plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
            plt.savefig(filename,dpi=args.dpi)
            sys.stdout.write("Plot saved: {}\n".format(filename))

    # Generate stress plots
    if args.stress_elements:
# TODO: check that stress_set is valid
        StressHistoryPlot().plot(stress_set, args.stress_elements)
        
    if args.num_sweeps:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().diagnostic_plot(args.event_file, "num_sweeps", min_year=args.min_year, max_year=args.max_year, min_mag=args.min_magnitude)
        for i, event_set in enumerate(events):
            DiagnosticPlot().plot_number_of_sweeps(fig, i, event_set, args.event_file[i].split("events_")[-1])
        plt.legend(loc='best', fontsize=8)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.event_shear_stress:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().diagnostic_plot(args.event_file, "shear_stress", min_year=args.min_year, max_year=args.max_year, min_mag=args.min_magnitude)
        for i, event_set in enumerate(events):
            DiagnosticPlot().plot_shear_stress_changes(fig, i, event_set, args.event_file[i].split("events_")[-1])
        plt.legend(loc='best', fontsize=8)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.event_normal_stress:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().diagnostic_plot(args.event_file, "normal_stress", min_year=args.min_year, max_year=args.max_year, min_mag=args.min_magnitude)
        for i, event_set in enumerate(events):
            DiagnosticPlot().plot_normal_stress_changes(fig, i, event_set, args.event_file[i].split("events_")[-1])
        plt.legend(loc='best', fontsize=8)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.event_shear_stress_vs_magnitude:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().diagnostic_plot(args.event_file, "shear_stress_vs_mag", min_year=args.min_year, max_year=args.max_year, min_mag=args.min_magnitude)
        for i, event_set in enumerate(events):
            DiagnosticPlot().plot_shear_stress_changes_vs_magnitude(fig, i, event_set, args.event_file[i].split("events_")[-1])
        plt.legend(loc='best', fontsize=8)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.event_mean_slip:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().diagnostic_plot(args.event_file, "mean_slip", min_year=args.min_year, max_year=args.max_year, min_mag=args.min_magnitude, combine=args.combine_file)
        for i, event_set in enumerate(events):
            DiagnosticPlot().plot_mean_slip(fig, i, event_set, args.event_file[i].split("events_")[-1])
        plt.legend(loc='best', fontsize=8)
        plt.savefig(filename,dpi=args.dpi)
        sys.stdout.write("Plot saved: {}\n".format(filename))
    if args.event_movie:
        if args.event_file is None or args.event_id is None or args.model_file is None:
            raise BaseException("\nMust specify event file, event id, and model file.")
        # If multiple event files are given, only use the first
        event_file = args.event_file[0]
        events = events[0]
        sim_sweeps = Sweeps(event_file, event_number=args.event_id)
        save_file = SaveFile().event_movie(event_file, args.event_id)
        sim_sweeps.event_movie(geometry, events, save_file)
    if args.spacetime:
        if args.use_faults is None:
            raise BaseException("\nMust specify a single fault for a spacetime plot, --use_faults #")
        elif len(args.use_faults) != 1:
            raise BaseException("\nMust specify a single fault for a spacetime plot, --use_faults #")
        # TODO: Make compatible with multiple event files
        fig = plt.figure()
        ax = fig.add_subplot(111)
        filename = SaveFile().event_plot(args.event_file, "spacetime", args.min_magnitude, args.min_year, args.max_year, args.combine_file)
        PLOT_TITLE = events[0].plot_str()
        if args.no_titles: PLOT_TITLE = " "
        stp = SpaceTimePlot(geometry=geometry, min_year=args.min_year, max_year=args.max_year, event_file=args.event_file[0], trigger_fault=args.use_faults[0], title=PLOT_TITLE)
        stp.plot(fig)
        plt.savefig(filename, dpi=args.dpi)
        sys.stdout.write("\nPlot saved: {}\n".format(filename))

    # Validate data if requested
    err = False
    if args.validate_slip_sum:
        events = events[0]
        mean_slip = sum(events.event_mean_slip())
        if abs(mean_slip-args.validate_slip_sum)/args.validate_slip_sum > 0.01: err = True
        sys.stdout.write("Calculated mean slip:", mean_slip, "vs. expected:", args.validate_slip_sum)

    if args.validate_mean_interevent:
        events = events[0]
        ie_times = events.interevent_times()
        mean_ie = sum(ie_times)/len(ie_times)
        if abs(mean_ie-args.mean_interevent)/args.mean_interevent > 0.02: err = True
        sys.stdout.write("Calculated mean interevent:", mean_interevent, "vs. expected:", args.mean_interevent)

    if err: exit(1)
