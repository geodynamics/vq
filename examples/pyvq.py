#!/usr/bin/env python

from __future__ import print_function

import math
import sys
import argparse
import quakelib

scipy_available = True
try:
    import scipy.stats
except ImportError:
    scipy_available = False

matplotlib_available = True
try:
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib_available = False

numpy_available = True
try:
    import numpy as np
except ImportError:
    numpy_available = False 

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
    def __init__(self, model, section_list):
        self._section_list = section_list
        self._elem_to_section_map = {elem_num: model.element(elem_num).section_id() for elem_num in range(model.num_elements())}

    def test_event(self, event):
        event_elements = event.getInvolvedElements()
        for elem_num in event_elements:
            elem_section = self._elem_to_section_map[elem_num]
            if not elem_section in self._section_list: return True

        return False

    def plot_str(self):
        return "my_string"

class Events:
    def __init__(self, event_file, event_file_type, sweep_file = None):
        self._events = quakelib.ModelEventSet()
        if event_file_type == "hdf5":
            self._events.read_file_hdf5(event_file)
        elif event_file_type == "text" and sweep_file != None:
            self._events.read_file_ascii(event_file, sweep_file)
        else:
            sys.exit("event_file_type must be hdf5 or text. If text, a sweep_file is required.")
            
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
        event_times = [self._events[evnum].getEventYear() for evnum in self._filtered_events]
        return [event_times[i+1]-event_times[i] for i in xrange(len(event_times)-1)]

    def event_years(self):
        return [self._events[evnum].getEventYear() for evnum in self._filtered_events]

    def event_rupture_areas(self):
        return [self._events[evnum].calcEventRuptureArea() for evnum in self._filtered_events]

    def event_magnitudes(self):
        return [self._events[evnum].getMagnitude() for evnum in self._filtered_events]

    def event_mean_slip(self):
        return [self._events[evnum].calcMeanSlip() for evnum in self._filtered_events]

class BasePlotter:
    def create_plot(self, plot_type, log_y, x_data, y_data, plot_title, x_label, y_label):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        if plot_type == "scatter":
            ax.scatter(x_data, y_data)
        elif plot_type == "line":
            ax.plot(x_data, y_data)
        plt.show()

class MagnitudeRuptureAreaPlot(BasePlotter):
    def plot(self, events):
        ra_list = events.event_rupture_areas()
        mag_list = events.event_magnitudes()
        ra_renorm_list = [quakelib.Conversion().sqm2sqkm(ra) for ra in ra_list]
        self.create_plot("scatter", True, mag_list, ra_renorm_list, events.plot_str(), "Magnitude", "Rupture Area (square km)")

class MagnitudeMeanSlipPlot(BasePlotter):
    def plot(self, events):
        slip_list = events.event_mean_slip()
        mag_list = events.event_magnitudes()
        self.create_plot("scatter", True, mag_list, slip_list, events.plot_str(), "Magnitude", "Mean Slip (meters)")

class FrequencyMagnitudePlot(BasePlotter):
    def plot(self, events):
        mag_list = events.event_magnitudes()
        min_mag = min(mag_list)
        max_mag = max(mag_list)
        mag_width = (max_mag - min_mag)/99
        mag_vals = [min_mag+i*mag_width for i in range(100)]
        mag_counts = [0]*100
        for mag in mag_list:
            mag_counts[int((mag-min_mag)/mag_width)] += 1
        mag_counts.reverse()
        for i in range(len(mag_counts)-1):
            mag_counts[i+1] += mag_counts[i]
        mag_counts.reverse()
        mag_norm = [m/float(len(mag_list)) for m in mag_counts]
        self.create_plot("scatter", True, mag_vals, mag_norm, events.plot_str(), "Magnitude", "Frequency")
        #print(mag_norm)

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
        
        
class ProbabilityPlot(BasePlotter):
    def plot_p_of_t(self, events):
        # Cumulative probability P(t) as a function of interevent time t
        intervals = np.array(events.interevent_times())
        prob = {}
        prob['x'] = np.sort(intervals)
        prob['y'] = np.arange(float(intervals.size))/float(intervals.size)
        self.create_plot("line", False, prob['x'], prob['y'], events.plot_str(), "P(t)","t [years]")
        
    def plot_conditional_fixed_dt(self, events, fixed_dt=30.0):
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
        self.create_plot("line", False, prob_dt['x'], prob_dt['y'], events.plot_str(), "P(t0 + dt, t0)","t0 [years]")
    """    
    def plot_p_of_t_multi(self, events):
        # Cumulative conditional probability P(t,t0) as a function of
        # interevent time t, computed for multiple t0
        intervals = np.array(events.interevent_times())
        conditional = {}
        t0_to_eval = np.arange(0.0,int(intervals.max())+.01,1.0)
        for t0 in t0_to_eval:
            int_t0_dt = intervals[np.where( intervals > t0+fixed_dt)]
        self.create_plot("line", False, prob['x'], prob['y'], events.plot_str(), "P(t)","t [years]")
    """   
         
class Distributions:
    def weibull(X, beta, tau):
        # Return the Weibull distribution for the parameters given, 
        # for the x point or for the array of x values given.
        if len(X) == 1:
            return 1-np.exp( -(X/float(tau))**beta)
        else:
            return np.array([1-np.exp( -(X/float(tau))**beta) for x in X])
            
    def cond_weibull(X, t0, beta, tau):
        # Return the conditional Weibull distribution at a single point
        # or for an array of points.
        if len(X) == 1:
            return 1-np.exp( (t0/float(tau))**beta - (X/float(tau))**beta)
        else:
            return np.array([1-np.exp( (t0/float(tau))**beta - (x/float(tau))**beta) for x in X])

if __name__ == "__main__":
    # Specify arguments
    parser = argparse.ArgumentParser(description="PyVQ.")

    # Event/model file arguments
    parser.add_argument('--event_file', required=True,
            help="Name of event file to analyze.")
    parser.add_argument('--event_file_type', required=True,
            help="Event file type, either hdf5 or text.")
    parser.add_argument('--sweep_file', required=False,
            help="Name of sweep file to analyze.")
    parser.add_argument('--model_file', required=False,
            help="Name of model (geometry) file to use in analysis.")
    parser.add_argument('--stress_index_file', required=False,
            help="Name of stress index file to use in analysis.")
    parser.add_argument('--stress_file', required=False,
            help="Name of stress file to use in analysis.")

    # Event filtering arguments
    parser.add_argument('--min_magnitude', type=float, required=False,
            help="Minimum magnitude of events to process.")
    parser.add_argument('--max_magnitude', type=float, required=False,
            help="Maximum magnitude of events to process.")
    parser.add_argument('--min_year', type=float, required=False,
            help="Minimum year of events to process.")
    parser.add_argument('--max_year', type=float, required=False,
            help="Maximum year of events to process.")
    parser.add_argument('--min_event_num', type=float, required=False,
            help="Minimum event number of events to process.")
    parser.add_argument('--max_event_num', type=float, required=False,
            help="Maximum event number of events to process.")
    parser.add_argument('--use_sections', type=int, nargs='+', required=False,
            help="List of model sections to use (all sections used if unspecified).")

    # Event plotting arguments
    parser.add_argument('--plot_freq_mag', required=False,
            help="Generate frequency magnitude plot.")
    parser.add_argument('--plot_mag_rupt_area', required=False,
            help="Generate magnitude vs rupture area plot.")        
    parser.add_argument('--plot_mag_mean_slip', required=False,
            help="Generate magnitude vs mean slip plot.")
    parser.add_argument('--plot_prob_vs_t', required=False,
            help="Generate earthquake recurrence probability at time t plot, --use_sections required.")
    parser.add_argument('--plot_prob_vs_t_fixed_dt', required=False,
            help="Generate earthquake recurrence probability at time t + dt vs t plot, --use_sections required.")

    # Stress plotting arguments
    parser.add_argument('--stress_elements', type=int, nargs='+', required=False,
            help="List of elements to plot stress history for.")

    # Validation/testing arguments
    parser.add_argument('--validate_slip_sum', required=False,
            help="Ensure the sum of mean slip for all events is within 1 percent of the specified value.")
    parser.add_argument('--validate_mean_interevent', required=False,
            help="Ensure the mean interevent time for all events is within 2 percent of the specified value.")

    args = parser.parse_args()

    # Read the event and sweeps files
    if args.event_file_type == "hdf5":
        events = Events(args.event_file, args.event_file_type)
    else:
        events = Events(args.event_file, args.event_file_type, sweep_file = args.sweep_file)

    # Read the geometry model if specified
    if args.model_file:
        model = quakelib.ModelWorld()
        model.read_file_ascii(args.model_file)
    else:
        model = None

    # Read the stress files if specified
    if args.stress_index_file and args.stress_file:
        stress_set = quakelib.ModelStressSet()
        stress_set.read_file_ascii(args.stress_index_file, args.stress_file)
    else:
        stress_set = None

    # Set up filters
    event_filters = []
    if args.min_magnitude or args.max_magnitude:
        event_filters.append(MagFilter(min_mag=args.min_magnitude, max_mag=args.max_magnitude))

    if args.min_year or args.max_year:
        event_filters.append(YearFilter(min_mag=args.min_year, max_mag=args.max_year))

    if args.min_event_num or args.max_event_num:
        event_filters.append(EventNumFilter(min_mag=args.min_event_num, max_mag=args.max_event_num))

    if args.use_sections:
        event_filters.append(SectionFilter(model, args.use_sections))

    events.set_filters(event_filters)

    # Generate plots
    if args.plot_freq_mag:
        FrequencyMagnitudePlot().plot(events)
    if args.plot_mag_rupt_area:
        MagnitudeRuptureAreaPlot().plot(events)
    if args.plot_mag_mean_slip:
        MagnitudeMeanSlipPlot().plot(events)
    if args.plot_p_of_t:
        ProbabilityPlot().plot_p_of_t(events)
    if args.plot_prob_vs_t_fixed_dt:
        ProbabilityPlot().plot_conditional_fixed_dt(events)`

    # Generate stress plots
    if args.stress_elements:
        # TODO: check that stress_set is valid
        StressHistoryPlot().plot(stress_set, args.stress_elements)

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