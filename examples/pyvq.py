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
    def __init__(self, event_file, sweep_file):
        self._events = quakelib.ModelEventSet()
        self._events.read_file_ascii(event_file, sweep_file)
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
    def scatter_plot(self, x_data, y_data, plot_title, x_label, y_label):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        ax.set_yscale('log')
        ax.scatter(x_data, y_data)
        plt.show()

class MagnitudeFrequencyPlot(BasePlotter):
    def plot(self, events):
        self.scatter_plot(x,y,xlabel,ylabel)

class MagnitudeRuptureAreaPlot(BasePlotter):
    def plot(self, events):
        ra_list = events.event_rupture_areas()
        mag_list = events.event_magnitudes()
        ra_renorm_list = [quakelib.Conversion().sqm2sqkm(ra) for ra in ra_list]
        self.scatter_plot(mag_list, ra_renorm_list, events.plot_str(), "Magnitude", "Rupture Area (square km)")

class MagnitudeMeanSlipPlot(BasePlotter):
    def plot(self, events):
        slip_list = events.event_mean_slip()
        mag_list = events.event_magnitudes()
        self.scatter_plot(mag_list, slip_list, events.plot_str(), "Magnitude", "Mean Slip (meters)")

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
        self.scatter_plot(mag_vals, mag_norm, events.plot_str(), "Magnitude", "Frequency")
        #print(mag_norm)

if __name__ == "__main__":
    # Specify arguments
    parser = argparse.ArgumentParser(description="PyVQ.")

    # Event/model file arguments
    parser.add_argument('--event_file', required=True,
            help="Name of event file to analyze.")
    parser.add_argument('--sweep_file', required=True,
            help="Name of sweep file to analyze.")
    parser.add_argument('--model_file', required=False,
            help="Name of model (geometry) file to use in analysis.")

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

    # Plotting arguments
    parser.add_argument('--plot_freq_mag', type=float, required=False,
            help="Generate frequency magnitude plot.")

    # Validation/testing arguments
    parser.add_argument('--validate_slip_sum', type=float, required=False,
            help="Ensure the sum of mean slip for all events is within 1 percent of the specified value.")
    parser.add_argument('--validate_mean_interevent', type=float, required=False,
            help="Ensure the mean interevent time for all events is within 2 percent of the specified value.")

    args = parser.parse_args()

    # Read the event and sweeps files
    events = Events(args.event_file, args.sweep_file)

    # Read the geometry model if specified
    if args.model_file:
        model = quakelib.ModelWorld()
        model.read_file_ascii(args.model_file)
    else:
        model = None

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

