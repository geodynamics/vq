#!/usr/bin/env python

from __future__ import print_function

import math
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
    def __init__(self, min_mag=-float("inf"), max_mag=float("inf")):
        self._min_mag = min_mag
        self._max_mag = max_mag

    def test_event(self, event):
        return (event.getMagnitude() >= self._min_mag and event.getMagnitude() <= self._max_mag)

    def plot_str(self):
        return "my_string"

class YearFilter:
    def __init__(self, min_year=0, max_year=10000):
        self._min_year = min_year
        self._max_year = max_year

    def test_event(self, event):
        return (event.getEventYear() >= self._min_year and event.getEventYear() <= self._max_year)

    def plot_str(self):
        return "my_string"

class EventNumFilter:
    def __init__(self, min_event_num=0, max_event_num=1000):
        self._min_event_num = min_event_num
        self._max_event_num = max_event_num

    def test_event(self, event):
        return (event.getEventNumber() >= self._min_event_num and event.getEventNumber() <= self._max_event_num)

    def plot_str(self):
        return "my_string"

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

    def set_filter(self, fltr):
        self._filtered_events = [evnum for evnum in range(len(self._events)) if fltr.test_event(self._events[evnum])]

    def interevent_times(self):
        event_times = [self._events[evnum].getEventYear() for evnum in self._filtered_events]
        return [event_times[i+1]-event_times[i] for i in xrange(len(event_times)-1)]

    def event_years(self):
        return {self._events[evnum].getEventNumber(): self._events[evnum].getEventYear() for evnum in self._filtered_events}

    def event_rupture_areas(self):
        return {self._events[evnum].getEventNumber(): self._events[evnum].calcEventRuptureArea() for evnum in self._filtered_events}

    def event_magnitudes(self):
        return {self._events[evnum].getEventNumber(): self._events[evnum].getMagnitude() for evnum in self._filtered_events}

    def event_mean_slip(self):
        return {self._events[evnum].getEventNumber(): self._events[evnum].calcMeanSlip() for evnum in self._filtered_events}

class Plotter:
    def scatter_plot(self, x_data, y_data, x_label, y_label):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.scatter(x_data, y_data)
        plt.show()

    def plot_magnitude_vs_frequency(self, events):
        self.scatter_plot(x,y,xlabel,ylabel)

    def plot_magnitude_vs_rupture_area(self, events):
        ra = events.event_rupture_areas()
        mag = events.event_magnitudes()
        ev_nums = ra.keys()
        ra_list = [quakelib.Conversion().sqm2sqkm(ra[event_num]) for event_num in ev_nums]
        mag_list = [mag[event_num] for event_num in ev_nums]
        self.scatter_plot(mag_list, ra_list, "Magnitude", "Rupture Area (square km)")

    def plot_magnitude_vs_mean_slip(self, events):
        slips = events.event_mean_slip()
        mag = events.event_magnitudes()
        ev_nums = slips.keys()
        slip_list = [slips[event_num] for event_num in ev_nums]
        mag_list = [mag[event_num] for event_num in ev_nums]
        self.scatter_plot(mag_list, slip_list, "Magnitude", "Mean Slip (meters)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PyVQ.")
    parser.add_argument('--event_file', nargs=1, required=True,
            help="Name of event file to analyze.")
    parser.add_argument('--sweep_file', nargs=1, required=True,
            help="Name of sweep file to analyze.")
    parser.add_argument('--model_file', nargs=1, required=False,
            help="Name of model (geometry) file to use in analysis.")

    # TODO: add arguments for min/max magnitude, min/max year, min/max event number, min/max section number

    args = parser.parse_args()

    event_file = args.event_file[0]
    sweep_file = args.sweep_file[0]
    if args.model_file: model_file = args.model_file[0]
    else: model_file = None

    model = quakelib.ModelWorld()
    if model_file: model.read_file_ascii(model_file)
    events = Events(event_file, sweep_file)
    #events.set_filter(SectionFilter(model, [0, 1]))
    #events.set_filter(YearFilter(max_year=300))
    events.set_filter(EventNumFilter(max_event_num=300))
    #events.set_filter(MagFilter(min_mag=3))

    plotter = Plotter()
    #plotter.plot_magnitude_vs_rupture_area(events)
    plotter.plot_magnitude_vs_mean_slip(events)

