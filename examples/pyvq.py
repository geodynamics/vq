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
    import matplotlib as mplt
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib_available = False

numpy_available = True
try:
    import numpy as np
except ImportError:
    numpy_available = False

class SaveFile:
    def __init__(self, event_file, plot_type):
        self.filename = plot_type+"_"+event_file.split(".")[0]+".png"

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
        return [self._events[evnum].getMagnitude() for evnum in self._filtered_events if self._events[evnum].getMagnitude() != float("-inf")]
# TODO: Handle -infinity magnitudes on the C++ side

    def event_mean_slip(self):
        return [self._events[evnum].calcMeanSlip() for evnum in self._filtered_events]

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
            ax.scatter(x_data, y_data)
        elif plot_type == "line":
            ax.plot(x_data, y_data)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

    def multi_line_plot(self, x_data, y_data, colors, labels, linewidths, plot_title, x_label, y_label, legend_str, filename):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if not (len(x_data) == len(y_data) and len(x_data) == len(colors) and len(colors) == len(labels) and len(linewidths) == len(colors)):
            sys.exit("These lists must be the same length: x_data, y_data, colors, labels, linewidths.")
        for i in range(len(x_data)):
            ax.plot(x_data[i], y_data[i], color=colors[i], label=labels[i], linewidth=linewidths[i])
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

    def scatter_and_errorbar(self, log_y, x_data, y_data, err_x, err_y, y_error, plot_title, x_label, y_label, filename, add_x = None, add_y = None, add_label = None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        ax.scatter(x_data, y_data)
        ax.errorbar(err_x, err_y, yerr = y_error, label="UCERF2", ecolor='r')
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
        ax.scatter(x_data, y_data)
        ax.plot(line_x, line_y, label = line_label, c = 'k')
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        ax.legend(loc = "best")
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

class MagnitudeRuptureAreaPlot(BasePlotter):
    def plot(self, events, filename):
        ra_list = events.event_rupture_areas()
        mag_list = events.event_magnitudes()
        ra_renorm_list = [quakelib.Conversion().sqm2sqkm(ra) for ra in ra_list]
        self.create_plot("scatter", True, mag_list, ra_renorm_list, events.plot_str(), "Magnitude", "Rupture Area (square km)", filename)

class MagnitudeMeanSlipPlot(BasePlotter):
    def plot(self, events, filename):
        slip_list = events.event_mean_slip()
        mag_list = events.event_magnitudes()
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
            add_y = 10**(math.log(freq_y[0],10)+freq_x[0]-add_x)
            add_label = "b==1"
        if UCERF2:
            self.scatter_and_errorbar(True, freq_x, freq_y, x_UCERF, y_UCERF, y_error_UCERF, events.plot_str(), "Magnitude (M)", "log(# events with mag > M /year)", filename, add_x=add_x, add_y=add_y, add_label=add_label)
        if b1 and not UCERF2:
            self.scatter_and_line(True, freq_x, freq_y, add_x, add_y, add_label, events.plot_str(), "Magnitude (M)", "log(# events with mag > M /year)", filename)
        if not UCERF2 and not b1:
            self.create_plot("scatter", True, freq_x, freq_y, events.plot_str(), "Magnitude (M)", "log(# events with mag > M /year)", filename)

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
        t0_to_eval = np.linspace(0, max_t0, num=numPoints)
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
        t0_to_plot  = np.linspace(0, int(max_t0/2.0), num=num_t0)
        match_inds  = [np.abs(np.array(t0_to_eval-t0_to_plot[i])).argmin() for i in range(len(t0_to_plot))]
        t0_to_plot  = np.array([t0_to_eval[k] for k in match_inds])
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
            labels = [round(t0_to_plot[k],2) for k in range(len(t0_to_plot))] + weib_labels
            linewidths = prob_lw + weib_lw
        else:
            colors = t0_colors
            x_data = x_data_prob
            y_data = y_data_prob
            labels = [round(t0_to_plot[k],1) for k in range(len(t0_to_plot))]
            linewidths = prob_lw
        legend_string = r't$_0$='
        y_lab         = r'P(t, t$_0$)'
        x_lab         = r't = t$_0$ + $\Delta$t [years]'
        plot_title    = ""
        self.multi_line_plot(x_data, y_data, colors, labels, linewidths, plot_title, x_lab, y_lab, legend_string, filename)

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
        t0_to_plot = np.linspace(0, int(max_t0/2.0), num=5)
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
    parser.add_argument('--plot_freq_mag', required=False, type=int,
            help="Generate frequency magnitude plot. 1: Only event data, 2: Plot b=1 Gutenberg-Richter relation, 3: Plot UCERF2 observed seismicity rates, 4: Plot UCERF2 and the b=1 line.")
    parser.add_argument('--plot_mag_rupt_area', required=False, action='store_true',
            help="Generate magnitude vs rupture area plot.")
    parser.add_argument('--plot_mag_mean_slip', required=False, action='store_true',
            help="Generate magnitude vs mean slip plot.")

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
    events = Events(args.event_file, args.event_file_type, args.sweep_file)

    # Read the geometry model if specified
    if args.model_file:
        model = quakelib.ModelWorld()
        model.read_file_ascii(args.model_file)
        # TODO: add HDF5 compatibility
    else:
        if args.use_sections:
            sys.exit("Model file required if specifying fault sections.")
        else:
            model = None

    # Check that if either beta or tau is given then the other is also given
    if (args.beta and not args.tau) or (args.tau and not args.beta):
        sys.exit("Must specify both beta and tau.")

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
        filename = SaveFile(args.event_file, "freq_mag").filename
        if args.plot_freq_mag == 1: UCERF2,b1 = False, False
        if args.plot_freq_mag == 2: UCERF2,b1 = False, True
        if args.plot_freq_mag == 3: UCERF2,b1 = True, False
        if args.plot_freq_mag == 4: UCERF2,b1 = True, True
        FrequencyMagnitudePlot().plot(events, filename, UCERF2=UCERF2, b1=b1)
    if args.plot_mag_rupt_area:
        filename = SaveFile(args.event_file, "mag_rupt_area").filename
        MagnitudeRuptureAreaPlot().plot(events, filename)
    if args.plot_mag_mean_slip:
        filename = SaveFile(args.event_file, "mag_mean_slip").filename
        MagnitudeMeanSlipPlot().plot(events, filename)
    if args.plot_prob_vs_t:
        filename = SaveFile(args.event_file, "prob_vs_time").filename
        ProbabilityPlot().plot_p_of_t(events, filename)
    if args.plot_prob_vs_t_fixed_dt:
        filename = SaveFile(args.event_file, "p_vs_t_fixed_dt").filename
        ProbabilityPlot().plot_conditional_fixed_dt(events, filename)
    if args.plot_cond_prob_vs_t:
        filename = SaveFile(args.event_file, "cond_prob_vs_t").filename
        if args.beta:
            ProbabilityPlot().plot_p_of_t_multi(events, filename, beta=args.beta, tau=args.tau)
        else:
            ProbabilityPlot().plot_p_of_t_multi(events, filename)
    if args.plot_waiting_times:
        filename = SaveFile(args.event_file, "waiting_times").filename
        ProbabilityPlot().plot_dt_vs_t0(events, filename)

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

