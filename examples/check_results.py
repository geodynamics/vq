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

def check_self_consistent(events):
    error = False
    for event in events:
        elements = event.getInvolvedElements()
        element_sweep_slip_sums = {}
        element_mu = {}
        element_area = {}
        for elem_id in elements: element_sweep_slip_sums[elem_id] = 0
        summed_moment = 0
        for sweep in event.getSweeps():
            element_sweep_slip_sums[sweep._element_id] += sweep._slip
            element_mu[sweep._element_id] = sweep._mu
            element_area[sweep._element_id] = sweep._area

        total_slips = {}
        for elem_num in elements:
            total_slips[elem_num] = event.getEventSlip(elem_num)
            summed_moment += element_sweep_slip_sums[elem_num]*element_area[elem_num]*element_mu[elem_num]

        # Confirm that the sum of sweep slips is equal to the total slip
        for elem_num in total_slips:
            if total_slips[elem_num] != element_sweep_slip_sums[elem_num]:
                print("ERROR: Total slip not equal to summed sweep slip for event", event.event_num, "element", elem_num)
                error = True

        # Confirm that the event magnitude is equal to the value determined from the sweeps
        # yoder: including the 1e7 term in the log argument can cause problems for really big numbers... which is likely indicative of a problem in and
        # of itself, but for now, let's just take it out so we can handle bigger numbers.
        #summed_mag = (2.0/3.0)*math.log10(1e7*summed_moment) - 10.7
        if (summed_moment <= 0):
            print("!!! Event {}, Moment {:.5f}, Mag {:.5f}".format(event.getEventNumber(), summed_moment, event.getMagnitude()))
        
        summed_mag = (2.0/3.0)*(7.0 + math.log10(summed_moment)) - 10.7
        #
        if abs(event.getMagnitude()-summed_mag) > 1e-5:
            print("ERROR: Recorded magnitude and summed sweep magnitude is not equal for event", event.event_num)
            error = True

    return error

def calc_mean_slip_sum(events):
    return sum([event.calcMeanSlip() for event in events])

def calc_mean_interevent(events):
    event_years = [event.getEventYear() for event in events]
    return sum([event_years[i+1] - event_years[i] for i in range(len(event_years)-1)])/(len(event_years)-1)

def calc_b_val(events):
    mags = [events.event_list[enum].magnitude for enum in events.event_list if events.event_list[enum].magnitude < 10 and events.event_list[enum].magnitude > 0]
    min_mag = reduce(lambda x,y: min(x,y), mags)
    max_mag = reduce(lambda x,y: max(x,y), mags)
    a_val = math.log10(len(mags))
    for i in range(10):
        cur_mag = min_mag + i*(max_mag-min_mag)/10.0
        n_above_mag = sum(1 for m in mags if m >= cur_mag)
        n_above_mag = math.log10(sum(1 for m in mags if m >= cur_mag)/float(len(mags)))
        print(cur_mag, n_above_mag)
    print(a_val, min_mag, max_mag)

def rupture_area_vs_mag(events):
    log_ra = []
    mag = []
    for event in events:
        rupture_area = event.calcEventRuptureArea()

        if not math.isnan(event.magnitude):
            log_ra.append(math.log10(rupture_area/1e6))
            mag.append(event.magnitude)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(log_ra, mag)
    print(slope, intercept)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze result file.")
    parser.add_argument('--event_file', nargs=1, required=True,
            help="Name of event file to analyze.")
    parser.add_argument('--sweep_file', nargs=1, required=True,
            help="Name of sweep file to analyze.")
    parser.add_argument('--check_consistent', action="store_true",
            help="Check internal self-consistency of result file.")
    parser.add_argument('--mean_slip', nargs=1, type=float,
            help="Perform mean slip analysis with specified expected value.")
    parser.add_argument('--mean_interevent', nargs=1, type=float,
            help="Perform mean interevent time analysis with specified expected value.")
    parser.add_argument('--gb_b_val', nargs=1, type=float,
            help="Calculate Gutenberg-Richter b value and compare with specified value.")
    if scipy_available:
        parser.add_argument('--rupture_area_vs_mag', action="store_true",
                help="Calculate rupture area vs magnitude, compare to Wells and Coppersmith.")

    args = parser.parse_args()

    event_file = args.event_file[0]
    sweep_file = args.sweep_file[0]
    mean_slip = args.mean_slip
    mean_interevent = args.mean_interevent

    # Check that inputs are valid
    if mean_slip and mean_slip[0] <= 0:
        print("ERROR: mean slip must be greater than 0")
        exit(1)

    if mean_interevent and mean_interevent[0] <= 0:
        print("ERROR: mean interevent must be greater than 0")
        exit(1)

    events = quakelib.ModelEventSet()
    events.read_file_ascii(event_file, sweep_file)

    err = False
    if args.check_consistent and check_self_consistent(events): err = True

    if args.mean_slip:
        expected_mean_slip = args.mean_slip[0]
        mean_slip = calc_mean_slip_sum(events)

        reldiff = abs(mean_slip-expected_mean_slip)/expected_mean_slip
        if reldiff > 0.01: err = True
        print("Calculated mean slip:", mean_slip, "vs. expected:", expected_mean_slip)

    if args.mean_interevent:
        expected_mean_interevent = args.mean_interevent[0]
        mean_interevent = calc_mean_interevent(events)

        reldiff = abs(mean_interevent-expected_mean_interevent)/expected_mean_interevent
        if reldiff > 0.02: err = True
        print("Calculated mean interevent:", mean_interevent, "vs. expected:", expected_mean_interevent)

    if args.gb_b_val:
        expected_gb_b_val = args.gb_b_val[0]
        b_val = calc_b_val(events)

    if scipy_available and args.rupture_area_vs_mag:
        rupture_area_vs_mag(events)

    if err: exit(1)

