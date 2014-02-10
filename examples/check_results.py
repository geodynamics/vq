#!/usr/bin/env python

from __future__ import print_function

import math
import sys
import argparse

class EventSweep:
    def __init__(self):
        self.sweep_num = -1
        self.block_id = self.slip = self.area = self.mu = -1
        self.shear_init = self.shear_final = -1
        self.normal_init = self.normal_final = -1

class Event:
    def __init__(self):
        self.event_num = self.year = self.trigger = self.magnitude = -1
        self.stress_init = self.normal_init = -1
        self.stress_final = self.normal_final = -1
        self.sweeps = []

    def get_blocks(self):
        block_set = set()
        for sweep in self.sweeps:
            block_set.add(sweep.block_id)
        return block_set

    def get_sweeps(self):
        return self.sweeps

    def get_total_slips(self, block_set):
        total_slips = {}
        for bid in block_set: total_slips[bid] = 0
        for sweep in self.sweeps:
            bid = sweep.block_id
            if bid in total_slips: total_slips[bid] += sweep.slip
        return total_slips

    def __str__(self):
        return str(self.year)+":M"+str(self.magnitude)

class EventFile:
    def __init__(self, event_file_name, sweep_file_name):
        event_file = open(event_file_name)
        self.event_list = {}
        while 1:
            line = event_file.readline()
            if not line: break
            data = line.split()
            new_event = Event()
            new_event.event_num = int(data[0])
            new_event.year = float(data[1])
            new_event.trigger = int(data[2])
            new_event.magnitude = float(data[3])
            new_event.stress_init = float(data[4])
            new_event.normal_init = float(data[5])
            new_event.stress_final = float(data[6])
            new_event.normal_final = float(data[7])
            self.event_list[new_event.event_num] = new_event

        event_file.close()

        sweep_file = open(sweep_file_name)
        while 1:
            line = sweep_file.readline()
            if not line: break
            sweep_data = line.split()
            new_sweep = EventSweep()
            event_num = int(sweep_data[0])
            new_sweep.sweep_num = int(sweep_data[1])
            new_sweep.block_id = int(sweep_data[2])
            new_sweep.slip = float(sweep_data[3])
            new_sweep.area = float(sweep_data[4])
            new_sweep.mu = float(sweep_data[5])
            new_sweep.shear_init = float(sweep_data[6])
            new_sweep.shear_final = float(sweep_data[7])
            new_sweep.normal_init = float(sweep_data[8])
            new_sweep.normal_final = float(sweep_data[9])
            self.event_list[event_num].sweeps.append(new_sweep)

        sweep_file.close()

def check_self_consistent(events):
    error = False
    for event_num in events.event_list:
        event = events.event_list[event_num]
        blocks = event.get_blocks()
        block_sweep_slip_sums = {}
        for block in blocks: block_sweep_slip_sums[block] = 0
        summed_moment = 0
        sweep_list = event.get_sweeps()
        for sweep in sweep_list:
            block_sweep_slip_sums[sweep.block_id] += sweep.slip
            summed_moment += sweep.slip*sweep.area*sweep.mu

        total_slips = event.get_total_slips(blocks)

        # Confirm that the sum of sweep slips is equal to the total slip
        for bnum in total_slips:
            if total_slips[bnum] != block_sweep_slip_sums[bnum]:
                print("ERROR: Total slip not equal to summed sweep slip for event", event.event_num, "block", bnum)
                error = True

        summed_mag = (2.0/3.0)*math.log10(1e7*summed_moment) - 10.7
        if abs(event.magnitude-summed_mag) > 1e-5:
            print("ERROR: Recorded magnitude and summed sweep magnitude is not equal for event", event.event_num)
            error = True

    return error

def calc_mean_slip(events):
    block_total_slips = {}
    for event_num in events.event_list:
        event = events.event_list[event_num]
        blocks = event.get_blocks()
        block_sweep_slip_sums = {}
        for block in blocks: block_sweep_slip_sums[block] = 0
        sweep_list = event.get_sweeps()
        for sweep in sweep_list:
            block_sweep_slip_sums[sweep.block_id] += sweep.slip
            if not block_total_slips.has_key(sweep.block_id): block_total_slips[sweep.block_id] = 0.0
            block_total_slips[sweep.block_id] += sweep.slip

        total_slips = event.get_total_slips(blocks)

    mean_total_slip = 0.0
    total_slips = [block_total_slips[bid] for bid in block_total_slips.keys()]

    return sum(total_slips)/float(len(total_slips))

def calc_mean_interevent(events):
    interevent_times = []
    last_year = -1
    for event_num in events.event_list:
        event = events.event_list[event_num]
        if last_year > 0: interevent_times.append(event.year - last_year)
        last_year = event.year

    return sum(interevent_times)/float(len(interevent_times))

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

    err = False
    events = EventFile(event_file, sweep_file)
    if args.check_consistent and check_self_consistent(events): err = True

    if args.mean_slip:
        expected_mean_slip = args.mean_slip[0]
        mean_slip = calc_mean_slip(events)

        reldiff = abs(mean_slip-expected_mean_slip)/expected_mean_slip
        if reldiff > 0.01: err = True
        print("Calculated mean slip:", mean_slip, "vs. expected:", expected_mean_slip)

    if args.mean_interevent:
        expected_mean_interevent = args.mean_interevent[0]
        mean_interevent = calc_mean_interevent(events)

        reldiff = abs(mean_interevent-expected_mean_interevent)/expected_mean_interevent
        if reldiff > 0.01: err = True
        print("Calculated mean interevent:", mean_interevent, "vs. expected:", expected_mean_interevent)

    if err: exit(1)

