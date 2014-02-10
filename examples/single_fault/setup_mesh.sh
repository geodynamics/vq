#!/bin/bash

ELEM_SIZE=$1
TAPER_METHOD=$2
DYNAMIC=$3

DIR_NAME=$TAPER_METHOD
mkdir -p $DIR_NAME
cd $DIR_NAME
BASE_NAME=single_fault_$ELEM_SIZE

../../src/mesher \
    --import_file=../../fault_traces/single_fault_trace.txt \
    --import_file_type=trace \
    --import_trace_element_size=$ELEM_SIZE \
    --taper_fault_method=$TAPER_METHOD \
    --export_file=$BASE_NAME.txt \
    --export_file_type=text \
    --export_file=$BASE_NAME.kml \
    --print_statistics=statistics_$ELEM_SIZE.txt \
    --export_file_type=kml &&
echo "sim.version               = 2.0
sim.time.end_year               = 10000
sim.greens.method               = standard
sim.greens.use_normal           = false
sim.friction.dynamic            = $DYNAMIC
sim.file.input                  = $BASE_NAME.txt
sim.file.input_type             = text
sim.file.output_event           = events_$ELEM_SIZE.txt
sim.file.output_sweep           = sweeps_$ELEM_SIZE.txt
sim.file.output_event_type      = text" > params_$ELEM_SIZE.d

exit $?

