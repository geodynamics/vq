#!/bin/bash

ELEM_SIZE=$1
TAPER_METHOD=$2

mkdir -p $TAPER_METHOD

cd $TAPER_METHOD
BASE_NAME=two_fault_$ELEM_SIZE

../../../build/src/mesher \
    --import_file=../../fault_traces/multiple_fault_trace.txt \
    --import_file_type=trace \
    --import_trace_element_size=$ELEM_SIZE \
    --taper_fault_method=$TAPER_METHOD \
    --export_file=$BASE_NAME.txt \
    --export_file_type=text \
    --export_file=$BASE_NAME.kml \
    --export_file_type=kml \
    --print_statistics=statistics_$ELEM_SIZE.txt &&

exit $?

