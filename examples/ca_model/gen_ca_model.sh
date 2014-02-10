#!/bin/sh

if [ $# -ne 1 ]
then
    echo "$0 <element size>"
    exit -1
fi

ELEM_SIZE=$1

# Create the command line argument based on the list of trace files in the directory
FILE_LIST=`ls ../fault_traces/ca_traces/trace_*.txt`

EDITOR_ARGS=
for FILE in $FILE_LIST
do
    EDITOR_ARGS="$EDITOR_ARGS--import_file=$FILE --import_file_type=trace --import_trace_element_size=$ELEM_SIZE "
done

../../src/vc_edit $EDITOR_ARGS \
    --export_eqsim_geometry=all_cal_eqsim_geom_$ELEM_SIZE.txt \
    --export_eqsim_friction=all_cal_eqsim_fric_$ELEM_SIZE.txt \
    --export_file=all_cal.txt \
    --export_file_type=text \
    --export_file=all_cal.h5 \
    --export_file_type=hdf5 \
    --export_file=all_cal.kml \
    --export_file_type=kml

