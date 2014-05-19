#!/bin/sh

if [ $# -ne 1 ]
then
    echo "$0 <element size>"
    exit -1
fi

ELEM_SIZE=$1

# Create the command line argument based on the list of trace files in the directory
#FILE_LIST=`ls ../fault_traces/ca_traces/trace_Greenvalley.txt ../fault_traces/ca_traces/trace_Hayward.txt`
FILE_LIST=`ls ../fault_traces/ca_traces/trace_*.txt`

EDITOR_ARGS=
for FILE in $FILE_LIST
do
    EDITOR_ARGS="$EDITOR_ARGS--import_file=$FILE --import_file_type=trace --import_trace_element_size=$ELEM_SIZE --taper_fault_method=none "
done

../../build/src/mesher $EDITOR_ARGS \
    --export_file=all_cal.txt \
    --export_file_type=text \
    --export_file=all_cal.h5 \
    --export_file_type=hdf5 \
    --export_file=all_cal.kml \
    --export_file_type=kml \
    --export_eqsim_geometry=all_cal_geom_${ELEM_SIZE}.dat \
    --export_eqsim_friction=all_cal_fric_${ELEM_SIZE}.dat \
    --merge_duplicate_verts

