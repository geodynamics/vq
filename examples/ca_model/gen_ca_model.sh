#!/bin/sh

if [ $# -ne 1 ]
then
    echo "$0 <element size>"
    exit -1
fi

ELEM_SIZE=$1

# Create the command line argument based on the list of trace files in the directory
#FILE_LIST=`ls ../fault_traces/ca_traces/trace_SAF-Parkfield.txt ../fault_traces/ca_traces/trace_Hayward.txt`

# The list below grabs all California fault sections into
#    one large fault model.
FILE_LIST=`ls ../fault_traces/ca_traces/trace_*.txt`

EDITOR_ARGS=
for FILE in $FILE_LIST
do
    EDITOR_ARGS="$EDITOR_ARGS--import_file=$FILE --import_file_type=trace --import_trace_element_size=$ELEM_SIZE --taper_fault_method=none "
done

../../build/src/mesher $EDITOR_ARGS \
    --export_file=all_cal_fault_${ELEM_SIZE}.txt \
    --export_file_type=text \
    --export_file=all_cal_fault_${ELEM_SIZE}.kml \
    --export_file_type=kml \
    --export_file=all_cal_fault_${ELEM_SIZE}.h5 \
    --export_file_type=hdf5 \
    --export_eqsim_geometry=all_cal_${ELEM_SIZE}.dat \
    --merge_duplicate_verts

