#!/bin/bash

ELEM_SIZE=$1
DYNAMIC=$2

echo ""
echo "dynamic triggering factor:   $DYNAMIC"
echo "element size [m]:            $ELEM_SIZE"
echo "****parameter file written:  params_$ELEM_SIZE.d"
echo ""

BASE_NAME=two_fault_$ELEM_SIZE

sed -e s/DYN/$DYNAMIC/ \
    -e s/RES/$ELEM_SIZE/g \
    -e s/OUTFILE/$BASE_NAME/ <../../../examples/diff_fault/params.d >params_two_$ELEM_SIZE.d

exit $?

