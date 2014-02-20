#!/bin/bash

ELEM_SIZE=$1
TAPER_METHOD=$2
DYNAMIC=$3

echo ""
echo "dynamic triggering factor:   $DYNAMIC"
echo "taper method:                $TAPER_METHOD"
echo "element size [m]:            $ELEM_SIZE"
echo "****parameter file written:  $TAPER_METHOD/params_$ELEM_SIZE.d"
echo ""

BASE_NAME=single_fault_$ELEM_SIZE

sed -e s/DYN/$DYNAMIC/ \
    -e s/RES/$ELEM_SIZE/g \
    -e s/OUTFILE/$BASE_NAME/ <params.d >$TAPER_METHOD/params_$ELEM_SIZE.d

exit $?

