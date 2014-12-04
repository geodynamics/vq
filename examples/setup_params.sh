#!/bin/bash

ELEM_SIZE=$1
DYNAMIC=$2
PREFIX=$3
PARAM_TEMPLATE=$4
PARAM_OUTPUT=$5

echo ""
echo "dynamic triggering factor:   ${DYNAMIC}"
echo "element size [m]:            ${ELEM_SIZE}"
echo "parameter file written:      params_${ELEM_SIZE}.d"
echo ""

BASE_NAME=${PREFIX}_${ELEM_SIZE}

sed -e s/DYNAMIC/${DYNAMIC}/g \
    -e s/ELEM_SIZE/${ELEM_SIZE}/g \
    -e s/INPUTFILE/${BASE_NAME}/g < ${PARAM_TEMPLATE} > ${PARAM_OUTPUT}

exit $?

