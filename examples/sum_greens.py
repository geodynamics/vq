#!/usr/bin/env python

import h5py
import sys

if len(sys.argv) != 4:
    print(sys.argv[0]+" file_name expected_normal_value expected_shear_value")
    exit(1)

file_name = sys.argv[1]
expected_normal = float(sys.argv[2])
expected_shear = float(sys.argv[3])

#fp = h5py.File(file_name, "r")
with h5py.File(file_name, "r") as fp:
	greens_normal = fp["greens_normal"][()] # copy directly to arrays using [()] syntax.
	greens_shear = fp["greens_shear"][()]
normal_sum = sum([sum(row) for row in greens_normal])
shear_sum = sum([sum(row) for row in greens_shear])
#fp.close()

normal_err = abs(expected_normal - normal_sum)/abs(max(expected_normal, normal_sum))
shear_err = abs(expected_shear - shear_sum)/abs(max(expected_shear, shear_sum))
print("Type", "Expected", "Actual", "Error")
print("Normal", expected_normal, normal_sum, normal_err)
print("Shear", expected_shear, shear_sum, shear_err)

if normal_err > 1e-5 or shear_err > 1e-5:
    exit(1)

exit(0)

