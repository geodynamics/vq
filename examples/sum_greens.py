#!/usr/bin/env python

import h5py
import sys

if len(sys.argv) != 4:
    print(sys.argv[0]+" file_name expected_normal_value expected_shear_value")
    exit(1)

# yoder:
# note that the default input shear/normal greens values (as per the specific greens tests built by Eric H.) are: 6.9056016275796917e-08 -91753588.690448046
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

# yoder: for small values, this error test is not very good. instead, let's look for an explicitly geometric error (aka, e=a/b).

normal_err = abs(1.0-(normal_sum/expected_normal))
shear_err  = abs(1.0-(shear_sum/expected_shear))
#
# yoder
normal_diff = expected_normal-normal_sum
shear_diff  = expected_shear-shear_sum
#
#normal_err = abs(expected_normal - normal_sum)/abs(max(expected_normal, normal_sum))
#shear_err = abs(expected_shear - shear_sum)/abs(max(expected_shear, shear_sum))
print("Type", "Expected", "Actual", "Error")
print("Normal", expected_normal, normal_sum, normal_err, expected_normal-normal_sum)
print("Shear", expected_shear, shear_sum, shear_err, expected_shear-shear_sum)

normal_ok = (normal_err<1e5 or normal_diff<1e-6)
shear_ok  = (shear_err<1e5  or shear_diff <1e-6)

#if normal_err > 1e-5 or shear_err > 1e-5:
#    exit(1)
if not (normal_ok and shear_ok):
    print("error.")
    exit(1)
print("ok.")
exit(0)

