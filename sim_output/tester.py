#!/usr/bin/env python

import matplotlib
import pylab as plt
import pyvq

if __name__ == "__main__":
	# for now, just use this tester script to run this one script. BUT, in the future, copy Eric/Kasey's argument parser
	# to run pyvq funcitons from the command line.
	#
	#
	my_hist = pyvq.betas.plot_greens_hist('all_cal_greens_3000.h5')
	plt.show()
