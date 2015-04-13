'''
# greens.py
# PyVQ greens functions analyzier, fixer, etc. functions.
# this module is intended to host scripts that analyze, modify, etc. greens function values. in particular,
# we are addressing issues like "Exploding California", in which a few greens function matrix values are non-physically
# large (in magnitude), and so lead to the (simulated) massive self destruction of the planet. 
#
# as usual, all software falls under OpenSource licensing and is available "as is", and may contain significant errors.
'''

import matplotlib
#
import numpy
import math
import pylab as plt
import h5py
import itertools

class PyGreens(object):
	# class to manage, analyze, modify greent functions. might eventually try to code this up as a quakelib class.
	#
	def __init__(self, greens_fname='model1_greens_3000.h5', gr_shear=True, gr_normal=False):
		'''
		# greens_fname: file name of greens functions to load. maybe later build in an option to load an array directly.
		# gr_shear {True/False}: load greens_shear array
		# gr_normal {True/False}: load greens_normal array   (for large simulations, it will likely be desirable to load only one of these.
		# 
		'''
		#
		self.greens_fname = greens_fname
		self.gr_shear = gr_shear
		self.gr_normal = gr_normal
		#
		# assign some functions as well:
		self.get_h5_greens_array = get_h5_greens_array
		#
		ary_shear  = None
		ary_normal = None
		#
		if gr_shear:  ary_shear  = get_h5_greens_array(greens_fname=greens_fname, shear_normal='shear')
		if gr_normal: ary_normal = get_h5_greens_array(greens_fname=greens_fname, shear_normal='normal')
		#
	#
	def get_shear(self):
		return ary_shear
	def get_normal(self):
		return ary_normal
	#
	def get_shear_h5(self, greens_fname=None):
		if greens_fname==None: greens_fname = self.greens_fname
		return get_h5_greens_array(greens_fname=greens_fname, shear_normal='shear')
	def get_normal_h5(self, greens_fname=None):
		if greens_fname==None: greens_fname = self.greens_fname
		return get_h5_greens_array(greens_fname=greens_fname, shear_normal='normal')
	#
	
	#
	def plot_greens(self, shear_normal='shear', greens_ary=None, fnum=0, do_sort=True, y_scale='log', x_scale='linear'):
		'''
		# default prams: (self, shear_normal='shear', greens_ary=None, fnum=0, do_sort=True, y_scale='log', x_scale='linear')
		#
		# note that this can memory-explode for large greens files. it might be better to use plot_green_hist with (cumulative=True) parameter
		# and bins= or n_bins={something reasonable} .
		'''
		#
		if greens_ary == None:
			shear_normal = shear_normal_aliases(shear_normal)
			if shear_normal=='greens_shear': g_data = self.get_shear()
			if shear_normal=='greens_normal': g_data = self.get_normal()
		#
		sh_0 = g_data.shape
		g_data.shape = (1, g_data.size)
		#
		plt.figure(fnum)
		plt.ion()
		plt.clf()
		ax = plt.gca()
		ax.set_xscale(x_scale)
		ax.set_yscale(y_scale)

		if not do_sort:
			ax.plot(xrange(g_data.size), g_data[0], '.')
			#plt.vlines(xrange(g_data.size), g_data[0], numpy.zeros(len(g_data[0])))
		#
		else:
			# and a distribution:
			# (but note, as stated above, for large arrays, this can be a problem; consider using plot_greens_hist() ).
			#	
			#print "lens: ", len(X), " ", len(Y)
			ax.plot([x+1 for x in xrange(len(g_data[0]))], sorted(g_data[0]), '.-')
		#
		del(g_data)
	#
	def plot_greens_hist(self, shear_normal='shear', greens_ary=None, fnum=0, do_clf=True, n_bins=1000, **hist_kwargs):
		'''
		# plot_greens_hist: plot a histogram of greens funciton values. besides giving a histogram, as oppposed to a cumulative type dist.,
		# this might be useful for really really big data sets that don't plot in memory.
		# prams (in order): (shear_normal='shear', greens_fname='model1_greens_3000.h5', fnum=0, n_bins=1000, **kwargs)
		# use kwargs to provide other hist() arguments (range, normed, weights, cumulative, bottom, histtype, align, log, ...)
		#   note that some hist_kwards will be allocated by default (bins=n_bins=1000, log=True, histtype='step')
		'''
		#
		if greens_ary == None:
			shear_normal = shear_normal_aliases(shear_normal)
			if shear_normal=='greens_shear':  g_data = self.get_shear()
			if shear_normal=='greens_normal': g_data = self.get_normal()
		#
		#n_bins = hist_kwargs.get('bins', n_bins)
		hist_kwargs['bins'] = hist_kwargs.get('bins', n_bins)
		hist_kwargs['log'] = hist_kwargs.get('log', True)
		hist_kwargs['histtype'] = hist_kwargs.get('histtype', 'step')
		#print hist_kwargs
		#
		g_data = greens_array(greens_fname=greens_fname, shear_normal=shear_normal)
		sh_0 = g_data.shape
		g_data.shape = (1, g_data.size)
		#
		#n_bins = min(n_bins, g_data.size/2)
		print "Some stats:"
		gr_val_mean   = numpy.mean(g_data[0])
		gr_val_stdev  = numpy.std(g_data[0])
		gr_val_median = numpy.median(g_data[0])
		gr_val_max, gr_val_min = max(g_data[0]), min(g_data[0])
		gr_val_max_abs, gr_val_min_abs = max(abs(g_data[0])), min(abs(g_data[0]))
		print "mean(greens): ", gr_val_mean
		print "median(greens): ", gr_val_median
		print "stdev(greens): ", gr_val_stdev
		print "max/min: ", gr_val_max, gr_val_min
		print "max/min abs: ", gr_val_max_abs, gr_val_min_abs
		#
		#
		plt.figure(fnum)
		#plt.ion()
		if do_clf: plt.clf()
		my_hist = plt.hist(g_data[0], **hist_kwargs)
		max_h = max(my_hist[0])
		plt.vlines(sorted([gr_val_mean + j*gr_val_stdev, gr_val_mean-j*gr_val_stdev] for j in xrange(4)), numpy.zeros(4), max_h*.9*numpy.ones(4), lw=1.5, alpha=.8, color='r')
		#
		#
#
def cap_greens(greens_fname='model1_greens_3000.h5', shear_normal='shear', fnum=0, n_bins=1000, top_n=.95, bottom_n=.95, **hist_kwargs):
	#
	#
	#n_bins = hist_kwargs.get('bins', n_bins)
	hist_kwargs['bins'] = hist_kwargs.get('bins', n_bins)
	hist_kwargs['log'] = hist_kwargs.get('log', True)
	hist_kwargs['histtype'] = hist_kwargs.get('histtype', 'step')
	#print hist_kwargs
	#
	g_data = greens_array(greens_fname=greens_fname, shear_normal=shear_normal)
	#
	if isinstance(top_n, float): top_n = min(1, int((1.0-top_n)*g_data.size))
	if isinstance(bottom_n, float): bottom_n = int((1.0-bottom_n)*g_data.size)
	#
	sh_0 = g_data.shape
	g_data.shape = (1, g_data.size)
	g_data[0].sort()
	#
	min_thresh = g_data[0][bottom_n]
	max_thresh = g_data[0][-(top_n+1)]
	print "min,max thresholds: %f, %f" % (min_thresh, max_thresh)
	#
	# plot first histogram:
	plt.figure(fnum)
	plt.clf()
	gr_hist_0 = plt.hist(g_data[0], **hist_kwargs)
	#
	del(g_data)
	#
	g_data = greens_array(greens_fname=greens_fname, shear_normal=shear_normal)
	#
	with h5py.File(greens_fname, 'r+') as gr_file:
		# spin through the data file; "correct" extreme values:
		gr_data = gr_file[shear_normal]
		#for j,k in itertools.product(xrange(gr_data.shape[0]), xrange(gr_data.shape[1])):
		for j,k in itertools.product(xrange(g_data.shape[0]), xrange(g_data.shape[1])):
			if g_data[j][k]<min_thresh:
				gr_data[j][k] = min_thresh
				g_data[j][k]  = min_thresh
				print "minning: %d, %d, %f" % (j,k,min_thresh)
			if g_data[j][k]>max_thresh:
				gr_data[j][k] = max_thresh
				g_data[j][k]  = max_thresh
				print "maxing: %d, %d, %f" % (j,k, max_thresh)
			#
		#
		#gr_data.flush()
		gr_file.flush()
		sh_0 = g_data.shape
		g_data.shape=(1, g_data.size)
		#
		gr_hist_1 = plt.hist(g_data[0], **hist_kwargs)
	#
	print "and hist one more time..."
	f_hist = plot_greens_hist(greens_fname=greens_fname, shear_normal=shear_normal, fnum=fnum+1, do_clf=True, n_bins=n_bins, **hist_kwargs)
	
	
#
def greens_consistency_check(greens_fname='model1_greens_3000.h5', shear_normal='shear', n_bins=1000, lowmem=False, fnum=0, **hist_kwargs):
	'''
	# one or more consistency check on greens functions. we're getting this "exploding california" problem, which is probably related
	# to some bogus greens function values... from some bogus, intersecting fault segments.
	#
	# first test, see that (g_ij^2/(g_ii)(g_ji) <<1? ~1? maybe use (g_ij * g_ji)/g_ii*g_jj.
	#
	# (generally not sure if this function is of much value).
	'''
	#
	hist_kwargs['bins'] = hist_kwargs.get('bins', n_bins)
	#
	#g_data = greens_array(greens_fname=greens_fname, shear_normal=shear_normal)
	if shear_normal.lower() in ('shear', 'shr', 'greensshear'):
		shear_normal = 'greens_shear'
	elif shear_normal.lower() in ('normal', 'nrml', 'norm', 'normalshear'):
		shear_normal = 'greens_normal'
	else:
		shear_normal = 'greens_shear'
	#
	with h5py.File(greens_fname) as gr_data:
		if lowmem:
			g_data = gr_data[shear_normal]
		else:
			g_data = gr_data[shear_normal][()]
		#
		sh = g_data.shape
		#g_data.shape = (1, g_data.size)
		print "greens array shape: ", sh
		#
		print "begin analyzing greens data:"
		#X = [[j,k, (g_data[j][k]**2.)/(g_data[j][j]*g_data[k][k]) ] for j,k in itertools.product(xrange(sh[0]), xrange(sh[1]))]
		X = [(g_data[j][k]**2.)/(g_data[j][j]*g_data[k][k]) for j,k in itertools.product(xrange(sh[0]), xrange(sh[1]))]
		print "finished analyzing greens data; now plot."
	#
	plt.hist(X, **hist_kwargs)
	#
	#del(X)
	del(g_data)
	return X

def shear_normal_aliases(shear_normal=None):
	if shear_normal==None: return None
	#
	if shear_normal.lower() in ('shear', 'shr', 'greensshear'):
		return = 'greens_shear'
	elif shear_normal.lower() in ('normal', 'nrml', 'norm', 'normalshear'):
		return = 'greens_normal'
	else:
		return = 'greens_shear'

def get_h5_greens_array(greens_fname='model1_greens_3000.h5', shear_normal='shear'):
	# return a greens arrray from a greens file.
	#
	if shear_normal.lower() in ('shear', 'shr', 'greensshear'):
		shear_normal = 'greens_shear'
	elif shear_normal.lower() in ('normal', 'nrml', 'norm', 'normalshear'):
		shear_normal = 'greens_normal'
	else:
		shear_normal = 'greens_shear'
	#
	with h5py.File(greens_fname) as gr_data:
		g_data = gr_data[shear_normal][()]
	#
	return g_data
	
