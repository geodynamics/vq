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
import scipy.optimize
from scipy.stats import norm


class PyGreens(object):
	# class to manage, analyze, modify greent functions. might eventually try to code this up as a quakelib class.
	#
	def __init__(self, greens_fname='model1_greens_3000.h5', do_shear=True, do_normal=False):
		'''
		# greens_fname: file name of greens functions to load. maybe later build in an option to load an array directly.
		# do_shear {True/False}: load greens_shear array
		# do_normal {True/False}: load greens_normal array   (for large simulations, it will likely be desirable to load only one of these.
		# 
		'''
		#
		self.greens_fname = greens_fname
		self.do_shear = do_shear
		self.do_normal = do_normal
		#
		# assign some functions as well:
		self.get_h5_greens_array = get_h5_greens_array
		#
		self.ary_shear  = None
		self.ary_normal = None
		#
		if do_shear:  self.ary_shear  = get_h5_greens_array(greens_fname=greens_fname, shear_normal='shear')
		if do_normal: self.ary_normal = get_h5_greens_array(greens_fname=greens_fname, shear_normal='normal')
		#
	#
	def get_shear(self):
		return self.ary_shear
	def get_normal(self):
		return self.ary_normal
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
	# note: another way to do this is to make a sub-class like greens_hist(object), and then child classes:
	# plot_normal_hist(greens_hist), in which shear_normal='normal' is defined in the child objects __init__() function.
	def plot_normal_hist(self, greens_ary=None, fnum=0, do_clf=True, n_bins=1000, **hist_kwargs):
		return self.plot_greens_hist(shear_normal='normal', greens_ary=greens_ary, fnum=fnum, do_clf=do_clf, n_bins=n_bins, **hist_kwargs)
	def plot_shear_hist(self, greens_ary=None, fnum=0, do_clf=True, n_bins=1000, **hist_kwargs):
		return self.plot_greens_hist(shear_normal='shear', greens_ary=greens_ary, fnum=fnum, do_clf=do_clf, n_bins=n_bins, **hist_kwargs)

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
			print "shear_normal (translated): ", shear_normal
			if shear_normal=='greens_shear':  greens_ary = self.get_shear()
			if shear_normal=='greens_normal': greens_ary = self.get_normal()
		#
		#n_bins = hist_kwargs.get('bins', n_bins)
		hist_kwargs['bins'] = hist_kwargs.get('bins', n_bins)
		hist_kwargs['log'] = hist_kwargs.get('log', True)
		hist_kwargs['histtype'] = hist_kwargs.get('histtype', 'step')
		hist_kwargs['normed'] = hist_kwargs.get('normed', False)
		#print hist_kwargs
		#
		sh_0 = greens_ary.shape
		greens_ary.shape = (1, greens_ary.size)
		#
		#n_bins = min(n_bins, greens_ary.size/2)
		print "Some stats:"
		gr_val_mean   = numpy.mean(greens_ary[0])
		gr_val_stdev  = numpy.std(greens_ary[0])
		gr_val_median = numpy.median(greens_ary[0])
		gr_val_max, gr_val_min = max(greens_ary[0]), min(greens_ary[0])
		gr_val_max_abs, gr_val_min_abs = max(abs(greens_ary[0])), min(abs(greens_ary[0]))
		print "mean(greens): ", gr_val_mean
		print "median(greens): ", gr_val_median
		print "stdev(greens): ", gr_val_stdev
		print "max/min: ", gr_val_max, gr_val_min
		print "max/min abs: ", gr_val_max_abs, gr_val_min_abs
		#
		#
		plt.figure(fnum)
		#plt.ion()
		ax = plt.gca()
		if do_clf: plt.clf()
		gr_hist = plt.hist(greens_ary[0], **hist_kwargs)
		bin_edges=gr_hist[1]		# contains the left edges + right edge of final entry.
		bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
		#
		max_h = max(gr_hist[0])
		n_v_lines = 4
		#plt.vlines(sorted([gr_val_mean + float(j)*gr_val_stdev, gr_val_mean-j*gr_val_stdev] for j in xrange(n_v_lines)), numpy.ones(n_v_lines), max_h*.9*numpy.ones(n_v_lines), lw=1.5, alpha=.8, color='r')
		plt.vlines([gr_val_mean + 4.*float(j-n_v_lines)*gr_val_stdev for j in xrange(2*n_v_lines)], .9*min([x for x in gr_hist[0] if x!=0]), max_h, lw=1.5, alpha=.8, color='r')
		#
		# now, get a gaussian fit. we can use this to strip away extraneous values...
		# ... but this is not as straight forward as it looks; there's some log-scaling, then of course negative numbers, and
		# generally the fit needs some coersion...
		'''
		#print "begin fitting to gauss model..."
		#x_hist, y_hist = zip(*[[x,math.log10(y)] for x,y in zip(bin_centers, gr_hist[0]) if y>0])
		x_hist, y_hist = [], []
		for j,x in enumerate(bin_centers):
			#
			if gr_hist[0][j]>0:
				y_hist += [gr_hist[0][j]]
			else:
				if len(y_hist)==0: continue
				y_hist += [y_hist[-1]]
			#
			x_hist += [x]
		
		plt.figure(1)
		plt.clf()
		plt.plot(x_hist, y_hist, '-')
		#for j in xrange(len(x_hist)): print "[%f, %f]" % (x_hist[j], y_hist[j])
		#return x_hist, y_hist
		#plt.figure(0)
		
		#gauss_p0 = [numpy.log10(max_h), 0., numpy.std(numpy.log10(gr_hist[0]))]
		std_0 = numpy.std(y_hist)
		gauss_p0 = [math.log10(max(y_hist)), 0., 10.*std_0]
		#gauss_p0 = [1., 0., 1.]
		#
		print "begin fit: A, mu, sigma = ", gauss_p0
		#coeff, var_matrix = scipy.optimize.curve_fit(gauss_pdf, numpy.array(bin_centers), numpy.array(gr_hist[0]), p0=gauss_p0)
		coeff, var_matrix = scipy.optimize.curve_fit(gauss_pdf, x_hist, y_hist, p0=gauss_p0)
		coeff[0]=10**coeff[0]
		#

		print "fit complete: A, mu, sigma = ", coeff, gauss_p0
		#
		#plt.figure(1)
		#plt.clf()
		hist_fit = gauss_pdf(numpy.array(x_hist), *coeff)
		#hist_fit = [gauss_pdf(float(x), *coeff) for x in x_hist]
		
		#plt.plot(bin_centers, hist_fit, '-', lw=1.5, alpha=.7, label='gauss fit: A=%f, mu=%f, sigma=%f' % (coeff[0], coeff[1], coeff[2]))
		ax = plt.gca()
		#ax.set_yscale('log')
		plt.figure(2)
		plt.clf()
		ax=plt.gca()
		ax.plot(x_hist, hist_fit, '-', lw=1.5, alpha=.7, label='gauss fit: A=%f, mu=%f, sigma=%f' % (coeff[0], coeff[1], coeff[2]))
		'''
		#
		# return to original shape.
		greens_ary.shape=sh_0
		#
		return gr_hist
		#, x_hist, hist_fit
#
def gauss_fit_mc(y,x,nits=1000, A=1.0, mu=0.0, sigma=1.0, dy=0.0, dA=None, dmu=None, dsigma=None, ddy=None):
	
	# note: allow a 4th parameter, the y-lift.
	#if len(p0)<4: p0+=[0.]
	#
	fitfunc  = lambda p, x: (10.**p[0])*exp(-0.5*((x-p[1])/p[2])**2)+p[3]
	errfunc  = lambda p, x, y: (y - fitfunc(p, x))
	Rs = [random.Random() for k in p0]
	#
	if dA  == None: dA=max()*2.0
	if dmu == None: dmu=.5*(max(x)-min(x))
	if dsigma == None: dsigma = max(y)
	if ddy == None: ddy = .25*numpy.mean(y)
	#
	# for now, just slop through this (little or no optimization):
	for n in xrange(nits):
		# get random guesses:
		this_A = A + dA*Rs[0].random()
		this_mu = mu + dmu*(.5 - Rs[1].random())
		this_sigma = sigma + dsigma
		
	
	
def err_gauss_pdf(y, x, *p):
	return (y-gauss_pdf(x, *p))
#
def gauss_pdf(x, *p):
	'''
	# gaussian pdf for line fitting (the CDF is probably better, but...)
	# *p should be like A,mu,sigma
	so calling is like y = gauss_pdf(x, my_A, my_mu, my_sigma)
	'''
	A,mu,sigma=p
	#print "A, mu, sigma: ", A, mu, sigma
	#
	return (10.**A)*numpy.exp(-((x-mu)**2.)/(2.*sigma**2.))
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
	
	
##########################
##########################
# Module level functions (helper functions, scripts, etc.)
#######

def plot_greens_hist(greens_fname='model1_greens_3000.h5', shear_normal='shear', greens_ary=None, fnum=0, do_clf=True, n_bins=1000, **hist_kwargs):
	str_shr_norm = shear_normal_aliases(shear_normal=shear_normal)
	print "plot hist for %s array" % str_shr_norm
	#
	obj_gr = PyGreens(greens_fname=greens_fname, do_shear=(str_shr_norm=='greens_shear'), do_normal=(str_shr_norm=='greens_normal') )
	#return obj_gr
	# change this to allow normal plots too:
	gr_hist = obj_gr.plot_shear_hist(greens_ary=greens_ary, fnum=fnum, do_clf=do_clf, n_bins=n_bins, **hist_kwargs)
	#
	return gr_hist


def shear_normal_aliases(shear_normal=None):
	if shear_normal==None: return None
	#
	if shear_normal.lower() in ('shear', 'shr', 'greensshear'):
		return 'greens_shear'
	elif shear_normal.lower() in ('normal', 'nrml', 'norm', 'normalshear'):
		return 'greens_normal'
	else:
		return 'greens_shear'

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
	#print "fetching h5 data with (%s)" % shear_normal
	with h5py.File(greens_fname) as gr_data:
		g_data = gr_data[shear_normal][()]
	#
	return g_data

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
	
