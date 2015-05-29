import matplotlib
#
import numpy
import math
import pylab as plt
import h5py
import itertools
#
#plt.ion()

default_events = 'vq_output_hattonsenvy_3k/events_3000_d.h5'
events_2 = 'ca_model_hattonsenvy_105yrs_3km/events_3000.hdf5'

def quick_figs(vc_data_file=default_events, fnum_0=0, events_start=0, events_end=None, m0=7.0):
	# make some quick figures for preliminary analysis.
	with h5py.File(vc_data_file) as vc_data:
		#
		events = vc_data['events']
		#
		if events_start==None: events_start=0
		if events_end==None:   events_end=len(events)-1
		events = events[events_start:events_end]
		#		
		print "get magnitudes and then sort..."
		mags = sorted(events['event_magnitude'].tolist())
		#
		print "get delta_ts..."
		T=events['event_year']
		#dts = [[t, t - f['events'][j]['event_year']] for j,t in enumerate(f['events']['event_year'])]
		dts = [[t, t - T[j]] for j,t in enumerate(T[1:])]
		#
		print "... and bigmags "
		big_mags = [[rw['event_year'], rw['event_magnitude']] for rw in events if rw['event_magnitude']>=m0]
		big_mag_dts = [[rw[0], rw[0]-big_mags[j][0]] for j, rw in enumerate(big_mags[1:])]
		#
		print "Some summary stats:"
		mean_dt_m0 = numpy.mean(zip(*big_mag_dts)[1])
		std_dt_m0 = numpy.std(zip(*big_mag_dts)[1])
		print "mean interval (N=%d) for m>%f: %f +/- %f" % (len(big_mags), m0, mean_dt_m0, std_dt_m0)
		#
		print "and now plot..."
		#
		figs=[]
		figs+=[plt.figure(len(figs)+fnum_0)]
		plt.clf()
		#
		# first: magnitude distributions
		f=figs[-1]
		ax = plt.gca()
		ax.set_yscale('log')
		#ax.plot(mags, reversed(xrange(1, len(mags)+1)), '.-')
		ax.plot(*zip(*[[m,len(mags)-j] for j,m in enumerate(mags)]), color='b', marker='.', ls='-', zorder=4, label='Cumulative $N(>m)$')
		# and the pdf...
		dolog=True
		ax.hist(mags,bins=200, range=[min(mags), max(mags)], log=dolog, histtype='step', label='Prob. Density')
		plt.legend(loc=0, numpoints=1)
		plt.title('Magnitudes')
		#
		# magnitudes PDF only.
		figs+=[plt.figure(len(figs)+fnum_0)]
		f=figs[-1]
		f.clf()
		ax=f.gca()
		dolog=True
		ax.hist(mags,bins=200, range=[min(mags), max(mags)], log=dolog)
		plt.title('Magnitudes (pdf)')
		#
		# intervals, magnitudes time series:
		figs+=[plt.figure(len(figs)+fnum_0)]
		f=figs[-1]
		f.clf()
		ax=f.gca()
		ldT = numpy.log10(zip(*dts)[1])
		ax.set_yscale('log')
		#ax.plot(T[1:], ldT, marker='.', ls='-', color='b', label='dt(t)')
		ax.plot(T[1:], zip(*dts)[1], marker='.', ls='-', color='b', label='dt(t)')
		# set up dt range:
		dts_sorted = sorted(zip(*dts)[1])
		#dt_max = 1.1*dts_sorted[int(.9995*len(dts_sorted))]
		#dt_max = dts_sorted[-4]
		#print "dt_max at: %f (%d)" % (dt_max, int(.9*len(dts_sorted)))
		ax.set_ylim(.9*min(zip(*dts)[1]), 1.1*max(zip(*dts)[1]))
		ax.set_ylabel('Intervals $\\Delta t$')
		#ax.draw()
		ax_mags = ax.twinx()
		#ax.vlines(*(zip(*big_mags)),[3.0 for x in big_mags], color='r')
		ax_mags.vlines(*(zip(*big_mags)), ymax=[3.0 for x in big_mags], color='r', lw=2, zorder=5, label='m>%.2f' % m0)
		ax_mags.vlines(T,[3.0 for m in mags], events['event_magnitude'], color='g', zorder=3, label='magnitudes')
		ax_mags.set_ylim(2.0, 9.5)
		ax_mags.set_ylabel('magnitude')
		plt.legend(loc=0, numpoints=1)
		#
		# big-mag intervals:
		# big_mag_dts
		print "... big-mag time-series:"
		figs+=[plt.figure(len(figs)+fnum_0)]
		f=figs[-1]
		f.clf()
		ax=f.gca()
		ax.set_yscale('log')
		ax.plot(zip(*big_mag_dts)[0], zip(*big_mag_dts)[1], '.-', label='$m>%.2f intervals')
		ax_mags = ax.twinx()
		#ax.vlines(*(zip(*big_mags)),[3.0 for x in big_mags], color='r')
		ax_mags.vlines(*(zip(*big_mags)), ymax=[3.0 for x in big_mags], color='r', lw=2, zorder=5, label='m>%.2f' % m0)
		plt.legend(loc=0, numpoints=1)
		#
		# interval distributions:
		print "... and interval distribuiton..."
		figs+=[plt.figure(len(figs)+fnum_0)]
		f=figs[-1]
		f.clf()
		ax=f.gca()
		ax.set_yscale('log')
		ax.set_xscale('log')
		N=len(dts_sorted)
		ax.plot(dts_sorted, [N-j for j,dt in enumerate(dts_sorted)], '.-')
		plt.title('intervals distribuiton')
		plt.xlabel('intervals $\Delta t$')
		plt.ylabel('N(<dt)')
		
		figs+=[plt.figure(len(figs)+fnum_0)]
		f=figs[-1]
		f.clf()
		ax=f.gca()
		dolog=True
		X = numpy.log10(dts_sorted)
		ax.hist(X, bins=200, range=[min(X), max(X)], log=dolog, histtype='stepfilled')
		plt.title('intervals distribuiton (hist)')
		plt.xlabel('log intervals $\\log \left( \\Delta t \\right)$')
		plt.ylabel('N(dt)')
		
		
	return dts
#
#def plot_recurrence(
class Sweep(object):
	def __init__(self, event_number=0, vc_data_file=default_events, block_id=None):
		self.sweep_sequences=sweep_sequence(event_number=event_number, block_id=block_id, vc_data_file=vc_data_file)
		self.shear_stress_sequences = shear_stress_sequence(sweepses=self.sweep_sequences, do_print=False)
		#
		b_id_list = self.sweep_sequences['block_id'].tolist()
		self.block_ids = {x:b_id_list.count(x) for x in b_id_list}
		#self.block_ids = list(set(self.sweep_sequences['block_id'].tolist()))
		#
		# we could also, at this point, parse out the individual block sequences, maybe make a class Block().
	#
	def plot_slips(self, block_ids=None, fignum=0):
		#if block_ids==None: block_ids=self.block_ids.keys()
		#if isinstance(block_ids, float): block_ids=[int(block_ids)]
		#if isinstance(block_ids, int): block_ids = [block_ids]
		block_ids = self.check_block_ids_list(block_ids)
		#
		plt.figure(fignum)
		plt.clf()
		#
		for block_id in block_ids:
			rws = numpy.core.records.fromarrays(zip(*filter(lambda x: x['block_id']==block_id, self.shear_stress_sequences)), dtype=self.shear_stress_sequences.dtype)
			plt.plot(rws['sweep_number'], rws['block_slip'], '.-', label='block_id: %d' % block_id)
		plt.legend(loc=0, numpoints=1)
		plt.title('Block slip sequences')
		plt.xlabel('sweep number')
		plt.ylabel('slip')
	#
	def plot_stress_drop(self, block_ids=None, fignum=0):
		block_ids = self.check_block_ids_list(block_ids)
		#
		plt.figure(fignum)
		plt.clf()
		#
		for block_id in block_ids:
			rws = numpy.core.records.fromarrays(zip(*filter(lambda x: x['block_id']==block_id, self.shear_stress_sequences)), dtype=self.shear_stress_sequences.dtype)
			plt.plot(rws['sweep_number'], rws['shear_diff'], '.-', label='block_id: %d' % block_id)
		plt.plot([min(self.shear_stress_sequences['sweep_number']), max(self.shear_stress_sequences['sweep_number'])], [0., 0.], 'k-')
		plt.legend(loc=0, numpoints=1)
		plt.title('Block shear_stress drop sequences')
		plt.xlabel('sweep number')
		plt.ylabel('shear stress drop')
	#
	def plot_stress(self, block_ids=None, fignum=0):
		block_ids = self.check_block_ids_list(block_ids)
		#
		plt.figure(fignum+1)
		ax1=plt.gca()
		plt.clf()
		plt.figure(fignum)
		plt.clf()
		ax0=plt.gca()

		#
		for block_id in block_ids:
			rws = numpy.core.records.fromarrays(zip(*filter(lambda x: x['block_id']==block_id, self.shear_stress_sequences)), dtype=self.shear_stress_sequences.dtype)
			stress_seq = []
			for rw in rws:
				stress_seq += [[rw['sweep_number'], rw['shear_init']]]
				stress_seq += [[rw['sweep_number'], rw['shear_final']]]
			X,Y = zip(*stress_seq)
			#
			ax0.plot(X,Y, '.-', label='block_id: %d' % block_id)
			#
			plt.figure(fignum+1)
			plt.plot(rws['sweep_number'], rws['shear_init'], '.-', label='block_id: %d' % block_id)
			plt.plot(rws['sweep_number'], rws['shear_final'], '.-', label='block_id: %d' % block_id)
			plt.figure(fignum)
		ax0.plot([min(self.shear_stress_sequences['sweep_number']), max(self.shear_stress_sequences['sweep_number'])], [0., 0.], 'k-')
		ax0.legend(loc=0, numpoints=1)
		plt.figure(fignum)
		plt.title('Block shear_stress sequences')
		plt.xlabel('sweep number')
		plt.ylabel('shear stress')
	#
	def check_block_ids_list(self, block_ids):
		if block_ids==None: block_ids=self.block_ids.keys()
		if isinstance(block_ids, float): block_ids=[int(block_ids)]
		if isinstance(block_ids, int): block_ids = [block_ids]
		#
		return block_ids

#
def shear_stress_sequence(block_id=None, event_number=0, vc_data_file=default_events, do_print=True, sweepses=None):
	if sweepses==None: sweepses = sweep_sequence(block_id=block_id, event_number=event_number, vc_data_file=vc_data_file)
	#
	outsies = [[rw['sweep_number'], rw['block_id'], rw['block_slip'], rw['shear_init'], rw['shear_final'], rw['shear_init']-rw['shear_final'], (rw['shear_init']-rw['shear_final'])/rw['shear_final']] for rw in sweepses]
	#
	if do_print:
		for rw in outsies: print rw
	#
	cols = ['sweep_number', 'block_id', 'block_slip', 'shear_init', 'shear_final', 'shear_diff', 'shear_diff_norm']
	#outsies = numpy.core.records.fromarrays(zip(*outsies), names=cols, formats = [type(x).__name__ for x in outsies[0]])
	#return outsies
	return numpy.core.records.fromarrays(zip(*outsies), names=cols, formats = [type(x).__name__ for x in outsies[0]])
#
def sweep_sequence(event_number=0, block_id=None, vc_data_file=default_events):
	# sweep sequence for a single block in a single event.
	#
	with h5py.File(vc_data_file) as vc_data:
		sweep_range = [vc_data['events'][event_number]['start_sweep_rec'], vc_data['events'][event_number]['end_sweep_rec']]
		sweeps = vc_data['sweeps'][sweep_range[0]:sweep_range[1]][()]
	#
	# so we could filter out all the blocks != block_id, but let's just assume that we (might) want all the blocks (for default None value).
	#if block_id==None or block_id not in (sweeps['block_id']): block_id=sweeps['block_id'][0]
	if block_id!=None:
		d_type = sweeps.dtype
		#sweeps = filter(lambda x: x['block_id']==block_id, sweeps)
		sweeps = numpy.core.records.fromarrays(zip(*filter(lambda x: x['block_id']==block_id, sweeps)), dtype=d_type)		
		
	#
	return sweeps
	
def get_h5_col(col_name, vc_data_file=default_events):
	#
	if isinstance(col_name, str): col_name=[col_name]
	if col_name[0] not in ('events', 'sweeps'): col_name.insert(0,'events')
	#
	with h5py.File(vc_data_file) as vc_data:
		vc1 = vc_data[col_name[0]]
		#
		col = vc_data
		for cl in col_name:
			#
			col=col[cl]
		#
		#
	#
	return col
		
