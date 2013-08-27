#!/usr/bin/env python

import quakelib
import unittest
import math
import os

# Set of unit tests for QuakeLib library

class TestQuakeLibEQSimEvent(unittest.TestCase):
	def testPrep(self):
		#event definitions
		event_file=quakelib.EQSimEventWriter()
		es1=quakelib.EQSimEventSet()
		s1=quakelib.EQSimEventSummary()
		s2=quakelib.EQSimEventSummary()
		sm1=quakelib.EQSimEventSlipMap()
		sm2=quakelib.EQSimEventSlipMap()
		sm3=quakelib.EQSimEventSlipMap()
		sm4=quakelib.EQSimEventSlipMap()
		element=[]
		for i in range(8):
			element.append(quakelib.EQSimEventSlipElement(i+1, -1))
			#Check for equality
			self.assertEqual(element[i].element_id(), i+1)
		summarylist = [s1, s2]
		sliplist= [sm1, sm2, sm3, sm4]
		for i in range(4):
			sliplist[i].add_slip_entry(element[2*i])
			sliplist[i].add_slip_entry(element[2*i+1])
		#structure
		es1.add_event_summary(s1)
		es1.add_event_summary(s2)
		s1.add_slip_map(sm1)
		s1.add_slip_map(sm2)
		s2.add_slip_map(sm3)
		s2.add_slip_map(sm4)
		#Population of EQSimEvent Summaries
		for i in range(2):
			j = i+1
			summarylist[i].set_event_id(j)
			summarylist[i].set_magnitude(j)
			summarylist[i].set_time(j)
			summarylist[i].set_duration(j)
			summarylist[i].set_sid(j)
			summarylist[i].set_depth_lo(-2*j)
			summarylist[i].set_depth_hi(-1*j)
			summarylist[i].set_das_lo(j)
			summarylist[i].set_das_hi(2*j)
			summarylist[i].set_hypo_depth(2*j)
			summarylist[i].set_hypo_das(3*j)
			summarylist[i].set_area(j)
			summarylist[i].set_mean_slip(4*j)
			summarylist[i].set_moment(j)
			summarylist[i].set_shear_before(j)
			summarylist[i].set_shear_after(6*j)
			summarylist[i].set_normal_before(j)
			summarylist[i].set_normal_after(2*j)
			# Check for equality
			self.assertEqual(summarylist[i].event_id(), j)
			self.assertEqual(summarylist[i].magnitude(), j)
			self.assertEqual(summarylist[i].time(), j)
			self.assertEqual(summarylist[i].duration(), j)
			self.assertEqual(summarylist[i].sid(), j)
			self.assertEqual(summarylist[i].depth_lo(), -2*j)
			self.assertEqual(summarylist[i].depth_hi(), -1*j)
			self.assertEqual(summarylist[i].das_lo(), j)
			self.assertEqual(summarylist[i].das_hi(), 2*j)
			self.assertEqual(summarylist[i].hypo_depth(), 2*j)
			self.assertEqual(summarylist[i].hypo_das(), 3*j)
			self.assertEqual(summarylist[i].area(), j)
			self.assertEqual(summarylist[i].mean_slip(), 4*j)
			self.assertEqual(summarylist[i].moment(), j)
			self.assertEqual(summarylist[i].shear_before(), j)
			self.assertEqual(summarylist[i].shear_after(), 6*j)
			self.assertEqual(summarylist[i].normal_before(), j)
			self.assertEqual(summarylist[i].normal_after(), 2*j)

		#Population of EQSimEvent Slip Map
		for i in range(4):
			j= i+1
			sliplist[i].set_depth_lo(-2*j)
			sliplist[i].set_depth_hi(-1*j)
			sliplist[i].set_das_lo(j)
			sliplist[i].set_das_hi(2*j)
			sliplist[i].set_area(j)
			sliplist[i].set_mean_slip(4*j)
			sliplist[i].set_moment(j)
			sliplist[i].set_shear_before(j)
			sliplist[i].set_shear_after(6*j)
			sliplist[i].set_normal_before(j)
			sliplist[i].set_normal_after(2*j)
		# Check for equality
			self.assertEqual(sliplist[i].depth_lo(), -2*j)
			self.assertEqual(sliplist[i].depth_hi(), -1*j)
			self.assertEqual(sliplist[i].das_lo(), j)
			self.assertEqual(sliplist[i].das_hi(), 2*j)
			self.assertEqual(sliplist[i].area(), j)
			self.assertEqual(sliplist[i].mean_slip(), 4*j)
			self.assertEqual(sliplist[i].moment(), j)
			self.assertEqual(sliplist[i].shear_before(), j)
			self.assertEqual(sliplist[i].shear_after(), 6*j)
			self.assertEqual(sliplist[i].normal_before(), j)
			self.assertEqual(sliplist[i].normal_after(), 2*j)

		err_list = quakelib.EQSimErrors()
		event_file.validate(err_list)

	def testAll(self):
		#quakelib.EQSimEventWriter().flush()
		event_file = quakelib.EQSimEventWriter()
		err_list = quakelib.EQSimErrors()
		event_file_name = "test_event.dat"
		event_file.open(event_file_name)
		#event definitions
		s1=quakelib.EQSimEventSummary()
		s2=quakelib.EQSimEventSummary()
		es1=quakelib.EQSimEventSet()
		sm1=quakelib.EQSimEventSlipMap()
		sm2=quakelib.EQSimEventSlipMap()
		sm3=quakelib.EQSimEventSlipMap()
		sm4=quakelib.EQSimEventSlipMap()
		element=[]
		for i in range(8):
			element.append(quakelib.EQSimEventSlipElement(i+1, -1))
		summarylist = [s1, s2]
		sliplist= [sm1, sm2, sm3, sm4]
		for i in range(4):
			sliplist[i].add_slip_entry(element[2*i])
			sliplist[i].add_slip_entry(element[2*i+1])
		es1.add_event_summary(s1)
		es1.add_event_summary(s2)
		s1.add_slip_map(sm1)
		s1.add_slip_map(sm2)
		s2.add_slip_map(sm3)
		s2.add_slip_map(sm4)
		for i in range(2):
			j = i+1
			summarylist[i].set_event_id(j)
			summarylist[i].set_magnitude(j)
			summarylist[i].set_time(j)
			summarylist[i].set_duration(j)
			summarylist[i].set_sid(j)
			summarylist[i].set_depth_lo(-2*j)
			summarylist[i].set_depth_hi(-1*j)
			summarylist[i].set_das_lo(j)
			summarylist[i].set_das_hi(2*j)
			summarylist[i].set_hypo_depth(2*j)
			summarylist[i].set_hypo_das(3*j)
			summarylist[i].set_area(j)
			summarylist[i].set_mean_slip(4*j)
			summarylist[i].set_moment(j)
			summarylist[i].set_shear_before(j)
			summarylist[i].set_shear_after(6*j)
			summarylist[i].set_normal_before(j)
			summarylist[i].set_normal_after(2*j)
		for i in range(4):
			j= i+1
			sliplist[i].set_depth_lo(-2*j)
			sliplist[i].set_depth_hi(-1*j)
			sliplist[i].set_das_lo(j)
			sliplist[i].set_das_hi(2*j)
			sliplist[i].set_area(j)
			sliplist[i].set_mean_slip(4*j)
			sliplist[i].set_moment(j)
			sliplist[i].set_shear_before(j)
			sliplist[i].set_shear_after(6*j)
			sliplist[i].set_normal_before(j)
			sliplist[i].set_normal_after(2*j)

		#test

		event_file.add_event_summary(s1)
		event_file.add_event_summary(s2)

		event_file.validate(err_list)
		event_file.write()
		self.assertEqual(err_list.count(), 0)
		event_file.close()
		event_file.flush()
		event_file_in=quakelib.EQSimEventReader()
		event_file_in.parse_file(event_file_name)
		ies=event_file_in.event_summaries
		for i in ies:
			j=i+1
			self.assertEqual(ies[i].event_id(), j)
			self.assertEqual(ies[i].magnitude(), j)
			self.assertEqual(ies[i].time(), j)
			self.assertEqual(ies[i].duration(), j)
			self.assertEqual(ies[i].sid(), j)
			self.assertEqual(ies[i].depth_lo(), -2*j)
			self.assertEqual(ies[i].depth_hi(), -1*j)
			self.assertEqual(ies[i].das_lo(), j)
			self.assertEqual(ies[i].das_hi(), 2*j)
			self.assertEqual(ies[i].hypo_depth(), 2*j)
			self.assertEqual(ies[i].hypo_das(), 3*j)
			self.assertEqual(ies[i].area(), j)
			self.assertEqual(ies[i].mean_slip(), 4*j)
			self.assertEqual(ies[i].moment(), j)
			self.assertEqual(ies[i].shear_before(), j)
			self.assertEqual(ies[i].shear_after(), 6*j)
			self.assertEqual(ies[i].normal_before(), j)
			self.assertEqual(ies[i].normal_after(), 2*j)


		eset = quakelib.EQSimEventSet()
		os.remove(event_file_name)

if __name__ == '__main__':
	unittest.main()
