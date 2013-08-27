#!/usr/bin/env python

import quakelib
import unittest
import math
import os

# Set of unit tests for QuakeLib library

class TestQuakeLibEQSimGeometry(unittest.TestCase):
	# Test the vertex value setting and retrieval functions and comparison functions
	def testGetSetCompareVertex(self):
		ind_val = 123
		latlonval = quakelib.LatLonDepth(30,40,1000)
		das_val = 456.78
		trace_val = quakelib.MIDDLE_TRACE
		err = quakelib.EQSimErrors()

		# Set up v1
		v1 = quakelib.EQSimGeometryVertex()
		v1.set_index(ind_val)
		v1.set_loc(latlonval)
		v1.set_das(das_val)
		v1.set_trace_flag(trace_val)

		# Set up v2 in the same way
		v2 = quakelib.EQSimGeometryVertex()
		v2.set_index(ind_val)
		v2.set_loc(latlonval)
		v2.set_das(das_val)
		v2.set_trace_flag(trace_val)

		# Confirm that v1 retains the values assigned to it
		self.assertEqual(v1.line_num(), -1)
		self.assertEqual(v1.index(), ind_val)
		self.assertEqual(v1.loc(), latlonval)
		self.assertEqual(v1.das(), das_val)
		self.assertEqual(v1.trace_flag(), trace_val)

		# Confirm that v1 is equal to v2
		self.assertEqual(v1, v2)

		# Change an attribute of v2 and confirm it is no longer equal to v1
		v2.set_loc(quakelib.LatLonDepth(40,40,1000))
		self.assertNotEqual(v1, v2)

		# Confirm that v1 and v2 pass correctness checks
		v1.validate(err)
		self.assertEqual(err.count(), 0)
		v2.validate(err)
		self.assertEqual(err.count(), 0)

		# Confirm that changing trace_flag to invalid value breaks correctness
		v1.set_trace_flag(quakelib.UNDEFINED_TRACE_STATUS)
		v1.validate(err)
		self.assertEqual(err.count(), 1)

	def testGetSetTriangle(self):
		ind_val = 123
		v_inds = [4, 5, 6]
		rake_val = 90.5
		slip_rate_val = 1.9
		aseis_factor_val = 0.5
		strike_val = 135.2
		dip_val = 14.9
		err = quakelib.EQSimErrors()

		# Set up the triangles
		t1 = quakelib.EQSimGeometryTriangle()
		t1.set_index(ind_val)
		for i, v in enumerate(v_inds): t1.set_vertex(i, v)
		t1.set_rake(rake_val)
		t1.set_slip_rate(slip_rate_val)
		t1.set_aseismic(aseis_factor_val)
		t1.set_strike(strike_val)
		t1.set_dip(dip_val)

		# Confirm that getting/setting a vertex out of bounds throws an exception
		self.assertRaises(IndexError, t1.vertex, 900)
		self.assertRaises(IndexError, t1.set_vertex, 900, 1)

		# Confirm that values read are those that were written
		self.assertEqual(t1.line_num(), -1)
		self.assertEqual(t1.index(), ind_val)
		for i, v in enumerate(v_inds): self.assertEqual(t1.vertex(i), v)
		self.assertEqual(t1.rake(), rake_val)
		self.assertEqual(t1.slip_rate(), slip_rate_val)
		self.assertEqual(t1.aseismic(), aseis_factor_val)
		self.assertEqual(t1.strike(), strike_val)
		self.assertEqual(t1.dip(), dip_val)

		# Confirm that t1 passes correctness checks
		t1.validate(err)
		self.assertEqual(err.count(), 0)

		# Confirm that changing aseismic to > 1 or < 0 generates errors
		t1.set_aseismic(-1)
		t1.validate(err)
		t1.set_aseismic(2)
		t1.validate(err)
		self.assertEqual(err.count(), 2)

		# TODO: write test of apply_remap

	def testAll(self):
		num_sec = 3
		num_tri = 3
		num_rect = 3
		err_list = quakelib.EQSimErrors()
		geom_file_name = "test_geom.dat"
		geom_file = quakelib.EQSimGeometryWriter()
		geom_file.open(geom_file_name)
		# Create 3 new sections
		for i in range(num_sec):
			rect_vert_ids = []
			new_sec = geom_file.new_section()
			# In each section create rectangles, triangles, and vertices
			for n in range(num_tri):
				tri = new_sec.new_triangle()
				for p in range(3):
					vert = new_sec.new_vertex()
					vert.set_trace_flag(quakelib.NOT_ON_TRACE)
					tri.set_vertex(p, vert.index())

			for n in range(num_rect):
				rect = new_sec.new_rectangle()
				rect.set_perfect_flag(0)
				for p in range(4):
					vert = new_sec.new_vertex()
					vert.set_trace_flag(quakelib.NOT_ON_TRACE)
					rect.set_vertex(p, vert.index())

		geom_file.validate(err_list)
#TODO: fix this
		#self.assertEqual(err_list.count(), 0)
		geom_file.write()
		geom_file.close()

		geom_in = quakelib.EQSimGeometryReader()
		geom_in.parse_file(geom_file_name)
		self.assertEqual(geom_in.num_sections(), num_sec)
		self.assertEqual(geom_in.num_vertices(), num_sec*(num_tri*3+num_rect*4))
		self.assertEqual(geom_in.num_triangles(), num_sec*num_tri)
		self.assertEqual(geom_in.num_rectangles(), num_sec*num_rect)

		os.remove(geom_file_name)

#os.remove(geom_file_name)

if __name__ == '__main__':
	unittest.main()

