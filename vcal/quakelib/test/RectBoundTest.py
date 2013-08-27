#!/usr/bin/env python

import quakelib
import unittest
import math

# Set of unit tests for QuakeLib library RectBound object

class TestRectBound(unittest.TestCase):
	def testInvalidRB(self):
		rb = quakelib.RectBound3()
		vec = quakelib.Vec3()
		self.assertTrue(math.isnan(rb.max_length()))
		self.assertTrue(math.isnan(rb.center()[0]))
		self.assertEqual(rb.get_child_subdivision(vec), 0)
		self.assertFalse(rb.get_child_bound(0).valid())
		self.assertFalse(rb.in_bound(vec))
		self.assertNotEqual(rb, quakelib.RectBound3())
		rb.extend_bound(vec)
		self.assertTrue(rb.valid())

	def testRBCreation(self):
		vec00 = quakelib.Vec3(0.0, 0.0, 0.0)
		vec05 = quakelib.Vec3(0.5, 0.5, 0.5)
		rb1 = quakelib.RectBound3(vec00, vec05)
		rb2 = quakelib.RectBound3(vec05, vec00)
		self.assertEqual(rb1, rb2)

	def testNormalRB(self):
		vec00 = quakelib.Vec3(0.0, 0.0, 0.0)
		vec05 = quakelib.Vec3(0.5, 0.5, 0.5)
		vec10 = quakelib.Vec3(1.0, 1.0, 1.0)
		vecneg = quakelib.Vec3(-1.0, -1.0, -1.0)
		rb = quakelib.RectBound3(vec00, vec10)
		rb2 = quakelib.RectBound3(vec00, vec05)
		self.assertEqual(rb.max_length(), 1.0)
		self.assertEqual(rb.center(), vec05)
		self.assertEqual(rb.get_child_subdivision(vec00), 0)
		self.assertEqual(rb.get_child_subdivision(vec05), 0)
		self.assertEqual(rb.get_child_subdivision(vec10), 7)
		self.assertEqual(rb.get_child_bound(0), rb2)
		self.assertTrue(rb.in_bound(vec00))
		self.assertTrue(rb.in_bound(vec05))
		self.assertFalse(rb.in_bound(vec10))
		self.assertFalse(rb.in_bound(vecneg))
		rb.extend_bound(vecneg)
		self.assertTrue(rb.in_bound(vecneg))
		rb.extend_bound(vec10)
		self.assertTrue(rb.in_bound(vec10))

if __name__ == '__main__':
	unittest.main()

