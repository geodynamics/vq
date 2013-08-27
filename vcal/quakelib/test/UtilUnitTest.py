#!/usr/bin/env python

import quakelib
import unittest
import math
import os
import random

# Set of unit tests for QuakeLib utilities

class TestVectors(unittest.TestCase):
    # Start with a set of basic test vectors for each test
    def setUp(self):
        self.x = quakelib.Vec3(1,2,3)
        self.y = quakelib.Vec3(2,4,6)
        self.z = quakelib.Vec3(-1,3,2)
        self.xax = quakelib.Vec3(1,0,0)
        self.yax = quakelib.Vec3(0,1,0)

    # Ensure that the [] operator works correctly in reading and assigning values
    # Also ensure that it correctly raises an exception
    def test_read_assign(self):
        x = quakelib.Vec3(1,3,5)
        self.assertEqual(x[0], 1)
        self.assertEqual(x[1], 3)
        self.assertEqual(x[2], 5)
        x[1] = -1
        self.assertEqual(x[1], -1)
        self.assertRaises(IndexError, x.__getitem__, 3)

    # Test basic vector arithmetic
    def test_arithmetic(self):
        self.assertEqual(self.x*2, self.y)
        self.assertEqual(self.x*2-self.y, quakelib.Vec3())
        # TODO: add division, other operations

    # Ensure dot products are correct and commutative
    def test_dot_prod(self):
        self.assertEqual(self.x.dot_product(self.y), 1*2+2*4+3*6)
        self.assertEqual(self.y.dot_product(self.x), self.x.dot_product(self.y))

    # Ensure cross product is correct and raises an exception for non-3D vectors
    def test_cross_prod(self):
        self.assertEqual(self.x.cross(self.y), quakelib.Vec3())
        self.assertEqual(self.y.cross(self.x), quakelib.Vec3())
        self.assertEqual(self.x.cross(self.z), quakelib.Vec3(-5,-5,5))
        self.assertRaises(ValueError, quakelib.Vec2().cross, (quakelib.Vec2()))

    # Ensure angles between vectors are calculated correctly
    def test_vector_angle(self):
        self.assertEqual(self.xax.vector_angle(self.yax), math.pi/2)
        self.assertAlmostEqual(self.xax.vector_angle(self.xax+self.yax), math.pi/4)

    # Ensure distance measurements between vectors are correct
    def test_dist(self):
        self.assertEqual(self.x.dist(self.y), math.sqrt(14))
        self.assertEqual(self.x.dist(self.z), math.sqrt(6))

    # Ensure magnitude calculations are correct
    def test_mag(self):
        self.assertEqual(self.x.mag(), math.sqrt(14))
        self.assertEqual(self.y.mag(), math.sqrt(56))
        self.assertEqual(self.z.mag(), math.sqrt(14))

    # Ensure unit vectors are generated correctly
    def test_unit_vector(self):
        self.assertEqual(self.x.unit_vector(), self.y.unit_vector())

    # Ensure axis based rotations work correctly
    def test_axis_rotate(self):
        rot_vec = self.x.rotate_around_axis(self.yax, math.pi/2)
        self.assertAlmostEqual(rot_vec[0], -3)
        self.assertAlmostEqual(rot_vec[1], 2)
        self.assertAlmostEqual(rot_vec[2], 1)

    # Ensure that the object representation can be evaluated to the original object
    def test_str_repr(self):
        x = quakelib.Vec3(1, 2, 3)
        y = quakelib.Vec2(1, 2)
        xp = eval(repr(x))
        yp = eval(repr(y))
        self.assertEqual(x, xp)
        self.assertEqual(y, yp)

class TestLatLonDepth(unittest.TestCase):
    # Ensure that out-of-bounds assignment and equality work correctly
    def test_assign(self):
        self.assertRaises(ValueError, quakelib.LatLonDepth, 91, 0, 0)
        self.assertRaises(ValueError, quakelib.LatLonDepth, -91, 0, 0)
        self.assertRaises(ValueError, quakelib.LatLonDepth, 0, 181, 0)
        self.assertRaises(ValueError, quakelib.LatLonDepth, 0, -181, 0)
        x = quakelib.LatLonDepth(1, 2, 3)
        y = quakelib.LatLonDepth(1, 2, 3)
        z = quakelib.LatLonDepth(3, 2, 3)
        self.assertRaises(ValueError, x.set_lat, -91)
        self.assertRaises(ValueError, x.set_lon, -181)
        self.assertEqual(x, y)
        self.assertNotEqual(x, z)

    # Ensure that the object representation can be evaluated to the original object
    def test_str_repr(self):
        x = quakelib.LatLonDepth(1, 2, 3)
        y = eval(repr(x))
        self.assertEqual(x, y)

class TestConversion(unittest.TestCase):
    # Check that conversions are symmetric
    def test_unit_conversion(self):
        c = quakelib.Conversion()
        self.assertEqual(c.deg2rad(c.rad2deg(1)), 1)
        self.assertEqual(c.year2sec(c.sec2year(1)), 1)
        self.assertEqual(c.m2km(c.km2m(1)), 1)
        self.assertEqual(c.sqkm2sqm(c.sqm2sqkm(1)), 1)
        self.assertEqual(c.pascal2bar(c.bar2pascal(1)), 1)

    # TODO: double check these conversions
    def test_deg_km_accuracy(self):
        c = quakelib.Conversion(quakelib.LatLonDepth(0,0))

        # Check that 360 * length of 1 longitude degree is equal to the circumference of the equator
        # Confirm accuracy is within 1 meter
        one_deg_len = c.convert2xyz(quakelib.LatLonDepth(0,1)).mag()
        self.assertAlmostEqual(one_deg_len*360.0/1000, 40075.016, 2)

        # Check that 4 * length of 90 degree vertical arc is equal to the polar circumference
        # Confirm accuracy is within 1 meter
        ninety_deg_len = c.convert2xyz(quakelib.LatLonDepth(90,0)).mag()
        self.assertAlmostEqual(ninety_deg_len*4.0/1000, 40007.860, 2)

        # Check that inverse of conversion results in the same value
        for base_lat in range(-90,91,5):
            for base_lon in range(-180, 180, 5):
                base_pt = quakelib.LatLonDepth(base_lat, base_lon)
                conv = quakelib.Conversion(base_pt)
                test_lat = math.fmod(base_lat+random.uniform(-45,45), 90)
                test_lon = math.fmod(base_lon+random.uniform(-45,45), 180)
                test_pt = quakelib.LatLonDepth(test_lat, test_lon)
                new_xyz = conv.convert2xyz(test_pt)
                rev_pt = conv.convert2LatLon(new_xyz)
                # Ensure accuracy to within 1e-7 degrees (~1 cm)
                self.assertAlmostEqual(test_lat, rev_pt.lat(), 7)
                self.assertAlmostEqual(test_lon, rev_pt.lon(), 7)

if __name__ == '__main__':
    unittest.main()

