#!/usr/bin/env python

import quakelib
import unittest
import random
import math
import sys

# Set of unit tests for QuakeLib library Octree object

class TestRectBound(unittest.TestCase):
    def testAddPoint(self):
        vec00 = quakelib.Vec3(0,0,0)
        vec10 = quakelib.Vec3(1,1,1)
        rb = quakelib.RectBound3(vec00, vec10)
        octree = quakelib.Octree3(rb)
        self.assertTrue(octree.add_point(vec00, 0))
        self.assertFalse(octree.add_point(vec10, 0))
        self.assertEqual(octree.num_descendents(), 1)

    def testAddMultiRegular(self):
        lg2_num_pts = 2
        num_dim_pts = 2 ** lg2_num_pts
        step = 1.0/float(num_dim_pts)
        half_step = 1.0/float(2*num_dim_pts)
        quarter_step = 1.0/float(2*2*num_dim_pts)
        quarter_step_vec = quakelib.Vec3(quarter_step, quarter_step, quarter_step)
        vec00 = quakelib.Vec3(0,0,0)
        vec10 = quakelib.Vec3(1,1,1)
        rb = quakelib.RectBound3(vec00, vec10)
        octree = quakelib.Octree3(rb)
        pt_list = [[quakelib.Vec3(x*step+half_step, y*step+half_step, z*step+half_step), x+y*num_dim_pts+z*num_dim_pts*num_dim_pts] for x in range(num_dim_pts) for y in range(num_dim_pts) for z in range(num_dim_pts)]
        for pt in pt_list: self.assertTrue(octree.add_point(pt[0], pt[1]))
        num_branches = 0
        for i in range(lg2_num_pts): num_branches += 8 ** i
        # Confirm that the number of branches is as expected
        self.assertEqual(octree.num_descendents() - octree.num_leaves(), num_branches)
        self.assertEqual(octree.num_leaves(), num_dim_pts**3)
        # Confirm that the tree is correctly balanced
        self.assertEqual(octree.max_depth(), lg2_num_pts)
        # Confirm that identical points return the same id
        for pt in pt_list: self.assertEqual(octree.get_leaf_containing_point(pt[0]).id(), pt[1])
        # And that slightly offset points return the same id
        for pt in pt_list: self.assertEqual(octree.get_leaf_containing_point(pt[0]+quarter_step_vec).id(), pt[1])
        for pt in pt_list: self.assertEqual(octree.get_leaf_containing_point(pt[0]-quarter_step_vec).id(), pt[1])

    def testAddMultiRandom(self):
        vec00 = quakelib.Vec3(-1,-1,-1)
        vec10 = quakelib.Vec3(1,1,1)
        rb = quakelib.RectBound3(vec00, vec10)
        octree = quakelib.Octree3(rb)
        pt_list = [[quakelib.Vec3(random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)), i] for i in range(100)]
        for pt in pt_list: self.assertTrue(octree.add_point(pt[0], pt[1]))
        self.assertEqual(octree.num_leaves(), 100)

    def testAddMultiBigRange(self):
        DBL_MAX = sys.float_info.max/2
        vec00 = quakelib.Vec3(-DBL_MAX,-DBL_MAX,-DBL_MAX)
        vec10 = quakelib.Vec3(DBL_MAX,DBL_MAX,DBL_MAX)
        rb = quakelib.RectBound3(vec00, vec10)
        octree = quakelib.Octree3(rb)
        pt_list = [[quakelib.Vec3(random.uniform(-DBL_MAX, DBL_MAX), random.uniform(-DBL_MAX, DBL_MAX), random.uniform(-DBL_MAX, DBL_MAX)), i] for i in range(1000)]
        for pt in pt_list: self.assertTrue(octree.add_point(pt[0], pt[1]))
        self.assertEqual(octree.num_leaves(), 1000)

    def testAddMultiIdentical(self):
        vec00 = quakelib.Vec3(0,0,0)
        vec10 = quakelib.Vec3(1,1,1)
        rb = quakelib.RectBound3(vec00, vec10)
        octree = quakelib.Octree3(rb)
        pt = quakelib.Vec3(0.5, 0.5, 0.5)
        self.assertTrue(octree.add_point(pt, 0))
        self.assertTrue(octree.add_point(pt, 0))
        self.assertFalse(octree.add_point(pt, 1))

if __name__ == '__main__':
    unittest.main()

