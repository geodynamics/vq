#!/usr/bin/env python

import quakelib
import unittest
import math

class TestGreenFunctionCalc(unittest.TestCase):
    def setUp(self):
        self.ok = quakelib.Okada()
        self.slip_rad = math.pi
        self.slip_vec = quakelib.Vec3(math.sin(self.slip_rad), math.cos(self.slip_rad), 0)
        self.dip_rad = math.pi/2.0
        # Tolerate at most a 1e-10 absolute difference in magnitude
        self.mag_tol = 1e-9

    # Check that displacement scales properly as the fault size grows
    def testDispCalcFaultSize(self):
        baseline = [0.102036860007, 0.152747014539, 0.178027759544, 0.182776615629, 0.183223384675, 0.183254816562, 0.183256842643, 0.183256970265, 0.183256978257, 0.183256978757]
        for i in range(10):
            fault_length = math.pow(2, i)
# Note that the location can't be too close or we trigger the boundary case and the results don't fit a curve
            loc = quakelib.Vec3(fault_length/2.0, 1, 0)
            source_dim = quakelib.Vec3(fault_length, 1, 0)
            # Ensure displacement is within acceptable bounds
            disp = self.ok.calc_displacement_vector(loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1],      self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
            rel_err = abs(baseline[i]-disp.mag())/baseline[i]
            self.assertTrue(rel_err < 1e-8)

    # Check that the Greens functions return symmetric results for a vertical strike slip fault
    # In other words, for a fault centered at (0,0), the displacement, dudx, etc
    # at (x, y) should be equivalent to that at (-x, -y) with a change in sign of Z
    def testGreenSymmetricDisplacement(self):
        source_dim = quakelib.Vec3(1, 1, 1)

        for x in range(-2, 8):
            for y in range(-2, 8):
                # Set up the test location and mirror location
                z = 0
                xloc = source_dim[0]/2.0+2**x
                yloc = 2**y
                zloc = -z
                orig_loc = quakelib.Vec3(xloc, yloc, zloc)
                xloc = source_dim[0]/2.0-2**x
                yloc = -2**y
                mirror_loc = quakelib.Vec3(xloc, yloc, zloc)

                # Calculate the displacements
                orig_disp = self.ok.calc_displacement_vector(orig_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                mirror_disp = self.ok.calc_displacement_vector(mirror_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                abs_err = abs(orig_disp.mag()-mirror_disp.mag())
                self.assertTrue(abs_err < self.mag_tol)

    def testGreenSymmetricDuDx(self):
        source_dim = quakelib.Vec3(1, 1, 1)

        for x in range(-2, 8):
            for y in range(-2, 8):
                # Set up the test location and mirror location
                z = 0
                xloc = source_dim[0]/2.0+2**x
                yloc = 2**y
                zloc = -z
                orig_loc = quakelib.Vec3(xloc, yloc, zloc)
                xloc = source_dim[0]/2.0-2**x
                yloc = -2**y
                mirror_loc = quakelib.Vec3(xloc, yloc, zloc)

                # Calculate dudx
                orig_dudx = self.ok.calc_dudx(orig_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                mirror_dudx = self.ok.calc_dudx(mirror_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                abs_err = abs(orig_dudx.mag()-mirror_dudx.mag())
                self.assertTrue(abs_err < self.mag_tol)

    def testGreenSymmetricDuDy(self):
        source_dim = quakelib.Vec3(1, 1, 1)

        for x in range(-2, 8):
            for y in range(-2, 8):
                # Set up the test location and mirror location
                z = 0
                xloc = source_dim[0]/2.0+2**x
                yloc = 2**y
                zloc = -z
                orig_loc = quakelib.Vec3(xloc, yloc, zloc)
                xloc = source_dim[0]/2.0-2**x
                yloc = -2**y
                mirror_loc = quakelib.Vec3(xloc, yloc, zloc)

                # Calculate dudy
                orig_dudy = self.ok.calc_dudy(orig_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                mirror_dudy = self.ok.calc_dudy(mirror_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                abs_err = abs(orig_dudy.mag()-mirror_dudy.mag())
                self.assertTrue(abs_err < self.mag_tol)

    def testGreenSymmetricDuDz(self):
        source_dim = quakelib.Vec3(1, 1, 1)

        for x in range(-2, 8):
            for y in range(-2, 8):
                # Set up the test location and mirror location
                z = 0
                xloc = source_dim[0]/2.0+2**x
                yloc = 2**y
                zloc = -z
                orig_loc = quakelib.Vec3(xloc, yloc, zloc)
                xloc = source_dim[0]/2.0-2**x
                yloc = -2**y
                mirror_loc = quakelib.Vec3(xloc, yloc, zloc)

                # Calculate dudz
                orig_dudz = self.ok.calc_dudz(orig_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                mirror_dudz = self.ok.calc_dudz(mirror_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                abs_err = abs(orig_dudz.mag()-mirror_dudz.mag())
                self.assertTrue(abs_err < self.mag_tol)

    def testGreenSymmetricTensor(self):
        source_dim = quakelib.Vec3(1, 1, 1)

        for x in range(-2, 8):
            for y in range(-2, 8):
                # Set up the test location and mirror location
                z = 0
                xloc = source_dim[0]/2.0+2**x
                yloc = 2**y
                zloc = -z
                orig_loc = quakelib.Vec3(xloc, yloc, zloc)
                xloc = source_dim[0]/2.0-2**x
                yloc = -2**y
                mirror_loc = quakelib.Vec3(xloc, yloc, zloc)

                # Calculate tensor and determine shear/normal stresses
                orig_tensor = self.ok.calc_stress_tensor(orig_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                mirror_tensor = self.ok.calc_stress_tensor(mirror_loc, source_dim[2], self.dip_rad, source_dim[0], source_dim[1], self.slip_vec[0], self.slip_vec[1], self.slip_vec[2], 1, 1)
                rake_vec = quakelib.Vec3(1, 0, 0)
                normal_vec = quakelib.Vec3(0, 1, 0)
                orig_stress_vec = orig_tensor*normal_vec
                mirror_stress_vec = mirror_tensor*normal_vec
                # Shear stresses should be exactly opposite in sign
                orig_shear_stress = orig_stress_vec.dot_product(rake_vec)
                mirror_shear_stress = mirror_stress_vec.dot_product(rake_vec)
                abs_err = abs(orig_shear_stress + mirror_shear_stress)
                self.assertTrue(abs_err < self.mag_tol)
                # Normal stresses
                orig_normal_stress = orig_stress_vec.dot_product(normal_vec)
                mirror_normal_stress = mirror_stress_vec.dot_product(normal_vec)
                abs_err = abs(orig_normal_stress + mirror_normal_stress)
                self.assertTrue(abs_err < self.mag_tol)

if __name__ == "__main__":
    unittest.main()

