#!/usr/bin/env python

import quakelib
import unittest
import math
import os

# Set of unit tests for QuakeLib library
# TODO: add test for non-existent file
# TODO: add test for get/set of friction attributes
# TODO: add test for exception handling

class TestEQSimFriction(unittest.TestCase):
	def testAll(self):
		err = quakelib.EQSimErrors()
		fric_file_name = "test_fric.dat"
		fric_file = quakelib.EQSimFrictionWriter()
		fric_file.open(fric_file_name)

		# Set up initial values in the file
		fric_file.set_lame_lambda_mu(1e3, 1e4)
		fric_file.set_strengths(1, 2.1, 3.1)
		fric_file.set_strengths(7, 8.1, 9.1)
		rs1 = quakelib.EQSimFrictionRateState(2, 3, 4, 5, 6)
		rs2 = quakelib.EQSimFrictionRateState(2.1, 3.1, 4.1, 5.1, 6.1)
		fric_file.set_rs_param(1, rs1)
		fric_file.set_rs_param(7, rs2)
		fric_file.validate(err)
		self.assertEqual(err.count(), 0)
		fric_file.write()
		fric_file.close()

		fric_file_in = quakelib.EQSimFrictionReader()
		fric_file_in.parse_file(fric_file_name)

		# Confirm that Lame parameters are the same
		self.assertEqual(fric_file_in.get_lame_lambda(), 1e3)
		self.assertEqual(fric_file_in.get_lame_mu(), 1e4)

		# Confirm that strengths are the same
		self.assertEqual(fric_file_in.get_static_strength(1), 2.1)
		self.assertEqual(fric_file_in.get_dynamic_strength(1), 3.1)
		self.assertEqual(fric_file_in.get_static_strength(7), 8.1)
		self.assertEqual(fric_file_in.get_dynamic_strength(7), 9.1)

		# Confirm that rate-state parameters are the same
		self.assertEqual(fric_file_in.get_rs_param(1).A(), 2)
		self.assertEqual(fric_file_in.get_rs_param(1).B(), 3)
		self.assertEqual(fric_file_in.get_rs_param(1).L(), 4)
		self.assertEqual(fric_file_in.get_rs_param(1).f0(), 5)
		self.assertEqual(fric_file_in.get_rs_param(1).V0(), 6)
		self.assertEqual(fric_file_in.get_rs_param(7).A(), 2.1)
		self.assertEqual(fric_file_in.get_rs_param(7).B(), 3.1)
		self.assertEqual(fric_file_in.get_rs_param(7).L(), 4.1)
		self.assertEqual(fric_file_in.get_rs_param(7).f0(), 5.1)
		self.assertEqual(fric_file_in.get_rs_param(7).V0(), 6.1)

		os.remove(fric_file_name)

if __name__ == '__main__':
	unittest.main()

