#!/usr/bin/env python

import quakelib
import unittest
import os
import tempfile

# Set of unit tests for QuakeLib library initial condition classes
# TODO: add test for non-existent file

class TestQuakeLibEQSimCondition(unittest.TestCase):
	def testGetSet(self):
		cond_file = quakelib.EQSimConditionWriter()
		cond_file.set_stresses(1, 2.1, 3.1)
		cond_file.set_stresses(7, 4.1, 5.1)

		self.assertEqual(cond_file.get_shear_stress(1), 2.1)
		self.assertEqual(cond_file.get_normal_stress(1), 3.1)
		self.assertEqual(cond_file.get_shear_stress(7), 4.1)
		self.assertEqual(cond_file.get_normal_stress(7), 5.1)

	def testExceptions(self):
		cond_file = quakelib.EQSimConditionWriter()
		self.assertRaises(IndexError, cond_file.get_shear_stress, 123)
		self.assertRaises(IndexError, cond_file.get_normal_stress, 123)
		self.assertRaises(IndexError, cond_file.get_rate_state, 123)

	def testFileReadWrite(self):
		cond_file_name = "test_cond.dat"
		cond_file = quakelib.EQSimConditionWriter()
		cond_file.open(cond_file_name)
		cond_file.set_stresses(1, 2, 3)
		cond_file.set_stresses(4, 5.1, 6.1)
		cond_file.set_rate_state(2, 3)
		cond_file.set_rate_state(7, 9.1)
		self.assertEqual(cond_file.num_elements(), 2)

		err = quakelib.EQSimErrors()
		cond_file.validate(err)
		self.assertEqual(err.count(), 0)

		cond_file.write()
		cond_file.close()

		cond_file_in = quakelib.EQSimConditionReader()
		cond_file_in.parse_file(cond_file_name)
		self.assertEqual(cond_file_in.get_shear_stress(1), 2)
		self.assertEqual(cond_file_in.get_normal_stress(1), 3)
		self.assertEqual(cond_file_in.get_shear_stress(4), 5.1)
		self.assertEqual(cond_file_in.get_normal_stress(4), 6.1)
		self.assertEqual(cond_file_in.get_rate_state(2), 3)
		self.assertEqual(cond_file_in.get_rate_state(7), 9.1)

		err = quakelib.EQSimErrors()
		cond_file_in.validate(err)
		self.assertEqual(err.count(), 0)

		os.remove(cond_file_name)

if __name__ == '__main__':
	unittest.main()

