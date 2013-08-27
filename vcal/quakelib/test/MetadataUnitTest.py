#!/usr/bin/env python

import quakelib
import unittest
import math
import os

# Set of unit tests for EqSim library
# TODO: add test for non-existent file

class TestQuakeLibEQSimMetadata(unittest.TestCase):
	# Ensure that passing invalid arguments will raise exceptions
	def testExceptionRaising(self):
		record_types = [quakelib.META_COMMENT, quakelib.META_SIGNATURE, quakelib.META_INFO,
						quakelib.META_TITLE, quakelib.META_AUTHOR, quakelib.META_DATE]
		bad_char_list = ["\r", "\n"]
		md_file = quakelib.EQSimMetadata()
		for rec_type in record_types:
			# Check that out of bounds indices raise an exception
			self.assertRaises(IndexError, md_file.meta_get_record, rec_type, 1)
			self.assertRaises(IndexError, md_file.meta_set_record, rec_type, 1, "empty")
			self.assertRaises(IndexError, md_file.meta_erase_record, rec_type, 1)

			# Check that bad characters raise exceptions
			for bad_char in bad_char_list:
				bad_string = "My"+bad_char+"metadata"
				self.assertRaises(ValueError, md_file.meta_add_record, rec_type, bad_string)

		# Confirm that bad rec types will cause an error
		self.assertRaises(ValueError, md_file.meta_num_records, -1)
		self.assertRaises(ValueError, md_file.meta_clear_record, -1)
		self.assertRaises(ValueError, md_file.meta_get_record, -1, 0)
		self.assertRaises(ValueError, md_file.meta_add_record, -1, "empty")
		self.assertRaises(ValueError, md_file.meta_set_record, -1, 0, "empty")
		self.assertRaises(ValueError, md_file.meta_erase_record, -1, 0)

	# Ensure that reading/writing works for metadata
	def testWriteRead(self):
		md_name = "./test_metadata.dat"

		# Test metadata record adding, setting, and counts
		md_file = quakelib.EQSimConditionWriter()
		md_file.open(md_name)
		md_file.meta_add_record(quakelib.META_COMMENT, "Comment 1")
		md_file.meta_add_record(quakelib.META_COMMENT, "Comment 2")
		md_file.meta_add_record(quakelib.META_COMMENT, "Comment 3")
		md_file.meta_set_record(quakelib.META_COMMENT, 1, "NewComment 2")
		self.assertEqual(md_file.meta_num_records(quakelib.META_COMMENT), 3)

		md_file.meta_add_record(quakelib.META_SIGNATURE, "Sig_1")
		md_file.meta_add_record(quakelib.META_SIGNATURE, "Sig_2")
		self.assertEqual(md_file.meta_get_record(quakelib.META_SIGNATURE, 0), "Sig_2")
		md_file.meta_set_record(quakelib.META_SIGNATURE, 0, "Sig_3")
		self.assertEqual(md_file.meta_num_records(quakelib.META_SIGNATURE), 1)

		md_file.meta_add_record(quakelib.META_INFO, "Info 1")
		md_file.meta_add_record(quakelib.META_INFO, "Info 2")
		md_file.meta_add_record(quakelib.META_INFO, "Info 3")
		md_file.meta_set_record(quakelib.META_INFO, 1, "NewInfo 2")
		self.assertEqual(md_file.meta_num_records(quakelib.META_INFO), 3)

		md_file.meta_add_record(quakelib.META_TITLE, "Title 1")
		md_file.meta_add_record(quakelib.META_TITLE, "Title 2")
		self.assertEqual(md_file.meta_get_record(quakelib.META_TITLE, 0), "Title 2")
		md_file.meta_set_record(quakelib.META_TITLE, 0, "Title 3")
		self.assertEqual(md_file.meta_num_records(quakelib.META_TITLE), 1)

		md_file.meta_add_record(quakelib.META_AUTHOR, "Author 1")
		md_file.meta_add_record(quakelib.META_AUTHOR, "Author 2")
		md_file.meta_add_record(quakelib.META_AUTHOR, "Author 3")
		md_file.meta_set_record(quakelib.META_AUTHOR, 1, "NewAuthor 2")
		self.assertEqual(md_file.meta_num_records(quakelib.META_AUTHOR), 3)

		md_file.meta_add_record(quakelib.META_DATE, "Date 1")
		md_file.meta_add_record(quakelib.META_DATE, "Date 2")
		self.assertEqual(md_file.meta_get_record(quakelib.META_DATE, 0), "Date 2")
		md_file.meta_set_record(quakelib.META_DATE, 0, "Date 3")
		self.assertEqual(md_file.meta_num_records(quakelib.META_DATE), 1)

		err = quakelib.EQSimErrors()
		md_file.validate(err)
		self.assertEqual(err.count(), 0)

		md_file.write()
		md_file.close()

		# Read data back in to ensure validity, test erase function
		md_file_in = quakelib.EQSimConditionReader()
		md_file_in.parse_file(md_name)
		# TODO: check that there are no parse_errors
		#self.assertEqual(err.count(), 0)

		self.assertEqual(md_file_in.meta_get_record(quakelib.META_COMMENT, 0), "Comment 1")
		self.assertEqual(md_file_in.meta_get_record(quakelib.META_COMMENT, 1), "NewComment 2")
		self.assertEqual(md_file_in.meta_get_record(quakelib.META_COMMENT, 2), "Comment 3")
		md_file_in.meta_erase_record(quakelib.META_COMMENT, 1)
		self.assertEqual(md_file_in.meta_get_record(quakelib.META_COMMENT, 1), "Comment 3")
		self.assertEqual(md_file_in.meta_num_records(quakelib.META_COMMENT), 2)

		self.assertEqual(md_file_in.meta_get_record(quakelib.META_SIGNATURE, 0), "Sig_3")
		md_file_in.meta_erase_record(quakelib.META_SIGNATURE, 0)
		self.assertRaises(IndexError, md_file_in.meta_get_record, quakelib.META_SIGNATURE, 0)
		self.assertEqual(md_file_in.meta_num_records(quakelib.META_SIGNATURE), 0)

		self.assertEqual(md_file_in.meta_get_record(quakelib.META_INFO, 0), "Info 1")
		self.assertEqual(md_file_in.meta_get_record(quakelib.META_INFO, 1), "NewInfo 2")
		self.assertEqual(md_file_in.meta_get_record(quakelib.META_INFO, 2), "Info 3")
		md_file_in.meta_erase_record(quakelib.META_INFO, 1)
		self.assertEqual(md_file_in.meta_get_record(quakelib.META_INFO, 1), "Info 3")
		self.assertEqual(md_file_in.meta_num_records(quakelib.META_INFO), 2)

		self.assertEqual(md_file_in.meta_get_record(quakelib.META_TITLE, 0), "Title 3")
		md_file_in.meta_erase_record(quakelib.META_TITLE, 0)
		self.assertRaises(IndexError, md_file_in.meta_get_record, quakelib.META_TITLE, 0)
		self.assertEqual(md_file_in.meta_num_records(quakelib.META_TITLE), 0)

		self.assertEqual(md_file_in.meta_get_record(quakelib.META_AUTHOR, 0), "Author 1")
		self.assertEqual(md_file_in.meta_get_record(quakelib.META_AUTHOR, 1), "NewAuthor 2")
		self.assertEqual(md_file_in.meta_get_record(quakelib.META_AUTHOR, 2), "Author 3")
		md_file_in.meta_erase_record(quakelib.META_AUTHOR, 1)
		self.assertEqual(md_file_in.meta_get_record(quakelib.META_AUTHOR, 1), "Author 3")
		self.assertEqual(md_file_in.meta_num_records(quakelib.META_AUTHOR), 2)

		self.assertEqual(md_file_in.meta_get_record(quakelib.META_DATE, 0), "Date 3")
		md_file_in.meta_erase_record(quakelib.META_DATE, 0)
		self.assertRaises(IndexError, md_file_in.meta_get_record, quakelib.META_DATE, 0)
		self.assertEqual(md_file_in.meta_num_records(quakelib.META_DATE), 0)

		err = quakelib.EQSimErrors()
		md_file_in.validate(err)
		self.assertEqual(err.count(), 0)

		os.remove(md_name)

if __name__ == '__main__':
	unittest.main()

