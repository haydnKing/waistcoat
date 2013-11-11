import tempfile, unittest, shutil, SequenceTest

from os.path import join as pjoin
from os.path import split as psplit
import os, shutil

from waistcoat import preprocess

preprocess.verbose = False

DATA_DIR = pjoin(psplit(__file__)[0], "data/")

class PreprocessTest(SequenceTest.SequenceTest):

	def setUp(self):
		self.tempdir = tempfile.mkdtemp(prefix="test")

	def tearDown(self):
		shutil.rmtree(self.tempdir)

	def test_split(self):
		reads = pjoin(DATA_DIR, 'test_reads.fq')
		code1 = pjoin(DATA_DIR, 'expected_1.fq')
		code2 = pjoin(DATA_DIR, 'expected_2.fq')
		unmapped = pjoin(DATA_DIR, 'expected_u.fq')

		barcodes = {"TCCA": "barcode_1", "TCTT": "barcode_2",}
		barcode_fmt = [0,1,2,6]

		expected_files = ["barcode_1.fq", "barcode_2.fq", "unmapped.fq",]

		preprocess.split_by_barcode(reads,
				barcode_fmt,
				barcodes,
				self.tempdir)

		#check that the correct files were produced
		self.assertEqual(sorted(os.listdir(self.tempdir)), expected_files)

		#check all of the files
		self.assertSequences(code1, pjoin(self.tempdir, 'barcode_1.fq'))
		self.assertSequences(code2, pjoin(self.tempdir, 'barcode_2.fq'))
		self.assertSequences(unmapped, pjoin(self.tempdir, 'unmapped.fq'))
		
	def test_clean(self):
		"""Test clean_distributions"""
		barcode_fmt = [0,1,2,6]
		input_file = pjoin(self.tempdir, 'clean_test.fq')
		output_file = pjoin(self.tempdir, 'clean_test_nopolyA.fq')
		test_input = pjoin(DATA_DIR, 'clean_test.fq')
		test_output = pjoin(DATA_DIR, 'clean_test_out.fq')

		#copy test file to tempdir
		shutil.copyfile(test_input, input_file)

		preprocess.clean_distributions([input_file,], barcode_fmt)

		#check that the file was deleted
		self.assertFalse(os.path.exists(input_file), 
			"Test file {} was not deleted".format(input_file))

		#check that the correct file was produced
		self.assertTrue(os.path.exists(output_file),
			"Output file {} was not produced".format(output_file))

		#check that the contents of the output file are correct
		self.assertSequences(test_output, output_file)


class SettingsTest(unittest.TestCase):
	"""Test loading of settings"""
	valid = "valid.json"
	invalid = ["nobarcode.json",  
						 "nofmt.json",  
						 "wrongchars1.json",
						 "wrongchars.json",
						 "wronglength.json",
						 "wrongtype1.json",
						 "wrongtype.json",] 

	def assertRaisesMsg(self, msg, etype, func, *args, **kwargs):
		try:
			func(*args, **kwargs)
			self.fail(msg)
		except etype as e:
			if not isinstance(e, etype):
				raise
			pass

	def test_valid(self):
		"""Test loading valid settings"""
		f = pjoin(DATA_DIR, 'settings/', self.valid)
		output = preprocess.read_settings(f)

		self.assertEqual(output, (
			((0,1,2,6,7), (3,4,5)), 
			['ATGCA', 'CTACT', 'CTACG',],))

	def test_invalid(self):
		for name in self.invalid:
			f = pjoin(DATA_DIR, 'settings/', name)
			self.assertRaisesMsg("Failed to raise ValueError in file {}".format(name),
					ValueError, preprocess.read_settings, f)

