import tempfile, unittest, shutil, SequenceTest

from os.path import join as pjoin
from os.path import split as psplit
import os, shutil

from waistcoat import preprocess, settings

preprocess.verbose = False

DATA_DIR = pjoin(psplit(__file__)[0], "data/preprocess/")


class PreprocessTest(SequenceTest.SequenceTest):

	def setUp(self):
		self.tempdir = tempfile.mkdtemp(prefix="test")

	def tearDown(self):
		shutil.rmtree(self.tempdir)

	def test_split(self):
		reads = pjoin(DATA_DIR, 'test_reads.fq')
		code1 = pjoin(DATA_DIR, 'expected_1.fq')
		code2 = pjoin(DATA_DIR, 'expected_2.fq')

		s = settings.Settings({
			'barcodes': {"barcode_1": 'TCCA', "barcode_2": 'TCTT',},
			'barcode_format': "BBBNNNB",
			'target': 'null',})
		
		files = preprocess.split_by_barcode(reads, s, self.tempdir)

		#check that the correct files were produced
		self.assertEqual([os.path.join(self.tempdir, f) for f in 
				sorted(os.listdir(self.tempdir))],
			sorted(files.values()))

		#check that the correct samples were produced
		self.assertEqual(sorted(files.keys()), sorted(s.barcodes.keys()))

		#check all of the files
		self.assertSequences(code1, files['barcode_1'])
		self.assertSequences(code2, files['barcode_2'])
		
	def test_clean(self):
		"""Test clean_distributions"""
		barcode_fmt = [0,1,2,6]
		input_file = pjoin(self.tempdir, 'clean_test.fq')
		output_file = pjoin(self.tempdir, 'clean_test_nopolyA.fq')
		test_input = pjoin(DATA_DIR, 'clean_test.fq')
		test_output = pjoin(DATA_DIR, 'clean_test_out.fq')

		#copy test file to tempdir
		shutil.copyfile(test_input, input_file)

		s = settings.Settings({
			'barcodes': {"barcode_1": 'TCCA', "barcode_2": 'TCTT',},
			'barcode_format': "BBBNNNB",
			'target': 'null',})
		preprocess.clean_distributions([input_file,], s)

		#check that the file was deleted
		self.assertFalse(os.path.exists(input_file), 
			"Test file {} was not deleted".format(input_file))

		#check that the correct file was produced
		self.assertTrue(os.path.exists(output_file),
			"Output file {} was not produced".format(output_file))

		#check that the contents of the output file are correct
		self.assertSequences(test_output, output_file)


