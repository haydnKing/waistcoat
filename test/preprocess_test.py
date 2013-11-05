import tempfile, unittest, shutil

from os.path import join as pjoin
from os.path import split as psplit
import os, shutil
from Bio import SeqIO

from waistcoat import preprocess

preprocess.verbose = False

DATA_DIR = pjoin(psplit(__file__)[0], "data/")

class PreprocessTest(unittest.TestCase):

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

	def assertSequences(self, expected_file, actual_file):
		"""assert that contain all the same sequences"""
		expected = sorted(list(SeqIO.parse(expected_file, 'fastq')),
				key=lambda x: str(x.seq))
		actual = sorted(list(SeqIO.parse(actual_file, 'fastq')),
				key=lambda x: str(x.seq))

		self.assertEqual(len(expected), len(actual), 
				"Length mismatch in file \'{}\'.\n{}".format(
					actual_file, print_seqs(expected,
						actual)))

		for i,(exp, act) in enumerate(zip(expected, actual)):
			self.assertEqual(str(exp.seq), str(act.seq),
				"Sequence mismatch in file \'{}\'\n{}".format(
					actual_file, print_seqs(expected,actual)))

def print_seqs(expected, actual):
	return ("Expected:\n" + "\n".join([str(x.seq) for x in expected]) +
			"\nActual:\n" + "\n".join([str(x.seq) for x in actual]))
