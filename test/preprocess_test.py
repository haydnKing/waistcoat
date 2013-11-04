import tempfile, unittest, shutil

from os.path import join as pjoin
from os.path import split as psplit
import os
from Bio import SeqIO

from waistcoat import preprocess


DATA_DIR = pjoin(psplit(__file__)[0], "data/")

class PreprocessTest(unittest.TestCase):
	reads = pjoin(DATA_DIR, 'test_reads.fq')
	code1 = pjoin(DATA_DIR, 'expected_1.fq')
	code2 = pjoin(DATA_DIR, 'expected_2.fq')
	unmapped = pjoin(DATA_DIR, 'expected_u.fq')

	barcodes = {"TCCA": "barcode_1", "TCTT": "barcode_2",}
	barcode_fmt = [0,1,2,6]

	expected_files = ["barcode_1.fq", "barcode_2.fq", "unmapped.fq",]

	def setUp(self):
		self.tempdir = tempfile.mkdtemp(prefix="test")

	def tearDown(self):
		shutil.rmtree(self.tempdir)

	def test_split(self):

		preprocess.split_by_barcode(self.reads,
				self.barcode_fmt,
				self.barcodes,
				self.tempdir)

		#check that the correct files were produced
		self.assertEqual(sorted(os.listdir(self.tempdir)), self.expected_files)

		#check all of the files
		self.assertSequences(self.code1, pjoin(self.tempdir, 'barcode_1.fq'))
		self.assertSequences(self.code2, pjoin(self.tempdir, 'barcode_2.fq'))
		self.assertSequences(self.unmapped, pjoin(self.tempdir, 'barcode_u.fq'))
		

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
				"Sequence mismatch in file \'{}\'\nExpected:\n{}\nActual:\n{}".format(
					actual_file, str(exp.seq), str(act.seq)))

def print_seqs(expected, actual):
	return ("Expected:\n" + "\n".join([str(x.seq) for x in expected]) +
			"\nActual:\n" + "\n".join([str(x.seq) for x in actual]))
