import os, os.path, tempfile, unittest, shutil

from Bio import SeqIO

from waistcoat import preprocess


DATA_DIR = os.path.join( os.path.split(__file__)[0], "data/")

def check_sequences(file1, file2):
	"""Return true if file1 and file2 contain all the same sequences"""
	f1 = sorted(list(SeqIO.parse(file1, 'fastq')),
			key=lambda x: x.seq)
	f2 = sorted(list(SeqIO.parse(file2, 'fastq')),
			key=lambda x: x.seq)

	if len(f1) != len(f2):
		return False

	for (a,b) in zip(f1,f2):
		if a != b:
			return False

	return True
	

class PreprocessTest(unittest.TestCase):
	reads = os.path.join(DATA_DIR, 'test_reads.fq')
	code1 = os.path.join(DATA_DIR, 'expected_1.fq')
	code2 = os.path.join(DATA_DIR, 'expected_2.fq')
	unmapped = os.path.join(DATA_DIR, 'expected_u.fq')

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
		self.assertTrue(self.code1, os.path.join(self.tempdir, 'barcode_1.fq'))
		self.assertTrue(self.code2, os.path.join(self.tempdir, 'barcode_2.fq'))
		self.assertTrue(self.unmapped, os.path.join(self.tempdir, 'barcode_u.fq'))
		



