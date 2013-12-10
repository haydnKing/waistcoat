import tempfile, unittest, shutil, testcases

from os.path import join as pjoin
from os.path import split as psplit
import os, shutil

from waistcoat import preprocess, settings, statistics

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#preprocess.verbose = False

DATA_DIR = pjoin(psplit(__file__)[0], "data/preprocess/")


class PreprocessTest(testcases.TestFastQ):

	def setUp(self):
		self.tempdir = tempfile.mkdtemp(prefix="test")
		statistics.recording = False

	def tearDown(self):
		shutil.rmtree(self.tempdir)
		statistics.recording = True

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
		self.assertFastQ(code1, files['barcode_1'])
		self.assertFastQ(code2, files['barcode_2'])

	def test_process_sample(self):
		"""Test process_sample"""
		input_file = pjoin(self.tempdir, 'process_test.fq')
		test_input = pjoin(DATA_DIR, 'process_test.fq')
		test_output = pjoin(DATA_DIR, 'process_test_out.fq')

		#copy test file to tempdir
		shutil.copyfile(test_input, input_file)

		s = settings.Settings({
			'barcodes': {"barcode_1": 'TCCA', "barcode_2": 'TCTT',},
			'barcode_format': "BBBNNNB",
			'target': 'null',})
		files = preprocess.process_sample(
				{'sample_1': input_file,}, s, self.tempdir)

		#check that the sample persisted
		self.assertEqual(files.keys(), ['sample_1',])
		output_file = files['sample_1']

		#check that the file was deleted
		self.assertFalse(os.path.exists(input_file), 
			"Test file {} was not deleted".format(input_file))

		#check that the correct file was produced
		self.assertTrue(os.path.exists(output_file),
			"Output file {} was not produced".format(output_file))

		#check that the contents of the output file are correct
		self.assertFastQ(test_output, output_file)

	def test_run(self):
		"""test preprocess.run"""
	
		reads = pjoin(DATA_DIR, 'run_in.fq')

		s = settings.Settings({
			'barcodes': {"barcode_1": 'TCCA', "barcode_2": 'TCTT',},
			'barcode_format': "BBBNNNB",
			'target': 'null',})
		
		files = preprocess.run(reads, s, self.tempdir)

		#check that the correct files were produced
		self.assertEqual(sorted(files.keys()), ['barcode_1','barcode_2',])
		self.assertEqual(sorted(files.values()), sorted(
			[pjoin(self.tempdir, p) for p in os.listdir(self.tempdir)]))

		#check that the sequences are correct
		self.assertFastQ(pjoin(DATA_DIR, 'run_out_bc1.fq'), files['barcode_1'])
		self.assertFastQ(pjoin(DATA_DIR, 'run_out_bc2.fq'), files['barcode_2'])
