from waistcoat import tophat
import unittest, os, os.path, tempfile, shutil 

DATA_DIR = os.path.join( os.path.split(__file__)[0], "data/")

class TophatTest(unittest.TestCase):
	"""Test the Tophat runner"""

	def setUp(self):
		self.tophat = tophat.TopHat()
		self.output = tempfile.mkdtemp()

	def tearDown(self):
		#Delete tmpdir
		shutil.rmtree(self.output)

	def test_options(self):
		"""Test that TopHat.getOptions works"""

		self.assertEqual(self.tophat.getOptions(), [])

		self.tophat.read_gap_length = 15
		self.assertEqual(self.tophat.getOptions(), 
				['--read-gap-length', '15',])

		self.tophat.read_gap_length = "/path/to/a File\'s Location.ext"
		self.assertEqual(self.tophat.getOptions(),
				['--read-gap-length', "/path/to/a File\'s Location.ext"])

		self.tophat.read_gap_length = True
		self.assertEqual(self.tophat.getOptions(), ['--read-gap-length',])

		self.tophat.read_gap_length = False
		self.assertEqual(self.tophat.getOptions(), [])

	def test_run(self):
		"""test that I can execute tophat on test data"""
		
		data = os.path.join(DATA_DIR, 'tophat_data/')
		
		expected_files = ['accepted_hits.sam',
					'deletions.bed',
					'insertions.bed',
					'junctions.bed',
					'logs',
					'prep_reads.info',
					'unmapped.bam']


		self.tophat.mate_inner_dist = 20
		self.tophat.no_convert_bam = True
		self.tophat.output_dir = self.output

		self.tophat.run(
				os.path.join(data, 'reads_1.fq'),
				os.path.join(data, 'reads_2.fq'),
				os.path.join(data, 'test_ref'))

		files = sorted(os.listdir(self.output))

		self.assertEqual(files, expected_files)


		
