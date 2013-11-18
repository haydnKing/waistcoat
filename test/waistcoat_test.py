import unittest, os.path, tempfile, shutil, pysam
from waistcoat import waistcoat

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data/")

class PipelineTest(unittest.TestCase):
	"""Test the pipeline"""
	settings_file = os.path.join(DATA_DIR, 'pipeline_data/settings.json')
	reads = os.path.join(DATA_DIR, 'pipeline_data/reads.fq')


	@classmethod
	def setUpClass(self):
		self.tempdir = tempfile.mkdtemp(prefix='test')
		waistcoat.run(self.settings_file, self.reads, self.tempdir)

	@classmethod
	def tearDownClass(self):
		shutil.rmtree(self.tempdir)
		self.tempdir = None

	def test_output(self):
		"""Test the pipeline outputs the correct files"""
		print "Output: {}".format(', '.join(os.listdir(self.tempdir)))
		self.assertTrue(os.path.exists(os.path.join(self.tempdir, 'sample_1.bam')), 
			"Sample 1 output file missing")
		self.assertTrue(os.path.exists(os.path.join(self.tempdir, 'sample_2.bam')), 
			"Sample 2 output file missing")

		raise AssertionError("Haven't written this yet")	

