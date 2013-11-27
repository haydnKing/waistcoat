import unittest, testcases, tempfile, shutil, os.path
from waistcoat import postprocess

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data/postprocess/")

class TestPostprocess(testcases.TestSamfile):

	def setUp(self):
		self.tempdir = tempfile.mkdtemp(prefix='postprocesstest')

	def tearDown(self):
		shutil.rmtree(self.tempdir)

	def test_accepted_hits(self):
		"""Test the postprocessing"""
		shutil.copy(os.path.join(DATA_DIR, 'accepted_hits.bam'), self.tempdir)
		shutil.copy(os.path.join(DATA_DIR, 'accepted_hits.bam.bai'), self.tempdir)

		postprocess.run(os.path.join(self.tempdir, 'accepted_hits.bam'),
				os.path.join(DATA_DIR, 'genome.fa'))
		
		
		output = [
				'GTCCGTAGTCCTAGTCGTCATCCCCGTA',
				'ACTGGACTATTTAGGACGATCGGACTGA',
			]

		self.assertBAM(os.path.join(self.tempdir, 'accepted_hits.bam'), 
				os.path.join(DATA_DIR, 'genome.fa'))

