import unittest, testcases, tempfile, shutil, os.path
from waistcoat import postprocess, statistics

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data/postprocess/")

class TestPostprocess(testcases.TestSamfile):

	def setUp(self):
		self.tempdir = tempfile.mkdtemp(prefix='postprocesstest')
		statistics.recording = False

	def tearDown(self):
		shutil.rmtree(self.tempdir)
		statistics.recording = True

	def test_accepted_hits(self):
		"""Test the postprocessing"""
		sample = "sampleName"
		tophatDir = os.path.join(self.tempdir, sample)
		os.mkdir(tophatDir)

		shutil.copy(os.path.join(DATA_DIR, 'accepted_hits.bam'), tophatDir)
		shutil.copy(os.path.join(DATA_DIR, 'accepted_hits.bam.bai'), tophatDir)

		postprocess.run(self.tempdir, sample,
				os.path.join(DATA_DIR, 'genome.fa'),
				extend=True)
		
		
		output = [
				'GTCCGTAGTCCTAGTCGTCATCCCCGTA',
				'ACTGGACTATTTAGGACGATCGGACTGA',
			]

		self.assertBAM(os.path.join(self.tempdir, '{}.bam'.format(sample)), 
				os.path.join(DATA_DIR, 'genome.fa'))

