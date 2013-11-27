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
		output = []
		for path, subdirs, files in os.walk(self.tempdir):
			for name in sorted(files):
				output.append(os.path.relpath(os.path.join(path,name), self.tempdir))

		self.assertEqual(output, [
				'sample 1/accepted_hits.bam',
				'sample 1/accepted_hits.bam.bai',
				'sample 1/deletions.bed',
				'sample 1/insertions.bed',
				'sample 1/junctions.bed',
				'sample 1/prep_reads.info',
				'sample 1/logs/bowtie.left_kept_reads.log',
				'sample 1/logs/prep_reads.log',
				'sample 1/logs/reports.log',
				'sample 1/logs/reports.samtools_sort.log0',
				'sample 1/logs/run.log',
				'sample 1/logs/tophat.log',
				'sample 2/accepted_hits.bam',
				'sample 2/accepted_hits.bam.bai',
				'sample 2/deletions.bed',
				'sample 2/insertions.bed',
				'sample 2/junctions.bed',
				'sample 2/prep_reads.info',
				'sample 2/logs/bowtie.left_kept_reads.log',
				'sample 2/logs/prep_reads.log',
				'sample 2/logs/reports.log',
				'sample 2/logs/reports.samtools_sort.log0',
				'sample 2/logs/run.log',
				'sample 2/logs/tophat.log'])

	def test_accepted_hits(self):
		"""Test the outputted hits"""
		
		output = {
				'sample 1': [
					'ACTACTATCTGACTAGACTGGAGGCGCT',
					'GCTCGACGCTCAGCCGTAGCGCCGCGCG',
					],
				'sample 2': [
					'ACTGGACTATTTAGGACGATCGGACTGA',
					'TACTGGACTATTTAGGACGATCGGACTG',
					],
				}

		track1 = pysam.Samfile(
				os.path.join(self.tempdir, 'sample 1/accepted_hits.bam'), 'rb')
		for alg in track1.fetch():
			print alg.name

		track2 = pysam.Samfile(
				os.path.join(self.tempdir, 'sample 2/accepted_hits.bam'), 'rb')
		for alg in track2.fetch():
			print alg.name

		raise AssertionError("Haven't written this yet")	

