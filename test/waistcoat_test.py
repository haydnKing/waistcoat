import unittest, os.path, tempfile, shutil, testcases
from waistcoat import waistcoat

waistcoat.verbose = False
waistcoat.check_output = False
waistcoat.preprocess.verbose = False

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data/pipeline_data/")

class PipelineTest(testcases.TestSamfile):
	"""Test the pipeline"""
	settings_file = os.path.join(DATA_DIR, 'settings.json')
	reads = os.path.join(DATA_DIR, 'reads.fq.gz')
	genome = os.path.join(DATA_DIR, 'genome.fa')

	@classmethod
	def setUpClass(self):
		self.tempdir = tempfile.mkdtemp(prefix='test_waistcoat.')
		try:
			waistcoat.run(self.settings_file, self.reads, self.tempdir, 
					temp_loc=self.tempdir, extend=True)
		except Exception:
			shutil.rmtree(self.tempdir)
			raise

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

		output = sorted(output)

		self.assertEqual(output, sorted([
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
				'sample 2/logs/tophat.log',
				'statistics/lengths.csv',
				'statistics/pipeline.csv',
				'sample 1.bam',
				'sample 1.bam.bai',
				'sample 2.bam',
				'sample 2.bam.bai']))

	def test_accepted_hits(self):
		"""Test the outputted hits"""
			
		self.assertBAM(os.path.join(self.tempdir, 'sample 1/accepted_hits.bam'), 
				self.genome)
		self.assertBAM(os.path.join(self.tempdir, 'sample 2/accepted_hits.bam'),
				self.genome)


