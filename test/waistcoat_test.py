import unittest, os.path, tempfile, shutil, pysam, re
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
					'GTCCGTAGTCCTAGTCGTCATCCCCGTA',
					'ACTGGACTATTTAGGACGATCGGACTGA',
					],
				}

		self.validate_samfile('sample 1/accepted_hits.bam', output['sample 1'])
		self.validate_samfile('sample 2/accepted_hits.bam', output['sample 2'])

	def validate_samfile(self, samfile, expected_seqs):
		"""Check whether the samfile is valid"""

		track = pysam.Samfile(
				os.path.join(self.tempdir, samfile), 'rb')
		for alg,exp in zip(track.fetch(), expected_seqs):
			(name, pos, length, direction) = re.match('(\w+)_(\d+)_(\d+)_([FR])', 
																									alg.qname).groups()
			self.assertEqual(name, 'genome', 
					"Alignment \'{}\' not expected in \'{}\'".format(alg.qname,samfile))
			if direction == 'F':
				self.assertFalse(alg.is_reverse, 
					"Alignment \'{}\' mapped to wrong strand in \'{}\'".format(
						alg.qname, samfile))
			else:
				self.assertTrue(alg.is_reverse, 
					"Alignment \'{}\' mapped to wrong strand in \'{}\'".format(
						alg.qname, samfile))
				exp = reverse_complement(exp)
				
			self.assertEqual(alg.seq.upper(), exp,
					"Sequence mismatch \'{}\' in \'{}\':\n\texpected: {}\n\t  actual: {}"
					.format(alg.qname, samfile, exp, alg.seq.upper()))
			self.assertEqual(alg.pos+1, int(pos),
					"Position mismatch: \'{}\' expected at {} but matched {} in \'{}\'"
					.format(alg.qname, int(pos), alg.pos+1, samfile))

def reverse_complement(seq):
	l = {
			'A': 'T',
			'T': 'A',
			'G': 'C',
			'C': 'G',
			}
	ret = []
	for b in reversed(seq.upper()):
		ret.append(l[b])
	return ''.join(ret)

