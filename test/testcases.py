"""Define a class which adds assertSequence"""

import unittest, pysam, re
from Bio import SeqIO

def print_seqs(expected, actual):
	return ("Expected:\n" + "\n".join([str(x.seq) for x in expected]) +
			"\nActual:\n" + "\n".join([str(x.seq) for x in actual]))

class TestFastQ(unittest.TestCase):
	"""A class for testing fastQ files"""

	def assertFastQ(self, expected_file, actual_file):
		"""assert that contain all the same sequences"""
		expected = sorted(list(SeqIO.parse(expected_file, 'fastq')),
				key=lambda x: str(x.seq))
		actual = sorted(list(SeqIO.parse(actual_file, 'fastq')),
				key=lambda x: str(x.seq))

		self.assertEqual(len(expected), len(actual), 
				"Length mismatch in \'{}\' vs \'{}\'\n{}".format(
					expected_file, actual_file, print_seqs(expected, actual)))

		for i,(exp, act) in enumerate(zip(expected, actual)):
			self.assertEqual(str(exp.seq), str(act.seq),
				"Sequence mismatch in \'{}\' vs \'{}\'\n{}".format(
					expected_file, actual_file, print_seqs(expected,actual)))
			self.assertEqual(exp.id, act.id, 
				"IDs mismatch in \'{}\' vs \'{}\': \'{}\' != \'{}\'".format(
					expected_file, actual_file, exp.id, act.id))
			self.assertEqual(exp.name, act.name,
				"Name mismatch in \'{}\' vs \'{}\': \'{}\' != \'{}\'".format(
					expected_file, actual_file, exp.name, act.name))
			self.assertEqual(exp.annotations, act.annotations,
				"annotation mismatch in \'{}\' vs \'{}\':\n\t" +
				"exp: {}\n\tact: {}".format(
					expected_file, actual_file, exp.annotations, act.annotations))
			self.assertEqual(exp.letter_annotations, act.letter_annotations,
				"letter_annotation mismatch in \'{}\' vs \'{}\':\n\t" +
				"exp: {}\n\tact: {}".format(
					expected_file, actual_file, exp.letter_annotations, 
					act.letter_annotations))

class TestSamfile(unittest.TestCase):
	"""A class for testing samfiles"""

	def assertBAM(self, samfile, target):
		"""Check whether the samfile is valid"""
		self.assertSamfile(samfile, target, 'rb')

	def assertSAM(self, samfile, target):
		"""check whether the samfile is valid"""
		self.assertSamfile(samfile, target, 'r')

	def assertSamfile(self, samfile, target, mode='rb'):
		
		#load the targets into RAM
		targets = {}
		for seq in SeqIO.parse(target, 'fasta'):
			targets[seq.name] = seq


		track = pysam.Samfile(samfile, mode)
		for alg in track.fetch():
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
				
			rname = track.getrname(alg.tid)
			seq = targets[rname][int(pos)-1:int(pos)-1 + int(length)]
			self.assertEqual(alg.seq.upper(), str(seq.seq).upper(),
					"Sequence mismatch \'{}\' in \'{}\':\n\texpected: {}\n\t  actual: {}"
					.format(alg.qname, samfile, str(seq.seq).upper(), alg.seq.upper()))
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

