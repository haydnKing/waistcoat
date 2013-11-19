"""Define a class which adds assertSequence"""

import unittest
from Bio import SeqIO

def print_seqs(expected, actual):
	return ("Expected:\n" + "\n".join([str(x.seq) for x in expected]) +
			"\nActual:\n" + "\n".join([str(x.seq) for x in actual]))

class SequenceTest(unittest.TestCase):


	def assertSequences(self, expected_file, actual_file, desc='unknown'):
		"""assert that contain all the same sequences"""
		expected = sorted(list(SeqIO.parse(expected_file, 'fastq')),
				key=lambda x: str(x.seq))
		actual = sorted(list(SeqIO.parse(actual_file, 'fastq')),
				key=lambda x: str(x.seq))

		self.assertEqual(len(expected), len(actual), 
				"Length mismatch in file \'{}\' (desc=\'{}\')\n{}".format(
					actual_file, desc, print_seqs(expected, actual)))

		for i,(exp, act) in enumerate(zip(expected, actual)):
			self.assertEqual(str(exp.seq), str(act.seq),
				"Sequence mismatch in file \'{}\' (desc=\'{}\')\n{}".format(
					actual_file, desc, print_seqs(expected,actual)))
