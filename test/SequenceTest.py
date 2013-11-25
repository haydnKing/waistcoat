"""Define a class which adds assertSequence"""

import unittest
from Bio import SeqIO

def print_seqs(expected, actual):
	return ("Expected:\n" + "\n".join([str(x.seq) for x in expected]) +
			"\nActual:\n" + "\n".join([str(x.seq) for x in actual]))

class SequenceTest(unittest.TestCase):


	def assertSequences(self, expected_file, actual_file):
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
