import unittest

from waistcoat import statistics

class StatisticsTest(unittest.TestCase):
	"""Test the statistics module"""

	def setUp(self):
		statistics.clear()

	def test_names(self):
		"""test that sample names must be correct"""
		statistics.setUp(['sample1', 'sample2'])
		self.assertRaises(KeyError, statistics.addValues, "test", {
			'sample1': 5,
			'sample2': 3,
			'unknownSample': 4,})

