from waistcoat import tophat
import unittest

class TophatTest(unittest.TestCase):
	"""Test the Tophat runner"""

	def setUp(self):
		self.tophat = tophat.TopHat()

	def test_options(self):
		"""Test that TopHat.getOptions works"""

		self.assertEqual(self.tophat.getOptions(), [])

		self.tophat.read_gap_length = 15
		self.assertEqual(self.tophat.getOptions(), 
				['--read-gap-length', '15',])

		self.tophat.read_gap_length = "/path /to /a \"File\""
		self.assertEqual(self.tophat.getOptions(),
				['--read-gap-length', '\"/path /to /a \\\"File\\\"\"'])

		self.tophat.read_gap_length = True
		self.assertEqual(self.tophat.getOptions(), ['--read-gap-length',])

		self.tophat.read_gap_length = False
		self.assertEqual(self.tophat.getOptions(), [])
