import unittest

from waistcoat import command

class TestAvailable(unittest.TestCase):
	should_have = "ls"
	should_not_have = "A_COMMAND_WHICH_DOESNT_EXIST"

	def test_exists(self):
		"""Test that I know when a function exists"""
		self.assertTrue(command.isAvailable(self.should_have))

	def test_no_exists(self):
		"""Test that I know when a function doesn't exist"""
		self.assertFalse(command.isAvailable(self.should_not_have))


