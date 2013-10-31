import unittest, os, os.path

from waistcoat import command

DATA_DIR = os.path.join( os.path.split(__file__)[0], "data/")

class TestAvailable(unittest.TestCase):
	should_have = "ls"
	should_not_have = "A_COMMAND_WHICH_DOESNT_EXIST"

	def test_exists(self):
		"""Test that I know when a function exists"""
		self.assertTrue(command.isAvailable(self.should_have))

	def test_no_exists(self):
		"""Test that I know when a function doesn't exist"""
		self.assertFalse(command.isAvailable(self.should_not_have))

class TestCommand(command.Command):
	def __init__(self):
		self.cmd = 'sh'
		self.default_args = [os.path.join(DATA_DIR, 'return_args.sh'),]

class TestCommandCall(unittest.TestCase):
	"""Test the Command.call function"""
	fmt = "My arguments were \"{}\"\n"
	cmd = TestCommand()
	
	def test_none(self):
		"""Test only default arguments"""
		self.assertEqual(self.cmd.call()[0], self.fmt.format(''))
		
	def test_args(self):
		"""Test argument passing"""
		self.assertEqual(self.cmd.call(['-a', '-b', '--aba', 'something'])[0], 
				self.fmt.format('-a -b --aba something'))	

