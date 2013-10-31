import unittest, os, os.path, subprocess

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

class ArgTester(command.Command):
	def __init__(self):
		self.cmd = 'sh'
		self.default_args = [os.path.join(DATA_DIR, 'return_args.sh'),]

class ReturnTester(command.Command):
	def __init__(self):
		self.cmd = 'sh'
		self.default_args = [os.path.join(DATA_DIR, 'exit_failure.sh'),]

class TestCommandCall(unittest.TestCase):
	"""Test the Command.call function"""
	fmt = "My arguments were \"{}\"\n"
	arg = ArgTester()
	ret = ReturnTester()
	
	def test_none(self):
		"""Test only default arguments"""
		self.assertEqual(self.arg.call()[0], self.fmt.format(''))
		
	def test_args(self):
		"""Test argument passing"""
		self.assertEqual(self.arg.call(['-a', '-b', '--aba', 'something'])[0], 
				self.fmt.format('-a -b --aba something'))	
	
	def test_failure(self):
		"""Test that a non-zero return value generates an exception"""
		with self.assertRaises(subprocess.CalledProcessError):
			self.ret.call()

