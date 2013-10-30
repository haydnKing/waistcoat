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
		self.cmd = os.path.join(DATA_DIR, 'return_args.sh')

class TestCommandCall(unittest.TestCase):
	"""Test the Command.call function"""
	fmt = "My arguments were \"{}\"\\n"
	cmd = TestCommand()
	
	def test_none(self):
		print self.cmd.cmd
		self.assertEqual(self.cmd.call()[0], fmt.format(''))
		
	def test_switches(self):
		self.assertEqual(self.cmd.call(['a', 'b', 'aba'])[0], 
				fmt.format('-a -b --aba'))	
		self.assertEqual(self.cmd.call(switches = ['a', 'b', 'aba'])[0], 
				fmt.format('-a -b --aba'))	
		self.assertEqual(self.cmd.call('a')[0], fmt.format('-a'))	

	def test_args(self):
		self.assertEqual(self.cmd.call(args = {'a': 'one', 'b': 'two', 'aba': 'three'})[0], 
				fmt.format('-a one -b two --aba three'))	
