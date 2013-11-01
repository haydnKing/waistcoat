"""Interface with the tophat program"""

import command

class TopHat(command.Command):
	"""Class to which interfaces with tophat"""
	def __init__(self):
		self.cmd = "tophat"


def test_programs():
	commands = {'tophat': False, 
		    'bowtie2': False, 
		    'bowtie2-inspect': False, 
		    'bowtie2-build': False, 
		    'samtools': False,}
	for cmd in commands.keys():
		commands[cmd] = command.isAvailable(cmd)

	errors = [cmd for cmd,available in commands.iteritems() if not available]

	if errors:
		raise ImportError("Required program(s) \'{}\' were not found.".format(
			"\', \'".join(errors)))

test_programs()
		
	
