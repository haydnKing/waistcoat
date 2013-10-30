"""Class to execute an external command after checking arguments etc"""

from subprocess import Popen, PIPE, STDOUT, call
import os

class Command:
	"""A class which calls an external command and checks return values"""

	#the command to call
	cmd = ""
	def __init__(self, command):
		self.cmd = command

	def call(self, switches=[], args={}, update_fn=None):
		"""Call the function with the argument and check return codes"""

		p = Popen([self.cmd, self.getArgs(switches, args)],
							stdout=PIPE, stdin=PIPE, stderr=PIPE)

		#wait for lines
		if hasattr(update_fn, '__call__'):
			while True:
				nextline = process.stdout.readline()
				if nextline == '' and p.poll() != None:
					break
				update_fn(nextline)

		#collect the output
		out = p.communicate()

	def getArgs(self, switches=[], args={}):
		"""Return the argument list"""
		def format(arg):
			if len(arg) > 1:
				return "--{}".format(arg)
			return "-{}".format(arg)

		ret = [format(s) for s in switches]

		for k, v in kwargs.iteritems():
			ret += [format(k), v]

		return ret
				
def isAvailable(command):
	try:
		call([command,])
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			# handle file not found error.
			return False
		else:
			# Something else went wrong
				raise(e)
	return True

