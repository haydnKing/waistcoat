"""Class to execute an external command after checking arguments etc"""

from subprocess import Popen, PIPE, STDOUT

class Command:

	#the command to call
	cmd = ""

	def call(self, switches=[], args={}, update_fn=None):
		"""Call the function with the argument and check return codes"""

		p = Popen([self.cmd, self.getArgs(switches, args)],
							stdout=PIPE, stdin=PIPE, stderr=PIPE)

		#wait for lines
		if hasattr(update_fn, '__call__'):
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
				

