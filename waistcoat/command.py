"""Class to execute an external command after checking arguments etc"""

import subprocess, os

class Command:
	"""A class which calls an external command and checks return values"""

	#the command to call
	cmd = ""

	def call(self, switches=[], args={}, update_fn=None):
		"""Call the function with the argument and check return codes"""

		print "Call: {}".format([self.cmd,] + self.getArgs(switches, args))
		p = subprocess.Popen([self.cmd,] + self.getArgs(switches, args),
							stdout = subprocess.PIPE, 
							stdin  = subprocess.PIPE,
							stderr = subprocess.PIPE)

		#wait for lines
		if hasattr(update_fn, '__call__'):
			while True:
				nextline = p.stdout.readline()
				if nextline == '' and p.poll() != None:
					break
				update_fn(nextline)

		#collect the output
		out = p.communicate()

		#check returncode
		if p.returncode < 0:
			#error
			raise subprocess.CalledProcessError(
				"Program \"{}\" returned error code {}".format(cmd, p.returncode),
				out[1])

		return out

	def getArgs(self, switches=[], args={}):
		"""Return the argument list"""
		def format(arg):
			if len(arg) > 1:
				return "--{}".format(arg)
			return "-{}".format(arg)

		if isinstance(switches, basestring):
			switches = [switches,]

		ret = [format(s) for s in switches]

		for k, v in args.iteritems():
			ret += [format(k), v]

		return ret
				
def isAvailable(command):
	try:
		subprocess.call([command,])
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			# handle file not found error.
			return False
		else:
			# Something else went wrong
				raise(e)
	return True

