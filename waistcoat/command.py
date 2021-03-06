"""Class to execute an external command after checking arguments etc"""

import subprocess, os, shlex

class Command:
	"""A class which calls an external command and checks return values"""

	#the command to call
	cmd = ""
	default_args = []

	def call(self, args=[], update_fn=None, stderr=False):
		"""Call the function with the argument and check return codes
		update_fn: function to update with each line of output. must not block.
		stderr: if true, listen to stderr instead of stdout
		"""

		if isinstance(args, basestring):
			args = shlex.split(args)

		p = subprocess.Popen([self.cmd,] + self.default_args + args,
							stdout = subprocess.PIPE, 
							stdin  = subprocess.PIPE,
							stderr = subprocess.PIPE)

		#wait for lines
		if update_fn:
			while True:
				if stderr:
					nextline = p.stderr.readline()
				else:
					nextline = p.stdout.readline()
				if nextline == '' and p.poll() != None:
					break
				update_fn(nextline.rstrip())

		#collect the output
		out = p.communicate()

		#check returncode
		if p.returncode != 0:
			#error
			raise subprocess.CalledProcessError(
				p.returncode,
				self.cmd,
				output = out[1])

		return out
		
def isAvailable(command):
	try:
		subprocess.call([command,],
							stdout = subprocess.PIPE, 
							stdin  = subprocess.PIPE,
							stderr = subprocess.PIPE)
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			# handle file not found error.
			return False
		else:
			# Something else went wrong
				raise(e)
	return True

