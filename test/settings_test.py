
import unittest, os.path
import simplejson as json
from waistcoat import settings

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data/")

class SettingsTest(unittest.TestCase):
	"""Test loading of settings"""
	data = os.path.join(DATA_DIR, 'settings/')
	valid = "valid.json"
	invalid = ["nobarcode.json",  
						 "nofmt.json",  
						 "wrongchars1.json",
						 "wrongchars.json",
						 "wronglength.json",
						 "wrongtype1.json",
						 "wrongtype.json",] 

	def assertRaisesMsg(self, msg, etype, func, *args, **kwargs):
		try:
			func(*args, **kwargs)
			self.fail(msg)
		except etype as e:
			if not isinstance(e, etype):
				raise
			pass

	def test_valid(self):
		"""Test loading valid settings"""
		output = settings.loadf(os.path.join(self.data, self.valid))

		self.assertEqual(output, (
			((0,1,2,6,7), (3,4,5)), 
			['ATGCA', 'CTACT', 'CTACG',],))

	def test_invalid(self):
		for name in self.invalid:
			f = os.path.join(DATA_DIR, 'settings/', name)
			self.assertRaisesMsg("Failed to raise ValueError in file {}".format(name),
					ValueError, settings.loadf, f)

