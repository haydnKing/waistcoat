
import unittest, os.path
import simplejson as json
from waistcoat import settings
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data/")

class SettingsTest(unittest.TestCase):
	"""Test loading of settings"""
	data = os.path.join(DATA_DIR, 'settings/')
	valid = "valid.json"
	invalid = ['badbarcode.json','badtophat1.json','badtophat3.json',
			'nobarcodes.json','noindex2.json','settingsnodiscard.json',
			'badbarcodes.json','badtophat2.json','nobarcode.json','noindex1.json',
			'notarget.json','unknownsetting.json','wronglength.json',] 

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
		mySettings = settings.loadf(os.path.join(self.data, self.valid))

		self.assertEqual(mySettings.barcode_format, "BBBNNNNBB")
		self.assertEqual(mySettings.barcodes, {
        "sample 1": "ACCTA",
        "sample 2": "GCGAT"
        })
		self.assertEqual(mySettings.discard, [
			('../tophat_data/test_ref'    , {'max_insertion_length':4,},),
			('../tophat_data/test_ref_two', {'max_insertion_length':5,},),
			])
		self.assertEqual(mySettings.target, 
				('../tophat_data/test_ref', {'max_insertion_length':3,}))

	def test_invalid(self):
		for name in self.invalid:
			f = os.path.join(DATA_DIR, 'settings/', name)
			self.assertRaisesMsg("Failed to raise ValueError in file {}".format(name),
					ValueError, settings.loadf, f)


class SettingsObjectTest(unittest.TestCase):
	"""Test settings object"""
	data = os.path.join(DATA_DIR, 'settings/valid.json')
	record = SeqRecord(Seq('ACTGCTATCTGATC'))

	def setUp(self):
		self.settings = settings.loadf(self.data)

	def test_barcode(self):
		"""Test the parsing of barcode/UMI"""
		self.assertEqual(self.settings.parse_barcode(self.record),
				('ACTTC', 'GCTA',))

	def test_strip(self):
		"""Test removal of barcode/UMI"""
		self.assertEqual(str(self.settings.strip_barcode(self.record).seq),
				str(self.record[9:].seq))

