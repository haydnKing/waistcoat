"""
Read and validate settings for waistcoat
"""

import simplejson as json
import os.path, tophat

def loadf(fname):
	"""Read settings from a JSON formatted file"""
	try:
		return loads(file(fname).read())
	except SettingsError, e:
		e.message = "File {}: {}".format(fname, e.message)
		raise

def loads(string):
	"""Read settings from a JSON formatted string"""
	return loadd(json.loads(string))

def loadd(data):
	"""Read settings from a dictionary"""
	validate(data)
	return Settings(data)

class Settings(object):
	"""Store the settings fora tophat run"""

	def __init__(self, valid_data):
		self.barcode_format = valid_data['barcode_format']
		self.barcodes = valid_data['barcodes']

		self.discard = []
		if valid_data.has_key('discard'):
			for index in valid_data['discard']:
				self.discard.append((index, 
					valid_data.get(os.path.basename(index) + '_settings', {})),)

		self.target = (valid_data['target'], valid_data.get('target_settings',{}),)














def validate(data):
	"""validate the settings in data"""
	#test required keys
	try:
		validate_barcode_format(data['barcode_format'])
		validate_barcodes(data['barcodes'], data['barcode_format'])
		validate_map(data['target'])
	except KeyError, e:
		raise SettingsError(
			'Required key \'{}\' not found'.format(
				e.message))

	#test discard
	if data.has_key('discard'):
		validate_discard(data['discard'])
		validate_discard_settings(data)

	#test map settings
	if data.has_key('map_settings'):
		validate_tophat_settings('target_settings', 
			data['target_settings'])

def validate_tophat_settings(desc, settings):
	try:
		t = tophat.tophat_from_settings(settings)
	except (TypeError, ValueError), e:
		raise(SettingsError('\'{}\': {}'.format(desc, e.message)))

def validate_discard_settings(data):
	#test individual settings
	discard = data['discard']

	for i, path in enumerate(discard):
		index_name = os.path.basename(path)
		if data.has_key(index_name + '_settings'):
			validate_tophat_settings(index_name + '_settings', 
					data[index_name + '_settings'])
		#test global settings
		if data.has_key('discard_settings'):
			validate_tophat_settings('discard_settings', 
				data['discard_settings'])

		
def validate_map(map_to):
	if not index_exists(map_to):
		raise(SettingsError('\'map_to\' index \'{}\' does not exist'
			.format(map_to)))

def validate_discard(discard):
	
	if not isinstance(discard, (list, basestring)):
		raise SettingsError(
			'\'discard\' should be a string or list of strings, not {}'
				.format(type(discard)))
	else:
		if isinstance(discard, basestring):
			discard = [discard,]
		for i,path in enumerate(discard):
			if not isinstance(path, basestring):
				raise SettingsError(
				('\'discard\' should be a string or a list of strings, ' +
					'item {} is a {}').format(i,type(path)))
			else:
				if not index_exists(path):
					raise SettingsError(
							('could not find index \'{}\'').format(path))

def index_exists(path):
	index_fmt = [
			'{}.1.bt2', 
			'{}.2.bt2', 
			'{}.3.bt2', 
			'{}.4.bt2',
			'{}.rev.1.bt2',
			'{}.rev.2.bt2',]

	def index_in_dir(path, index):
		for fmt in index_fmt:
			if not os.path.exists(os.path.join(path,fmt.format(index))):
				return False
		return True

	(head, tail) = os.path.split(path)
	
	if head == '':
		head = '.'

	if os.path.exists(head):
		if index_in_dir(head, tail):
			return True
		if index_in_dir(os.path.join(head,'index/'), tail):
			return True
	
	p = os.environ.get('BOWTIE2_INDEXES')
	if os.path.dirname(path) == '' and p:
		if index_in_dir(p, tail):
			return True
	
	return False


def validate_barcode_format(barcode_fmt):
	"""Check that the given format is valid"""

	#format must be a string
	if not isinstance(barcode_fmt, basestring):
		#cannot recover from this
		raise SettingsError(
				'\'barcode_fmt\' should be a string, not {}'.format(
					type(barcode_fmt)))

	#format can only contain Bs and Ns
	for char in barcode_fmt.upper():
		if char not in "BN":
			raise SettingsError('\'barcode_fmt\' characters must be B or N')
	
def validate_barcodes(barcodes, barcode_fmt):
	length = sum([1 for b in barcode_fmt.upper() if b == 'B'])

	#barcodes must be a dictionary of strings
	if not isinstance(barcodes, dict):
		raise SettingsError(
				'\'barcodes\' should be a dictionary mapping samples to barcodes')
	
	for sample,barcode in barcodes.iteritems():
		if not isinstance(sample, basestring):
			raise SettingsError('\'barcodes\' sample name should be a string, not{}'
					.format(type(sample)))
		if not isinstance(barcode, basestring):
			raise SettingsError('\'barcodes\' barcode should be a string, not {}'
					.format(type(barcode)))

		#barcode strings can only contain ATGC
		for char in barcode.upper():
			if char not in "ATGC":
				raise SettingsError(
						'barcodes can only contain A,T,C or Gs, not {}'.format(char))
				break

		#check whether the barcode are the right length
		if len(barcode) != length:
			raise SettingsError(
					('barcodes must be the same length as the format' + 
					', barcode = {}, format	= {}'.format(len(barcode), length)))
			break

example_settings = """{
"_comment": "This is an example settings file for waistcoat. Comments begin with underscores",

"_barcode_format": "REQUIRED. Format of the barcode - B = Barcode character N = Random (UMI) character",
"barcode_format" : "BBBNNNNBB",

"_barcodes": "REQUIRED. A dictionary mapping sample names to barcodes",
"barcodes" : {
	"sample 1": "ACCTA",
	"sample 2": "GCGAT"
	},

"_discard": "OPTIONAL. Sequences which map to indexes listed here will be discarded",
"discard": [
	"path/to/discard_index_base1",
	"discard_index_base2"
	],

"_discard_settings": "OPTIONAL. Default settings for tophat when discarding",
"discard_settings": {
	"_comment": "settings go here as key value pairs, e.g. this sets --max-insertion-length 5",
	"max_insertion_length": 5
},

"_discard_index_base1_settings": "OPTIONAL. Override default discard settings for a specific index by setting [index_name]_settings",
"discard_index_base1_settings": {
	"_comment": "override settings here"
	},

"_target": "REQUIRED. Index to perform final map against", 
"target": "final_ref",

"_target_settings": "OPTIONAL. Settings for final mapping",
"target_map_settings": {
	"_comment": "Settings go here"
	}
}
"""

class SettingsError(ValueError):
	pass
