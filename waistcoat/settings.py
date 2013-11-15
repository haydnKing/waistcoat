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

	#test required keys
	try:
		validate_barcode_format(data['barcode_format'])
		validate_barcodes(data['barcodes'], data['barcode_format'])
		validate_map(data['tophat_map'])
	except KeyError, e:
		raise SettingsError(
			'Required key \'{}\' not found'.format(
				e.message))

	#test discard
	if data.has_key('tophat_discard'):
		validate_discard(data['tophat_discard'])
		validate_discard_settings(data)

	#test map settings
	if data.has_key('tophat_map_settings'):
		validate_tophat_settings('tophat_map_settings', 
			data['tophat_map_settings'])

def validate_tophat_settings(desc, settings):
	try:
		t = tophat.tophat_from_settings(settings)
	except (TypeError, ValueError), e:
		raise(SettingsError('\'{}\': {}'.format(desc, e.message)))

def validate_discard_settings(data):
	#test individual settings
	discard = data['tophat_discard']

	for i, path in enumerate(discard):
		index_name = os.path.basename(path)
		if data.has_key(index_name + '_settings'):
			validate_tophat_settings(index_name + '_settings', 
					data[index_name + '_settings'])
		#test global settings
		if data.has_key('tophat_discard_settings'):
			validate_tophat_settings('tophat_discard_settings', 
				data['tophat_discard_settings'])

		
def validate_map(map_to):
	if not index_exists(map_to):
		raise(SettingsError('\'map_to\' index \'{}\' does not exist'
			.format(map_to)))

def validate_discard(discard):
	
	if not isinstance(discard, (list, basestring)):
		raise SettingsError(
			'\'tophat_discard\' should be a string or list of strings, not {}'
				.format(type(discard)))
	else:
		if isinstance(discard, basestring):
			discard = [discard,]
		for i,path in enumerate(discard):
			if not isinstance(path, basestring):
				raise SettingsError(
				('\'tophat_discart\' should be a string or a list of strings, ' +
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

"_tophat_discard": "OPTIONAL. Sequences which map to indexes listed here will be discarded",
"tophat_discard": [
	"path/to/discard_index_base1",
	"discard_index_base2"
	],

"_tophat_discard_settings": "OPTIONAL. Default settings for tophat when discarding",
"tophat_discard_settings": {
	"_comment": "settings go here as key value pairs, e.g. this sets --max-insertion-length 5",
	"max_insertion_length": 5
},

"_discard_index_base1_settings": "OPTIONAL. Override default discard settings for a specific index by setting [index_name]_settings",
"discard_index_base1_settings": {
	"_comment": "override settings here"
	},

"_tophat_map": "REQUIRED. Index to perform final map against", 
"tophat_map": "final_ref",

"_tophat_map_settings": "OPTIONAL. Settings for final mapping",
"tophat_map_settings": {
	"_comment": "Settings go here"
	}
}
"""



class SettingsError(ValueError):
	pass
