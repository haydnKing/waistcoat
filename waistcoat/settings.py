"""
Read and validate settings for waistcoat
"""

import simplejson as json
import logging

def read(fname):
	"""Read settings from a JSON formatted file"""

	data = json.loads(file(fname).read())

	#test required keys
	try:
		validate_barcode_format(data['barcode_format'])
		validate_barcodes(data['barcodes'], data['barcode_format'])
		validate_map(data['tophat_map'])
	except KeyError, e:
		raise SettingsError(
			'Required key \'{}\' not found in settings file \"{}\"'.format(
				e.message, fname))

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
					data[index_name + '_settings')
		#test global settings
		if data.has_key('tophat_discard_settings'):
			validate_tophat_settings('tophat_discard_settings', 
				data['tophat_discard_settings'])

		

def validate_discard(discard):
	
	if not isinstance(discard, (list, basestring)):
		raise SettingsError(
			'\'tophat_discard\' should be a string or list of strings, not {}'
				.format(type(discard)))
	else:
		if isinstace(discard, basestring):
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
	
	if head = '':
		head = '.'

	if os.path.exists(head):
		if index_in_dir(head, tail):
			return True
		if index_in_dir(os.path.join(head,'index/'), tail):
			return True
	
	p = os.environ.get('BOWTIE2_INDEXES')
	if os.path.dirname(path) = '' and p:
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

	#barcodes must be a list of strings
	if isinstance(barcodes, basestring):
		raise SettingsError('\'barcodes\' should be a list of strings')
		
	for i,s in enumerate(barcodes):
		if not isinstance(s, basestring):
			raise SettingsError(
					'\'barcodes\'[{}] should be a string, not {}'.format(type(s)))
		#barcode strings can only contain ATGC
		for char in s.upper():
			if char not in "ATGC":
				raise SettingsError(
						'barcodes can only contain A,T,C or Gs, not {}'.format(char))
				break

		#check whether the barcodes are the right length
		for barcode in barcodes:
			if len(barcode) != length:
				raise SettingsError(
						'barcodes must be the same length as the format')
				break

example_settings = """{
"_comment": "This is an example settings file for waistcoat. Comments begin
with underscores",

"_barcode_fmt": "REQUIRED. Format of the barcode - 
												B = Barcode character
												N = Random (UMI) character",
"barcode_fmt" : "BBBNNNNBB",

"_barcodes": "REQUIRED. A dictionary mapping sample names to barcodes",
"barcodes" : {
	"sample 1": "ACCTA",
	"sample 2": "GCGAT"
	}

"_tophat_discard": "OPTIONAL. Sequences which map to indexes listed here will be discarded",
"tophat_discard": [
	"path/to/discard_index_base1",
	"discard_index_base2"
	],

"_tophat_discard_settings": "OPTIONAL. Default settings for tophat when
dicarding"
"tophat_discard_settings": {
	"_comment": "settings go here as key value pairs, e.g. this sets 
		--max-insertion-length 5",
	"max_insertion_length": 5
},

"_discard_index_base1_settings": "OPTIONAL. Override default discard settings
for a specific index by setting [index_name]_settings",
"discard_index_base1_settings": {
	"_comment": "override settings here"
	}

"_tophat_map": "REQUIRED. Index to perform final map against", 
"tophat_map": "final_index",

"_tophat_map_settings": "OPTIONAL. Settings for final mapping",
"tophat_map_settings": {
	"_comment": "Settings go here",
	}
}
"""



class SettingsError(ValueError):
	pass
