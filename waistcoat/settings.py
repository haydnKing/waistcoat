"""
Read and validate settings for waistcoat
"""

import simplejson as json
import os.path, tophat

def loadf(fname):
	"""Read settings from a JSON formatted file"""
	try:
		data = json.loads(file(fname).read())
		validate(data, fname)
		return Settings(data)
	except SettingsError, e:
		e.message = "File {}: {}".format(fname, e.message)
		raise

def loads(string, fname=None):
	"""Read settings from a JSON formatted string"""
	return loadd(json.loads(string), fname)

def loadd(data, fname=None):
	"""Read settings from a dictionary"""
	validate(data)
	return Settings(data)

class Settings(object):
	"""Store the settings fora tophat run"""

	def __init__(self, valid_data):
		fname = valid_data.get('_filename', './settings.json')

		self.barcode_format = valid_data['barcode_format'].upper()
		self.barcodes = valid_data['barcodes']

		self.discard = []
		if valid_data.has_key('discard'):
			#for each discard index
			for index in valid_data['discard']:
				#find the correct settings
				discard_settings = {}
				if valid_data.has_key(os.path.basename(index) + '_settings'):
					discard_settings = valid_data[os.path.basename(index) + '_settings']
				elif valid_data.has_key('discard_settings'):
					discard_settings = valid_data['discard_settings']

				#add to list
				self.discard.append((get_index(index, fname), discard_settings,))

		self.target = (get_index(valid_data['target'], fname), 
				valid_data.get('target_settings',{}),)


	def strip_header(self, record):
		"""Remove the barcode from the record"""
		return record[len(self.barcode_format):]


	def parse_header(self, record):
		"""Return a tuple of (barcode, UMI) as strings"""
		return (self.parse_barcode(record), self.parse_UMI(record))

	def parse_barcode(self, record):
		"""return the sequence barcode"""
		return self._parse(record, 'B')

	def parse_UMI(self, record):
		"""return the sequence UMI"""
		return self._parse(record, 'N')

	def UMI_len(self):
		"""length of the UMI"""
		return len([x for x in self.barcode_format if x == 'N'])

	def barcode_len(self):
		"""length of the UMI"""
		return len([x for x in self.barcode_format if x == 'B'])

	def _parse(self, record, the_code='B'):
		ret = []
		for base,code in zip(str(record.seq).upper(), self.barcode_format):
			if code == the_code:
				ret.append(base)
		return ''.join(ret).upper()




def validate(data, fname='./settings.json'):
	"""validate the settings in data, fname = settings file for relative paths"""
	data['_filename'] = fname

	#test required keys
	try:
		validate_barcode_format(data['barcode_format'])
		validate_barcodes(data['barcodes'], data['barcode_format'])
		validate_map(data['target'], fname)
	except KeyError, e:
		raise SettingsError(
			'Required key \'{}\' not found'.format(
				e.message))

	#test discard
	if data.has_key('discard'):
		validate_discard(data['discard'], fname)
		validate_discard_settings(data)
	elif data.has_key('discard_settings'):
		raise SettingsError('\'discard_settings\' given but no \'discard\'')

	#test map settings
	if data.has_key('target_settings'):
		validate_tophat_settings('target_settings', 
			data['target_settings'])

	#check for unused settings
	known_settings = ['barcode_format', 'barcodes', 'target', 'target_settings',
			'discard', 'discard_settings',]
	known_settings += [os.path.basename(index) + '_settings' for index in 
			data.get('discard', [])]
	for setting in data.iterkeys():
		if not isinstance(setting, basestring):
			raise SettingsError('settings must be strings, not {}'.format(
				type(setting)))
		if setting[0] == '_' or setting in known_settings:
			continue
		else:
			raise SettingsError('unknown key \'{}\''.format(setting))

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

		
def validate_map(map_to, fname):
	if not index_exists(map_to, fname):
		raise(SettingsError('\'map_to\' index \'{}\' does not exist'
			.format(map_to)))

def validate_discard(discard, fname):
	
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
				if not index_exists(path, fname):
					raise SettingsError(
							('could not find index \'{}\'').format(path))

def get_index(path, settings_file= './settings.json'):
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

	(spath, index) = os.path.split(path)
	#get the path relative to the CWD
	path = os.path.join(os.path.dirname(settings_file), spath)
	#remove loops
	path = os.path.relpath(os.path.abspath(path))

	#if the index exists in the dir
	if index_in_dir(path, index):
		return os.path.join(path, index)
	#if the path is in index/
	if spath == '' and index_in_dir(os.path.join(path, 'indexes/'), index):
		return os.path.join(path, 'indexes/', index)
	#if the path is in the BOWTIE2_INDEXES variable
	p = os.environ.get('BOWTIE2_INDEXES')
	if spath == '' and p and index_in_dir(p, index):
		return os.path.join(p, index)
	
	return None

def index_exists(path, settings_file='./settings.json'):
	return bool(get_index(path, settings_file))

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
