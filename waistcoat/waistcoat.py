"""Main pipeline file"""

import argparse, simplejson as json

def main():
	
	#parse command line
	my_args = get_arguments()

	#Read and validate settings for waistcoat
	print "Reading settings from \'{}\'".format(my_args.settings)
	settings = read_settings(my_args.settings)

	#split by barcode

	#clean files

	#filter unmapped

	#map to genome

	#index output SAM file


def get_arguments():
	parser = argparse.ArgumentParser(
		description="Process RNA-seq reads and map them to a genome")

	parser.add_argument('settings', help='Path to the waistcoat settings file')
	parser.add_argument('reads.fq.gz', help='GZIPed fastq file containing reads')
	parser.add_argument('output', nargs='?', default='waistcoat_out/',
			help='Directory to store output files')

	return parser.parse_args()

def read_settings(fname):
	"""Read settings from a JSON formatted file"""

	data = json.loads(file(fname).read())

	errors = []

	#there must be a format key
	if not data.has_key('format'):
		errors.append('No format key found in settings file \"{}\"'.format(fname))

	#format must be a string
	if not isinstance(data.get('format'), basestring):
		errors.append('format key should be a string, not {}'.format(
			type(data.get('format'))))
	else:
		#format can only contain Bs and Ns
		for char in data.get('format', ''):
			if char not in "BN":
				errors.append('format characters must be B or N')
				break
	
	#There must be a barcodes key
	if not data.has_key('barcodes'):
		errors.append('No barcodes key found in settings file \"{}\"'.format(fname))
	else:
		#barcodes must be a list of strings
		if isinstance(data.get('barcodes', []), basestring):
			errors.append('barcodes key should be a list of strings')
		
		for i,s in enumerate(data.get('barcodes', [])):
			if not isinstance(s, basestring):
				errors.append('barcodes[{}] should be a string, not {}'.format(type(s)))
			#barcode strings can only contain ATGC
			for char in s:
				if char.upper() not in "ATGC":
					errors.append('barcodes can only contain A,T,C or Gs')
					break

		#check whether the barcodes are the right length
		if data.has_key('format') and isinstance(data['format'], basestring):
			length = data['format'].count('B')
			for barcode in data['barcodes']:
				if len(barcode) != length:
					errors.append('barcodes must be the same length as the format')
					break

	if errors:
		raise ValueError("Failed to find settings file:\n\t" + "\n\t".join(errors))

	barcode_fmt = (tuple((i for i,j in enumerate(data['format']) if j=='B')),
			tuple((i for i,j in enumerate(data['format']) if j=='N')),)
	barcodes = [s.upper() for s in data['barcodes']]

	return (barcode_fmt, barcodes)

if __name__ == "__main__":
	main()
