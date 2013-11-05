"""Preprocess the sequence data"""

from Bio import SeqIO

import tempfile, os, os.path

import simplejson as json

verbose = True

def split_by_barcode(in_file, barcode_fmt, barcodes, outdir=None):
	"""
	Split all the sequences found in in_file into seperate files depending on
	their barcode, and saves unique sequences to files

	arguments:
		in_file: path to a FASTQ file containing the reads
		barcode_fmt: list of integers containting positions of the barcode
		barcodes: dictionary mapping (uppercase) barcode to sample name
		outdir (optional): output directory - 
			a temporary one is created if None

	return:
		A dictionary mapping barcodes to files, 
			including an extra file "unmapped"
		for reads which did not appear to contain a barcode
		outputs file to outdir
	"""
	#create a temp dir if necessary
	if not outdir:
		outdir = tempfile.mkdtemp(prefix='barcodes')

	files = {'unmapped': os.path.join(outdir, 'unmapped.fq')}
	count = {'unmapped': 0,}
	for barcode, sample in barcodes.iteritems():
		files[barcode] = os.path.join(outdir, '{}_nonunique.fq'.format(sample))
		count[barcode] = 0

	#open the files
	for barcode, name in files.iteritems():
		files[barcode] = (name, open(name, 'w'))

	if verbose:
		print "Splitting sequences by barcode..."

	for seq in SeqIO.parse(in_file, 'fastq'):
		#get the barcode
		barcode = str(''.join([seq.seq[x] for x in barcode_fmt]).upper())

		if barcode in barcodes.iterkeys():
			SeqIO.write(seq, files[barcode][1], 'fastq')
			count[barcode] += 1
		else:
			SeqIO.write(seq, files['unmapped'][1], 'fastq')
			count['unmapped'] += 1

	#close the files
	for barcode, (fname, out) in files.iteritems():
		out.close()
		files[barcode] = fname

	if verbose:
		print "Found {} sequences".format(sum(count.values()))
		print str_count(count)
		print "Removing duplicates..."

	removed = 0
	for barcode, sample in barcodes.iteritems():
		out_file = os.path.join(outdir, '{}.fq'.format(sample))
		out = remove_duplicates(files[barcode], out_file)
		files[barcode] = out_file
		rm = out[0] - out[1]
		count[barcode] -= rm
		removed += rm

	if verbose:
		print "Removed {} Sequences, {} remaining:".format(removed,
				sum(count.values()))
		print str_count(count)

	return files

def str_count(count):
	"""return a string representation of the count"""
	ret = []
	for barcode, number in count.iteritems():
		ret.append("\t{}: {}".format(barcode, number))

	return "\n".join(ret)

def str_dist(dist):
	low = 0
	high = len(dist)
	for d in dist:
		if d != 0:
			break
		low += 1
	for d in reverse(dist):
		if d!= 0:
			break
		high -= 1

	ret = []
	scale = 60.0 / float(max(dist))
	for i in range(low, high):
		ret.append("\t{:3d}: ({:03d}) |{}".format(
			i, dist[i], "*" * floor(dist[i] * scale)))
	
	ret.append("scale = {}".format(scale))

	return "\n".join(ret)


def remove_duplicates(input_file, output_file):
	"""
	Removes duplicate sequences from the file

	Arguments:
		input_file: name of the file to parse
		output_file: name of the file to output too

	Returns:
		a tuple of (before, after) giving the number of records found in 
		the file before and after the function is called
	"""
	before = after = 0

	sequences = []
	out = open(output_file, 'w')
	for seq in SeqIO.parse(input_file, 'fastq'):
		before += 1
		if not str(seq.seq) in sequences:
			SeqIO.write(seq, out, 'fastq')
			sequences.append(str(seq.seq))
			after += 1

	out.close()

	try:
		os.remove(input_file)
	except TypeError:
		pass

	return (before, after)

def clean_distributions(in_files, barcode_fmt, min_length = 15, 
		suffix = '_nopolyA', remove_input=True):
	"""
	Removes terminal As
	Reatains only unique reads
	Removes Bar Codes
	Discards reads shorter than min_length

	arguments:
		in_files: list of files to parse
		barcode_fmt: list of integers containing positions of the barcode
		min_length: the minimum length of read to retain
		suffix: suffix to add to output files
		remove_input: whether or not to delete the input file

	output:
		A list of outputted files
	"""
	files = []
	lengths = [0] * 1024

	for f in in_files:
		if verbose:
			print "Cleaning sequences in {}...".format(os.path.split(f)[1])

		temp_file = tempfile.NamedTemporaryFile()
		for seq in SeqIO.parse(f, 'fastq'):
			#Remove all terminal As and annotations
			newseq = seq.seq.rstrip('aA')
			seq = seq[0:len(newseq)]

			#remove barcode
			barcode_length = max(barcode_fmt)
			seq = seq[max(barcode_fmt)+1:]
			
			if len(seq.seq) < 1024:
				lengths[len(seq.seq)] += 1

			if len(seq.seq) > min_length:
				SeqIO.write(seq, temp_file, 'fastq')

		if remove_input:
			os.remove(f)

		out_file = "{}{}.fq".format(f[0:f.rfind('.')], suffix)
		temp_file.seek(0)
		remove_duplicates(temp_file, out_file)
		files.append(out_file)

	if verbose:
		print "Cleaned all sequences, length distribution:"
		print str_dist(lengths)

	return files

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

