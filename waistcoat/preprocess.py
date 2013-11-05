"""Preprocess the sequence data"""

from Bio import SeqIO

import tempfile, os, os.path

verbose = True

def split_by_barcode(in_file, barcode_fmt, barcodes, outdir=None):
	"""
	Split all the sequences found in in_file into seperate files depending on
	their barcode, and saves unique sequences to files

	arguments:
		in_file: path to a FASTQ file containing the reads
		barcode_fmt: list of integers containting positions of the barcode
		barcodes: dictionary mapping (uppercase) barcode to sample name
		outdir (optional): output directory - a temporary one is created if None

	return:
		A dictionary mapping barcodes to files, including an extra file "unmapped"
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

def remove_duplicates(input_file, output_file):
	"""
	Removes duplicate sequences from the file

	Arguments:
		input_file: name of the file to parse
		output_file: name of the file to output too

	Returns:
		a tuple of (before, after) giving the number of records found in the file
		before and after the function is called
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

	os.remove(input_file)

	return (before, after)

def clean_distributions(in_files, barcode_fmt, min_length = 15, outdir=None,
		remove_input=False):
	"""
	Removes terminal As
	Reatains only unique reads
	Removes Bar Codes
	Discards reads shorter than min_length

	arguments:
		in_files: list of files to parse
		barcode_fmt: list of integers containing positions of the barcode
		min_length: the minimum length of read to retain
		outdir: directory to output to - a temporary one is created if None
		remove_input: whether or not to delete the input file

	output:
		A list of outputted files
		file names are chosed as:
			original_name.ext -> original_name_nopolyA.ext
	"""
	for f in in_files:
		pass


