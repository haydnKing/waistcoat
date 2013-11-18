"""Preprocess the sequence data"""

from Bio import SeqIO

import tempfile, os, os.path, settings

import simplejson as json

verbose = True

def split_by_barcode(in_file, my_settings, outdir=None):
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

	#initialte count and open files
	files = {}
	count = {}
	total = 0
	for sample,barcode in my_settings.barcodes.iteritems():
		name = os.path.join(outdir, '{}_nonunique.fq'.format(sample))
		files[barcode] = (name, open(name, 'w'))
		count[barcode] = 0

	if verbose:
		print "Splitting sequences by barcode..."

	for seq in SeqIO.parse(in_file, 'fastq'):
		total += 1
		#get the barcode
		barcode = my_settings.parse_barcode(seq)[0]

		if barcode in files.iterkeys():
			SeqIO.write(seq, files[barcode][1], 'fastq')
			count[barcode] += 1

	#close the files
	for barcode, (fname, out) in files.iteritems():
		out.close()
		files[barcode] = fname

	if verbose:
		print "Found {} sequences".format(total)
		for sample,barcode in my_settings.barcodes.iteritems():
			print "{}: {}".format(sample, count[sample])
		print "unmapped: {}".format(total - sum(count.itervalues()))
		print "Removing duplicates..."

	removed = {}
	for sample,barcode in my_settings.barcodes.iteritems():
		out_file = os.path.join(outdir, '{}.fq'.format(sample))
		out = remove_duplicates(files[barcode], out_file)
		files[barcode] = out_file
		removed[barcode] = out[0] - out[1]

	if verbose:
		for sample,barcode in my_settings.barcodes.iteritems():
			print "{}: {} -> {} | {}% unique".format(sample, count[sample],
					count[sample]-removed[sample], 100.0 * count[sample] / (float) (
						count[sample]-removed[sample]))

	return files.values()


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

def clean_distributions(in_files, my_settings, min_length = 15, 
		suffix = '_nopolyA', remove_input=True):
	"""
	Removes terminal As
	Reatains only unique reads
	Removes Bar Codes
	Discards reads shorter than min_length

	arguments:
		in_files: list of files to parse
		settings: settings object
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
			seq = my_settings.strip_barcode(seq)
			
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


