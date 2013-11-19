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
		my_settings: Settings object
		outdir (optional): output directory - outputs to the same directory as
			in_file if None

	return:
		A dictionary mapping sample names to file names
	"""
	#create a temp dir if necessary
	if not outdir:
		outdir = os.path.dirname(in_file)

	#initialte count and open files
	files = {}
	count = {}
	total = 0
	for sample,barcode in my_settings.barcodes.iteritems():
		name = os.path.join(outdir, '{}_nonunique.fq'.format(sample))
		files[sample] = (name, open(name, 'w'))
		count[barcode] = 0

	if verbose:
		print "Splitting sequences by barcode..."

	for seq in SeqIO.parse(in_file, 'fastq'):
		total += 1
		#get the barcode
		barcode = my_settings.parse_barcode(seq)[0]
		sample = my_settings.barcodes.get(barcode, '')
		if sample in files.iterkeys():
			SeqIO.write(seq, files[sample][1], 'fastq')
			count[barcode] += 1

	#close the files
	for sample, (fname, out) in files.iteritems():
		out.close()
		files[sample] = fname

	if verbose:
		print "Found {} sequences".format(total)
		for sample,barcode in my_settings.barcodes.iteritems():
			print "{}: {}".format(sample, count[sample])
		print "unmapped: {}".format(total - sum(count.itervalues()))

	return files	

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


def remove_duplicates(in_files, my_settings, outdir=None):
	"""
	Removes duplicate sequences from the file

	Arguments:
		input_file: name of the files to parse
		output_file: name of the file to output too

	Returns:
		a dictionary mapping samples to files	
	"""
	files = {}

	stats = {}

	for sample, f in in_files.iteritems():
		if verbose:
			print "Remove Duplicates in {}...".format(sample)

		out_file = tempfile.mkstemp(prefix='no_duplicates', dir=outdir)
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
		os.remove(f)
		stats[sample] = (before, after)

	if verbose:
		print "sample, before, after, % unique"
		for sample, (before, after) in stats.iteritems():
			print "{}, {}, {}, {:.1f}%".format(sample, before, after, 
					100.0 * after / float(before))

	return files

def clean_distributions(in_files, my_settings, outdir=None, min_length = 15, 
		remove_input=True):
	"""
	Removes terminal As
	Reatains only unique reads
	Removes Bar Codes
	Discards reads shorter than min_length

	arguments:
		in_files: a dictionary mapping samples to files
		settings: settings object
		min_length: the minimum length of read to retain
		remove_input: whether or not to delete the input file

	output:
		A dictionary mapping samples to files
	"""
	files = {}
	lengths = [0] * 1024

	for sample,f in in_files.iteritems():
		if verbose:
			print "Cleaning {}...".format(sample)

		(out_file,out_file_name) = tempfile.mkstemp(dir=outdir)
		files[sample] = out_file_name

		for seq in SeqIO.parse(f, 'fastq'):
			#Remove all terminal As and annotations
			newseq = seq.seq.rstrip('aA')
			seq = seq[0:len(newseq)]

			#remove barcode
			seq = my_settings.strip_barcode(seq)
			
			if len(seq.seq) < 1024:
				lengths[len(seq.seq)] += 1

			if len(seq.seq) > min_length:
				SeqIO.write(seq, out_file, 'fastq')

		if remove_input:
			os.remove(f)
		out_file.close()


	if verbose:
		print "Cleaned all sequences, length distribution:"
		print str_dist(lengths)

	return files


