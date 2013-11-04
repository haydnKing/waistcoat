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
	for barcode, sample in barcodes.iteritems():
		files[barcode] = os.path.join(outdir, '{}.fq'.format(sample))

	barcodes = barcodes.items()

	#open the files
	for barcode, name in files.iteritems():
		files[barcode] = (name, open(name, 'w'))

	for seq in SeqIO.parse(in_file, 'fastq'):
		#get the barcode
		barcode = ''.join([seq.seq[x] for x in barcode_fmt]).upper()

		if barcode in barcodes:
			SeqIO.write(seq, files[barcode][1], 'fastq')

	#close the files
	for barcode, (fname, out) in files.iteritems():
		out.close()
		files[barcode] = fname

	return files

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


