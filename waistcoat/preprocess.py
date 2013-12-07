"""Preprocess the sequence data"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import tempfile, os, os.path, settings, statistics, itertools, shutil

import simplejson as json

verbose = True

DIFF_THRESH = 0.04
LEN_THRESH = 15

def run(in_file, my_settings, outdir=None):
	"""Run the preprocessing pipeline"""

	files = split_by_barcode(in_file, my_settings, outdir)
	files = process_sample(files, my_settings, outdir)

	return files

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
		files[sample] = tempfile.mkstemp(prefix=sample+'.', 
				suffix='.barcode', dir=outdir)
		files[sample] = (os.fdopen(files[sample][0], 'w'), files[sample][1])
		count[sample] = 0

	if verbose:
		print "Splitting sequences by barcode..."

	barcode2sample = {}
	for sample,barcode in my_settings.barcodes.iteritems():
		barcode2sample[barcode] = sample

	for seq in SeqIO.parse(in_file, 'fastq'):
		total += 1
		#get the barcode
		barcode = my_settings.parse_barcode(seq)
		sample = barcode2sample.get(barcode, '')
		if sample in files.iterkeys():
			SeqIO.write(seq, files[sample][0], 'fastq')
			count[sample] += 1

	#close the files
	for sample, (out, fname) in files.iteritems():
		out.close()
		files[sample] = fname

	#drop any samples which have no reads
	drop_empty(files, count)

	if verbose:
		print "\tFound {} sequences".format(total)
		for sample,barcode in my_settings.barcodes.iteritems():
			print "\t\t{}: {} ({:.1f}%)".format(sample, count[sample], 
					100.0*float(count[sample]) / float(total))
		um = total - sum(count.itervalues())
		print "\t\tunmapped: {:.1f} ({}%)".format(um, 100.0 * um / total)

	statistics.addValues('start_seqs', count)

	return files	

def str_dist(dist):
	low = 0
	high = len(dist)
	for d in dist:
		if d != 0:
			break
		low += 1
	for d in reversed(dist):
		if d!= 0:
			break
		high -= 1

	ret = []
	scale = 60.0 / float(max(dist))
	for i in range(low, high):
		ret.append("\t{:3d}: ({:03d}) |{}".format(
			i, dist[i], "*" * int(dist[i] * scale)))
	
	ret.append("\t\t* = {:.3f} reads".format(1.0 / scale))

	return "\n".join(ret)

def all_seqs(length):
	arg = ['ATGC',]*length
	return itertools.product(*arg)

def process_sample(in_files, my_settings, outdir=None):
	"""Clean reads and remove non-unique reads"""
	out_files = {}
	count = {}

	lengths = [0,]*1024


	#for each sample
	for i,(sample,in_file) in enumerate(in_files.iteritems()):
		if verbose: 
			print "Cleaning sample \'{}\' ({}/{})".format(sample, i+1, len(in_files))
		count[sample] = 0
		#prepare output files
		tempdir = tempfile.mkdtemp(dir=outdir, prefix=sample+'.', suffix='.by_umi')
		UMIs = {}
		UMI_usage = {}
		for c in all_seqs(my_settings.UMI_len()):
			umi = ''.join(c)
			(f,fname) = tempfile.mkstemp(
					dir=tempdir, prefix=sample+'.', suffix='.'+umi)
			UMIs[umi] = (os.fdopen(f, 'w'), fname)
			UMI_usage[umi] = 0

		#open the input file
		for seq in SeqIO.parse(in_file, 'fastq'):
			umi = my_settings.parse_UMI(seq)
			seq = clean_read(seq, my_settings)
			UMI_usage[umi] += 1
			if len(seq) > LEN_THRESH:
				SeqIO.write(seq, UMIs[umi][0], 'fastq')
		
		#delete input file and close all the output files
		os.remove(in_file)
		for umi, f in UMIs.iteritems():
			f[0].close()
			UMIs[umi] = f[1]

		if verbose:
			print "\tUMI usage (min, avg, max) = ({} {} {})".format(
					min(UMI_usage.itervalues()),
					sum(UMI_usage.itervalues()) / float(len(UMI_usage)),
					max(UMI_usage.itervalues()))
		#GC???

		#open final output file
		(fd, fname) = tempfile.mkstemp(
				dir=outdir, prefix=sample+'.', suffix='.clean')
		out_files[sample] = fname
		out_file = os.fdopen(fd, 'w')
		#for each UMI
		for umi, umi_file in UMIs.iteritems():
			#load reads into memory, group similar
			reads = []
			for seq in SeqIO.parse(umi_file, 'fastq'):
				if not len(seq):
					continue
				#try to find a match. We should do K-means or something
				for conflict in reads:
					if diff_seqs(seq, conflict[0]) <= DIFF_THRESH:
						conflict.append(seq)
						seq = None
						break
				#if this is unique
				if seq:
					reads.append([seq,])

			#resolve conflicts and write
			for conflict in reads:
				seq = resolve_conflict(conflict)
				if len(seq) < 1024:
					lengths[len(seq)] += 1
				SeqIO.write(seq, out_file, 'fastq')
				count[sample] += 1

			#clear reads
			del reads

		out_file.close()
		shutil.rmtree(tempdir)

	statistics.addValues('unique_seqs', count)

	if verbose:
		print "Length distribution:"
		print str_dist(lengths)

	return out_files


def diff_seqs(lhs, rhs):
	d = sum(1 for l,r in zip(str(lhs.seq).upper(), str(rhs.seq).upper())
				if l != r)
	return float(d) / float(min(len(rhs), len(lhs)))


def resolve_conflict(seqs):
	"""resolve conflicts between sequences with similar sequence"""
	if len(seqs) == 1:
		return seqs[0]

	seqs = [(s, sum(a for a in s.letter_annotations['phred_quality'])) 
			for s in seqs]
	max_q = max(s[1] for s in seqs)
	seqs = [s for s,qs in seqs if qs == max_q]

	return seqs[0]

def clean_read(seq_in, my_settings):
	"""
	Remove the poly(A) tail and barcode from seq_in and return the result
	"""
	#Remove all terminal As and annotations and barcode
	newseq = seq_in.seq.rstrip('aA')
	return my_settings.strip_header(seq_in[0:len(newseq)])

def drop_empty(files, count):
	"""Drop samples which are empty"""
	for sample in files.iterkeys():
		if count[sample] == 0:
			if verbose: print "Dropping \'{}\' as there are no reads"
			os.remove(files[sample])
			del files[sample]
