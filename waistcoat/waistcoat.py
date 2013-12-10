"""Main pipeline file"""

import argparse, sys, statistics, gzip

import settings, tophat, preprocess, postprocess, tempfile, shutil, os,os.path, pysam

tophat.verbose = False

verbose = True
check_output = True

def main():
	
	#parse command line
	my_args = get_arguments()
	run(my_args.settings, my_args.reads, my_args.output)

def run(settings_file, reads, outdir, temp_loc=None):

	if os.path.exists(outdir):
		if (check_output and not 
				query_yes_no("Output path \"{}\" already exists, overwrite?"
					.format(outdir))):
			print "Cannot continue as output already exists"
			sys.exit(1)
		
		shutil.rmtree(outdir)

	os.mkdir(outdir)

	tempdir = tempfile.mkdtemp(prefix='waistcoat', dir=temp_loc)

	#Read and validate settings for waistcoat
	if verbose:
		print "Reading settings from \'{}\'".format(settings_file)
	my_settings = settings.loadf(settings_file)

	statistics.setUp(my_settings.barcodes.keys())

	#run the preprocessing pipeline
	if verbose:
		print "\n========== Preprocessing =========="
	remove_input = False
	if reads.endswith('.gz'):
		print "Inflating..."
		gzfile = gzip.GzipFile(reads, 'r')
		(out, reads) = tempfile.mkstemp(dir=tempdir, prefix='input.', 
											suffix='.inflated')
		out = os.fdopen(out, 'w')
		out.writelines(gzfile)
		gzfile.close()
		out.close()
		remove_input = True

	files = preprocess.run(reads, my_settings, tempdir, remove_input)

	#discard those which map to discard
	if verbose:
		print "\n========== Discard =========="
	for i,(index, dcs) in enumerate(my_settings.discard):
		new_files = {}
		count = {}
		if verbose:
			print "Removing reads which map to \'{}\' ({}/{})...".format(index,
					i+1, len(my_settings.discard))
		for sample,f in files.iteritems():
			if verbose:
				print "\tScanning \'{}\'".format(sample)
			(new_files[sample], count[sample]) = (
					tophat.discard_mapped(f, index, tophat_settings = dcs))
		files = new_files
		statistics.addValues('discard_' + os.path.basename(index), count)

	#map to genome
	(target, target_settings) = my_settings.target
	if verbose:
		print "\n========== Map to {} ==========".format(os.path.basename(target))
	th = tophat.tophat_from_settings(target_settings)
	for i,(sample,f) in enumerate(files.iteritems()):
		th.output_dir = os.path.join(outdir, sample)
		os.mkdir(th.output_dir)
		if verbose: print "Mapping {} ({}/{})...".format(sample, i+1, len(files))
		th.run(f, index_base = target)
		os.remove(f)
		
	if verbose: print "\n========== Postprocess =========="
	count = {}
	for i,(sample,f) in enumerate(files.iteritems()):
		if verbose: print "{} ({}/{})...".format(sample, i+1, len(files))
		
		ah = os.path.join(outdir, sample, 'accepted_hits.bam')
		count[sample] = postprocess.run(ah,	"{}.fa".format(target))
		statistics.collectFinalStats(sample, ah)
		
	statistics.addValues('final_seqs', count)

	statistics.write(os.path.join(outdir, 'statistics'))
	
	shutil.rmtree(tempdir)

	if verbose:
		print "\n__________ Pipeline Statistics __________"
		print statistics.prettyString()

def get_arguments():
	parser = argparse.ArgumentParser(
		description="Process RNA-seq reads and map them to a genome")

	parser.add_argument('settings', help='Path to the waistcoat settings file')
	parser.add_argument('reads', help='GZIPed fastq file containing reads')
	parser.add_argument('output', nargs='?', default='waistcoat_out/',
			help='Directory to store output files')

	return parser.parse_args()

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")





if __name__ == "__main__":
	main()

