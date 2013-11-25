"""Main pipeline file"""

import argparse

import settings, tophat, preprocess, tempfile, shutil, os,os.path

def main():
	
	#parse command line
	my_args = get_arguments()
	run(my_args.settings, my_args.reads, my_args.output)

def run(settings_file, reads, outdir):
	tempdir = tempfile.mkdtemp()

	#Read and validate settings for waistcoat
	print "Reading settings from \'{}\'".format(settings_file)
	my_settings = settings.loadf(settings_file)

	#run the preprocessing pipeline
	files = preprocess.run(reads, my_settings, tempdir)

	#discard those which map to discard
	for (index, dcs) in my_settings.discard:
		new_files = {}
		for sample,f in files.iteritems():
			print "tophat.discard_mapped(\'{}\', \'{}\', {})".format(
					sample, index, dcs)
			new_files[sample] = tophat.discard_mapped(f, index, 
																												tophat_settings = dcs)
		files = new_files

	#map to genome
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	(target, target_settings) = my_settings.target
	th = tophat.tophat_from_settings(target_settings)
	for sample,f in files.iteritems():
		th.output_dir = os.path.join(outdir, sample)
		os.mkdir(th.output_dir)
		th.run(f, index_base = target)
		os.remove(f)
	
	shutil.rmtree(tempdir)

	#index output SAM file


def get_arguments():
	parser = argparse.ArgumentParser(
		description="Process RNA-seq reads and map them to a genome")

	parser.add_argument('settings', help='Path to the waistcoat settings file')
	parser.add_argument('reads', help='GZIPed fastq file containing reads')
	parser.add_argument('output', nargs='?', default='waistcoat_out/',
			help='Directory to store output files')

	return parser.parse_args()

if __name__ == "__main__":
	main()
