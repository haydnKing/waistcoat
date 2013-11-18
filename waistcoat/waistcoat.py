"""Main pipeline file"""

import argparse

import settings, tophat, preprocess

def main():
	
	#parse command line
	my_args = get_arguments()
	run(my_args.settings, myargs.reads.fq.gz, output)

def run(settings_file, reads, outdir):
	#Read and validate settings for waistcoat
	print "Reading settings from \'{}\'".format(settings_file)
	s = settings.loadf(settings_file)

	#split by barcode
	files = preprocess.split_by_barcode(reads, s)

	#clean files
	files = preprocess.clean_distributions(files, s)

	#discard those which map to discard
	for (index, dcs) in s.discard:
		new_files = []
		for f in files:
			new_files.append(tophat.discard_mapped(f, index, settings = dcs))
		files = new_files

	#map to genome
	th = tophat.tophat_from_settings(s.map_settings)
	for f in files:
		th.output_dir = os.path.join(outdir, os.path.basename(f).split('.')[0])
		os.mkdir(th.output_dir)
		th.run(f, index_base = s.map)
		os.remove(f)
	
	#index output SAM file


def get_arguments():
	parser = argparse.ArgumentParser(
		description="Process RNA-seq reads and map them to a genome")

	parser.add_argument('settings', help='Path to the waistcoat settings file')
	parser.add_argument('reads.fq.gz', help='GZIPed fastq file containing reads')
	parser.add_argument('output', nargs='?', default='waistcoat_out/',
			help='Directory to store output files')

	return parser.parse_args()

if __name__ == "__main__":
	main()
