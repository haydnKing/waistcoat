"""Main pipeline file"""

import argparse

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

if __name__ == "__main__":
	main()
