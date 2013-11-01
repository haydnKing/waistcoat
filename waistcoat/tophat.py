"""Interface with the tophat program"""

import command

class TopHat(command.Command):
	"""Class to interface with tophat
			See http://tophat.cbcb.umd.edu/manual.shtml

			To set a particular command line option, e.g. "--read-gap-length 5"
			set self.read_gap_length = 5
	"""

	options = ['read_mismatches',
						'read_gap_length',
						'read_edit_dist',
						'read_realign_edit_dist',
						'bowtie1',
						'output_dir',
						'mate_inner_dist',
						'mate_std_dev',
						'min_anchor_length',
						'splice_mismatches',
						'min_intron_length',
						'max_intron_length',
						'max_insertion_length',
						'max_deletion_length',
						'solexa_quals',
						'solexa1',
						'quals',
						'integer_quals',
						'color',
						'num_threads',
						'max_multihits',
						'report_secondary_alignments',
						'no_discordant',
						'no_mixed',
						'no_coverage_search',
						'coverage_search',
						'microexon_search',
						'library_type',
						'bowtie_n',
						'segment_mismatches',
						'segment_length',
						'min_segment_intron',
						'max_segment_intron',
						'min_coverage_intron',
						'max_coverage_intron',
						'keep_tmp',
						'keep_fasta_order',
						'no_sort_bam',
						'no_convert_bam',
						'resume',
						'zpacker',
						'b2_very_fast',
						'b2_fast',
						'b2_sensitive',
						'b2_very_sensitive',
						'b2_N',
						'b2_L',
						'b2_i',
						'b2_n_ceil',
						'b2_gbar',
						'b2_mp',
						'b2_np',
						'b2_rdg',
						'b2_rfg',
						'b2_score_min',
						'b2_D',
						'b2_R',
						'fusion_search',
						'fusion_anchor_length',
						'fusion_min_dist',
						'fusion_read_mismatches',
						'fusion_multireads',
						'fusion_multipairs',
						'fusion_ignore_chromosomes',
						'raw_juncs',
						'no_novel_juncs',
						'GTF',
						'transcriptome_index',
						'transcriptome_index',
						'transcriptome_index',
						'transcriptome_only',
						'transcriptome_max_hits',
						'prefilter_multihits',
						'insertions',
						'deletions',
						'no_novel_indels',]

	def __init__(self):
		self.cmd = "tophat"

		self.ebwt_base = ""

		for o in self.options:
			setattr(self, o, None)

	def getOptions(self):
		opts = []
		for o in self.options:
			opt = getattr(self, o)
			o_ = '--{}'.format(o.replace('_', '-'))
			#if the option is set
			if opt is not None:
				#flags
				if isinstance(opt, bool):
					if opt: opts.append(o_)
				#strings
				elif isinstance(opt, basestring):
					opts += [o_, shellquote(opt),]
				#numbers
				elif isinstance(opt, (int, long, float)):
					opts += [o_, str(opt),]

		return opts
		
def shellquote(s):
	return "\"{}\"".format(s.replace("\"", "\\\""))

def test_programs():
	commands = {'tophat': False, 
		    'bowtie2': False, 
		    'bowtie2-inspect': False, 
		    'bowtie2-build': False, 
		    'samtools': False,}
	for cmd in commands.keys():
		commands[cmd] = command.isAvailable(cmd)

	errors = [cmd for cmd,available in commands.iteritems() if not available]

	if errors:
		raise ImportError("Required program(s) \'{}\' were not found.".format(
			"\', \'".join(errors)))

test_programs()
		
	
