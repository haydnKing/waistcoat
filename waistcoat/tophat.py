"""Interface with the tophat program"""

import command

import simplejson as json

class TopHat(command.Command):
	"""Class to interface with tophat
			See http://tophat.cbcb.umd.edu/manual.shtml

			To set a particular command line option, e.g. "--read-gap-length 5"
			set self.read_gap_length = 5
	"""

	options = {
			'read_mismatches' : int,
			'read_gap_length' : int,
			'read_edit_dist' : int,
			'read_realign_edit_dist' : int,
			'bowtie1' : bool,
			'output_dir' : basestring,
			'mate_inner_dist' : int,
			'mate_std_dev' : int,
			'min_anchor_length' : int,
			'splice_mismatches' : int,
			'min_intron_length' : int,
			'max_intron_length' : int,
			'max_insertion_length' : int,
			'max_deletion_length' : int,
			'solexa_quals' : bool,
			'solexa1' : bool,
			'quals' : bool,
			'integer_quals' : bool,
			'color' : bool,
			'num_threads' : int,
			'max_multihits' : int,
			'report_secondary_alignments' : bool,
			'no_discordant' : bool,
			'no_mixed' : bool,
			'no_coverage_search' : bool,
			'coverage_search' : bool,
			'microexon_search' : bool,
			'library_type' : basestring,
			'bowtie_n' : bool,
			'segment_mismatches' : int,
			'segment_length' : int,
			'min_segment_intron' : int,
			'max_segment_intron' : int,
			'min_coverage_intron' : int,
			'max_coverage_intron' : int,
			'keep_tmp' : bool,
			'keep_fasta_order' : bool,
			'no_sort_bam' : bool,
			'no_convert_bam' : bool,
			'resume' : basestring,
			'zpacker' : basestring,
			'b2_very_fast' : bool,
			'b2_fast' : bool,
			'b2_sensitive' : bool,
			'b2_very_sensitive' : bool,
			'b2_N' : int,
			'b2_L' : int,
			'b2_i' : int,
			'b2_n_ceil' : int,
			'b2_gbar' : int,
			'b2_mp' : (int,float),
			'b2_np' : (int,float),
			'b2_rdg' : (int,float),
			'b2_rfg' : (int,float),
			'b2_score_min' : basestring,
			'b2_D' : int,
			'b2_R' : int,
			'fusion_search' : bool,
			'fusion_anchor_length' : int,
			'fusion_min_dist' : int,
			'fusion_read_mismatches' : int,
			'fusion_multireads' : int,
			'fusion_multipairs' : int,
			'fusion_ignore_chromosomes' : bool,
			'raw_juncs' : basestring,
			'no_novel_juncs' : bool,
			'GTF' : basestring,
			'transcriptome_index' : basestring,
			'transcriptome_only' : bool,
			'transcriptome_max_hits' : int,
			'prefilter_multihits' : bool,
			'insertions' : basestring,
			'deletions' : basestring,
			'no_novel_indels' : bool,}

	def __init__(self):
		self.cmd = "tophat"

		self.index_base = ""

		for o in self.options.iterkeys():
			setattr(self, o, None)

	def run(self, reads_1, reads_2 = '', index_base = ''):
		"""Run tophat"""

		#make sure we allways have a list of files
		if isinstance(reads_1, basestring):
			reads_1 = [reads_1,]
		if isinstance(reads_2, basestring):
			reads_2 = [reads_2]

		index_base = index_base if index_base else self.index_base
		if not index_base:
			raise ValueError("Location of index_base not set")

		self.call( self.getOptions() + [index_base,
			','.join([read for read in reads_1]),
			','.join([read for read in reads_2])],
			update_fn = self.__the_callback)

	def __the_callback(self, line):
		print "tophat: {}".format(line)

	def getOptions(self):
		opts = []
		for o in self.options.keys():
			opt = getattr(self, o)
			o_ = '--{}'.format(o.replace('_', '-'))
			#if the option is set
			if opt is not None:
				#flags
				if isinstance(opt, bool):
					if opt: opts.append(o_)
				#strings
				elif isinstance(opt, basestring):
					opts += [o_, opt,]
				#numbers
				elif isinstance(opt, (int, long, float)):
					opts += [o_, str(opt),]

		return opts


def load_settings(file_name):
	"""Read a JSON settings file and return a TopHat object with those options
	set
	Raises ValueError if the file is invalid"""
	pass

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
		
	
