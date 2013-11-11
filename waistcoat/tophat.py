"""Interface with the tophat program"""

import command, tempfile, pysam, os.path, shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import simplejson as json

verbose = True

class TopHat(command.Command):
	"""Class to interface with tophat
			See http://tophat.cbcb.umd.edu/manual.shtml

			To set a particular command line option, e.g. "--read-gap-length 5"
			set self.read_gap_length = 5
	"""

	options = {
			'index_base': basestring,
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
			'b2_i' : basestring,
			'b2_n_ceil' : basestring,
			'b2_gbar' : int,
			'b2_mp' : basestring,
			'b2_np' : int,
			'b2_rdg' : basestring,
			'b2_rfg' : basestring,
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

	cmd = "tophat"

	def __init__(self):

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
			update_fn = self.__the_callback,
			stderr = True)

	def __the_callback(self, line):
		if verbose:
			print ">Tophat: {}".format(line)

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

	def __setattr__(self, name, value):
		if self.options.has_key(name):
			if name == 'library_type':
				if value in ('fr-unstranded', 'fr-firststrand', 'fr-secondstrand',
						None):
					self.__dict__[name] = value
					return
				else:
					raise ValueError(("Option \'library_type\' supports values " + 
							"\'fr-unstranded\', \'fr-firststrand\', and " + 
							"\'fr-secondstrand\' only, not \'{}\'").format(value))
			#have to handle int seperately to make sure we aren't given a bool to
			# store as an int
			if not self._option_is_int(name) and ( 
					isinstance(value, self.options[name]) or value is None):
				self.__dict__[name] = value
			elif self._option_is_int(name) and (
					isinstance(value, self.options[name]) or value is None) and (
					not isinstance(value, bool)):
				self.__dict__[name] = value
			else:
				raise TypeError(
					"Setting \'{}\' must be of type {} or None, not {}".format(
						name, str(self.options[name]), str(type(value))))
		else:
			raise ValueError("Tophat has no option \'{}\'".format(name))

	def _option_is_int(self, name):
		if hasattr(self.options[name], '__iter__'):
			return int in self.options[name]
		else:
			return self.options[name] == int

def load_settings(file_name):
	"""Read a JSON settings file and return a TopHat object with those options
	set
	Raises ValueError if the file is invalid"""
	
	data = json.loads(open(file_name).read())

	#object to return
	t = TopHat()

	for option,value in data['options'].iteritems():
		setattr(t, option, value)

	return t

def discard_mapped(reads_file, index_base, tophat_settings = None, 
			suffix="_nomapping"):
	"""Map reads to the index and discard all the reads which map successfully"""
	tempd = tempfile.mkdtemp(prefix='waistcoat')

	outfile_name = reads_file[0:reads_file.rfind('.')] + suffix + ".fq"
	
	if tophat_settings:
		th = load_settings(tophat_settings)
	else:
		th = TopHat()
	
	th.keep_fasta_order = True
	th.output_dir = tempd

	#run tophat
	th.run(reads_file, index_base= index_base)

	#open the hits
	samfile = pysam.Samfile(os.path.join(tempd, "unmapped.bam"), "rb")

	outfile = open(outfile_name, "wb")
	#write out FASTQ records
	for read in samfile:
		SeqIO.write(build_fastq(read), outfile, 'fastq')


	shutil.rmtree(tempd)
	os.remove(reads_file)

	outfile.close()

	return outfile_name


# Precompute conversion table
SANGER_SCORE_OFFSET = ord("!")
q_mapping = dict()
for letter in range(0, 255):
	q_mapping[chr(letter)] = letter - SANGER_SCORE_OFFSET

def build_fastq(read):
	'''Build FASTQ SeqRecords out of BAM reads
	
	Note: chunks copied from Bio.SeqIO.QualityIO.py

	This code adapted from issue 137 from pysam:
	https://code.google.com/p/pysam/issues/attachmentText?id=137&aid=1370000000&name=build_fastq_seqrecords.py
	'''

	# Get the sequence first
	descr = read.qname
	id = read.qname
	name = id
	from Bio.Alphabet import IUPAC
	record = SeqRecord(Seq(read.seq, IUPAC.ambiguous_dna),
										 id=id, name=name, description=descr)

	# Get the qualities second
	qualities = [q_mapping[letter] for letter in read.qual]
	if qualities and (min(qualities) < 0 or max(qualities) > 93):
			raise ValueError("Invalid character in quality string")

	#For speed, will now use a dirty trick to speed up assigning the
	#qualities. We do this to bypass the length check imposed by the
	#per-letter-annotations restricted dict (as this has already been
	#checked by FastqGeneralIterator). This is equivalent to:
	#record.letter_annotations["phred_quality"] = qualities
	dict.__setitem__(record._per_letter_annotations,
									 "phred_quality", qualities)
	
	return record



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
		
	
