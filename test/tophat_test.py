from waistcoat import tophat
tophat.verbose = False
import unittest, os, os.path, tempfile, shutil, testcases
import simplejson as json

DATA_DIR = os.path.join( os.path.split(__file__)[0], "data/")

general_options = {
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

specific_options = {
		'library_type': ('fr-unstranded', 'fr-firststrand', 'fr-secondstrand'),
		}

class TophatTestOptions(unittest.TestCase):
	"""Test the Tophat class options"""

	def setUp(self):
		self.tophat = tophat.TopHat()

	def tearDown(self):
		pass

	def test_getOptions(self):
		"""Test that TopHat.getOptions works"""

		self.assertEqual(self.tophat.getOptions(), [])

		self.tophat.read_gap_length = 15
		self.assertEqual(self.tophat.getOptions(), 
				['--read-gap-length', '15',])
		self.tophat.read_gap_length = None

		self.tophat.output_dir = "/path/to/a File\'s Location"
		self.assertEqual(self.tophat.getOptions(),
				['--output-dir', "/path/to/a File\'s Location"])
		self.tophat.output_dir = None

		self.tophat.bowtie1 = True
		self.assertEqual(self.tophat.getOptions(), ['--bowtie1',])

		self.tophat.bowtie1 = False
		self.assertEqual(self.tophat.getOptions(), [])

	def test_options_type(self):
		"""Test that options can be set to the right type"""
		values = {
				basestring: "A string",
				int: 5,
				bool: True,
				(int,float): 5.4,
				}

		#shouldn't fail
		for option,o_type in general_options.iteritems():
			setattr(self.tophat, option, values[o_type])
			self.assertEqual(getattr(self.tophat, option), values[o_type])

	def test_options_invalid(self):
		"""Test that trying to set an option which doesn't exist raises an
		exception"""
		self.assertRaises(ValueError, setattr,
				self.tophat, 'not_an_option', None)

	def test_options_none(self):
		"""Test that options can be set to none"""
		#also shouldn't fail
		for option in general_options.iterkeys():
			setattr(self.tophat, option, None)
			self.assertEqual(getattr(self.tophat, option), None)

	def test_options_wrong_type(self):
		#test one of each type
		#String
		self.assertRaises(TypeError, setattr, 
				self.tophat, 'index_base', 57)
		self.assertRaises(TypeError, setattr, 
				self.tophat, 'index_base', True)
		self.assertRaises(TypeError, setattr, 
				self.tophat, 'index_base', 57.4)

		#int
		self.assertRaises(TypeError, setattr,
				self.tophat, 'read_mismatches', 'a string')
		self.assertRaises(TypeError, setattr,
				self.tophat, 'read_mismatches', False)
		self.assertRaises(TypeError, setattr,
				self.tophat, 'read_mismatches', 2.5)
		
		#Bool
		self.assertRaises(TypeError, setattr, 
				self.tophat, 'integer_quals', 'a string')
		self.assertRaises(TypeError, setattr, 
				self.tophat, 'integer_quals', 5)
		self.assertRaises(TypeError, setattr, 
				self.tophat, 'integer_quals', 5.5)



	def test_specific_options(self):
		"""Test options which only accept specific strings"""

		for option,possible_values in specific_options.iteritems():
			for value in possible_values:
				setattr(self.tophat, option, value)
			#test that setting an unsupported value generates an exception
			self.assertRaises(ValueError, setattr, self.tophat, option, 
				"invalid option value")

def jread(filename):
	return json.loads(open(filename).read())['options']

class TopHatSettingsFile(unittest.TestCase):
	"""Test the parsing of the settings file"""
	valid = jread(os.path.join(DATA_DIR, "tophat_settings/valid.json"))
	invalid = jread(os.path.join(DATA_DIR, "tophat_settings/invalid.json"))

	def test_load_valid(self):
		"""Test loading a valid file"""
		th = tophat.tophat_from_settings(self.valid)

		for option,value in self.valid.iteritems():
			self.assertEqual(getattr(th, option), value, 
					"option \'{}\' should be \'{}\', not \'{}\'".format(option,
						value, getattr(th,option)))

	def test_load_invalid(self):
		"""test that loading an invalid file generates an exception"""
		self.assertRaises(BaseException, tophat.tophat_from_settings, 
				self.invalid)
	
	def test_is_valid(self):
		"""test that tophat.is_valid returns correctly"""
		self.assertTrue(tophat.is_valid(self.valid))
		self.assertFalse(tophat.is_valid(self.invalid))

class TophatTestRun(testcases.TestFastQ):
	"""Test that the tophat class can actually run tophat"""

	def setUp(self):
		self.tophat = tophat.TopHat()
		self.output = tempfile.mkdtemp(prefix='test')

	def tearDown(self):
		#Delete tmpdir
		shutil.rmtree(self.output)


	def test_run(self):
		"""test that I can execute tophat on test data"""
		
		data = os.path.join(DATA_DIR, 'tophat_data/')
		
		expected_files = ['accepted_hits.sam',
					'deletions.bed',
					'insertions.bed',
					'junctions.bed',
					'logs',
					'prep_reads.info',
					'unmapped.bam']


		self.tophat.mate_inner_dist = 20
		self.tophat.no_convert_bam = True
		self.tophat.output_dir = self.output

		self.tophat.run(
				os.path.join(data, 'reads_1.fq'),
				os.path.join(data, 'reads_2.fq'),
				os.path.join(data, 'test_ref'))

		files = sorted(os.listdir(self.output))

		self.assertEqual(files, expected_files)

	def test_discard_mapped(self):
		data = os.path.join(DATA_DIR, 'tophat_data/')
		files = ['test_ref.1.bt2',
						 'test_ref.2.bt2',
						 'test_ref.3.bt2',
						 'test_ref.4.bt2',
						 'test_ref.rev.1.bt2',
						 'test_ref.rev.2.bt2',
						 'test_ref.fa',
						 'reads_unmapped.fq']

		expected = os.path.join(DATA_DIR, 'tophat_unmapped_out.fq')
		expected_name = os.path.join(self.output, 
				'reads_unmapped_nomapping.fq')

		for f in files:
			shutil.copyfile(os.path.join(data, f), os.path.join(self.output,f))

		(fname, length) = tophat.discard_mapped(
							os.path.join(self.output, 'reads_unmapped.fq'),
							os.path.join(self.output, 'test_ref'))

		self.assertEqual(length, 2)

		self.assertEqual(fname, expected_name,
				"discard_mapped expected \'{}\', produced \'{}\'".format(
					expected_name, fname))
		self.assertTrue(os.path.exists(expected_name), 
			"Expected output file \'{}\' not produced".format(expected_name))

		
		self.assertFastQ(expected, expected_name)


		
