
import pysam, os.path, os, shutil, tempfile
from Bio import SeqIO

def run(samfile, genome, target_length = 28):
	"""extend hits by adding As up until target length where possible"""

	#load the targets into RAM
	targets = {}
	for seq in SeqIO.parse(genome, 'fasta'):
		targets[seq.name] = seq

	#touch the temporary file
	(temp, tempname) = tempfile.mkstemp('w', prefix='postprocess')
	os.close(temp)

	#open the files
	instream = pysam.Samfile(samfile, 'rb')
	outstream = pysam.Samfile(tempname, 'wb', template = instream)

	for alg in instream.fetch():
		rname = instream.getrname(alg.tid)
		if targets.has_key(rname):
			extend_limit = count_a(alg, targets[rname])
			extend = max(extend_limit, target_length - alg.alen)
			if not alg.is_reverse:
				qual = alg.qual
				alg.seq = alg.seq + ('A' * extend)
				alg.qual = qual + ('!' * extend)
			else:
				qual = alg.qual
				alg.seq = ('T' * extend) + alg.seq
				alg.qual = ('T' * extend) + qual
				alg.pos = alg.pos - extend
		outstream.write(alg)
	
	#close up
	instream.close()
	outstream.close()

	#move to original location
	shutil.move(tempname, samfile)

	#rebuild index
	pysam.index(samfile)
			

def count_a(alg, target):
	"""Count the number of As after the read"""
	
	if not alg.is_reverse:
		p = int(alg.aend)
		while str(target[p+1]).upper() == 'A':
			p += 1
		return p - alg.aend
	else:
		p = int(alg.pos)
		while str(target[p-1]).upper() == 'A':
			p -= 1
		return alg.pos - p

