
import pysam, os.path, os, shutil, tempfile, statistics
from Bio import SeqIO

def run(outdir, sample, genome, target_length = 28, extend=False):
	"""extend hits by adding As up until target length where possible"""

	samfile = os.path.join(outdir, sample, 'accepted_hits.bam')
	pysam.index(samfile)
	if extend:
		count = extend_short_reads(samfile, genome, target_length)
	else:
		count = count_reads(samfile)


	# sort and index
	outsam = os.path.join(outdir, '{}'.format(sample))
	pysam.sort(samfile, outsam)
	outsam = "{}.bam".format(outsam)
	pysam.index(outsam)

	return count
			
def count_reads(samfile):
	instream = pysam.Samfile(samfile, 'rb')
	count = instream.mapped
	instream.close()
	return count

def extend_short_reads(samfile, genome, target_length):
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

	count = 0
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
		count += 1
	
	#close up
	instream.close()
	outstream.close()

	#move to original location
	shutil.move(tempname, samfile)

	#rebuild index
	os.remove(samfile + ".bai")
	pysam.index(samfile)

	return count

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

