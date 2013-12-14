from Bio import SeqIO, Seq, SeqRecord

import random

barcode = "TCCA" 

seqs = []

length = 60;

out = open("process_seq_out.fq", "wb")

#generate 10 seqs
for i in  range(10):
	s = ''.join(random.choice("ATGC") for k in range(random.randint(25,60))) +'G'
	seq = SeqRecord.SeqRecord(Seq.Seq(s), "seq{}".format(i),
			name="seq{}".format(i))
	seq.letter_annotations = {
			'phred_quality': [20, ] * len(s),
			}
	SeqIO.write(seq, out, "fastq")
	

	s = ("TCC" + random.choice("ATCG") + "A" + s + 
			('A'*random.randint(0,30)))
	seq = SeqRecord.SeqRecord(Seq.Seq(s), "seq{}".format(i),
			name="seq{}".format(i))
	seq.letter_annotations = {
			'phred_quality': [20, ] * len(s),
			}
	seqs.append(seq)

out.close()
in_file = open("process_seq_in.fq", "wb")

for i in range(100):
	SeqIO.write(random.choice(seqs), in_file, "fastq")

in_file.close()

