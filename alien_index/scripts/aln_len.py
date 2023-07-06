import Bio
from Bio import Seq
from Bio import SeqIO
import sys 
import os  

###############################################
# aln_len.py 
# Jen Wisecaver
# 20201018 
# input: an alignment in fasta format 
# output: length of alignment
################################################

infile = sys.argv[1]

seq_record = next(SeqIO.parse(infile, "fasta"))

sequence = str(seq_record.seq)
seqlen = len(sequence)
print(seqlen)
