#!/depot/jwisecav/apps/envs/env.genomics/bin/python
from Bio import SeqIO
import sys
import os

###############################################
# parse_braker.py
# Jen Wisecaver
# 20220414 
# input: Braker directory containing predicted 
#			proteins and coding sequences 
# output: protein and coding sequence fasta
#			all and longest per gene
################################################


indir = sys.argv[1]
indir = indir.rstrip("/")
#print(indir)
assm = indir.split('/')[-1]
outdir = os.path.dirname(os.path.realpath(indir))
#print(outdir)


# Parse protein sequences
infile = indir + '/augustus.hints.aa'
#print(infile)
outfile1 = outdir + '/' + assm + '-proteins.fa'
print(outfile1)

fo1 = open(outfile1, 'w')
seqDict = {}
for record in SeqIO.parse(infile, "fasta"):
    gene = record.id.split(".")[0]

    fo1.write('>' + record.id + '\n')
    fo1.write(str(record.seq) + '\n')        

    if gene not in seqDict:
        seqDict[gene] = [str(record.seq), record.id]
        continue

    if len(seqDict[gene][0]) < len(record.seq):
        seqDict[gene] = [str(record.seq), record.id]
        
fo1.close()

        
outfile2 = outdir + '/' + assm + '-proteins-longest.fa'
with open(outfile2, 'w') as fo2:
    for gene in seqDict:

        fo2.write('>' + seqDict[gene][1] + '\n')
        fo2.write(seqDict[gene][0] + '\n')        



# Parse coding sequences
infile = indir + '/augustus.hints.codingseq'
outfile1 = outdir + '/' + assm + '-codingseq.fa'

fo1 = open(outfile1, 'w')
seqDict = {}
for record in SeqIO.parse(infile, "fasta"):
    gene = record.id.split(".")[0]

    fo1.write('>' + record.id + '\n')
    fo1.write(str(record.seq) + '\n')        

    if gene not in seqDict:
        seqDict[gene] = [str(record.seq), record.id]
        continue

    if len(seqDict[gene][0]) < len(record.seq):
        seqDict[gene] = [str(record.seq), record.id]

fo1.close()

        
outfile2 = outdir + '/' + assm + '-codingseq-longest.fa'
with open(outfile2, 'w') as fo2:
    for gene in seqDict:

        fo2.write('>' + seqDict[gene][1] + '\n')
        fo2.write(seqDict[gene][0] + '\n')        


