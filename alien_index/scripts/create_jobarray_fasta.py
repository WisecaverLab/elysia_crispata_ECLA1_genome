from Bio import SeqIO
import sys 
import os  

###############################################
# create_jobarray_fasta.py 
# Jen Wisecaver
# 20201018 
# input: a fasta file 
# output: a directory with fasta files in jobarray
#         format 
################################################

infile = sys.argv[1]
outdir = sys.argv[2]
jobs = int(sys.argv[3]) #number of fasta files to create
filename = 'query' # name for each subfile

if not os.path.exists(outdir):
    print('Creating new directory for split fastas:', outdir)
    os.makedirs(outdir)
else:
    print('Writing split fastas to existing directory:', outdir)
    
totalSeqs = 0
for seq_record in SeqIO.parse(infile, "fasta"):
    totalSeqs += 1
    
fileCount = 0

seqDict = {}

for seq_record in SeqIO.parse(infile, "fasta"):
    fileCount += 1
    if fileCount > jobs:
        fileCount = 1

    outfile = outdir + '/' + filename + '.' + str(fileCount)
    if outfile not in seqDict:
        seqDict[outfile] = {}
    
    header = seq_record.description
    sequence = seq_record.seq
    
    seqDict[outfile][header] = sequence
    
for outfile in seqDict:
    fo = open(outfile, 'w')
    
    for header in seqDict[outfile]:
        sequence = seqDict[outfile][header]
        fo.write('>' + header + '\n' + str(sequence) + '\n')
    
    fo.close()
    
    