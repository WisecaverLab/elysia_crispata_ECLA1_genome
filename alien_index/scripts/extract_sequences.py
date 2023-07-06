import glob
import subprocess
import sys
import re
import getopt

# queryfile = '/depot/jwisecav/data/jwisecav/endophytes_Feb23/phylopipe_Xylaria_intraflava_YMJ_725_-_Xylint1/query/Xylaria_intraflava_YMJ_725_-_Xylint1.aa.fasta'
# qseqid = 'Xylint1_518656'
# indir = '/depot/jwisecav/data/jwisecav/endophytes_Feb23/phylopipe_Xylaria_intraflava_YMJ_725_-_Xylint1/normbit-out'
# database = '/depot/jwisecav/data/dbs/custom-refseq-current-release/refseq_w_supp_cd95_special_xylariaceae.fa'
# outdir = '/depot/jwisecav/data/jwisecav'
# minseq = 1
# maxseq = 100
# maxevalue = 1e-20

###############################################
# Usage: python extract_sequences.py -i [indir] -s [query sequence] -q [query file] -d [database file]-o [outdir] [options]
#
# Jen Wisecaver
# 20201018
################################################

def usage():
    print('\n  Usage: '+sys.argv[0]+' -i <indir> -s <query sequence id> -q <query fasta file> -d <database fasta file> -o <outdir> [options]')
    print("    -i|--indir <FILENAME>  path to directory containing diamon/phmmer output with normalized bitscores")
    print("    -s|--qseqid <FILENAME>  target sequence id")
    print("    -q|--qfasta <FILENAME> path to query fasta file (must be pre indexed with esl-sfetch)")
    print("    -d|--dbfasta <FILENAME> path to database fasta file (must be pre indexed with esl-sfetch)")
    print("    -o|--outdir <FILENAME> path to output directory")
    print("\n    OPTIONAL:")
    print("    -n|--minhits <INTEGER> minimum number of total sequences to print to fasta file (default = 1)")
    print("    -x|--maxhits <INTEGER> maximum number of database homologs to print to fasta file (default = 100)")
    print("    -e|--evalue <FLOAT> maximum evalue for sequence to be included in downstream analysis (default = 1e-10)\n\n")

# Read in command line arguments
try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hi:s:q:d:o:n:x:e:', ['help', 'indir=', 'qseqid=', 'qfasta=', 'dbfasta=', 'outdir=', 'minhits=', 'maxhits=', 'evalue='])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err)  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

minseq = 1
maxseq = 100
maxevalue = 1e-10

for opt, arg in options:
    if opt in ('-h', '--help'):
        usage()
        sys.exit()
    elif opt in ('-o', '--outdir'):
        outdir = arg
    elif opt in ('-i', '--indir'):
        indir = arg
    elif opt in ('-s', '--qseqid'):
        qseqid = arg
    elif opt in ('-q', '--qfasta'):
        queryfile = arg
    elif opt in ('-d', '--dbfasta'):
        database = arg
    elif opt in ('-x', '--maxhits'):
        maxseq = int(arg)
    elif opt in ('-n', '--minhits'):
        minseq = int(arg)
    elif opt in ('-e', '--evalue'):
        maxevalue = float(arg)
command = " ".join(sys.argv)

try:
    indir
    outdir
    qseqid
    queryfile
    database
except NameError:
    usage()
    sys.exit(2)
else:
    sys.stdout.write('\nCOMMAND: python ' + command + '\n\n')

#############################
command = "grep -P '^" + qseqid + "' " + indir + "/normbitout* | head -n 1 | cut -f 1"
sys.stdout.write('Selecting infile:\n' + command + '\n')
infile = subprocess.getoutput(command).split(':')[0]
sys.stdout.write('\nInput file: ' + infile + '\n')

command = 'esl-sfetch ' + queryfile + ' ' + qseqid

sys.stdout.write('\nExtracting query sequence:\n' + command + '\n')
qfasta = subprocess.getoutput(command)
#print(qfasta)
qfasta = re.sub('>', '>QUERY-', qfasta)
qfasta = re.sub('\*$', '', qfasta)
qseq = re.sub('^>.*', '', qfasta)
qseq = re.sub('\n', '', qseq)

sys.stdout.write('\n' + qfasta + '\n')

hitDict = {}
queryDict = {}

# Parse diamond/phmmer
fi = open(infile)

for line in fi:
    seqCounter = len(hitDict)
    if seqCounter >= maxseq:
        break
    
    col = line.rstrip().split('\t')
    if col[0] != qseqid:
        if seqCounter > 1:
            break
        else:
            continue
    
    sseqid = col[1]
    evalue = float(col[3])
    
    if qseqid == sseqid:
        continue
        
    if evalue > maxevalue:
        sys.stdout.write('\tSkipping ' + sseqid + ' : poor score (evalue' + str(evalue) + ' > max ' + str(maxevalue) + ')\n')
        continue

    command = 'esl-sfetch ' + database + ' ' + sseqid
    fasta = subprocess.getoutput(command)

	# if searching for the hit sequence in the database does not work
	# The sequence must come from query file
    if fasta[0] != '>':
        command2 = 'esl-sfetch ' + queryfile + ' ' + sseqid
        fasta = subprocess.getoutput(command2)
        fasta = re.sub('\*$', '', fasta)
        queryDict[sseqid] = fasta
                    
    else:
        fasta = re.sub('\*$', '', fasta)
        hitDict[sseqid] = fasta
        

fi.close()

if len(hitDict) < minseq:
    sys.stdout.write('\nSkipping Query ' + qseqid + ' : retained hits ' + str(len(hitDict)) +  ' < minimum (' + str(minseq) + ')\n')
else:
    outfile = outdir + '/' + qseqid + '.fa'
    seqCounter = len(queryDict) + len(hitDict)
    sys.stdout.write('\nprinting ' + str(seqCounter+1) + ' fulllength sequences to ' + outfile + '..........\n')

    fo = open(outfile, 'w')
    fo.write(qfasta + '\n')

    for sseqid in queryDict:
        fo.write(queryDict[sseqid] + '\n')

    for sseqid in hitDict:
        fo.write(hitDict[sseqid] + '\n')
        
        
    fo.close()
    