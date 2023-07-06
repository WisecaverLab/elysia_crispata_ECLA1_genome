import os
import os.path
import subprocess
import sys
import getopt
import time
import glob

def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print ('command failed')
        print (cmd)
        sys.exit(1)

###############################################
# Usage: python parse_diamond.py -i [diamondfile] -s [self_diamond] -o [outfile] [options]
#
# Jen Wisecaver
# 20201018
################################################

def usage():
    print('\n  Usage: '+sys.argv[0]+' -i <diamondfile> -s <self_diamond file> -o <outfile> [options]')
    print("    -i|--infile <FILENAME>  path to diamond file")
    print("    -s|--selfsearch <FILENAME>  path to self diamond file")
    print("    -o|--outfile <FILENAME> path to output file")
    print("\n    OPTIONAL:")
    print("    -n|--maxhits_per_query <INTEGER> max hits to retain per query sequence (default = 0 for no limit)")
    print("    -m|--maxhits_per_subject <INTEGER> max hits to retain per unique taxonomy id (default = 0 for no limit)")
    print("    -e|--evalue_threshold <FLOAT> max evalue to retain (default = 1)\n\n")

# Read in command line arguments
try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hi:s:o:n:m:e:', ['help', 'infile=', 'selfsearch=', 'outfile=', 'maxhits_per_query=', 'maxhits_per_subject=', 'evalue_threshold='])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err)  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

maxhits_per_query = 0 
maxhits_per_subject = 0 
evalue_threshold = 1 

for opt, arg in options:
    if opt in ('-h', '--help'):
        usage()
        sys.exit()
    elif opt in ('-o', '--outfile'):
        outfile = arg
    elif opt in ('-i', '--infile'):
        infile = arg
    elif opt in ('-s', '--selfsearch'):
        selfsearch = arg
    elif opt in ('-n', '--maxhits_per_query'):
        maxhits_per_query = int(arg)
    elif opt in ('-m', '--maxhits_per_subject'):
        maxhits_per_subject = int(arg)
    elif opt in ('-e', '--evalue_threshold'):
        evalue_threshold = float(arg)
command = " ".join(sys.argv)

try:
    infile
    selfsearch
    outfile
except NameError:
    usage()
    sys.exit(2)
else:
    sys.stdout.write('\nCOMMAND: ' + command + '\n\n')

################################################
# Store max possible bitscore per query sequence
################################################
sys.stdout.write('Saving max bitscores from ' + selfsearch + '\n')

queryDict = {}
sortResults = {}

fi = open(selfsearch)

for line in fi:
    if line[0] == '#':
        continue
    
    col = line.rstrip().split('\t')
    qseqid = col[0]
    sseqid = col[1]
    bitscore = float(col[11])
    
    if qseqid == sseqid:
        if qseqid not in queryDict:
            queryDict[qseqid] = bitscore
            
        else:
            if queryDict[qseqid] < bitscore:
                queryDict[qseqid] = bitscore

fi.close()

selfCount = len(queryDict)
sys.stdout.write('\tstored max bitscores for ' + str(selfCount) + ' queries \n\n')

################################################
# Calculate normalized bitscore for all phmmer hits to database
################################################

queryCounter = {}
hitCounter = {}

fi = open(infile)

t = time.time()
tmpdir = '/tmp/' + 'PHYLO' +str(t)

if os.path.exists(tmpdir) == False:
    cmd = 'mkdir ' + tmpdir
    sys.stdout.write('Creating output directory:' + tmpdir + '...\n')
    runCMD(cmd)
else:
    sys.stdout.write('Writing to existing output directory:' + tmpdir + '...\n')

current = 'n/a'

sys.stdout.write('Parsing phmmer hits from ' + infile + '\n')

for line in fi:
    if line[0] == '#':
        continue

    col = line.rstrip().split('\t')
    rest = '\t'.join(col[2:])

    qseqid = col[0]
    
    if qseqid not in queryCounter: 
        queryCounter[qseqid] = 1
        hitCounter[qseqid] = {}
    else:
        queryCounter[qseqid] += 1
    
    if qseqid != current: 
        if queryCounter[qseqid] > 1:
            fo.close()

        current = qseqid
        tmpfile = tmpdir + '/' + qseqid
        fo = open(tmpfile, 'w')

    sseqid = col[1]
    staxid = sseqid.split('-')[1]
    #print(staxid)

    if staxid not in hitCounter[qseqid]:
        hitCounter[qseqid][staxid] = 1
    else:
        hitCounter[qseqid][staxid] += 1
        
    bitscore = float(col[11])
    evalue = float(col[10])
    normbit = 0
    
    if evalue > evalue_threshold:
        sys.stdout.write('\tSkipping ' + qseqid + ' - ' + sseqid + ' : poor score (evalue ' + str(evalue) + ' > max ' + str(evalue_threshold) + ')\n')
        #print(line)
        continue

    if maxhits_per_query > 0:
        if queryCounter[qseqid] >= maxhits_per_query:
            sys.stdout.write('\tSkipping ' + qseqid + ' - ' + sseqid + ' : reached max hits (' + maxhits_per_query + ') for query\n')
            continue    
    
    if maxhits_per_subject > 0:
        if hitCounter[qseqid][staxid] >= maxhits_per_subject:
            sys.stdout.write('\tSkipping ' + qseqid + ' - ' + sseqid + ' : reached max hits (' + maxhits_per_subject + ') for hit taxid\n')
            continue    

    if qseqid in queryDict:
        normbit = bitscore / queryDict[qseqid]
        
        if normbit > 1 and normbit < 1.1:
            normbit = 1
        elif normbit > 1.1:
            sys.stdout.write('\tWARNING: normalized bitscore(' + str(normbit) + ') greater than 1 : ' + qseqid + ' - ' + sseqid + '\n')
            normbit = 1
    else:
        sys.stdout.write('\tWARNING: no max bitscore stored for' + qseqid + '\n')
        
    fo.write(sseqid + '\t' + str(normbit) + '\t' + str(evalue) + '\t' + rest + '\n')
    

fi.close()
fo.close()



sys.stdout.write('\nParsing phmmer hits from self search file' + selfsearch + '\n')
fi = open(selfsearch)
current = 'n/a'

for line in fi:
    if line[0] == '#':
        continue
    
    col = line.rstrip().split('\t')
    rest = '\t'.join(col[2:])
    
    qseqid = col[0]
    if qseqid not in queryCounter:
        continue
        
    if qseqid != current: 
        if current != 'n/a':
            fo.close()

        current = qseqid
        tmpfile = tmpdir + '/' + qseqid
        fo = open(tmpfile, 'a')
    
    sseqid = col[1]
    if qseqid == sseqid:
        continue
        
    bitscore = float(col[11])
    evalue = float(col[10])
    normbit = 0
    
    if evalue > evalue_threshold:
        sys.stdout.write('\tSkipping ' + qseqid + ' - ' + sseqid + ' : poor score (evalue ' + str(evalue) + ' > max ' + str(evalue_threshold) + ')\n')
        continue

    if qseqid in queryDict:
        normbit = bitscore / queryDict[qseqid]
        
        if normbit > 1 and normbit < 1.1:
            normbit = 1
        elif normbit > 1.1:
            sys.stdout.write('\tWARNING: normalized bitscore(' + str(normbit) + ') greater than 1 : ' + qseqid + ' - ' + sseqid + '\n')
            normbit = 1
    else:
        sys.stdout.write('\tWARNING: no max bitscore stored for' + qseqid + '\n')
        
    fo.write(sseqid + '\t' + str(normbit) + '\t' + str(evalue) + '\t' + rest + '\n')

fi.close()
fo.close()


fo = open(outfile, 'w')

for infile in glob.glob(tmpdir + '/*'):
    qseqid = infile.split('/')[-1]
    
    bitDict = {}
    lineDict = {}
    
    fi = open(infile)
    
    for line in fi:
        sseqid = line.split('\t')[0]
        normbit = line.split('\t')[1]
        if sseqid not in bitDict:
            bitDict[sseqid] = normbit
            lineDict[sseqid] = line
        else:
            if normbit > bitDict[sseqid]:
                bitDict[sseqid] = normbit
                lineDict[sseqid] = line
    fi.close()
    
    sort_normbits = sorted(bitDict.items(), key=lambda x: x[1], reverse=True)

    for sseqid in sort_normbits:
        #print(sseqid)
        line = lineDict[sseqid[0]]
        fo.write(qseqid + '\t' + line)
        
fo.close()

