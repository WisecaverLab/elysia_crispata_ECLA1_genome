import sys
import getopt
import os.path
import glob

###############################################
# Usage: python calculate_alien_index.py -i [indir] -o [outname] -t [taxdir] -g [group] -k [recipient] [options]
#
# Jen Wisecaver
# 20201018
################################################

def usage():
    print('\n  Usage: '+sys.argv[0]+' -i [indir] -o [outname] -t [taxdir] -g [group] -k [recipient] [options]')
    print("    -i|--indir <DIRNAME> path to directory containing parsed diamond/phmmer output")
    print("    -o|--outname <FULL_FILENAME> base for outfile; recipient and ancestral lineages will be appended ")
    print("    -l|--linfile <FULL_FILENAME> path to custome 'lineages.dmp' ")
    print("    -t|--taxdir <DIRNAME> path to directory containing the ncbi taxonomy files 'nodes.dmp' 'merged.dmp' and custom 'lineages.dmp'")
    print("    -a|--ancestralid <INTEGER> ncbi taxonomy id of clade that vertically inherited query sequences should group within")
    print("    -r|--recipientid <INTEGER> ncbi taxonomy id of clade which encompases the query sequences")
    print("\n    OPTIONAL:")
    print("    -n|--num_top_hits <INTEGER> number of top unique taxonomy ids to summarize (default = 100)")
    print("    -k|--skipid <LIST> comma-delimited list of taxonomy ids to skip")
    print("\n    NOTE: ")
    print("    ANCESTRALID and RECIPIENTID may be supplied as a common-delimited list of non-overlapping NCBI taxonomy ids")
    print("    The leftmost item in the list must be a name for this user-definied group (I typically write this name in CAPS)")
    print("         For example: -a ARCHAEPLASTIDA,33090,38254,2763")
    print("    You can also use this feature to give the ANCESTRAL and RECIPIENT lineages for human interpretable names")
    print("         EUKARYOTA,2759 \n\n")

################################
# Read in command line arguments
#
try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hi:o:l:t:a:r:n:k:', ['help', 'indir=', 'outname=', 'linfile=', 'taxdir=', 'ancestralid=', 'recipientid=', 'num_top_hits=', 'skipid='])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err)  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

num_top_hits = 100 
skipid = 'na'

for opt, arg in options:
    if opt in ('-h', '--help'):
        usage()
        sys.exit()
    elif opt in ('-o', '--outname'):
        outname = arg
    elif opt in ('-i', '--indir'):
        indir = arg
    elif opt in ('-l', '--linfile='):
        linfile = arg
    elif opt in ('-t', '--taxdir'):
        taxdir = arg
    elif opt in ('-a', '--ancestralid'):
        ancestralid = arg
    elif opt in ('-r', '--recipientid'):
        recipientid = arg
    elif opt in ('-n', '--num_top_hits'):
        num_top_hits = int(arg)
    elif opt in ('-k', '--skipid'):
        skipid = arg
command = " ".join(sys.argv)

try:
    outname
    indir
    linfile
    taxdir
    ancestralid
    recipientid
except NameError:
    usage()
    sys.exit(2)
else:
    sys.stdout.write("\nCOMMAND: " + command + '\n\n')


################################
# Parse skip, ancestralid and recipientid strings into lists of ncbi taxids
#
ancestral_list = ancestralid.split(',')
ancestral_name = ancestralid
if len(ancestral_list) > 1:
	ancestral_name = ancestral_list.pop(0)

recipient_list = recipientid.split(',')
if len(recipient_list) > 1:
	recipient_name = recipient_list.pop(0)

skiplist = []
if os.path.isfile(skipid) :
    fi = open(skipid)
    for line in fi:
    	line = line.rstrip()
    	skiplist.append(line)
    fi.close()
else:
    skiplist = skipid.split(',')
skipid = ','.join(skiplist)

# Specify outfile
outfile = outname + '_recipient' + recipient_name + '_ancestral' + ancestral_name + '.txt';


################################
# Initialize taxonomy lookup dictionarys
#
parentDict = {}
lin = {}

nodesfile = taxdir + '/nodes.dmp'
mergedfile = taxdir + '/merged.dmp'

nfi = open(nodesfile)
for line in nfi:
    #print(line)
    col = line.rstrip().split('\t|\t')
    node = col[0]
    parent = col[1]
    #print(node,parent)
    parentDict[node] = parent
    
nfi.close()

mfi = open(mergedfile)
for line in mfi:
    col = line.rstrip().split('\t|\t')
    node = col[0]
    parent = col[1].split('\t')[0]  
    #print(node,parent)
    parentDict[node] = parent
mfi.close()

lfi = open(linfile)
for line in lfi:
    col = line.rstrip().split('\t')
    lin[col[0]] = col[1]
lfi.close()


def taxdump(taxid):
    root = 'no'
    taxlist = []
    
    while root == 'no':
        #print(taxid)
        if taxid == '':
            break
        if taxid == '1':
            root = 'yes'
        taxlist.append(taxid)
        
        if taxid not in parentDict:
        	break
        	
        taxid = parentDict[taxid]
    
    return(taxlist)


################################
# Read in normbitout files and store relevant hits
#
scores = {} # dictionary for storing the best score, best ancestral score, and best other score for each query
count = {} 
lineages = {}

infiles = glob.glob(indir + '/normbitout.*')
if len(infiles) < 1:
    sys.stdout.write('NO INPUT: normbitout.* in ' + indir + '\n\n')
    sys.exit(2)    
    
for infile in infiles:
    sys.stdout.write('parsing hits: ' + infile + '\n')
    
    fi = open(infile)
    for line in fi:
        if line[0] =='#':
            continue
        col = line.rstrip().split('\t')
        qseqid = col[0]
        
        # initialize variables
        if qseqid not in scores:
            scores[qseqid] = {}
            scores[qseqid]['top'] = ['na', 'na', 'na', 0]            
            scores[qseqid]['ancestral'] = ['na', 'na', 'na', 0]            
            scores[qseqid]['other'] = ['na', 'na', 'na', 0]            
            
            count[qseqid] = {}
            count[qseqid]['total'] = set()
            count[qseqid]['recipient'] = set()
            count[qseqid]['ancestral'] = set()
            count[qseqid]['other'] = set()
            
            lineages[qseqid] = set()
                        
        sseqid = col[1]
        # skip if self hit
        if qseqid == sseqid.split('-')[0]:
            continue
        
        evalue = col[3]
        norm_bitscore = col[2]
        
        # Get taxonomic information for subject seq
        staxid = ''
        if len(sseqid.split('-')) == 4:
            staxid = sseqid.split('-')[1]
        else:
            continue
        
        # Skip if subject is in skiplist
        if staxid in skiplist:
            continue

        # Skip line if already stored information on better hit from same taxid
        if staxid in count[qseqid]['total']:
            continue        

        count[qseqid]['total'].add(staxid)
        running_count = len(count[qseqid]['total'])

        # skip line if already stored best ancestral hit and best other hit and info for top hits
        if scores[qseqid]['ancestral'][0] is not 'na' and scores[qseqid]['other'][0] is not 'na' and running_count > num_top_hits:
            continue
        
        # Get subject seq ncbi lineage
        slin_name = 'Not_specified'
        ncbi_lineage = taxdump(staxid)
        ncbi_lineage_list = []
        is_recipient = 'no'
        is_lineage = 'no'
        genomeid = ''
        
        for taxon_id in ncbi_lineage:
            ncbi_lineage_list.append(taxon_id)
            if taxon_id in lin and slin_name == 'Not_specified':
                slin_name = lin[taxon_id]
            
            if taxon_id in ancestral_list:
                is_lineage = 'yes'
                
            if taxon_id in recipient_list:
                is_recipient = 'yes'

                        
        # Bring it all home...                
        if slin_name == 'Not_specified':
            sys.stdout.write('Subject lineage is not specified: ' + qseqid + ', ' + sseqid + ', ' + staxid + '\n')
            continue         
                
        if is_recipient == 'yes':
            if running_count <= num_top_hits:
                count[qseqid]['recipient'].add(staxid)

                
        if running_count <= num_top_hits:
            lineages[qseqid].add(slin_name)


        if scores[qseqid]['top'][0] is 'na':
            scores[qseqid]['top'] = [sseqid, slin_name, evalue, norm_bitscore]

            
        if is_lineage == 'yes' and is_recipient == 'no':
            if running_count < num_top_hits:
                count[qseqid]['ancestral'].add(staxid)
                
            if scores[qseqid]['ancestral'][0] is 'na':
                scores[qseqid]['ancestral'] = [sseqid, slin_name, evalue, norm_bitscore]

                
        if is_lineage == 'no' and is_recipient == 'no':
            if running_count <= num_top_hits:
                count[qseqid]['other'].add(staxid)
                
            if scores[qseqid]['other'][0] is 'na':
                scores[qseqid]['other'] = [sseqid, slin_name, evalue, norm_bitscore]
        
    fi.close()

sys.stdout.write('\nFinished reading diamond output in: ' + indir + '\n\n')


fo = open(outfile, 'w')
fo.write('# RECIPIENT: ' + recipientid + '\n')
fo.write('# ANCESTRAL LINEAGE: ' + ancestralid + '\n')
fo.write('# NUMBER OF HITS IN \'TOP HITS\': ' + str(num_top_hits) + '\n')
fo.write('# SKIPPING TAXIDS: ' + skipid + '\n')
fo.write('# query\tbest hit\tsubclade of best hit\tbest hit evalue\tbest hit norm bitscore\t')
fo.write('best hit within ancestral lineage\tsubclade of best hit within ancestral lineage\tbest hit within ancestral lineage evalue\tbest hit within ancestral lineage norm bitscore\t')
fo.write('best hit outside ancestral lineage\tsubclade of best hit outside ancestral lineage\tbest hit outside ancestral lineage evalue\tbest hit outside ancestral lineage norm bitscore\t')
fo.write('alien index\tunique lineage(s) in top hits\tno. of unique lineage(s) in top hits\tno. recipient seqs in top hits\tno. ancestral lineage seqs in top hits\tno. other seqs in top hits\tno. total hits\n')

sys.stdout.write('# RECIPIENT: ' + recipientid + '\n')
sys.stdout.write('# ANCESTRAL LINEAGE: ' + ancestralid + '\n')
sys.stdout.write('# NUMBER OF HITS IN \'TOP HITS\': ' + str(num_top_hits) + '\n')
sys.stdout.write('# SKIPPING TAXIDS: ' + skipid + '\n')
sys.stdout.write('\nWriting results to ' + outfile + '\n')

for qseqid in scores:
    lin_count = len(lineages[qseqid])
    lin_list = ', '.join(lineages[qseqid])
    
    
    total_count = 0
    if len(count[qseqid]['total']) <= num_top_hits:
        total_count = len(count[qseqid]['total'])
    elif len(count[qseqid]['total']) > num_top_hits: 
    	total_count = '>' + str(num_top_hits)
                
    if scores[qseqid]['top'][3] == 0:
        continue
    
    ancestral_bitscore = scores[qseqid]['ancestral'][3]
    other_bitscore = scores[qseqid]['other'][3]
    
    ai = float(other_bitscore) - float(ancestral_bitscore)
    #print(qseqid, ai)
    
    line = [qseqid] + scores[qseqid]['top'] + scores[qseqid]['ancestral'] + scores[qseqid]['other'] + [ai, lin_list, lin_count, len(count[qseqid]['recipient']), len(count[qseqid]['ancestral']), len(count[qseqid]['other']), total_count]    
    line = [str(i) for i in line]
    line = '\t'.join(line)
    fo.write(line + '\n')

fo.close()        

