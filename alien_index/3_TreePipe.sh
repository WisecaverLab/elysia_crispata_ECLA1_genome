#! /usr/local/bin/bash

### USER SPECIFIED VARIABLES ###
QUERY='query/Elysia_final_proteins.fa'
GENELIST='query/hgt_candidates.txt'
JOBS=`wc -l $GENELIST | cut -d ' ' -f 1`
MAXJOBS=7
DB="dbs/refseq_release207_w_supp_cd95.fa"				# fasta database to search against (must be preindexed using esl-sfetch)
ITOLDIR='Batch'
ITOLID=$ITOLID								# https://itol.embl.de/help.cgi#batch
TAXDIR="../../dbs"
COLORSCHEME='microbial_euk'

### SLURM variables ###
USER=$EMAIL				# user email; either enter manually of add email to path, eg: export EMAIL='youremail@purdue.edu'
QUEUE='jwisecav'			# user queue
THREADS=24				
TIME='24:00:00'			# max walltime

# tree building parameters
SEARCHPROGRAM='diamond'			# can be 'phmmer' or 'diamond'
MAXEVALUE='1e-10'
MINSEQ=5						# minimum number of homologs to print to fasta file 
MAXSEQ=200						# maximum number of homologs to print to fasta file
MINLEN=50																				# minimum length of alignment. If shorter after trimming, tree will not be built.
MAFFTOPT="--reorder --bl 30 --op 1.0 --maxiterate 1000 --retree 1 --genafpair --quiet"	# mafft options
TRIMOPT="-gappyout"																		# trimal options
TREETYPE='iqtree'																	    # type of tree:'raxml', 'iqtree', or 'both'
RAXMLOPT='-f a -m PROTGAMMAAUTO -x 12345 -p 12345 -N 100'							    # raxml options
IQTREEOPT='-alrt 1000 -bb 1000'
RXSUPPORT=80
IQSUPPORT=95

# Environmental variables. You can specify the path to the executables if not in your $PATH
MAFFT='mafft'				# path to mafft executable 
TRIMAL='trimal'				# path to trimal executable 
RAXML='raxmlHPC-PTHREADS-AVX'		# path to raxml executable
IQTREE='iqtree'

### Directory information ###
SLURMDIR="$PWD/slurm-out"
SCRIPTS="$PWD/scripts"
PARSEDIR="$PWD/diamond-out"
if [[ "$SEARCHPROGRAM" == 'phmmer' ]]; then
	PARSEDIR="$PWD/phmmer-out"
fi

echo "Creating working directories"
if [ ! -d "$SLURMDIR" ]; then mkdir  $SLURMDIR; fi
if [ ! -d "tree-out" ]; then mkdir  tree-out; fi
if [ ! -d "tree-out/fasta" ]; then mkdir  tree-out/fasta; fi
if [ ! -d "tree-out/tree" ]; then mkdir  tree-out/tree; fi
if [ ! -d "tree-out/mafft" ]; then mkdir  tree-out/mafft; fi
if [ ! -d "tree-out/trimal" ]; then mkdir  tree-out/trimal; fi

export QUERY GENELIST DB ITOLDIR SCRIPTS PARSEDIR ITOLID COLORSCHEME TAXDIR THREADS
export MAXEVALUE MINSEQ MAXSEQ MINLEN MAFFTOPT TRIMOPT TREETYPE RAXMLOPT IQTREEOPT RXSUPPORT IQSUPPORT
export MAFFT TRIMAL RAXML IQTREE

##  2 - Submit jobs to align and trim sequences and build phylogeny
echo 'Submitting jobs for tree building'
echo "QUERY: $QUERY"
echo "GENELIST: $GENELIST"
echo "ITOL COLORSCHEME: $COLORSCHEME"

echo sbatch --array=1-$JOBS%$MAXJOBS --account=$QUEUE --job-name=buildtrees --mail-user=$USER --nodes=1 --ntasks=$THREADS --time=$TIME --output=$SLURMDIR/%x-%j.out $SCRIPTS/run_buildtree.sub
sbatch --array=1-$JOBS%$MAXJOBS --account=$QUEUE --job-name=buildtrees --mail-user=$USER --nodes=1 --ntasks=$THREADS --time=$TIME --output=$SLURMDIR/%x-%j.out $SCRIPTS/run_buildtree.sub


