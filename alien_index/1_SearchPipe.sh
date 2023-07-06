#! /usr/local/bin/bash

### USER SPECIFIED VARIABLES ###
QUERY='query/Elysia_final_proteins.fa'				# path to input protein fasta file
JOBS=3					# Split query file into subjobs that have no more than 20k sequences each. EG if your query has 28,000 sequences, set JOBS equal to 2 (or greater); if query has 63,000 sequences, set JOBS equal to 4 (or greater).
MAXJOBS=3 				# max number of jobs to run at a time

# SLURM variables
USER=$EMAIL				# user email; either enter manually or add $EMAIL to your .bashrc
QUEUE='jwisecav'		# user queue
THREADS=16
TIME='72:00:00'			# max walltime

# HMMER/DIAMOND variables
PROGRAM='diamond'								# can be 'phmmer' or 'diamond'
DB="dbs/refseq_release207_w_supp_cd95.dmnd"		# path to database 
HMMEROPTS='--F1 1e-5 --F2 1e-7 --F3 1e-10 --noali --notextw -o /dev/null' 	# HMMER Options
DMNDOPTS='blastp --sensitive --evalue 1e-3 --max-target-seqs 0'				# Diamond blast options

# Parsing variables
MAX_HITS_PER_QUERY=0 		# maximum number of hits to keep per query sequence (use 0 to retain all hits)
MAX_HITS_PER_SUBJECT=0 		# maximum number of hits to keep per uniq taxnomy id (use 0 to retain all hits)
EVALUE_THRESHOLD=1			# upper bound for evalue to include in parsed output 

### Directory information ###
QUERYDIR="$PWD/query"
SPLITDIR="$PWD/split"
SCRIPTS="$PWD/scripts"
HMMRDIR="$PWD/phmmer-out"		# output directory for phmmer
DMNDIR="$PWD/diamond-out"		# output directory for diamond or diamond-parsed phmmer results
SLURMDIR="$PWD/slurm-out"

echo "Creating working directories"
if [ ! -d "$SPLITDIR" ]; then mkdir $SPLITDIR; fi
if [ $PROGRAM == 'phmmer' ]; then if [ ! -d "$HMMRDIR" ]; then mkdir $HMMRDIR; fi; fi
if [ $PROGRAM == 'diamond' ]; then if [ ! -d "$DMNDIR" ]; then mkdir $DMNDIR; fi; fi
if [ ! -d "$SLURMDIR" ]; then mkdir  $SLURMDIR; fi

export QUERY THREADS HMMRDIR DMNDIR
export PROGRAM DB HMMEROPTS DMNDOPTS JOBS PWD
export MAX_HITS_PER_QUERY MAX_HITS_PER_SUBJECT EVALUE_THRESHOLD

##  1 - Create job array input files
echo "Creating job array input files"
module load conda-env/env.genomics
python $SCRIPTS/create_jobarray_fasta.py $QUERY $SPLITDIR $JOBS

##  2 - Submit query self-search
echo "Submitting query self-search job"
CMD=`sbatch --account=$QUEUE --job-name=selfsearch --mail-user=$USER --nodes=1 --ntasks=128 --time=$TIME --output=$SLURMDIR/%x-%j.out $SCRIPTS/run_selfsearch.sub`
FIRST=`echo $CMD | sed "s/Submitted batch job //"`

##  3 - Submit phmmer/diamond jobs against database
echo "Submitting phmmer/diamond jobs: 1-$JOBS max running: $MAXJOBS $SCRIPTS/run_search.sub"
SECOND=`sbatch --dependency=afterok:$FIRST --array=1-$JOBS%$MAXJOBS --account=$QUEUE --job-name=search --mail-user=$USER --nodes=1 --ntasks=128 --time=$TIME --output=$SLURMDIR/%x-%A-%a.out $SCRIPTS/run_search.sub`
