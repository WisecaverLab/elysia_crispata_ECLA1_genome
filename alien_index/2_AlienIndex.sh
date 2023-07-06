#! /usr/local/bin/bash

### USER SPECIFIED VARIABLES ###
BASE="query/alien_index"				#base name for outputfile
SEARCHPROGRAM='diamond'			# can be 'phmmer' or 'diamond'
TAXDIR="dbs/"
LINFILE="dbs/lineages.dmp"
RECIPIENT_LINEAGE='PLACOBRANCHOIDEA,71491' 					# ncbi taxonomy id of clade which encompases the query sequences
ANCESTRAL_LINEAGE='METAZOA,33208'					# ncbi taxonomy id of clade that vertically inherited query sequences should group within
NUM_TOP_HITS=200						# number of top hits to summarize (default = 100)
SKIPLIST="scripts/skiplist.txt"								# optional comma-delimited list of taxonomy ids to skip (e.g. genomes you suspect of having contamination)

### Directory information ###
SCRIPTS="$PWD/scripts"
SLURMDIR="$PWD/slurm-out"
PARSEDIR="$PWD/diamond-out"
if [[ "$SEARCHPROGRAM" == 'phmmer' ]]; then
        PARSEDIR="$PWD/phmmer-out"
fi

# SLURM variables
USER=$EMAIL				# user email; either enter manually of add email to path, eg: export EMAIL='youremail@purdue.edu'
QUEUE='jwisecav'			# user queue
TIME='24:00:00'			# max walltime

echo "Creating working directories"
if [ ! -d "$SLURMDIR" ]; then mkdir  $SLURMDIR; fi

export BASE RECIPIENT_LINEAGE ANCESTRAL_LINEAGE NUM_TOP_HITS SKIPLIST SCRIPTS PARSEDIR TAXDIR LINFILE

##  1 - Submit alien index job
echo "Submitting alien index job"
sbatch --account=$QUEUE --job-name=alienindex --mail-user=$USER --nodes=1 --ntasks=128 --time=$TIME --output=$SLURMDIR/%x-%j.out $SCRIPTS/run_alienindex.sub
