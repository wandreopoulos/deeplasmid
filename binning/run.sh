#!/bin/bash

set -e

export PYTHONPATH=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/assemblyqc/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/readqc/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/tools/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc/

FASTA=`cat ./run.config | sed -n 1p`

LIBFQ=`cat ./run.config | sed -n $(($SGE_TASK_ID+1))p`

OUTDIR=`pwd`

/global/projectb/sandbox/rqc/andreopo/src/bitbucket/rnaseq_microbe/bin/rnaseq_microbe -o  $OUTDIR/$LIBFQ  -r  $FASTA   -g  /global/projectb/scratch/andreopo/test_rnaseq_microbe/SEQQC-1807/Transcriptome_Analysis/2517572208/2517572208.gff  -n

