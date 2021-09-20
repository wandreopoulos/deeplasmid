#!/bin/bash -l

echo "Running chrsketch"

FASTA=$1

#To sketch a new db:  sketch.sh in=x.faa out=x.sketch amino persequence

shifter --image=bryce911/bbtools comparesketch.sh   -Xmx1000m -threads=1  in=$FASTA translate ref=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_clean/classifier/dl/asafl_plasmidPred/chrProt.faa.sketch persequence

ID=$!

pkill -P $$

kill -9 $ID

#shifter --image=bryce911/bbtools stats.sh  $FASTA

trap 'kill $(jobs -p)' EXIT

