#!/bin/bash -l

echo "Running prodigal"

FASTA=$1

shifter  --image=registry.services.nersc.gov/jgi/prodigal:latest  prodigal  -a  $FASTA.gene.faa  -d  $FASTA.gene.fasta  -i  $FASTA  -o  $FASTA.prodigal.out  -p meta


ID=$!

pkill -P $$

kill -9 $ID

#shifter --image=bryce911/bbtools stats.sh  $FASTA

trap 'kill $(jobs -p)' EXIT

