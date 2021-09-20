#!/bin/bash -l

echo "Running prodigal"

FASTA=$1

#####shifter  --image=registry.services.nersc.gov/jgi/prodigal:latest  prodigal  -a  $FASTA.gene.faa  -d  $FASTA.gene.fasta  -i  $FASTA  -o  $FASTA.prodigal.out  -p meta

shifter --image=ramanik/hmmer:3.2.1 hmmsearch --noali --cut_nc -o  $FASTA.out_pfam --domtblout $FASTA.domtblout --cpu 16  Pfam-A.hmm  $FASTA.gene.faa


ID=$!

pkill -P $$

kill -9 $ID

#shifter --image=bryce911/bbtools stats.sh  $FASTA

trap 'kill $(jobs -p)' EXIT

