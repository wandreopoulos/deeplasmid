#!/bin/bash -l

echo "Running pentamer.sh"

FASTA=$1

shifter  --image=bryce911/bbtools  commonkmers.sh  -Xmx1000m -threads=1  in=$FASTA  out=stdout  k=5  display=3  count

ID=$!

pkill -P $$

kill -9 $ID

#shifter --image=bryce911/bbtools stats.sh  $FASTA

trap 'kill $(jobs -p)' EXIT
