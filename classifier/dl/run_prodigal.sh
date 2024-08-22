#!/bin/bash -l

#set -m

trap 'pkill -P $$' EXIT

#echo "Running prodigal"

FASTA=$1

#shifter  --image=registry.services.nersc.gov/jgi/prodigal:latest  prodigal  -a  $FASTA.gene.faa  -d  $FASTA.gene.fasta  -i  $FASTA  -o  $FASTA.prodigal.out  -p meta

#SEQ=$2
#rm -f run_prodigal.sh.$FASTA
##mkfifo $FASTA.run_prodigal.sh
#echo ">$FASTA" > $FASTA.run_prodigal.sh
#echo $SEQ >> $FASTA.run_prodigal.sh &
#cp /dev/stdin run_prodigal.sh.$FASTA


./Prodigal/prodigal  -a  $FASTA.gene.faa  -d  $FASTA.gene.fasta  -i  $FASTA  -o  $FASTA.prodigal.out  -p meta   

#rm -f run_prodigal.sh.$FASTA

#& ID=$! ; fg

pkill -P $$

#kill -9 $ID

#shifter --image=bryce911/bbtools stats.sh  $FASTA

#trap 'kill $(jobs -p)' EXIT

exit 0

