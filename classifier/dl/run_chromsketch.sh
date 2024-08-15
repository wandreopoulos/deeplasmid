#!/bin/bash -l
#set -m

trap 'pkill -P $$' EXIT


#echo "Running chrsketch"

FASTA=$1

#To sketch a new db:  sketch.sh in=x.faa out=x.sketch amino persequence

#shifter --image=bryce911/bbtools comparesketch.sh   -Xmx1000m -threads=1  in=$FASTA translate ref=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_clean/classifier/dl/asafl_plasmidPred/chrProt.faa.sketch persequence

#SEQ=$2
#rm -f run_chromsketch.sh.$FASTA
##mkfifo $FASTA.run_chromsketch.sh
#echo ">$FASTA" > $FASTA.run_chromsketch.sh
#echo "$SEQ" >> $FASTA.run_chromsketch.sh &
#cp /dev/stdin run_chromsketch.sh.$FASTA

./bbmap/comparesketch.sh   -Xmx1000m -threads=1  in=/dev/shm/$FASTA translate ref=./asafl_plasmidPred/chrProt.faa.sketch persequence   

#rm -f run_chromsketch.sh.$FASTA
#& ID=$! ; fg

pkill -P $$

#kill -9 $ID

#shifter --image=bryce911/bbtools stats.sh  $FASTA

#trap 'kill $(jobs -p)' EXIT

exit 0


