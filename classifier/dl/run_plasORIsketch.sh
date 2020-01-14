#!/bin/bash -l
#set -m

trap 'pkill -P $$' EXIT


echo "Running plassketch"

FASTA=$1

#To sketch a new db:  sketch.sh in=x.faa out=x.sketch persequence

#shifter --image=bryce911/bbtools comparesketch.sh   -Xmx1000m -threads=1  in=$FASTA ref=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_clean/classifier/dl/asafl_plasmidPred/plasmid_originOfReplication.nr.fasta.sketch persequence

/srv/jgi-ml/classifier/dl/bbmap/comparesketch.sh   -Xmx1000m -threads=1  in=$FASTA ref=/srv/jgi-ml/classifier/dl/asafl_plasmidPred/plasmid_originOfReplication.nr.fasta.sketch persequence   

#& ID=$! ; fg

pkill -P $$

#kill -9 $ID

#shifter --image=bryce911/bbtools stats.sh  $FASTA

#trap 'kill $(jobs -p)' EXIT

exit 0


