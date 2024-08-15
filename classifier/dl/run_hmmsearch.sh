#!/bin/bash -l
#set -m

trap 'pkill -P $$' EXIT


#echo "Running hmmsearch"

FASTA=$1

if [ -z "$FASTA" ]; then
    echo "Error: No FASTA file provided."
    exit 1
fi

#####shifter  --image=registry.services.nersc.gov/jgi/prodigal:latest  prodigal  -a  $FASTA.gene.faa  -d  $FASTA.gene.fasta  -i  $FASTA  -o  $FASTA.prodigal.out  -p meta

if `./hmmer-3.3.2/src/hmmsearch --noali --cut_nc -o $FASTA.out_pfam --domtblout $FASTA.domtblout --cpu 16 Pfam-A.TMP2.hmm $FASTA.gene.faa`; then
    echo "hmmsearch completed successfully."
else
    echo "Error: hmmsearch failed."
    exit 1
fi

#pkill -P $$

#kill -9 $ID

#shifter --image=bryce911/bbtools stats.sh  $FASTA

#trap 'kill $(jobs -p)' EXIT

exit 0


