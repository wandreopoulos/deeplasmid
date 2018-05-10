#!/bin/bash -l

FASTA=$1

OUT=$2

DATETIME=`date "+%Y-%m-%d_%H:%M:%S"`

read_fasta2_plasmids.py  -i $FASTA -o $OUT/$DATETIME.$FASTA  > $OUT/$DATETIME.$FASTA.outerr 2>&1

mkdir $OUT/dlDataPath
chmod 777 $OUT/dlDataPath

#pass fasta and yml directory to format
format_Test.py  --inputfasta $FASTA --inputyml $OUT/$DATETIME.$FASTA/yml --dataPath $OUT/dlDataPath

predict_Plasmid.py --dataPath $OUT/dlDataPath  --given test  --kModelList  20-29  --dataSegment -1


