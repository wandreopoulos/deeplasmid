#!/bin/bash -l

echo "Running feature_DL_plasmid_predict.sh
Plasmid finder
Author: Bill Andreopoulos
Last maintained: May 10, 2018"

FASTA=$1

OUT=$2

DATETIME=`date "+%Y-%m-%d_%H:%M:%S"`

read_fasta2_plasmids.py  -i $FASTA -o $OUT/$DATETIME.$FASTA  > $OUT/$DATETIME.$FASTA.outerr 2>&1

if [ $? -ne 0 ];
then echo "read_fasta2_plasmids.py failed"
exit 1
fi

mkdir $OUT/dlDataPath
chmod 777 $OUT/dlDataPath

#pass fasta and yml directory to format
format_Test.py  --inputfasta $FASTA --inputyml $OUT/$DATETIME.$FASTA/yml --dataPath $OUT/dlDataPath

if [ $? -ne 0 ];
then echo "format_Test.py failed"
exit 1
fi

predict_Plasmid.py --dataPath $OUT/dlDataPath  --given test  --kModelList  20-29  --dataSegment -1

if [ $? -ne 0 ];
then echo "predict_Plasmid.py failed"
exit 1
fi

