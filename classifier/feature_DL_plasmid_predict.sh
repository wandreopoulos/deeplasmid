#!/bin/bash -l

echo "Running feature_DL_plasmid_predict.sh
Plasmid finder
Author: Bill Andreopoulos
Last maintained: July 17, 2018"

FASTA=$1

OUT=$2

mkdir $OUT
chmod 777 $OUT

DATETIME=`date "+%Y-%m-%d_%H:%M:%S"`

python2.7 /srv/jgi-ml/classifier/read_fasta2_plasmids.py  -i $FASTA -o $OUT/$DATETIME

if [ $? -ne 0 ];
then echo "read_fasta2_plasmids.py failed"
exit 1
fi

mkdir $OUT/dlDataPath
chmod 777 $OUT/dlDataPath
mkdir out
mkdir outPR
mkdir outKF
mkdir plotKF

#pass fasta and yml directory to format
python3 /srv/jgi-ml/classifier/format_Test.py  --inputfasta $FASTA --inputyml $OUT/$DATETIME/yml --dataPath $OUT/dlDataPath

if [ $? -ne 0 ];
then echo "format_Test.py failed"
exit 1
fi

python3 /srv/jgi-ml/classifier/predict_Plasmid.py --dataPath $OUT/dlDataPath  --given test  --kModelList  18-29  --dataSegment -1  --seedModel  /srv/jgi-ml/classifier/models/plasmid4k-

if [ $? -ne 0 ];
then echo "predict_Plasmid.py failed"
exit 1
fi

