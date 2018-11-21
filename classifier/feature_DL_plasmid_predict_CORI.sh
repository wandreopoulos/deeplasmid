#!/bin/bash -l

echo "Running feature_DL_plasmid_predict.sh
Plasmid finder
Author: Bill Andreopoulos
Last maintained: November 16, 2018"

#module load deeplearning

FASTA=$1

OUT=$2

mkdir $OUT
chmod 777 $OUT

DATETIME=`date "+%Y-%m-%d_%H:%M:%S"`

module unload python
module load python/2.7-anaconda-4.4

python2.7 read_fasta2_plasmids.py  -i $FASTA -o $OUT/$DATETIME

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

module unload python
module load python/3.6-anaconda-4.4

#pass fasta and yml directory to format
python3 format_Test.py  --inputfasta $FASTA --inputyml $OUT/$DATETIME/yml --dataPath $OUT/dlDataPath

if [ $? -ne 0 ];
then echo "format_Test.py failed"
exit 1
fi

python3 predict_Plasmid.py --dataPath $OUT/dlDataPath  --given test  --kModelList  18-29  --dataSegment -1  --seedModel  models/plasmid4k-

if [ $? -ne 0 ];
then echo "predict_Plasmid.py failed"
exit 1
fi

