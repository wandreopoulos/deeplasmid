#!/bin/bash -l

echo "Running feature_DL_plasmid_train.sh . This version is meant for Cori.
DelPlasmid - Plasmid finder
Author: Bill Andreopoulos
Last maintained: May 2, 2019"





#module load deeplearning
####################################
#First the plasmids (class 1)
####################################

FASTA=$1

OUT=$2

mkdir $OUT
chmod 777 $OUT

DATETIME=`date "+%Y-%m-%d_%H:%M:%S"`

module unload python
module load python/2.7-anaconda-4.4

PARENT=`dirname $0`


python2.7 $PARENT/read_fasta2_plasmids.py  -i $FASTA -o $OUT/dlFeatures.$DATETIME

if [ $? -ne 0 ];
then echo "read_fasta2_plasmids.py failed"
exit 1
fi





####################################
#Second the chromosomes (class 2)
####################################

FASTA2=$3

OUT2=$4

mkdir $OUT2
chmod 777 $OUT2

module unload python
module load python/2.7-anaconda-4.4

PARENT=`dirname $0`


python2.7 $PARENT/read_fasta2_plasmids.py  -i $FASTA2 -o $OUT2/dlFeatures.$DATETIME

if [ $? -ne 0 ];
then echo "read_fasta2_plasmids.py failed"
exit 1
fi





#TODO this should not be a hard-coded path
mv delplasmid_data4train delplasmid_data4train.$DATETIME
mkdir delplasmid_data4train
chmod 777 delplasmid_data4train

mkdir out
mkdir outPR
mkdir outKF
mkdir plotKF

module unload python
module load python/3.6-anaconda-5.2

export HDF5_USE_FILE_LOCKING=FALSE

#pass fasta and yml directory to format
#python3 $PARENT/format_Test.py  --inputfasta $FASTA --inputyml $OUT/$DATETIME/yml --dataPath $OUT/dlDataPath
python3 $PARENT/format_train.py --dataPath  delplasmid_data4train  -X   --inputfasta $FASTA --inputyml $OUT/dlFeatures.$DATETIME/yml   --inputfasta2 $FASTA2 --inputyml2 $OUT2/dlFeatures.$DATETIME/yml


if [ $? -ne 0 ];
then echo "format_train.py failed"
exit 1
fi




sbatch batchTrain.slr

