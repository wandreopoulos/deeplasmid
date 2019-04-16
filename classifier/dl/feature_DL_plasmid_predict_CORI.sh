#!/bin/bash -l

echo "Running feature_DL_plasmid_predict.sh . This version is meant for Cori.
DelPlasmid - Plasmid finder
Author: Bill Andreopoulos
Last maintained: April 1, 2019"

#module load deeplearning

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

mkdir $OUT/dlDataFormattedPred.$DATETIME
chmod 777 $OUT/dlDataFormattedPred.$DATETIME
#outPR.$DATETIME stored a pdf file with prediction data and yml files for each sequence, the plots with histograms and the predictions.txt file
mkdir $OUT/outPR.$DATETIME
chmod 777 $OUT/outPR.$DATETIME
#These dirs are meant for storing plots, but unused atm.
#mkdir out
#mkdir outKF
#mkdir plotKF

module unload python
module load python/3.6-anaconda-5.2

#NERSC requires this env var!
export HDF5_USE_FILE_LOCKING=FALSE

#pass fasta and yml directory to format
python3 $PARENT/format_predict.py  --inputfasta $FASTA --inputyml $OUT/dlFeatures.$DATETIME/yml --dataPath $OUT/dlDataFormattedPred.$DATETIME  --outPath  $OUT/outPR.$DATETIME  -X

if [ $? -ne 0 ];
then echo "format_predict.py failed"
exit 1
fi

python3 $PARENT/predict_Plasmid.py --dataPath $OUT/dlDataFormattedPred.$DATETIME   --outPath  $OUT/outPR.$DATETIME  --given test  --kModelList  18-29  --seedModel  /global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_rqcprod/classifier/models/plasmid4k-
### /global/cscratch1/sd/andreopo/plasmid4o-  
###$PARENT/models/plasmid4k-

if [ $? -ne 0 ];
then echo "predict_Plasmid.py failed"
exit 1
fi

#mv predictions.txt $OUT

