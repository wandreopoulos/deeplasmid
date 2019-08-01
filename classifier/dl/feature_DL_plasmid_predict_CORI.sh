#!/bin/bash -l

echo "Running feature_DL_plasmid_predict.sh . This version is meant for Cori.
DelPlasmid - Plasmid finder
Author: Bill Andreopoulos
Last maintained: July 30, 2019"

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

MODEL=/global/dna/shared/data/functests/Assembly/models/plasmid4k-
MODEL=/global/cscratch1/sd/andreopo/plasmid4t-
MODEL=/global/cscratch1/sd/andreopo/plasmid4u-
MODEL=/global/cscratch1/sd/andreopo/plasmid4v-
MODEL=/global/dna/shared/data/functests/Assembly/models/plasmid4v-
MODEL=/global/cscratch1/sd/andreopo/plasmid4w-
MODEL=/global/cscratch1/sd/andreopo/plasmid4x-newfeat2-
###/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml2/classifier/models/plasmid4k-
MODEL=/global/cscratch1/sd/andreopo/plasmid4x-newfeat3-
MODEL=/global/cscratch1/sd/andreopo/plasmidml_tests/plasmid4x-newfeat4-
MODEL=/global/cscratch1/sd/andreopo/plasmidml_tests/plasmid4z-newfeat5-
MODEL=/global/dna/shared/data/functests/Assembly/models/plasmid4z-newfeat5-
echo "Model used for prediction: $MODEL" 
echo "Model used for prediction: $MODEL" > $OUT/outPR.$DATETIME/model_path.txt
echo "Command used: $0 $1 $2"
echo "Command used: $0 $1 $2" >> $OUT/outPR.$DATETIME/model_path.txt


python3 $PARENT/predict_Plasmid.py --dataPath $OUT/dlDataFormattedPred.$DATETIME   --outPath  $OUT/outPR.$DATETIME   --dataSegment -1   --given test  --kModelList  18-29  --seedModel  $MODEL 
# /global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml2/classifier/models/plasmid4k-
### /global/cscratch1/sd/andreopo/plasmid4p-
### /global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_rqcprod/classifier/models/plasmid4k-
### /global/cscratch1/sd/andreopo/plasmid4o-  
###$PARENT/models/plasmid4k-

if [ $? -ne 0 ];
then echo "predict_Plasmid.py failed"
exit 1
fi

#mv predictions.txt $OUT

