#!/bin/bash -l
set -e

echo "Deeplasmid - Plasmid finder for microbial genome assemblies. 
Running feature_DL_plasmid_predict_CORI.sh .
This .sh script is meant for running Deeplasmid native on Cori.
Usage: please specify 2 arguments - the input fasta file and output directory - as follows:
     feature_DL_plasmid_predict_CORI.sh  in.fasta outdir
Contact person: Bill Andreopoulos, wandreopoulos@lbl.gov
Last maintained: September 19, 2021"

#module load deeplearning

if [ $# -ne 2 ]; then
    echo "Usage of Deeplasmid on the command line. Please specify 2 arguments (input fasta file and output directory):   feature_DL_plasmid_predict_CORI.sh  in.fasta outdir"
    exit 1
fi

trap 'pkill -P $$' EXIT

FASTA=$1

OUT=$2

mkdir -p $OUT
if [ $? -ne 0 ];
then echo "mkdir $OUT failed"
exit 1
fi

chmod 777 $OUT
if [ $? -ne 0 ];
then echo "chmod 777 $OUT failed"
exit 1
fi

DATETIME=`date "+%Y%m%d_%H%M%S"`
echo "Using $DATETIME for outdir suffix"

module unload python
module load python/2.7-anaconda-4.4

PARENT=`dirname $0`

python2.7 $PARENT/read_fasta2_plasmids.py  -i $FASTA -o $OUT/dlFeatures.$DATETIME

if [ $? -ne 0 ];
then echo "read_fasta2_plasmids.py failed"
exit 1
fi

mkdir -p $OUT/dlDataFormattedPred.$DATETIME
chmod 777 $OUT/dlDataFormattedPred.$DATETIME
#outPR.$DATETIME stored a pdf file with prediction data and yml files for each sequence, the plots with histograms and the predictions.txt file
mkdir -p $OUT/outPR.$DATETIME
chmod 777 $OUT/outPR.$DATETIME

module unload python
module load python/3.6-anaconda-5.2

#NERSC requires this env var!
export HDF5_USE_FILE_LOCKING=FALSE

#pass fasta and yml directory to format
python3 $PARENT/format_predict.py  --inputfasta $FASTA --inputyml $OUT/dlFeatures.$DATETIME/yml --dataPath $OUT/dlDataFormattedPred.$DATETIME  --outPath  $OUT/outPR.$DATETIME  --no-Xterm

if [ $? -ne 0 ];
then echo "format_predict.py failed"
exit 1
fi

#The model used for prediction
MODEL=/global/dna/shared/data/functests/Assembly/Plasmid_Models/plasmid4z-newfeat5-
MODEL=/global/cscratch1/sd/andreopo/plasmidml_tests/jgi-ml_paper/classifier/dl/plasmid4z-newfeat12-
echo "Model used for prediction: $MODEL" 
echo "Model used for prediction: $MODEL" > $OUT/outPR.$DATETIME/model_path.txt
echo "Command used: $0 $1 $2"
echo "Command used: $0 $1 $2" >> $OUT/outPR.$DATETIME/model_path.txt


python3 $PARENT/predict_Plasmid.py --dataPath $OUT/dlDataFormattedPred.$DATETIME   --outPath  $OUT/outPR.$DATETIME   --dataSegment -1   --given test  --kModelList  18-29  --seedModel  $MODEL  

if [ $? -ne 0 ];
then echo "predict_Plasmid.py failed"
exit 1
fi

#mv predictions.txt $OUT

