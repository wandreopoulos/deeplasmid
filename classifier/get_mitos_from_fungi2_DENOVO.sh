#!/bin/bash -l
#get_mitos_from_fungi2.sh <in fasta file>

usage(){
echo "
Written by Bill Andreopoulos, from February 2015 - present
Last modified January 28, 2016

Description:  This is a tool for finding mitos in fungi.
The run is based on a Naive Bayes classifier.
The classifier has been trained to separate fungi vs. mitochondrial
sequences on the basis of a set of predetermined features, including:
- GC % in the entire sequence.
- Min and Max GC% in any window of 100b.
- The longest homopolymer for each of A,C,G,T.
- The total nucleotides in long (>5n) homopolymers.
- Most frequent di-, tri-, tetranucleotide up to dekamers (python khmer package).

These features became obsolete (June 2015):
- Length of the sequence.
- Longest alignment to refseq.mito.
- The repeat and inverse repeat content.

Note: the model has been trained only for mitochondrial vs. fungi separation.

To run:   get_mitos_from_fungi2.sh <in fasta file>

The parameter is a fasta file of contigs. The run will produce an
output directory with a .fa file of the contigs that are predicted
to be mitochondrial.

Please contact Bill Andreopoulos at wandreopoulos@lbl.gov if you encounter
any problems
"
}

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
        usage
        exit
fi

#MODELS:
#Fungi model for refseq.fungi no mito/no plasmids is under:
#/global/scratch2/sd/andreopo/GAA-1290_plasmids/Kurt/refseq.fungi
#
#Microbial model for refseq.microbial no mito/no plasmids is under:
#/global/scratch2/sd/andreopo/GAA-1290_plasmids/Dingl/microbial/NEW_ALL_FEATURES/
#
#***Fungi model for jgi released non-mitos in fungal projects is under:
#/global/scratch2/sd/andreopo/GAA-1330_fungal/main_fungal_jgi_releases_ALLFEATURES3
#
#***Mito model for jgi released mitos in fungal projects is under:
#/global/scratch2/sd/andreopo/GAA-1330_fungal/mito_fungal_jgi_releases_ALLFEATURES3
#
#Plasmid model for refseq.plasmid is under:
# /global/scratch2/sd/andreopo/GAA-1290_plasmids/Dingl/plasmids/NEW2_ALL_FEATURES/

# read contigs.fa as command line option

# run read_fasta.py against contigs.fa
export CLASSPATH=/global/projectb/sandbox/rqc/andreopo/GAA-1330_fungal/weka-3-6-12:$CLASSPATH

export PYTHONPATH=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/assemblyqc/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/readqc/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/tools/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc/:$PYTHONPATH

export PATH=$PATH:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml/classifier/

export PATH=~qc_user/miniconda/miniconda2/bin:$PATH

FASTA=`realpath $1`
if [ -f $FASTA ];
then
   echo "File $FASTA exists and will be used as the input fasta of contigs."
else
   echo "File $FASTA does not exist."
   exit
fi
###`cat ./run.config | sed -n ${SGE_TASK_ID}p`

FASTA_FILE=`echo $FASTA | sed 's/^.*scratch//' | sed 's/\//_/g' | sed 's/$/_NEW/'`
###echo $(basename $FASTA) | sed 's/^.*organelleAssembly\/test\/NEW\///' | sed 's/\/firstVelvet\/contigs.fa$//' | 
###FASTA_FILE=$(basename $FASTA)

DIR=`pwd`

if [ -d "$DIR/$FASTA_FILE" ]; then
  echo "The DIRECTORY exists $DIR/$FASTA_FILE"
  NOW=$(date +%Y-%m-%d:%H:%M:%S)
  mv $DIR/$FASTA_FILE $DIR/$FASTA_FILE.$NOW
  echo "The DIRECTORY was moved to $DIR/$FASTA_FILE.$NOW"
fi

mkdir $DIR/$FASTA_FILE/

chmod -R 775 $DIR/$FASTA_FILE/

ln -s $FASTA $DIR/$FASTA_FILE/contigs.fa

echo "Output directory: $DIR/$FASTA_FILE/"
echo "Started computing prediction features....."

###module load khmer

###/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/
DIR2="$( dirname ${0})"

$DIR2/read_fasta2_DENOVO_test.py  -i  $FASTA  -o  $DIR/$FASTA_FILE  >  $DIR/$FASTA_FILE/organelle.out  2>  $DIR/$FASTA_FILE/organelle.err
###Filter: If the result did not align at all to refseq.mito then consider there to be no mito and exit, else continue.

echo "Finished computing prediction features"

#cp the features to a new file
for i in `find $DIR/$FASTA_FILE -name features.txt` ; do cp $i $i-NEW.arff ; done


#replace column 1 with T
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff` ; do sed -i 's/^-1 /T /g' $i ; done

# replace end columns
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff` ; do sed -i 's/ DATETIME.*$//g' $i ; done

#replace space with comma
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff` ; do sed -i 's/ /,/g' $i ; done

#replace header with text blob
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff` ; do sed -i '/^id.*$/d' $i ; done

for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff` ; do sed -i 's/,,/,?,/g' $i ; done

for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff` ; do cat /global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml/classifier/header_mito_fungi2 $i > $i.tmp && mv $i.tmp $i ; done
# /global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/header_mito_fungi2


###To produce the model use:
###andreopo@gpint108:/global/scratch2/sd/andreopo/GAA-1330_fungal/weka-3-6-12$ java weka.classifiers.bayes.NaiveBayes -c 1 -d  ./MODEL/mainfungal_mito_jgi_releases_ALLFEATURES_BALANCED.NB.model -t ./MODEL/mainfungal_mito_jgi_releases_ALLFEATURES_BALANCED.arff
###For the full sequence of commands to make the model see: /global/scratch2/sd/andreopo/GAA-1330_fungal/weka-3-6-12/MODEL/cmds.log
ML=weka.classifiers.bayes.NaiveBayes

MODEL=/global/projectb/sandbox/rqc/andreopo/GAA-1330_fungal/weka-3-6-12/MODEL/mainfungal_mito_jgi_releases_ALLFEATURESmanymers_BALANCED.NB.model

###Do the prediction
#run weka NB against saved model
#andreopo@gpint108:/global/scratch2/sd/andreopo/GAA-1330_fungal/weka-3-6-12$ for i in `find ../organelleAssemblyResults/*_NEW -name features.txt-NEW.arff` ; do  java weka.classifiers.bayes.NaiveBayes -c 1 -p 0 -l ./MODEL/mainfungal_mito_jgi_releases_ALLFEATURES_BALANCED.NB.model -T $i > $i.NB.pred ; done
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff` ; do  java $ML -c 1 -p 0 -l $MODEL -T $i > $i.NB.pred ; done

#parse weka output for 2:M lines
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff.NB.pred` ; do egrep "2:M" $i | sed "s/ 3:T.*$//g" > $i.MITOLINES ; done

#extract from contigs.fa the mito contigs and save into mito_contigs.fa
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff.NB.pred.MITOLINES` ; do iDIR=$(dirname $i); for j in `cat $i` ; do  sed -n $(($j+1))p  $iDIR/features.txt  | sed 's/^.* >//g' | sed 's/COV/cov/g' | sed 's/LENGTH/length/g' >> $i.MITOCONTIGS ; done ; done

###module load jigsaw

###for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff.NB.pred.MITOLINES.MITOCONTIGS` ; do iDIR=$(dirname $i); ~jfroula/Tools/Jazz/screen_list.pl $i $iDIR/contigs.fa keep > $i.fa ; done
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff.NB.pred.MITOLINES.MITOCONTIGS` ; do iDIR=$(dirname $i); shifter --image=bryce911/bbtools  filterbyname.sh names=$i in=$iDIR/contigs.fa  out=$i.fa  include=t ; done

NUMCONTIGS=`wc -l $DIR/$FASTA_FILE/features.txt-NEW.arff.NB.pred.MITOLINES | sed 's/ .*$//'`

echo "Prediction run complete"
echo "Number of mitochondrial contigs: $NUMCONTIGS"
echo "Prediction results: features.txt-NEW.arff.NB.pred"
echo "Mitochondrial fasta file: features.txt-NEW.arff.NB.pred.MITOLINES.MITOCONTIGS.fa"

#optionally: run quast against released mito fasta

if ! [ -z "$2" ]; then
    cp $DIR/$FASTA_FILE/features.txt-NEW.arff.NB.pred.MITOLINES.MITOCONTIGS.fa  $DIR/$2
    echo "Mitochondrial fasta file was copied to $DIR/$2 "
else
    echo "Mitochondrial fasta file is features.txt-NEW.arff.NB.pred.MITOLINES.MITOCONTIGS.fa"
fi

