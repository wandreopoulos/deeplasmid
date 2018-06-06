#!/bin/bash -l
#get_plasmids_from_fungi.sh <in fasta file>

usage(){
echo "
Written by Bill Andreopoulos, from February 2015 - present
Last modified June 2, 2015

Description:  This is a tool for finding plasmids in fungi.
The run is based on a J48 decision tree.
The decision tree has been trained to separate fungi vs. plasmid
sequences on the basis of a set of predetermined features, including:
- GC % in the entire sequence.
- Min and Max GC% in any window of 100b.
- The longest homopolymer for each of A,C,G,T.
- The total nucleotides in long (>5n) homopolymers.
- Most frequent di-, tri-, tetranucleotide (python khmer package).
- Length of the sequence.
- Longest alignment to refseq.plasmid.

Note: the model has been trained only for plasmid vs. fungi separation.

To run:   get_plasmids_from_fungi.sh <in fasta file>

The parameter is a fasta file of contigs. The run will produce an
output directory with a .fa file of the contigs that are predicted
to be plasmids.

Please contact Bill Andreopoulos at wandreopoulos@lbl.gov if you encounter
any problems
"
}

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
        usage
        exit
fi

#MODELS:
#***Fungi model for refseq.fungi no mito/no plasmids is under:
#/global/scratch2/sd/andreopo/GAA-1290_plasmids/Kurt/refseq.fungi
#
#Microbial model for refseq.microbial no mito/no plasmids is under:
#/global/scratch2/sd/andreopo/GAA-1290_plasmids/Dingl/microbial/NEW2_ALL_FEATURES/
#
#Fungi model for jgi released non-mitos in fungal projects is under:
#/global/scratch2/sd/andreopo/GAA-1330_fungal/main_fungal_jgi_releases_ALLFEATURES2
#
#Mito model for jgi released mitos in fungal projects is under:
#/global/scratch2/sd/andreopo/GAA-1330_fungal/mito_fungal_jgi_releases_ALLFEATURES2
#
#***Plasmid model for refseq.plasmid is under:
# /global/scratch2/sd/andreopo/GAA-1290_plasmids/Dingl/plasmids/NEW4_ALL_FEATURES/

# read contigs.fa as command line option

# run read_fasta.py against contigs.fa
export CLASSPATH=/global/scratch2/sd/andreopo/GAA-1330_fungal/weka-3-6-12:$CLASSPATH

export PYTHONPATH=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/assemblyqc/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/readqc/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/tools/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc/:$PYTHONPATH

FASTA=`realpath $1`
if [ -f $FASTA ];
then
   echo "File $FASTA exists and will be used as the input fasta of contigs."
else
   echo "File $FASTA does not exist."
   exit
fi
###`cat ./run.config | sed -n ${SGE_TASK_ID}p`

FASTA_FILE=`echo $(basename $FASTA) | sed 's/^.*organelleAssembly\/test\/NEW\///' | sed 's/\/firstVelvet\/contigs.fa$//' | sed 's/\//_/g' | sed 's/$/_NEW/'`

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

module load khmer

###/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/
read_fasta.py -b p  -i  $FASTA  -o  $DIR/$FASTA_FILE  >  $DIR/$FASTA_FILE/organelle.out  2>  $DIR/$FASTA_FILE/organelle.err

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

for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff` ; do cat /global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/header_plasmids_fungi $i > $i.tmp && mv $i.tmp $i ; done

###To produce the model use:
###andreopo@gpint108:/global/scratch2/sd/andreopo/GAA-1330_fungal/weka-3-6-12$ java weka.classifiers.trees.J48 -c 1 -d  ./MODEL/mainfungal_plasmid_jgi_releases_ALLFEATURES_BALANCED.J48.model -t ./MODEL/mainfungal_plasmid_jgi_releases_ALLFEATURES_BALANCED.arff
###For the full sequence of commands to make the model see: /global/scratch2/sd/andreopo/GAA-1330_fungal/weka-3-6-12/MODEL/cmds2.log
ML=weka.classifiers.trees.J48

MODEL=/global/scratch2/sd/andreopo/GAA-1330_fungal/weka-3-6-12/MODEL/mainfungal_plasmid_jgi_releases_ALLFEATURES_BALANCED.J48.model

###Do the prediction
#run weka J48 against saved model
#andreopo@gpint108:/global/scratch2/sd/andreopo/GAA-1330_fungal/weka-3-6-12$ for i in `find ../organelleAssemblyResults/*_NEW -name features.txt-NEW.arff` ; do  java weka.classifiers.trees.J48 -c 1 -p 0 -l ./MODEL/mainfungal_plasmid_jgi_releases_ALLFEATURES_BALANCED.J48.model -T $i > $i.J48.pred ; done
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff` ; do  java $ML -c 1 -p 0 -l $MODEL -T $i > $i.J48.pred ; done

#parse weka output for 2:M lines
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff.J48.pred` ; do egrep "2:P" $i | sed "s/ 3:T.*$//g" > $i.PLASMIDLINES ; done

#extract from contigs.fa the mito contigs and save into mito_contigs.fa
for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff.J48.pred.PLASMIDLINES` ; do iDIR=$(dirname $i); for j in `cat $i` ; do  sed -n $(($j+1))p  $iDIR/features.txt  | sed 's/^.* >//g' | sed 's/COV/cov/g' | sed 's/LENGTH/length/g' >> $i.PLASMIDCONTIGS ; done ; done

for i in `find $DIR/$FASTA_FILE -name features.txt-NEW.arff.J48.pred.PLASMIDLINES.PLASMIDCONTIGS` ; do iDIR=$(dirname $i); ~jfroula/Tools/Jazz/screen_list.pl $i $iDIR/contigs.fa keep > $i.fa ; done

NUMCONTIGS=`wc -l $DIR/$FASTA_FILE/features.txt-NEW.arff.J48.pred.PLASMIDLINES | sed 's/ .*$//'`

echo "Prediction run complete"
echo "Number of plasmid contigs: $NUMCONTIGS"
echo "Prediction results: features.txt-NEW.arff.J48.pred"
echo "Plasmids fasta file: features.txt-NEW.arff.J48.pred.PLASMIDLINES.PLASMIDCONTIGS.fa"

#optionally: run quast against released mito fasta
