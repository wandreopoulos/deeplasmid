
README file
August 1, 2019
Maintainer: Bill Andreopoulos, wandreopoulos@lbl.gov



DelPlasmid is a tool based on machine learning that separates plasmids from chromosomal sequences. The input sequences are in the form of contigs and could have been produced from any sequencing technology or assembly algorithm. The deep learning model was trained on a corpus of:
1) plasmids from ACLAME, and 
2) chromosomal sequences from refseq.microbial from which plasmids and mito were removed.

______________________

To run delplasmid training
Training on Cori is done as follows:
. mySetup_nersc.source

Training command for the entire model:
andreopo@cori21:/global/cscratch1/sd/andreopo/plasmidml_tests> SRC=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_clean/classifier/dl
andreopo@cori21:/global/cscratch1/sd/andreopo/plasmidml_tests> $SRC/feature_DL_plasmid_train_CORI.sh $SRC/../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta ./aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.OUT19 $SRC/../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta ./refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.OUT19

Old way to run out of the src dir:
 andreopo@cori12:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml2/classifier/dl> ./feature_DL_plasmid_train_CORI.sh  ../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta   ../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.yml ../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta ../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.yml

______________________
To run delplasmid prediction example
Predict plasmids in the IMG dataset:
andreopo@cori21:/global/cscratch1/sd/andreopo/plasmidml_tests> $SRC/feature_DL_plasmid_predict_CORI.sh  $SRC/../DATA/IMG/genome_list.fasta.MIN1kMAX330k.fasta

Old way to run out of the src dir: 
andreopo@nid00401:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml2/classifier/dl> ./feature_DL_plasmid_predict_CORI.sh  ../DATA/IMG/genome_list.fasta.MIN1kMAX330k.fasta  ../DATA/IMG/genome_list.fasta.MIN1kMAX330k.fasta.OUT2

 ./feature_DL_plasmid_predict_CORI.sh  ../DATA/ACLAME.REFSEQMICROB.testseg6/assayer4.plasm_main-scaff-split.yml.fasta.fasta   ../DATA/ACLAME.REFSEQMICROB.testseg6/assayer4.plasm_main-scaff-split.yml.fasta.fasta.OUT4

______________________

Plasmid Finding GoogleDoc:
https://docs.google.com/document/d/1W7B-O5xKuXbWA_CwHjCl9K0-qKLBLpTnvuxkslcOH30/edit#

LucidChart with the software design:
https://www.lucidchart.com/documents/edit/97c70c55-cdf2-483e-a095-dbf8bed3537c/0

Code:
/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_clean/classifier/dl
/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml2/classifier/dl
Was cloned from:
git clone git@bitbucket.org:berkeleylab/jgi-ml.git  jgi-ml2

Excel sheet with data results:  https://docs.google.com/spreadsheets/d/1TDPn9uOAnZOBS95dJzUfVtb6mhteTuRFrFcv-_IZ4z4/edit#gid=1266368803


______________________

Testing methodology on smaller datasets:
Test the training script - for training use a subset of the data, though training will still take some time to complete (30 epochs):
    cd $CSCRATCH
    mkdir plasmidml_tests
    cd plasmidml_tests/
    SRC=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_clean/classifier/dl
    $SRC/feature_DL_plasmid_train_CORI.sh $SRC/../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.SUB.fasta ./aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.SUB.OUT $SRC/../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.SUB.fasta ./refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.SUB.OUT

Test the prediction script - use a small dataset with 10 contigs:
$SRC/feature_DL_plasmid_predict_CORI.sh $SRC/../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.SUB.fasta ./aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.SUB.OUTPRED
The predictions will be under the file ./aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.SUB.OUTPRED/outPR.*/predictions.txt
All predictions should be PLASMID since the test dataset comes from the ACLAME dataset.

The path to the model used and cmd line is under: outPR.*/model_path.txt

____________________

Docker:
There is the Dockerfile and Dockerfile.dev under the dl directory, 
which can be used to create Docker images with the tool.
The latest Docker image can be pulled from:


Logging:

License:

References:

Publication:
