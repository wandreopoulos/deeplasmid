
Ran delplasmid training
Training on Cori is done as follows:
. mySetup_nersc.source
andreopo@cori12:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml2/classifier/dl> ./feature_DL_plasmid_train_CORI.sh  ../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta   ../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.yml ../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta ../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.yml


Ran delplasmid prediction example:
andreopo@nid00401:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml2/classifier/dl> ./feature_DL_plasmid_predict_CORI.sh  ../DATA/IMG/genome_list.fasta.MIN1kMAX330k.fasta  ../DATA/IMG/genome_list.fasta.MIN1kMAX330k.fasta.OUT2

 ./feature_DL_plasmid_predict_CORI.sh  ../DATA/ACLAME.REFSEQMICROB.testseg6/assayer4.plasm_main-scaff-split.yml.fasta.fasta   ../DATA/ACLAME.REFSEQMICROB.testseg6/assayer4.plasm_main-scaff-split.yml.fasta.fasta.OUT4

