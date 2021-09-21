
*********

README file
September 20, 2021
Maintainer: Bill Andreopoulos, wandreopoulos@lbl.gov
Codebase: master branch


Deeplasmid is a tool based on machine learning that separates plasmids from chromosomal sequences. The input sequences are in the form of contigs and could have been produced from any sequencing technology or assembly algorithm. The deep learning model was trained on a corpus of:
1) plasmids from ACLAME, and 
2) chromosomal sequences from refseq.microbial from which plasmids and mito were removed.

______________________

Deeplasmid Docker container for identifying plasmids in microbial isolated and metagenome assemblies
Input: a .fasta file
Tested on: Mac, Ubuntu

How to run the Docker container with deeplasmid:
Install Docker on your system. Register on dockerhub.
Pull the deeplasmid image from dockerhub as follows:

docker login
docker pull billandreo/deeplasmid

Please run deeplasmid prediction as follows. Substitute the /path/to/input/fasta and /path/to/output/directory below with your inut file and output dir full paths:
docker run -it -v /path/to/input/fasta:/srv/jgi-ml/classifier/dl/in.fasta -v /path/to/output/directory:/srv/jgi-ml/classifier/dl/outdir billandreo/deeplasmid feature_DL_plasmid_predict.sh in.fasta outdir

The file predictions.txt file is the file of plasmid prediction results for the contigs that we are interested in.
Each contig name will indicate if it was:
- Plasmid with score near 1.0
- Chromosome (non-plasmid) with score near 0.0
- Ambiguous have scores around 0.5 (gray zone)
- Longer than 330k bases (possibly a chromosome or megaplasmid)
- Shorter than 1k bases (inconclusive)

These 2 score files are also output from the tool; they are histograms that show the scores for samples and scaffolds that had or didn't have "plasmid" in the header. In reality headers won't have "plasmid" in the header, obviously, but these figures are useful when testing the tool on a dataset where the classes are known.
- samplescore_hist.png
- scaffscore_hist.png
These files are not output by default. The final plotting step is skipped by setting the -noXterm(-X) flag in the .sh script to true by default.

The public Docker repository is available here:
https://hub.docker.com/repository/docker/billandreo/deeplasmid

The present sourceforge repo has the branch "docker" codebase used for building the Docker image, which has been built for prediction purposes on any platform where Docker is installed.
Building the Docker container was done as follows:
docker build -t billandreo/deeplasmid -f Dockerfile.v2 .

Please see the Supplementary Information from the publication for 3 issues when building the Docker file: Prodigal and bbtools/sketch need to be built, and the model .h5 files from training are needed.

Testing

The 649989979.fna is a testing file, and is available for download from IMG:
/global/cscratch1/sd/andreopo/plasmidml_tests/jgi-ml_paper/classifier/dl/testing/649989979/649989979.fna
This way you can test to verify if your installation of the deeplasmid tool gives the same results as expected, which are shown below.

You can run on this file with:

andreopo@nid00245:/global/cscratch1/sd/andreopo/plasmidml_tests/jgi-ml_paper/classifier/dl/testing> ../feature_DL_plasmid_predict_CORI.sh 649989979/649989979.fna 649989979d.native.OUT
andreopo@nid00245:/global/cscratch1/sd/andreopo/plasmidml_tests/jgi-ml_paper/classifier/dl/testing> cat 649989979d.native.OUT/outPR.20200113_172443/predictions.txt
name,pred,conf
NZ_ADHJ01000001,LONGER_330000.0,1
NZ_ADHJ01000014,LONGER_330000.0,1
NZ_ADHJ01000017,LONGER_330000.0,1
NZ_ADHJ01000025,LONGER_330000.0,1
NZ_ADHJ01000027,SHORTER_1000.0,1
NZ_ADHJ01000030,SHORTER_1000.0,1
NZ_ADHJ01000032,SHORTER_1000.0,1
NZ_ADHJ01000037,LONGER_330000.0,1
NZ_ADHJ01000051,SHORTER_1000.0,1
nz_adhj01000009 paenibacillus vortex v453 cnt_pvor1000009, whole genome shotgun sequence.,GENOME,0.013 +/- 0.001
nz_adhj01000010 paenibacillus vortex v453 cnt_pvor1000010, whole genome shotgun sequence.,PLASMID,0.830 +/- 0.005
nz_adhj01000028 paenibacillus vortex v453 cnt_pvor1000028, whole genome shotgun sequence.,GENOME,0.018 +/- 0.002
nz_adhj01000033 paenibacillus vortex v453 cnt_pvor1000033, whole genome shotgun sequence.,GENOME,0.015 +/- 0.002
nz_adhj01000040 paenibacillus vortex v453 cnt_pvor1000040, whole genome shotgun sequence.,GENOME,0.007 +/- 0.001
nz_adhj01000043 paenibacillus vortex v453 cnt_pvor1000043, whole genome shotgun sequence.,GENOME,0.005 +/- 0.000
nz_adhj01000044 paenibacillus vortex v453 cnt_pvor1000044, whole genome shotgun sequence.,GENOME,0.057 +/- 0.003
nz_adhj01000004 paenibacillus vortex v453 cnt_pvor1000004, whole genome shotgun sequence.,PLASMID,0.798 +/- 0.007
nz_adhj01000039 paenibacillus vortex v453 cnt_pvor1000039, whole genome shotgun sequence.,PLASMID,0.739 +/- 0.006
nz_adhj01000042 paenibacillus vortex v453 cnt_pvor1000042, whole genome shotgun sequence.,GENOME,0.002 +/- 0.000
nz_adhj01000046 paenibacillus vortex v453 cnt_pvor1000046, whole genome shotgun sequence.,GENOME,0.315 +/- 0.011
nz_adhj01000048 paenibacillus vortex v453 cnt_pvor1000048, whole genome shotgun sequence.,PLASMID,0.847 +/- 0.004
nz_adhj01000052 paenibacillus vortex v453 cnt_pvor1000052, whole genome shotgun sequence.,GENOME,0.055 +/- 0.004
nz_adhj01000054 paenibacillus vortex v453 cnt_pvor1000054, whole genome shotgun sequence.,PLASMID,0.816 +/- 0.006
nz_adhj01000055 paenibacillus vortex v453 cnt_pvor1000055, whole genome shotgun sequence.,GENOME,0.425 +/- 0.006
nz_adhj01000056 paenibacillus vortex v453 cnt_pvor1000056, whole genome shotgun sequence.,PLASMID,0.571 +/- 0.009
nz_adhj01000003 paenibacillus vortex v453 cnt_pvor1000003, whole genome shotgun sequence.,GENOME,0.210 +/- 0.005
nz_adhj01000006 paenibacillus vortex v453 cnt_pvor1000006, whole genome shotgun sequence.,GENOME,0.018 +/- 0.001
nz_adhj01000007 paenibacillus vortex v453 cnt_pvor1000007, whole genome shotgun sequence.,GENOME,0.004 +/- 0.000
nz_adhj01000012 paenibacillus vortex v453 cnt_pvor1000012, whole genome shotgun sequence.,PLASMID,0.799 +/- 0.005
nz_adhj01000013 paenibacillus vortex v453 cnt_pvor1000013, whole genome shotgun sequence.,GENOME,0.010 +/- 0.000
nz_adhj01000015 paenibacillus vortex v453 cnt_pvor1000015, whole genome shotgun sequence.,GENOME,0.106 +/- 0.013
nz_adhj01000019 paenibacillus vortex v453 cnt_pvor1000019, whole genome shotgun sequence.,GENOME,0.034 +/- 0.004
nz_adhj01000021 paenibacillus vortex v453 cnt_pvor1000021, whole genome shotgun sequence.,PLASMID,0.686 +/- 0.008
nz_adhj01000022 paenibacillus vortex v453 cnt_pvor1000022, whole genome shotgun sequence.,GENOME,0.025 +/- 0.003
nz_adhj01000024 paenibacillus vortex v453 cnt_pvor1000024, whole genome shotgun sequence.,PLASMID,0.731 +/- 0.008
nz_adhj01000036 paenibacillus vortex v453 cnt_pvor1000036, whole genome shotgun sequence.,PLASMID,0.893 +/- 0.003
nz_adhj01000049 paenibacillus vortex v453 cnt_pvor1000049, whole genome shotgun sequence.,PLASMID,0.890 +/- 0.003
nz_adhj01000053 paenibacillus vortex v453 cnt_pvor1000053, whole genome shotgun sequence.,PLASMID,0.703 +/- 0.004
nz_adhj01000002 paenibacillus vortex v453 cnt_pvor1000002, whole genome shotgun sequence.,GENOME,0.186 +/- 0.004
nz_adhj01000008 paenibacillus vortex v453 cnt_pvor1000008, whole genome shotgun sequence.,PLASMID,0.878 +/- 0.004
nz_adhj01000016 paenibacillus vortex v453 cnt_pvor1000016, whole genome shotgun sequence.,GENOME,0.008 +/- 0.001
nz_adhj01000018 paenibacillus vortex v453 cnt_pvor1000018, whole genome shotgun sequence.,GENOME,0.454 +/- 0.006
nz_adhj01000045 paenibacillus vortex v453 cnt_pvor1000045, whole genome shotgun sequence.,PLASMID,0.759 +/- 0.007
nz_adhj01000011 paenibacillus vortex v453 cnt_pvor1000011, whole genome shotgun sequence.,GENOME,0.146 +/- 0.010
nz_adhj01000020 paenibacillus vortex v453 cnt_pvor1000020, whole genome shotgun sequence.,GENOME,0.002 +/- 0.000
nz_adhj01000026 paenibacillus vortex v453 cnt_pvor1000026, whole genome shotgun sequence.,GENOME,0.086 +/- 0.006
nz_adhj01000029 paenibacillus vortex v453 cnt_pvor1000029, whole genome shotgun sequence.,GENOME,0.000 +/- 0.000
nz_adhj01000050 paenibacillus vortex v453 cnt_pvor1000050, whole genome shotgun sequence.,PLASMID,0.889 +/- 0.003
nz_adhj01000005 paenibacillus vortex v453 cnt_pvor1000005, whole genome shotgun sequence.,GENOME,0.006 +/- 0.001
nz_adhj01000023 paenibacillus vortex v453 cnt_pvor1000023, whole genome shotgun sequence.,GENOME,0.450 +/- 0.006
nz_adhj01000031 paenibacillus vortex v453 cnt_pvor1000031, whole genome shotgun sequence.,GENOME,0.044 +/- 0.002
nz_adhj01000034 paenibacillus vortex v453 cnt_pvor1000034, whole genome shotgun sequence.,GENOME,0.131 +/- 0.006
nz_adhj01000035 paenibacillus vortex v453 cnt_pvor1000035, whole genome shotgun sequence.,GENOME,0.400 +/- 0.013
nz_adhj01000038 paenibacillus vortex v453 cnt_pvor1000038, whole genome shotgun sequence.,GENOME,0.041 +/- 0.003
nz_adhj01000041 paenibacillus vortex v453 cnt_pvor1000041, whole genome shotgun sequence.,PLASMID,0.923 +/- 0.001
nz_adhj01000047 paenibacillus vortex v453 cnt_pvor1000047, whole genome shotgun sequence.,PLASMID,0.691 +/- 0.010
andreopo@nid00245:/global/cscratch1/sd/andreopo/plasmidml_tests/jgi-ml_paper/classifier/dl/testing>

Training

The training of the deeplasmid deep learning model was done on Cori at NERSC. These are the training steps:

Get an salloc session with 48 hours allocation:
salloc: sal with 48 hours

Then run as follows:
feature_DL_plasmid_train_CORI.sh plasmid.fasta plasmid.fasta.OUTDIR nonplasmid.fasta nonplasmid.fasta.OUTDIR

For example:
andreopo@cori02:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_paper/classifier/dl> ./feature_DL_plasmid_train_CORI.sh  ../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta    ../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.OUT18/    ../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.UNION.lists_shortq.archaea.txt.fasta.MIN1kMAX330k.fasta       ../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.UNION.lists_shortq.archaea.txt.fasta.MIN1kMAX330k.fasta.OUT18/

Running feature_DL_plasmid_train.sh . This version is meant for Cori.

.....

After training completes you will have a trained model under the two directories you specified on the command line. You may use the trained model in the predictions file feature_DL_plasmid_predict_CORI.sh

The codebase used for training on Cori (as well as running prediction on Cori) is under https://bitbucket.org/berkeleylab/jgi-ml/src/master/

A single training data element consists of the label and two input words: xseq - a 300bp contiguous subsequence sampled randomly from the full original contig sequence and xf - a vector containing 16 features extracted from the full sequence, as described in the Table below. The number (m) of 300bp subsequences sampled from each contig is proportional to the square root of the contig length, such that longer contigs contribute more samples, but do not overwhelm the training.

Header	Header	Header
Name	Definition	Type
gc_content	GC content of contig	Float [0-1]
A(C/G/T)_longest_homopolymer	Length of longest homopolymer	Integer
A(C/G/T)_total_homopolymer	Total number of homopolymers of length >5	Integer
hit_chromosome_proteins	Hit to chromosome proteins	Boolean 0/1
hit_plasmid_proteins	Hit to plasmid proteins	Boolean 0/1
hit_plasmid_ORIs	Hit to plasmid ORI	Boolean 0/1
gene_count	Number of genes in scaffold (Prodigal)	Integer
gene_percent	Coding percent of scaffold (Prodigal)	Float [0-1]
polypeptide_aa_avg_len	Average length of aa sequence (Prodigal)	Integer
len_sequence	Scaffold seq length	Integer
Table 1. Definition of 16 features per sequence. These are the 16 features extracted from each sequence used in training. Some of these features are extracted by .sh scripts that are called by the predict.sh script. See the helper scripts: run_pentamer.sh, run_prodigal.sh, run_plassketch.sh, run_plasORIsketch.sh, run_chromsketch.sh, comparesketch.sh.


____________________

Docker:
The Dockerfile is available under the dl directory of the "docker" branch,
which can be used to create Docker images with the tool.
The latest Docker image can be pulled from dockerhub. For details please see the QuickStart page.
https://sourceforge.net/p/deeplasmid/wiki/deeplasmid%20Quick%20Start/


********************
DeePlasmid Copyright (c) 2019, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov referring to " DelPlasmid" (LBNL Ref 2019-037)."

NOTICE. This software was developed under funding from the U.S. Department of Energy. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, and perform publicly and display publicly. The U.S. Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.


