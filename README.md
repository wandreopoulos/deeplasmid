

## README file
### Codebase: deeplasmid master branch 
#### Note: the master branch was mainly used for native code runs on Cori-NERSC. Most users will prefer to use the docker branch, which is the default branch, unless you want to setup a native code run.
### December 22, 2022
### Maintainer: Bill Andreopoulos, wandreopoulos@lbl.gov


Deeplasmid is a tool based on machine learning that separates plasmids from chromosomal sequences. It can identify plasmids in microbial isolate or metagenome assemblies. The input sequences are in the form of contigs and could have been produced from any sequencing technology or assembly algorithm. The deep learning model was trained on a corpus of:
1) plasmids from ACLAME, 2) PLSDB, and 
3) chromosomal sequences from refseq.microbial (bacteria and archaea) from which plasmids and mito were removed.

Below are instructions for using both the GPU and CPU Docker images for deeplasmid. It is recommended to use GPU (Docker image built December 22, 2022). If you need help with rebuilding the Docker image please feel free to contact me (a reply may be delayed, but I will get back). 
The train and test data used can be found under: https://portal.nersc.gov/dna/microbial/assembly/deeplasmid/


The output file predictions.txt contains plasmid predictions for all contigs in the input fasta file for each run.
Each contig name will indicate if it was:
- Plasmid with score near 1.0
- Chromosome (non-plasmid) with score near 0.0
- Ambiguous have scores around 0.5 (gray zone)
- Longer than 330k bases (possibly a chromosome or megaplasmid)
- Shorter than 1k bases (inconclusive)

The following 2 score files are also output; they are histograms that show the scores for samples and scaffolds that had or didn't have "plasmid" in the header. In reality headers are unlikely to have "plasmid" in the header, but these figures may be useful if testing the tool on a dataset where the classes are known.
- samplescore_hist.png
- scaffscore_hist.png
These files are not output by default. The final plotting step is skipped by setting the -noXterm(-X) flag in the .sh script to true by default.


The Git branch "docker" codebase was used for building the Docker images for deeplasmid.
The public Docker repository is available under:
https://hub.docker.com/repository/docker/billandreo/deeplasmid



## Deeplasmid Docker image for GPU

Input: a .fasta file

Output: a directory of results

Built and tested on: Ubuntu 20.04.3 with NVIDIA GEFORCE RTX 3090

To run the deeplasmid Docker container first
install Docker on your system, then register on dockerhub.
Pull the deeplasmid image from dockerhub as follows:

```
docker login
docker pull billandreo/deeplasmid.tf.gpu2
```

Run deeplasmid for plasmid identification on GPU as follows (note you may need to run docker with sudo on your system):

```
~/Downloads/deeplasmid/classifier/dl$ sudo /usr/bin/docker run -it       -v `pwd`/testing/649989979/649989979.fna:/srv/jgi-ml/classifier/dl/in.fasta  -v  `pwd`/testing/649989979/649989979.fna.OUT:/srv/jgi-ml/classifier/dl/outdir   billandreo/deeplasmid.tf.gpu2   deeplasmid.sh  in.fasta outdir
```

GPU result:
```
~/Downloads/deeplasmid/classifier/dl/testing/649989979$ grep PLASMID 649989979.fna.OUT/outPR.20220209*/predictions.txt
nz_adhj01000031 paenibacillus vortex v453 cnt_pvor1000031, whole genome shotgun sequence.,PLASMID,0.999 +/- 0.000
nz_adhj01000041 paenibacillus vortex v453 cnt_pvor1000041, whole genome shotgun sequence.,PLASMID,0.899 +/- 0.000
nz_adhj01000046 paenibacillus vortex v453 cnt_pvor1000046, whole genome shotgun sequence.,PLASMID,0.509 +/- 0.002
```

### Building the Docker image for GPU

The Dockerfile.GPU2 Docker image was built on tensorflow/tensorflow:latest-gpu ( https://hub.docker.com/r/tensorflow/tensorflow/ https://github.com/tensorflow/tensorflow ).

To build the Docker image, stop or remove your unused containers and images to make space on your drive and use docker build with Dockerfile.GPU2:
```
    sudo docker rm $(sudo docker ps --filter status=exited -q)
    sudo docker images
    sudo docker rmi ...ids....
    sudo docker build -t billandreo/deeplasmid.tf.gpu2 -f Dockerfile.GPU2 .
```

If needed, make space on /var (https://askubuntu.com/questions/1219555/low-disk-space-on-var https://askubuntu.com/questions/178909/not-enough-space-in-var-cache-apt-archives). The second method is quick but may remove images you are using:
```
du -ks /var/* | sort -nr | more
docker system prune -a -f
```

Please see the Supplementary Information from the publication for things to consider when building the Docker image: Prodigal and bbtools/sketch need to be built, and the model .h5 files from training are needed, as well as several sketch files and Pfam-A.TMP2.hmm that can be downloaded from https://portal.nersc.gov/dna/microbial/assembly/deeplasmid/ .


### Troubleshooting a GPU run

Docker is the easiest way to run TensorFlow on a GPU since the host machine only requires the NVIDIA® driver (the NVIDIA® CUDA® Toolkit is not required).
Install the Nvidia Container Toolkit to add NVIDIA® GPU support to Docker. nvidia-container-runtime is only available for Linux. See the nvidia-container-runtime platform support FAQ for details.
https://www.tensorflow.org/install/docker
If you are getting an error when running Docker with the Nvidia Container Toolkit, see:
https://github.com/NVIDIA/nvidia-docker/blob/master/README.md#quickstart
https://github.com/NVIDIA/nvidia-docker/issues/1243 (follow the commands under the Sep 18, 2020 post)

In case your execution is slow it might not be using the GPU.
If you have DL packages like Pytorch on your host computer, check if those can detect the GPU:
```
>>> import torch
>>> torch.cuda.is_available()
True
>>> torch.cuda.device_count()
1
>>> torch.cuda.current_device()
0
>>> torch.cuda.get_device_name(0)
'NVIDIA GeForce RTX 3090'
```


Check if your GPUs are detected on your host with these commands: 
```
    nvidia-smi -q
    nvidia-smi -L
    nvidia-smi   
    sudo apt install nvidia-cuda-toolkit
    nvcc --version
    whereis cuda
    cd /usr/local/cuda*/samples/ && make
    ./bin/x86_64/linux/release/deviceQuery
    ./bin/x86_64/linux/release/bandwidthTest
    cat /usr/local/cuda/version.txt
```

For details on these commands see: https://www.geeksforgeeks.org/how-to-check-if-tensorflow-is-using-gpu/  https://varhowto.com/check-cuda-version-ubuntu-18-04/
In case the commands above don't work on your system here are suggestions and ideas: 
https://askubuntu.com/questions/902636/nvidia-smi-command-not-found-ubuntu-16-04 
https://stackoverflow.com/questions/43022843/nvidia-nvml-driver-library-version-mismatch
Sometimes rebooting your system may work too.

Cuda Toolkit libraries aren't needed for TensorFlow, just the NVIDIA Driver and NVIDIA Container Toolkit are (see above). 
If it can not see the GPU it may be an issue with your Cuda libraries or NVIDIA driver. For Cuda installation, see https://developer.nvidia.com/cuda-downloads https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=Ubuntu&target_version=20.04&target_type=deb_network https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#ubuntu-installation  . You can install NVIDIA drivers from your system settings, or to download new NVIDIA drivers see https://www.nvidia.com/Download/index.aspx .

If you want to setup tensorflow and run the code natively: https://www.tensorflow.org/install/gpu

To read more on Nvidia Docker images for deeplearning: https://docs.nvidia.com/deeplearning/frameworks/user-guide/index.html

#### CNTK (old):

To setup CNTK on your host (outside of the Docker image) see https://docs.microsoft.com/en-us/cognitive-toolkit/setup-cntk-on-your-machine
The old Dockerfile.GPU was extended from the keras-gpu container (https://github.com/gw0/docker-keras); however I downgraded to Cntk 2.3.1 since 2.4 appears to have a compatibility issue with Ryzen processors (see https://github.com/microsoft/CNTK/issues/2908). 
The docker-keras github repo above also provides keras-tensorflow-gpu dockerfiles, but I couldn't get them to build.
In case your execution is slow and you suspect it might not be using the GPU, you can verify if Keras and cntk uses a GPU inside the running Docker container, as follows:
 
```
$ sudo docker run -it  --rm   $(ls /dev/nvidia* | xargs -I{} echo '--device={}') $(ls /usr/lib/*-linux-gnu/{libcuda,libnvidia}* | xargs -I{} echo '-v {}:{}:ro')   -v `pwd`:/srv  a8d777f61654  /bin/sh
#python3 
Python 3.5.3 (default, Nov  4 2021, 15:29:10) 
[GCC 6.3.0 20170516] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import cntk as C
>>> C.device.all_devices()
(GPU[0] NVIDIA GeForce RTX 3090, CPU)
>>> from keras import backend as K
Using CNTK backend
Selected GPU[0] NVIDIA GeForce RTX 3090 as the process wide default device.>>> 
```




## Deeplasmid Docker image for CPU-only on Ubuntu 20.04.3 with Ryzen processor

Input: a .fasta file

Output: a directory of results

Built and tested on: Ubuntu 20.04.3 with Ryzen processor

To run the deeplasmid Docker container first
install Docker on your system, then register on dockerhub.
Pull the deeplasmid image from dockerhub as follows:

```
docker login
docker pull billandreo/deeplasmid-cpu-ubuntu2
```

You can run deeplasmid for plasmid identification as follows. Substitute the `/path/to/input/fasta` and `/path/to/output/directory` below with the full paths to your input file and output dir (note you may need to run docker with sudo):
```
docker run -it -v /path/to/input/fasta:/srv/jgi-ml/classifier/dl/in.fasta -v /path/to/output/directory:/srv/jgi-ml/classifier/dl/outdir billandreo/deeplasmid-cpu-ubuntu2 deeplasmid.sh in.fasta outdir
```

### Building the Docker image for CPU-only on Ubuntu with Ryzen processor

Dec 2022 update: I rebuilt a different image on Ubuntu using Dockerfile.CPU-UbuntuBuild2 with Tensorflow.                                
The Dockerfile.CPU-UbuntuBuild2 Docker image was built on tensorflow/tensorflow ( https://hub.docker.com/r/tensorflow/tensorflow/ https://github.com/tensorflow/tensorflow ).
An older deeplasmid-cpu-ubuntu build was with Cntk and I had downgraded to Cntk 2.3.1, since 2.4 appears to have a compatibility issue with Ryzen processors: 
```
~/Downloads/deeplasmid/classifier/dl$ sudo docker build -t billandreo/deeplasmid-cpu-ubuntu2 -f Dockerfile.CPU-UbuntuBuild2 .

docker pull billandreo/deeplasmid-cpu-ubuntu2
```



## Deeplasmid Docker image for CPU-only on MacBook Pro 2019 with Intel Processor (the Docker image may also work on other Unix-like systems)

Input: a .fasta file

Output: a directory of results

Built and tested on: MacBook Pro. Disclaimer: It's been a while since I run this one on MacBook.

To run the deeplasmid Docker container first
install Docker on your system, then register on dockerhub.
Pull the deeplasmid image from dockerhub as follows:

```
docker login
docker pull billandreo/deeplasmid
```

You can run deeplasmid for plasmid identification as follows. Substitute the `/path/to/input/fasta` and `/path/to/output/directory` below with the full paths to your input file and output dir (note you may need to run docker with sudo):
```
docker run -it -v /path/to/input/fasta:/srv/jgi-ml/classifier/dl/in.fasta -v /path/to/output/directory:/srv/jgi-ml/classifier/dl/outdir billandreo/deeplasmid feature_DL_plasmid_predict.sh in.fasta outdir
```


### Building the Docker image for CPU-only on MacBook

Building the Docker image (CPU-only) was done with Dockerfile.v2 on a MacBook Pro. This command is from an older build that I did with CNTK:
```
docker build -t billandreo/deeplasmid -f Dockerfile.v2 .
```


One of these Docker builds should work if you are on a Linux-based system, else you can try to rebuild the Docker image. I haven't tested on Windows.


## Testing

The 649989979.fna is a testing file, which was downloaded from IMG taxonoid 649989979.
This file can be found under the master branch in order to reduce space in the docker branch.
This way you can test to verify if your installation of the deeplasmid tool gives the same results as expected, which are shown below.

You can run deeplasmid on this input file as follows (note you may need to run docker with sudo):

```
~/Downloads/deeplasmid/classifier/dl$ sudo docker run -it -v `pwd`/testing/649989979/649989979.fna:/srv/jgi-ml/classifier/dl/in.fasta -v `pwd`/testing/649989979/649989979.fna.OUT:/srv/jgi-ml/classifier/dl/outdir billandreo/deeplasmid.tf.gpu2  deeplasmid.sh in.fasta outdir
Counts: Plasm=3  Ambig=0  Main=44  nCount=47

~/Downloads/deeplasmid/classifier/dl$ sudo docker run -it -v `pwd`/testing/649989979/649989979.fna:/srv/jgi-ml/classifier/dl/in.fasta -v `pwd`/testing/649989979/649989979.fna.OUT:/srv/jgi-ml/classifier/dl/outdir billandreo/deeplasmid-cpu-ubuntu2  deeplasmid.sh in.fasta outdir
Counts: Plasm=3  Ambig=0  Main=44  nCount=47

~/Downloads/deeplasmid/classifier/dl$ docker run -it -v `pwd`/testing/649989979/649989979.fna:/srv/jgi-ml/classifier/dl/in.fasta -v `pwd`/testing/649989979/649989979.fna.OUT:/srv/jgi-ml/classifier/dl/outdir billandreo/deeplasmid feature_DL_plasmid_predict.sh in.fasta outdir
Counts: Plasm=3  Ambig=0  Main=44  nCount=47
```

Then you can check the plasmid identified contigs as follows:

```
~/Downloads/deeplasmid/classifier/dl$ grep PLASM testing/649989979/649989979.fna.OUT/outPR.20220210_085843/predictions.txt
nz_adhj01000041 paenibacillus vortex v453 cnt_pvor1000041, whole genome shotgun sequence.,PLASMID,0.899 +/- 0.000
nz_adhj01000031 paenibacillus vortex v453 cnt_pvor1000031, whole genome shotgun sequence.,PLASMID,0.999 +/- 0.000
nz_adhj01000046 paenibacillus vortex v453 cnt_pvor1000046, whole genome shotgun sequence.,PLASMID,0.510 +/- 0.002
```





## Training

The training data and some testing data can be downloaded from https://portal.nersc.gov/dna/microbial/assembly/deeplasmid/ .
The training of the deeplasmid deep learning model was done on Cori at NERSC with Tensorflow in 2021. 
These are the training steps:

Get an salloc session with 48 hours allocation:
salloc: sal with 48 hours

Then run as follows:
```
./feature_DL_plasmid_train_CORI.sh plasmid.fasta plasmid.fasta.OUTDIR nonplasmid.fasta nonplasmid.fasta.OUTDIR
```

For example:
```
./feature_DL_plasmid_train_CORI.sh  ../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta    ../DATA/ACLAME.REFSEQMICROB/aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.OUT18/    ../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.UNION.lists_shortq.archaea.txt.fasta.MIN1kMAX330k.fasta       ../DATA/ACLAME.REFSEQMICROB/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.UNION.lists_shortq.archaea.txt.fasta.MIN1kMAX330k.fasta.OUT18/
Running feature_DL_plasmid_train.sh . This version is meant for Cori.
```

.....

After training completes you will have a trained model under the two directories you specified on the command line. You may use the trained model in the predictions file feature_DL_plasmid_predict_CORI.sh

The codebase used for training on Cori (as well as running prediction on Cori) is under the master branch of the code repository.

A single training data element consists of the label and two input words: xseq - a 300bp contiguous subsequence sampled randomly from the full original contig sequence and xf - a vector containing 16 features extracted from the full sequence, as described in the Table below. The number (m) of 300bp subsequences sampled from each contig is proportional to the square root of the contig length, such that longer contigs contribute more samples, but do not overwhelm the training.

```
Feature Name	Definition	Type
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
```



## Citation:

William B Andreopoulos, Alexander M Geller, Miriam Lucke, Jan Balewski, Alicia Clum, Natalia N Ivanova, Asaf Levy, Deeplasmid: deep learning accurately separates plasmids from bacterial chromosomes, Nucleic Acids Research, 2021;, gkab1115, https://doi.org/10.1093/nar/gkab1115 


## Docker:

The Dockerfile is available under the dl directory of the "docker" branch,
which can be used to create Docker images with the tool.
The latest Docker images can be pulled from dockerhub. 


********************
DeePlasmid Copyright (c) 2019, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov referring to " DelPlasmid" (LBNL Ref 2019-037)."

NOTICE. This software was developed under funding from the U.S. Department of Energy. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, and perform publicly and display publicly. The U.S. Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.


