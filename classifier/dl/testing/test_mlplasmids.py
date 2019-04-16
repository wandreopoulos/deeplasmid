#!/usr/bin/env python


import sys
import os
import hashlib
from nose import *
import glob


dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'../'))       ## rqc-pipeline/assemblyqc/lib
#sys.path.append(os.path.join(dir,'../lib'))    ## rqc-pipeline/lib
#sys.path.append(os.path.join(dir,'../tools'))  ## rqc-pipeline/tools
#sys.path.append(os.path.join(dir,'../../lib'))    ## rqc-pipeline/lib
#sys.path.append(os.path.join(dir,'../../tools'))  ## rqc-pipeline/tools

#from common import *

#from assemblyqc_utils import *
#from assemblyqc_constants import *
#from assemblyqc_report import *
#from assemblyqc import *

from time import time, strftime
import shutil

def setup_module(module):
    print ("") # this is to get a newline after the dots
    print ("setup_module before anything in this file")

def teardown_module(module):
    print ("teardown_module after everything in this file")

def my_setup_function():
    print ("my_setup_function")

def my_teardown_function():
    print ("my_teardown_function")


###For each of the cases below, run assemblyqc on the specific sequence and write to the output directory the results.
###Then open the rqc_stats file in the output directory for each:
###compare the contents of the rqc_stats file to a golden standard file that is under test/
###Then open the rqc_files file in the output directory for each:
###compare the contents of the rqc_files file to a golden standard file that is under test/ but also check if the file exists with non-zero size.


###Fungal


###Metagenomic


###Microbial
###fungal, microbial single cell, microbial isolate
###test_fastqs = ['/global/dna/dm_archive/sdm/illumina/00/79/53/7953.2.88059.GCCAAT.fastq.gz',
###'/global/dna/dm_archive/sdm/illumina/00/80/69/8069.1.89378.GGCTAC.fastq.gz',
###'/global/dna/dm_archive/sdm/illumina/00/79/45/7945.1.87371.GCCAA.fastq.gz',
test_fastas = [
#'/global/projectb/scratch/andreopo/AsafPlasmids/MockFASTAs/ALL.fasta'
'ALL.fasta' , 'ALL3.fasta' , 'final.assembly.fasta'
]



#GOLD_features_yml = ['/global/projectb/scratch/andreopo/AsafPlasmids/MockFASTAs/yml']
GOLD_formatting = ['ALL.fasta.test-scaff-split.yml' , 'ALL3.fasta.test-scaff-split.yml' , 'final.assembly.fasta.test-scaff-split.yml']
#'/global/projectb/sandbox/rqc/andreopo/src/bitbucket/plasmidoraclev6/dataTest_MockFASTAs/assayer4.test-scaff-split.yml']
#GOLD_DL = ['/global/projectb/sandbox/rqc/andreopo/src/bitbucket/plasmidoraclev6/typescript.data_MockFASTAs']
GOLD_predictions = ['ALL.fasta.predictions.txt' , 'ALL3.fasta.predictions.txt' , 'final.assembly.fasta.predictions.txt']
#'/global/projectb/sandbox/rqc/andreopo/src/bitbucket/plasmidoraclev6/dataTest_MockFASTAs.predictions.txt']
 
 
def compare_predictions(txt1, txt2):
    a = {}
    b = {}
    with open(txt1, 'r') as f1:
       for l1 in f1:
          n = l1.split(',')
          a[n[0]] = n[1]

    with open(txt2, 'r') as f2:
       for l2 in f2:
          n = l2.split(',')
          b[n[0]] = n[1]

    for i in a:
       assert i in b
       assert a[i] == b[i]


def compare_yaml(yaml1, yaml2):
    bulk1=read_yaml(yaml1)
    bulk2=read_yaml(yaml2)

    scaffD1={}
    for seg in range(6):
        if seg in bulk1:
            #print('ww',seg,type(bulk[seg]),len(scaffD))
            scaffD1.update(bulk1[seg])

    scaffD2={}
    for seg in range(6):
        if seg in bulk2:
            #print('ww',seg,type(bulk[seg]),len(scaffD))
            scaffD2.update(bulk2[seg])

    #features to compare
    #import numpy as np
    globFeatureL = ['gc_content','len_sequence']
    for x in 'ACTG' :
        globFeatureL.append(x+'_longestHomopol')
        globFeatureL.append(x+'_totalLongHomopol')

    for scaffN in scaffD1:  
        seqStr=scaffD1[scaffN]['seq']
        txt = scaffD1[scaffN]['text'].lower()
        assert scaffN in scaffD2
        assert seqStr == scaffD2[scaffN]['seq']
        assert txt == scaffD2[scaffN]['text'].lower()
    
        floatD1=scaffD1[scaffN]['features']
        floatL1=[ floatD1[x] for x in globFeatureL ]
        floatD2=scaffD2[scaffN]['features']
        floatL2=[ floatD2[x] for x in globFeatureL ]
        assert len(floatL1) == len(floatL2)
        for x in globFeatureL:
           assert int(floatD1[x]) == int(floatD2[x])
        #for i in range(len(floatL1)):
        #   assert floatL1[i] == floatL2[i]




def compare_csv(csv1,csv2):
    assert 1

@with_setup(my_setup_function, my_teardown_function)
def test_assembly():
    timestamp = strftime("%m%d%Y-%H%M%S")
    i=0
    for fasta in test_fastas:
    
        ###my_name = "TEST " + fasta
        base = os.path.basename(fasta)
        output_path = os.path.join( os.getcwd() , base+"."+timestamp )
        
        ###Run feature extraction
        ###Read yaml file
        ###Read gold yaml file
        ###Compare values . Assert equal or within +-1%
        cm="feature_DL_plasmid_predict_CORI.sh %s %s" % (fasta, output_path)
        os.system(cm)

        formatted_yaml = glob.glob(output_path + "/dlDataFormattedPred/*.yml")[0]

        #Ensure all yaml file values are close
        compare_yaml(GOLD_formatting[i], formatted_yaml)

        #Ensure the file exists and the predictions agree
        compare_predictions(GOLD_predictions[i], output_path + "/outPR/predictions.txt")

        i += 1
        '''
        ###Read the golden standard files
        golden_path = os.path.join( os.getcwd() , "gold_" + base )
        assert os.path.exists(golden_path)
        assert os.path.isdir(golden_path)
    
        ###Compare golden_path/yml with output_path/yml
        for file1 in 
           file2=
           assert file2 exists
           for tag in tags
                assert file1.tag in file2.tag +- 1%

        #open true plasmid and true genome files
        #check the DL output
        #assert precision > 90%, recall > 75%
        
        log_level = "INFO"
        max_subsampl_reads = "5000000"
        
        if os.path.isdir(output_path) and os.path.exists(output_path):
            shutil.move(output_path, output_path + "_" + timestamp)
            
        # create output_directory if it doesn't exist
        if not os.path.isdir(output_path):
            os.makedirs(output_path)	
    
        # initialize my logger
        log_file = os.path.join(output_path, "assemblyqc.log")
        #print "Started phix_error pipeline, writing log to: %s" % (log_file)
        
        # log = logging object
        log = get_logger("assemblyqc", log_file, log_level)
    
    
        log.info("Starting %s: %s", my_name, fastq)
    
        
        ## Validate output dir and fastq
        
        # create output_directory if it doesn't exist
        if not os.path.isdir(output_path):
            log.info("Cannot find %s, creating as new", output_path)
            os.makedirs(output_path)
        
        assert os.path.isdir(output_path)
        
        assert os.path.isfile(fastq)
            
        # Create the status log file location.
        status_log = os.path.join(output_path, "status.log")
        
        # create the file log and stats log locations.
        rqc_file_log = os.path.join(output_path, "rqc-files.txt")
        rqc_stats_log = os.path.join(output_path, "rqc-stats.txt")
        
        # cmd = "mv %s %s" % (rqc_file_log, rqc_new_file_log)
        # log.info("- cmd: %s", cmd)
        
        
    
        log.info( "The fastq file is " + str(fastq) )
        kmer = 0
        read_length = 0
        read_length_1 = 0
        read_length_2 = 0
        is_pe = True
    
        (read_length, read_length_1, read_length_2, is_pe) = get_working_read_length(fastq, log);
        assert read_length > 0
    
        #Save the read_lengths and pairedness to the statistics
        if is_pe == True:
            append_rqc_stats(rqc_stats_log, AssemblyStats.ILLUMINA_READ_LENGTH_1, read_length_1)
            append_rqc_stats(rqc_stats_log, AssemblyStats.ILLUMINA_READ_LENGTH_2, read_length_2)
        else:
            append_rqc_stats(rqc_stats_log, AssemblyStats.ILLUMINA_READ_LENGTH_1, read_length)
    
        kmer = 63
        maxkmer = 63
        kmer = read_length_to_kmer(read_length, maxkmer); #KMER_TEST
    
        log.info("kmer " + str(kmer))
    
        status, dest_fastq = illumina_subsampling(fastq, output_path, status_log, max_subsampl_reads, log)
        assert status.endswith("start")
        # The last True is for skipping megan because in testing it may cause problems.
        status = illumina_assemble_velvet(kmer, fastq, output_path, status_log, rqc_stats_log, rqc_file_log, log, 16, True)
        assert status.endswith("success") or status.endswith("complete")
        status = illumina_assembly_level_report(output_path, status_log, rqc_stats_log, rqc_file_log, kmer, read_length_1, log, True)
        assert status.endswith("success") or status.endswith("complete")
    
        compare(rqc_stats_log, gold_rqc_stats_log)
        
        compare(rqc_file_log, gold_rqc_file_log)
    
        print "Status: " + status
        
        checkpoint_step(status_log, status)
    
        if status.find("failed") > -1:
            checkpoint_step(status_log, "failed")
    
        log.info("Completed %s: %s", my_name, fastq)
        '''
