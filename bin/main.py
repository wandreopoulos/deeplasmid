#!/usr/bin/env python

'''
This is the file that contains the logic for binning, given a fasta file of contigs and a set of BAM files
'''


import MySQLdb

import os
import sys
import logging
# from optparse import OptionParser
import argparse
from collections import defaultdict
from subprocess import Popen, call, PIPE
from multiprocessing import Process, Queue
import operator
import re
import shutil

dir = os.path.dirname(__file__)
from common import get_logger, get_status, run_command, checkpoint_step, append_rqc_file, append_rqc_stats, get_run_path, post_mortem_cmd
from os_utility import runCommand
from db_access import db_connect

#from common import *
#from rqc_utility import *
#from rqc_constants import *
#from rqc_fastq import *

from time import time, strftime
import multiprocessing
import glob
import pysam

from metagenome_binning.rqc_metagenome_binning import run_unsup_binning_test_datasets, run_sup_binning_test_datasets, binning_cmd_sup_blastnt_megan, binning_cmd_sup_blastfungal_megan, binning_cmd_sup_blastmicrob_megan, binning_cmd_sup_blastfungal_taxmapper, binning_cmd_sup_blastmicrob_taxmapper, binning_cmd_unsup_metabat


#
# Create a timestamp with the format YYYYMMDDHHMMSS.
#
def create_timestamp():
    # 14 digits YYYYMMDDHHMMSS
    year = datetime.datetime.now().year
    month = datetime.datetime.now().month
    day = datetime.datetime.now().day
    hour = datetime.datetime.now().hour
    minute = datetime.datetime.now().minute
    second = datetime.datetime.now().second
    if (month < 10):
        month = "0" + str(month);
    if (day < 10):
        day = "0" + str(day);
    if (hour < 10):
        hour = "0" + str(hour);
    if (minute < 10):
        minute = "0" + str(minute);
    if (second < 10):
        second = "0" + str(second);
    res = str(month) + "/" + str(day) + "/" + str(year) + " " + str(hour) + ":" + str(minute) + ":" + str(second);
    return res;



#
# Use one-liner instead.
#
def get_timestamp():
    return strftime("%m/%d/%Y %H:%M:%S")


def delete_temp_files(nocleanup, output_path, log):
    return 0


pipeline_name = "JGI Microbial Metagenome Binning"
pipeline_version = "1.0"    
readme_headers = {}
readme_bodies = {}

##This constant variable determines if only read1 will be used throughtout the analysis.
##If true, only read1 will be used.
##If false, both reads will be used but ONLY if fastq is pair-ended to start with (else read1 will be used if fastq is single-ended).


def main(output_path, fasta, bam_files, nocleanup=None):

    startruntime = time()

    if nocleanup == True:
        print "The files will not be cleaned up after run"

    fasta_name = os.path.basename(fasta)
    ##TODO remove srf files, allow fastq(cat) and fastq.gz(zcat)
    ###TODO change the exit codes to 2 for file-related problems.
    #if fastq_name.endswith(".srf"):
    #        print "fastq_name should either end with fastq.gz or fastq . fastq_name given: " + str(fastq_name)
    #        exit(1)
    if not fasta_name.endswith('.fasta') and not fasta_name.endswith('.fa') and fasta_name.find('.fa') < 0:
            print "fasta_name should either end with fasta or fa or contain .fa . fasta_name given: " + str(fasta_name)
            return 2
    ##TODO double check if fastq is an abs path
    if not os.path.isabs(fasta):
            ###TODO Convert a relative path to abs path
            fasta = os.path.join(os.getcwd(), fasta)
    ###TODO Check the fastq file exists and is readable os.path.isfile
    if not os.path.exists(fasta):
            print("fasta path does not exist: " + fasta)
            return 2
    if not os.path.isfile(fasta):
            print("fasta file does not exist, is missing or is not readable: " + fasta)
            return 2
   
    if output_path.find("homes") > -1:
            print("The output path can not be under homes. The output_path given: " + output_path)
            return 2

    ## Validate output dir. Create output_directory if it doesn't exist
    if not os.path.isdir(output_path):
        print "Cannot find %s, creating as new" % output_path
        os.makedirs(output_path)    
        os.chmod(output_path, 0775)

    timestamp = strftime("%m%d%Y-%H%M%S")

    cat_cmd = "cat"
    if (fasta.endswith("gz")):
        cat_cmd = "zcat"
    
    # initialize my logger
    # log_file = os.path.join(output_path, "rnaseq_microbe_qc.log")
    # log = logging object
    # log = get_logger("rnaseq_microbe_qc", log_file, log_level)

    # initialize my logger
    log_level = "INFO"
    log_file = os.path.join(output_path, "binning_metagenome_cmd.log")
    #print "Started phix_error pipeline, writing log to: %s" % (log_file)
    
    # log = logging object
    log = get_logger("binning_metagenome_cmd", log_file, log_level)
    log.info("Starting %s: %s", pipeline_name, fasta)



    ###TODO strip off the fastq.gz or fastq from fastq_name. We assume the input file ends with fastq or fastq.gz, this condition was checked.
    fasta_name_prefix = fasta_name.replace(".fasta.gz", "").replace(".fasta", "")
    
    #Make the readme_file writer, which will be passed to the functions individually
    readme_file = os.path.join(output_path, fasta_name_prefix + ".readme")

    # Create the status log file location.
    status_log = os.path.join(output_path, "status.log")
    
    last_status = get_status(status_log, log)

    # create the file log and stats log locations.
    ###rqc_file_log = os.path.join(output_path, "rqc-files.tmp")
    ###rqc_stats_log = os.path.join(output_path, "rqc-stats.tmp")
    
    # cmd = "mv %s %s" % (rqc_file_log, rqc_new_file_log)
    # log.info("- cmd: %s", cmd)

    done = False
    cycle = 0
    
    
    
    if last_status == "complete":
        log.info("Status is complete, not processing.")
        done = True
        return 0

    status = "start"

    if last_status.startswith("resubmit more memory"):
        ###TODO in case we exited with exit code 1 because of memory issues,
        ###set the status intelligently based on the last failed status
        ###so we do not repeat the successful steps at every run.
        status = "start"
        log.info(status + " . Resubmitting with more memory.")
        checkpoint_step(status_log, status)


    TAXONOMIC = binning_cmd_sup_blastmicrob_megan ###binning_cmd_sup_blastnt_megan ###binning_cmd_sup_blastfungal_megan ###binning_cmd_sup_blastmicrob_megan  ###binning_cmd_sup_blastfungal_taxmapper ###binning_cmd_sup_blastmicrob_taxmapper
    UNSUPERVISED = binning_cmd_unsup_metabat
    
    ###TODO consider doing a flush
    readme_headers[0]  = pipeline_name
    
    ###TODO consider putting readme file writing in a separate function
    readme_headers[len(readme_headers)] = "Project Information"
    readme_bodies[len(readme_headers)-1] = ""


    readme_headers[len(readme_headers)] = "Tools and Input Files"
    
    readme_bodies[len(readme_headers)-1] = ""
    ###readme_bodies[len(readme_headers)-1] += "Taxonomic binning: " + TAXONOMIC + "\n"
    readme_bodies[len(readme_headers)-1] += "Unsupervised binning: " + UNSUPERVISED + "\n"
    readme_bodies[len(readme_headers)-1] += "Fasta " + fasta_name + "\n"
    #readme_bodies[len(readme_headers)-1] += "Reference Genome: " + ref_genome + "\n"
    #readme_bodies[len(readme_headers)-1] += "Reference GFF: " + gff_annots_file + "\n"
    #readme_bodies[len(readme_headers)-1] += "IMG Taxon OID: " + os.path.basename(ref_genome).split('.', 1)[0] + "\n"
    #readme_bodies[len(readme_headers)-1] += "Output Path: " + output_path + "\n"
    #readme_bodies[len(readme_headers)-1] += "Output BAM: " + bam_file + "\n"
    #readme_bodies[len(readme_headers)-1] += "Output README: " + readme_file
    ###if ref_genome != ref_genome_output_path:
    ###    readme_file_writer.write("Reference Genome was copied to location: " + ref_genome_output_path + "\n\n")
    
    readme_headers[len(readme_headers)] = "Commands"
    readme_bodies[len(readme_headers)-1] = ""
    readme_bodies[len(readme_headers)-1] += "\n_____________\npipeline: %s -f %s -o %s -b %s \n" % ( str(sys.argv[0]) , fasta, output_path, bam_files)

    unclass = []
    
    if status == "start":
        '''
        binning_cmd = TAXONOMIC ###="/global/projectb/scratch/andreopo/binning/taxonomy/run_blastnt_megan.sh "
        binning_cmd_params = fasta
        output_path2 = os.path.join( output_path , "TAXONOMIC" )
        
        if os.path.isdir(output_path2) and os.path.exists(output_path2):
            shutil.move(output_path2, output_path2 + "_" + timestamp)
            
        # create output_directory if it doesn't exist
        if not os.path.isdir(output_path2):
            os.makedirs(output_path2)
        
        (exit_code, unclass, bins) = run_sup_binning_test_datasets(binning_cmd, output_path2, binning_cmd_params, log)
        if exit_code != 0:
            status = "run_sup_binning_test_datasets failed"
            checkpoint_step(status_log, status)
            done = True
        else:
            status = "run_sup_binning_test_datasets success"
            checkpoint_step(status_log, status)

        if status == "run_sup_binning_test_datasets success" and len(unclass) > 0:
        '''
        binning_cmd = UNSUPERVISED ###="module load metabat; runMetaBat.sh "
        binning_cmd_params = "%s %s" % (fasta, bam_files) ###unclass[0]
        output_path3 = os.path.join( output_path , "UNSUPERVISED" )
        
        if os.path.isdir(output_path3) and os.path.exists(output_path3):
            shutil.move(output_path3, output_path3 + "_" + timestamp)
            
        # create output_directory if it doesn't exist
        if not os.path.isdir(output_path3):
            os.makedirs(output_path3)

        (exit_code, bins) = run_unsup_binning_test_datasets(binning_cmd, output_path3, binning_cmd_params, log)
        if exit_code != 0:
            status = "run_unsup_binning_test_datasets failed"
            checkpoint_step(status_log, status)
            done = True
            ###In this case we will exit with status 1, so we will resubmit with more memory
            status = "resubmit more memory"
            checkpoint_step(status_log, status)
            ###exit with 5 in case of memory problems for resubmission.
            return 5
        else:
            status = "run_unsup_binning_test_datasets success"
            checkpoint_step(status_log, status)


    readme_headers[len(readme_headers)] = "Authors"
    
    readme_bodies[len(readme_headers)-1] = "Bill Andreopoulos - wandreopoulos@lbl.gov"
    
    readme_headers[len(readme_headers)] = "Release Date: " + get_timestamp()

    readme_bodies[len(readme_headers)-1] = "This file was automatically generated by the " + pipeline_name + " pipeline, version " + pipeline_version + "\n\n"


    readme_file_writer = open(readme_file, 'w')
    for i in readme_headers.keys():
        readme_file_writer.write("\n\n" + str(i) + ") " + readme_headers.get(i))
        if i in readme_bodies:
            readme_file_writer.write("\n\n" + readme_bodies.get(i))    
    readme_file_writer.close()
    
    log.info( "Status: " + status)
    print "Status: " + status

    ###checkpoint_step(status_log, status)
    
    exit_code = delete_temp_files(nocleanup, output_path, log)
    if exit_code != 0:
        log.info("delete_temp_files failed")
        print "\ndelete_temp_files failed\n\n"

    log.info("Completed %s: %s" % (pipeline_name, fasta))
    print "Completed %s: %s" % (pipeline_name, fasta)

    totalruntime = str(time() - startruntime)
    log.info( "Total wallclock runtime: " + totalruntime )
    print "Total wallclock runtime: " + totalruntime

    ###If we reached this point, we intend to return exit code 0 so we will not rerun (e.g. no rerun with more memory)
    if status.find("success") > -1:
        checkpoint_step(status_log, "complete")
    if status.find("failed") > -1:
        checkpoint_step(status_log, "failed")
        return 1

    done = True
        
    return 0
