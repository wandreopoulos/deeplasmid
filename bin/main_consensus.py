#!/usr/bin/env python

'''
This is the file that contains the logic for consensus clustering or binning
with taxonomic information, given a fasta file of contigs and a set of BAM files
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

from rqc_metagenome_binning import run_unsup_binning_test_datasets, run_sup_binning_test_datasets, binning_cmd_sup_blastnt_megan, binning_cmd_sup_blastfungal_megan, binning_cmd_sup_blastmicrob_megan, binning_cmd_sup_blastfungal_taxmapper, binning_cmd_sup_blastmicrob_taxmapper, binning_cmd_unsup_metabat


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


pipeline_name = "JGI Binning"
pipeline_version = "1.0"    
readme_headers = {}
readme_bodies = {}

##This constant variable determines if only read1 will be used throughtout the analysis.
##If true, only read1 will be used.
##If false, both reads will be used but ONLY if fastq is pair-ended to start with (else read1 will be used if fastq is single-ended).


def main(output_path, fasta, bam_files, reftype="m", nocleanup=None):

    startruntime = time()
    
    print "Binning pipeline version %s %s. " % (pipeline_name, pipeline_version)
    print "Author: Bill Andreopoulos. Date: January 16, 2016."
    print "This binning pipeline is aimed for decontamination of datasets and will output 3 sets of clusters using different approaches: "
    print "  1) TAXONOMIC/BINS contains the bins as fasta files from a taxonomic analysis using Megan."
    print "  2) UNSUPERVISED/BINS contains the bins as fasta files from an unsupervised (cluster based) analysis using Metabat."
    print "  3) The files Level*.fa contain the bins as fasta files from a consensus clustering of the Megan and Metabat bins."
    print "     In this case, the Level1* bins are based on the taxonomic analysis and the Level2* bins are based on Metabat."
    print "Each contig can belong in only one bin/cluster."
    print "While the bins/clusters will represent the organisms in the dataset as reliably as possible,"
    print "interpretation of the clusters is always in the eye of the beholder."

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
    log_file = os.path.join(output_path, "binning_cmd.log")
    #print "Started phix_error pipeline, writing log to: %s" % (log_file)
    
    # log = logging object
    log = get_logger("binning_cmd", log_file, log_level)
    log.info("Starting %s: %s", pipeline_name, fasta)



    ###TODO strip off the fastq.gz or fastq from fastq_name. We assume the input file ends with fastq or fastq.gz, this condition was checked.
    fasta_name_prefix = fasta_name.replace(".fasta.gz", "").replace(".fasta", "")
    
    #Make the readme_file writer, which will be passed to the functions individually
    readme_file = os.path.join(output_path, fasta_name_prefix + ".readme")

    # Create the status log file location.
    status_log = os.path.join(output_path, "binning_status.log")
    
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
    if reftype == "m":
        TAXONOMIC = binning_cmd_sup_blastmicrob_megan
        print "The refseq.microbial database will be used for taxonomic binning. This option is suitable for decontamination and binning of microbial datasets"
    elif reftype == "f":
        TAXONOMIC = binning_cmd_sup_blastfungal_megan
        print "The refseq.fungi database will be used for taxonomic binning. This option is suitable for decontamination and binning of fungal datasets"
    elif reftype == "n":
        TAXONOMIC = binning_cmd_sup_blastnt_megan
        print "The nt database will be used for taxonomic binning. This option is suitable for decontamination and binning of all datasets"
    UNSUPERVISED = binning_cmd_unsup_metabat
    
    ###TODO consider doing a flush
    readme_headers[0]  = pipeline_name
    
    ###TODO consider putting readme file writing in a separate function
    readme_headers[len(readme_headers)] = "Project Information"
    readme_bodies[len(readme_headers)-1] = ""


    readme_headers[len(readme_headers)] = "Tools and Input Files"
    
    readme_bodies[len(readme_headers)-1] = ""
    readme_bodies[len(readme_headers)-1] += "Taxonomic binning: " + TAXONOMIC + "\n"
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
        binning_cmd = TAXONOMIC ###="/global/projectb/scratch/andreopo/binning/taxonomy/run_blastnt_megan.sh "
        binning_cmd_params = fasta
        output_path2 = os.path.join( output_path , "TAXONOMIC" )
        
        if os.path.isdir(output_path2) and os.path.exists(output_path2):
            shutil.move(output_path2, output_path2 + "_" + timestamp)
            
        # create output_directory if it doesn't exist
        if not os.path.isdir(output_path2):
            os.makedirs(output_path2)
        
        (exit_code, unclass, clusters1) = run_sup_binning_test_datasets(binning_cmd, output_path2, binning_cmd_params, log)
        if exit_code != 0:
            status = "run_sup_binning_test_datasets failed"
            checkpoint_step(status_log, status)
            done = True
        else:
            status = "run_sup_binning_test_datasets success"
            checkpoint_step(status_log, status)

    if status == "run_sup_binning_test_datasets success" and len(unclass) > 0:
        binning_cmd = UNSUPERVISED ###="module load metabat; runMetaBat.sh "
        binning_cmd_params = "%s %s" % (fasta, bam_files)
        output_path3 = os.path.join( output_path , "UNSUPERVISED" )
        
        if os.path.isdir(output_path3) and os.path.exists(output_path3):
            shutil.move(output_path3, output_path3 + "_" + timestamp)
            
        # create output_directory if it doesn't exist
        if not os.path.isdir(output_path3):
            os.makedirs(output_path3)

        (exit_code, clusters2) = run_unsup_binning_test_datasets(binning_cmd, output_path3, binning_cmd_params, log)
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

        dondepth = os.path.join( output_path3 , "depth.txt" )
        doncmd = "module load metabat; jgi_summarize_bam_contig_depths --outputDepth %s  %s" % (dondepth, bam_files)
        print(doncmd)
        std_out, std_err, exit_code = run_command(doncmd, True, log)
        if exit_code != 0:
            print "CMD %s failed with OUT %s ERR %s" % (doncmd, std_out, std_err)
            exit(2)
            
        donfile = os.path.join( output_path3 , "bin.dist" )
        doncmd2 = "/global/homes/d/ddkang/program/metabat_dist  -i %s -a %s -o %s" % (fasta, dondepth, donfile)
        print(doncmd2)
        std_out, std_err, exit_code = run_command(doncmd2, True, log)
        if exit_code != 0:
            print "CMD %s failed with OUT %s ERR %s" % (doncmd2, std_out, std_err)
            exit(2)

        ###Bilcom implementation
        ###Interface the python with C++ code
        ###Use UPC or Hadoop
        
        ###Cluster the contig vs. gene id matrix OR
        ###Run taxator-tk OR
        ###Run RapSearch2 to get the hits for the taxlist file, for example:
        ###unknown.fasta.metabat-bins-.14.fa.gene.faa.RAPSearch2.out.m8.parsed.tophit.taxlist
        #Input1: one taxlist file
        ###clusters1 = megan_fastas() ###FA under TAXONOMIC from megan
        
        ###Cluster all contigs that had any hits by the similarity of their taxlist hits and the confidence of the hits
        #Set of clusters1
        
        ###Run unsupervised to produce clusters based on all contigs
        #Input2: set of fasta files
        ###clusters2, donfile = metabat_fastas() ###FA under UNSUPERVISED from metabat
        
        ###Assign seeds in each cluster based on which contigs were placed in clusters1
        #(takes Input1 and Input2)
        #Set of clusters2
        seeds = []
        level1seedsMemberships = {}
        level1seedsMemberships_rev = {}
        level2seedsMemberships = {}
        level2seedsMemberships_rev = {}
        final_binMemberships = {}
        
        ###A count for the clusters, used as cluster ID (for both levels 1 and 2)
        count = 0
        file_lev1_labels = {}
        
        for cluster_fasta_file in clusters1:
            count += 1
            file_lev1_labels[count] = os.path.basename(cluster_fasta_file)
            level1seedsMemberships_rev[count] = []
            seeds_i = get_contigs_list(cluster_fasta_file)
            seeds += seeds_i
            for j in seeds_i:
                level1seedsMemberships[j] = count
                ###print "count %s level1seedsMemberships_rev.get_count=%s" % (count, level1seedsMemberships_rev.get(count))
                level1seedsMemberships_rev[count].append(j)
                final_binMemberships[j] = count
                
        seedset = set(seeds)
        
        sim1 = {}
        sim2 = {}
        
        subclustersLevel2 = {}
        
        separator_lev12 = count
        
        ###cluster1seed2seedSim()
        ###Refine the clusters:
        ###Compute for each seed the unsup pairwise similarity to all seeds in the same cluster1.
        ###def cluster1seed2seedSim():
        #read don's file linear and put in named 2d array cluster1_seed2seedSims
        donfilef = open(donfile, 'r')
        for i in donfilef: ###file with all sims by don
            sims = i.split("\t")
            sim1[frozenset((sims[0], sims[1]))] = float(sims[2])
            sim2[frozenset((sims[0], sims[1]))] = float(sims[3])
        ###for c in clusters1:
        ###upc(seed,)
        ###barrier
        
        ###cluster2seed2allcontigsSimSubclusters()
        ###Compute for each seed the unsup pairwise similarity to all contigs (including seeds) in the same unsupervised cluster2:
        ###Split cluster into subclusters:
        ###Form subclusters by assigning each non-seed contig to the closest seed in the same unsupervised cluster2.
        ###def cluster2seed2allcontigsSimSubclusters():
        #read don's file linear and put in named 2d array cluster2_seed2allcontigsSims
        #pragma omp parallel for
        for cluster_fasta_file in clusters2:
            print "----> Level2 cluster_fasta_file %s" % (cluster_fasta_file)
            count += 1
            level2seedsMemberships_rev[count] = []
            contigs_i = get_contigs_list(cluster_fasta_file)
            longest_contig_found = False
            for ci in contigs_i:
                ###print "       contigi %s" % (ci)
                seed_tmp = []
                nonseed_tmp = []
                if ci in seedset:
                    ###print "seed"
                    seed_tmp.append(ci)
                    subclustersLevel2[ci] = []
                    level2seedsMemberships[ci] = count
                    level2seedsMemberships_rev[count].append(ci)
                elif not longest_contig_found:
                    ###We set the longest contig as additional seed, so in the closest function the non-seeds
                    ###have an option where to go in case there is an unknown organism not detected in step 1.
                    longest_contig_found = True
                    seed_tmp.append(ci)
                    subclustersLevel2[ci] = []
                    level2seedsMemberships[ci] = count
                    level2seedsMemberships_rev[count].append(ci)
                else:
                    print "nonseed contigi %s" % (ci)
                    nonseed_tmp.append(ci)
                final_binMemberships[ci] = count
                if len(seed_tmp) > 0:
                    for n in nonseed_tmp:
                        ###for n in nonseed get closest seed
                        cj = closest(seed_tmp, n, sim1, sim2)
                        print "cj %s seed_tmp %s n %s" % (cj, seed_tmp, n)
                        subclustersLevel2[cj].append(n)
        ###upc(contig,)
        ###barrier
        
        ###compareEachSeedtoAllSeedsSimsinClusters1and2()
        ###Move the subclusters based on the clusters1 results (if their seeds were clustered together in clusters1),
        ###UNLESS the unsup pairwise similarities of the seeds are high, in which case leave in clusters2:
        ###For each seed assign to the seed with which it was in the same cluster1 [and to which it has
        ###the highest unsup pairwise similarity(UPS), as long as the UPS is at least as
        ###high as the highest unsup pairwise similarity to other seeds in clusters2], as long as
        ###the UPS to the closest seed in level 2 is not significant (>0.02).
        ###def compareEachSeedtoAllSeedsSimsinClusters1and2():
        #pragma omp parallel for
        for seed in subclustersLevel2.keys():
            ###(sim1_c1, seed1_c1) = nearest_sim_neighbor_seed(seed, sim1, seedset)
            ###(sim2_c1, seed2_c1) = nearest_sim_neighbor_seed(seed, sim2, seedset)
            #(sim1_c2, seed1_c2) = nearest_sim_neighbor_seed(seed, sim1, set(level2seedsMemberships_rev[level2seedsMemberships[seed]]))
            #(sim2_c2, seed2_c2) = nearest_sim_neighbor_seed(seed, sim2, set(level2seedsMemberships_rev[level2seedsMemberships[seed]]))
            #if sim1_c2 > 0.02 and sim2_c2 > 0.02 and seed1_c2 != "" and seed2_c2 != "":  ###sim1_c1
            if level1seedsMemberships.has_key(seed):
                print "------> switch for seed %s from level2 %s to level1 %s" % (seed, level2seedsMemberships[seed], level1seedsMemberships[seed])
                final_binMemberships[seed] = level1seedsMemberships[seed] ###seed1_c1]
                for s in subclustersLevel2[seed]:
                    print "     --------------> switch for member %s to level1 %s" % (s, level1seedsMemberships[seed])
                    final_binMemberships[s] = level1seedsMemberships[seed]
            ###if sim2_c2 > 0.1 :  ###sim2_c1
            ###    final_binMemberships[seed] = level1seedsMemberships[seed] ###seed2_c1]
            ###decide if subcluster is assigned based on cluster1_seed2seedSims or cluster2_seed2allcontigsSims  upc(seed, cluster1highestSim, cluster2highestSim)
        
        ###print str(final_binMemberships)
        inv_map = {}
        for k, v in final_binMemberships.iteritems():
            inv_map[v] = inv_map.get(v, [])
            inv_map[v].append(k)
        ###print str(inv_map)

        for k in inv_map:
            print "_________________"
            ###print str(inv_map.get(k))
            if k <= separator_lev12:
                filename1 = "Level1." + str(k) + ".____." + file_lev1_labels[k] + ".fa"
                print "Cluster level1 bin %s written to file %s" % (k, filename1)
                f = open( os.path.join( output_path, filename1 ), 'a' )
                for s in inv_map.get(k):
                    f.write(s + "\n")
                f.close()
            else:
                filename2 = "Level2." + str(k) + ".fa"
                print "Cluster level2 bin %s written to file %s" % (k, filename2)
                f = open( os.path.join( output_path, filename2 ), 'a' )
                for s in inv_map.get(k):
                    f.write(s + "\n")
                f.close()

        exit(0)

    exit(1)

def get_contigs_list(filename):
    contigs = {}
    s = ""
    prev_contig = ""
    f = open(filename, 'r')
    for i in f:
        if i.startswith(">"):
            #If L=length of s is >0:
            #Assign L to previous contig
            #Set string to ""
            if len(s) > 0 and len(prev_contig) > 0:
                if not len(s) in contigs:
                    contigs[len(s)] = []
                contigs.get(len(s)).append(prev_contig)
            ###contigs.append(i[1:].strip())
            prev_contig = i[1:].strip()
            s = ""
        else:
            s += i.strip()
    r = [value for (key, value) in sorted(contigs.items())]
    return list(itertools.chain(*r))

'''
Only attach nonseed to closest seed if the sim significance is <0.02,
since else it might be an unknown organism contig that we are searching for and does not match the seed.
'''
def closest(seed_array, nonseed_item, sim1, sim2):
    simholder1 = 0.02
    simholder2 = 0.02
    res = ""
    for s in seed_array:
        fset = frozenset((nonseed_item,s))
        if fset in sim1 and sim1[fset] <= simholder1 and fset in sim2 and sim2[fset] <= simholder2:
            simholder1 = sim1[fset]
            simholder2 = sim2[fset]
            res = s
    return res
    
def nearest_sim_neighbor_seed(seed, sim, seedset):    
    simholder = 10000.0
    res = ""
    for s in seedset:
        fset = frozenset((seed,s))
        if s != seed and fset in sim and sim[fset] < simholder:
            simholder = sim[fset]
            res = s
    return simholder, res

###TODO sims lower or higher is more significant?
###TODO assign subcluster nonseeds to level 1 membership as well.
