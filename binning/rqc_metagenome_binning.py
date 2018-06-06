#!/usr/bin/env python

'''
To run as test (not live):
cd /global/projectb/scratch/andreopo/test_rqc_metagenome_binning
export PYTHONPATH=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/assemblyqc/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/readqc/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/lib/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-pipeline/tools/:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc/
export PATH=/global/projectb/scratch/andreopo/binning/metabat:$PATH
source /global/dna/projectdirs/PI/rqc/prod/versions/jgi-rqc-pipeline-20140430/tools/external_tools/pmgw/bin/../config/megan_settings.sh
/global/projectb/sandbox/rqc/andreopo/src/bitbucket/rqc_metagenome_binning/metagenome_binning/rqc_metagenome_binning.py -o  /global/projectb/scratch/andreopo/test_rqc_metagenome_binning  -t  Albertsen2013
To run live create bams first. For this purpose, you can check run.sh under previous runs of rnaseq pipeline
'''
import sys
import os
root = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(root)

dir = os.path.dirname(__file__)
from common import get_logger, get_status, run_command, checkpoint_step, append_rqc_file, append_rqc_stats, get_run_path, post_mortem_cmd

import glob
import logging
# from optparse import OptionParser
import argparse
from collections import defaultdict
from subprocess import Popen, call, PIPE
import operator
import re
from time import time, strftime
import shutil

megablast_output_holder = None


#metagenome_binning_clustering_quality
#compute against all

###CONSTANTS
binning_cmd_unsup_metabat = "binning_cmd_unsup_metabat"

binning_cmd_sup_nt_bwa = "binning_cmd_sup_nt_bwa"
binning_cmd_sup_nt_mb = "binning_cmd_sup_nt_mb"
binning_cmd_sup_nt_bbmap = "binning_cmd_sup_nt_bbmap"
binning_cmd_sup_nt_blat = "binning_cmd_sup_nt_blat"
binning_cmd_sup_nt_bowtie2 = "binning_cmd_sup_nt_bowtie2"
binning_cmd_sup_nt_tophat = "binning_cmd_sup_nt_tophat"
binning_cmd_sup_nt_megan = "binning_cmd_sup_nt_megan"

binning_cmd_sup_microb_bwa = "binning_cmd_sup_microb_bwa"
binning_cmd_sup_microb_mb = "binning_cmd_sup_microb_mb"
binning_cmd_sup_microb_bbmap = "binning_cmd_sup_microb_bbmap"
binning_cmd_sup_microb_blat = "binning_cmd_sup_microb_blat"
binning_cmd_sup_microb_bowtie2 = "binning_cmd_sup_microb_bowtie2"
binning_cmd_sup_microb_tophat = "binning_cmd_sup_microb_tophat"
binning_cmd_sup_microb_megan = "binning_cmd_sup_microb_megan"

binning_cmd_sup_nt16s_bwa = "binning_cmd_sup_nt16s_bwa"
binning_cmd_sup_nt16s_mb = "binning_cmd_sup_nt16s_mb"
binning_cmd_sup_nt16s_bbmap = "binning_cmd_sup_nt16s_bbmap"
binning_cmd_sup_nt16s_blat = "binning_cmd_sup_nt16s_blat"
binning_cmd_sup_nt16s_bowtie2 = "binning_cmd_sup_nt16s_bowtie2"
binning_cmd_sup_nt16s_tophat = "binning_cmd_sup_nt16s_tophat"
binning_cmd_sup_nt16s_megan = "binning_cmd_sup_nt16s_megan"

binning_cmd_sup_nr_last = "binning_cmd_sup_nr_last"

binning_cmd_sup_blastnt_megan = "binning_cmd_sup_blastnt_megan"
binning_cmd_sup_blastfungal_megan = "binning_cmd_sup_blastfungal_megan"
binning_cmd_sup_blastmicrob_megan = "binning_cmd_sup_blastmicrob_megan"
binning_cmd_sup_blastfungal_taxmapper = "binning_cmd_sup_blastfungal_taxmapper"
binning_cmd_sup_blastmicrob_taxmapper = "binning_cmd_sup_blastmicrob_taxmapper"

###These tools will be tested in this order.
###Changing the order may affect megan runs, which depend on blast output.
holder_of_runs = [
    binning_cmd_sup_microb_megan   ###binning_cmd_sup_nt_megan   ###binning_cmd_sup_nt16s_megan, 
    ]

'''
more datasets: decontamination fungal, decontamination microbial, CAMI low/med/high, MetaHit.
more options: 2 binning methods classic and consensus.
'''

'''
    binning_cmd_unsup_metabat, 
                  binning_cmd_sup_nt_mb, binning_cmd_sup_microb_mb, binning_cmd_sup_nt16s_mb,
                  binning_cmd_sup_nt_megan, binning_cmd_sup_microb_megan, binning_cmd_sup_nt16s_megan
'''

'''
                  binning_cmd_sup_nt_bbmap, binning_cmd_sup_microb_bbmap, binning_cmd_sup_nt16s_bbmap,
                  binning_cmd_sup_nt_blat, binning_cmd_sup_microb_blat, binning_cmd_sup_nt16s_blat,
                  binning_cmd_sup_nr_last
'''


###MetaBat
binning_cmd_unsup = {
    binning_cmd_unsup_metabat : "module load metabat; runMetaBat.sh "
    ###binning_cmd_unsup_metabat : "/global/projectb/scratch/andreopo/binning/metabat/runMetaBat.sh "
}
binning_cmd_sup = {
    binning_cmd_sup_nt_mb : "/global/projectb/scratch/andreopo/binning/taxonomy/nt_megablast.sh  ",
    #binning_cmd_sup_nt_bbmap : "/global/projectb/scratch/andreopo/binning/taxonomy/nt_bbmap.sh  ",
    #binning_cmd_sup_nt_blat : "/global/projectb/scratch/andreopo/binning/taxonomy/nt_blat.sh  ",
    binning_cmd_sup_nt_megan : "/global/projectb/scratch/andreopo/binning/taxonomy/run_megan.sh  ",
    #binning_cmd_sup_nt_bwa : "/global/projectb/scratch/andreopo/binning/taxonomy/nt_bwa.sh -t 90  ",
    #binning_cmd_sup_nt_bowtie2 : "/global/projectb/scratch/andreopo/binning/taxonomy/nt_bowtie2.sh -t 90  ",
    #binning_cmd_sup_nt_tophat : "/global/projectb/scratch/andreopo/binning/taxonomy/nt_tophat.sh -t 90  ",
    binning_cmd_sup_microb_mb : "/global/projectb/scratch/andreopo/binning/taxonomy/microb_megablast.sh  ",
    #binning_cmd_sup_microb_bbmap : "/global/projectb/scratch/andreopo/binning/taxonomy/microb_bbmap.sh  ",
    #binning_cmd_sup_microb_blat : "/global/projectb/scratch/andreopo/binning/taxonomy/microb_blat.sh  ",
    binning_cmd_sup_microb_megan : "/global/projectb/scratch/andreopo/binning/taxonomy/run_megan.sh  ",
    #binning_cmd_sup_microb_bwa : "/global/projectb/scratch/andreopo/binning/taxonomy/microb_bwa.sh -t 90  ",
    #binning_cmd_sup_microb_bowtie2 : "/global/projectb/scratch/andreopo/binning/taxonomy/microb_bowtie2.sh -t 90  ",
    #binning_cmd_sup_microb_tophat : "/global/projectb/scratch/andreopo/binning/taxonomy/microb_tophat.sh -t 90  ",
    binning_cmd_sup_nt16s_mb : "/global/projectb/scratch/andreopo/binning/taxonomy/nt16s_megablast.sh  ",
    #binning_cmd_sup_nt16s_bbmap : "/global/projectb/scratch/andreopo/binning/taxonomy/nt16s_bbmap.sh  ",
    #binning_cmd_sup_nt16s_blat : "/global/projectb/scratch/andreopo/binning/taxonomy/nt16s_blat.sh  ",
    binning_cmd_sup_nt16s_megan : "/global/projectb/scratch/andreopo/binning/taxonomy/run_megan.sh  ",
    #binning_cmd_sup_nt16s_bwa : "/global/projectb/scratch/andreopo/binning/taxonomy/nt16s_bwa.sh -t 90  ",
    #binning_cmd_sup_nt16s_bowtie2 : "/global/projectb/scratch/andreopo/binning/taxonomy/nt16s_bowtie2.sh -t 90  ",
    #binning_cmd_sup_nt16s_tophat : "/global/projectb/scratch/andreopo/binning/taxonomy/nt16s_tophat.sh -t 90  ",
    #binning_cmd_sup_nr_last : "/global/projectb/scratch/andreopo/binning/taxonomy/nr_last.sh  "
    binning_cmd_sup_blastnt_megan : "/global/projectb/scratch/andreopo/binning/taxonomy/run_blastnt_megan.sh ",
    binning_cmd_sup_blastfungal_megan : "/global/projectb/scratch/andreopo/binning/taxonomy/run_blastfungal_megan.sh ",
    binning_cmd_sup_blastmicrob_megan : "/global/projectb/scratch/andreopo/binning/taxonomy/run_blastmicrob_megan.sh ",
    binning_cmd_sup_blastfungal_taxmapper : "/global/projectb/scratch/andreopo/binning/taxonomy/run_blastfungal_taxmapper.sh",
    binning_cmd_sup_blastmicrob_taxmapper : "/global/projectb/scratch/andreopo/binning/taxonomy/run_blastmicrob_taxmapper.sh"
}

###TODO evaluate metabat addons: prior taxonomy-aware, posterior plus GC and TNS.
###TODO compare other methods outside metabat on several test datasets and justify why not used them.
#binning_cmd2 = "concoct"
#binning_cmd3 = "maxbin"
#binning_cmd4 = "metacluster"
#binning_cmd5 = "multi-metagenome"


###TODO add more cases of test_datasets, gold_standard_bins and output_dirs (check Alicia's emails, Human gut microbiome data, Metabat paper dataset, Multi-metagenome paper dataset)

###The libraries of each project are placed in an array.
binning_cmd_params = {
binning_cmd_unsup_metabat : {
"CAMI_high" : "cami.fa bam bam",
"mockup_microbial" : "/global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/1030534.scaffolds.fasta   /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/clostridium_perfringensatcc_13124/clostridium_perfringensatcc_13124-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/clostridium_thermocellumvpi_7372_atcc_27405/clostridium_thermocellumvpi_7372_atcc_27405-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/coraliomargarita_akajimensis_dsm_45221/coraliomargarita_akajimensis_dsm_45221-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/corynebacterium_glutamicum_atcc_13032/corynebacterium_glutamicum_atcc_13032-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/desulfosporosinus_acidiphilus_sj4_dsm_22704/desulfosporosinus_acidiphilus_sj4_dsm_22704-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/desulfosporosinus_meridiei_dsm_13257/desulfosporosinus_meridiei_dsm_13257-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/desulfotomaculum_gibsoniae_dsm_7213/desulfotomaculum_gibsoniae_dsm_7213-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/e.coli_k12_atcc_700926/e.coli_k12_atcc_700926-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/fervidobacterium_pennivorans_dsm_9078/fervidobacterium_pennivorans_dsm_9078-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/echinicola_vietnamensis_dsm_17526/echinicola_vietnamensis_dsm_17526-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/frateuria_aurantia_dsm_6220/frateuria_aurantia_dsm_6220-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/halovivax_ruber_xh-70/halovivax_ruber_xh-70-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/hirschia_baltica_atcc_49814/hirschia_baltica_atcc_49814-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/meiothermus_silvanus_dsm_9946/meiothermus_silvanus_dsm_9946-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/natronobacterium_gregoryi_sp2/natronobacterium_gregoryi_sp2-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/natronococcus_occultus_dsm_3396/natronococcus_occultus_dsm_3396-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/nocardiopsis_dassonvillei_dsm_43111/nocardiopsis_dassonvillei_dsm_43111-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/olsenella_uli_dsm_7084/olsenella_uli_dsm_7084-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/pseudomonas_stutzeri_rch2/pseudomonas_stutzeri_rch2-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/salmonella_bongori_nctc_12419/salmonella_bongori_nctc_12419-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/salmonella_enterica_subsp._arizonae_serovar_62_z4_z23_-_strain_rsk2980/salmonella_enterica_subsp._arizonae_serovar_62_z4_z23_-_strain_rsk2980-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/segniliparus_rotundus_dsm_44985/segniliparus_rotundus_dsm_44985-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/spirochaeta_smaragdinae_dsm_11293/spirochaeta_smaragdinae_dsm_11293-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/streptococcus_pyogenes_m1_gas/streptococcus_pyogenes_m1_gas-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/terriglobus_roseus_dsm_18391/terriglobus_roseus_dsm_18391-genome_sorted.bam /global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/thermobacillus_composti_kwc4/thermobacillus_composti_kwc4-genome_sorted.bam ",
###Albertsen2013 - assembly of all libs followed by bam file for each lib
###"/global/dna/shared/data/functests/assembly/Meta/Albertsen2013/Albertsen2013.APMI.contigs.fsa   /global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/SRR770019/SRR770019-genome_sorted.bam    /global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/SRR770018/SRR770018-genome_sorted.bam    /global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/SRR770574/SRR770574-genome_sorted.bam ",
"Albertsen2013" : "  -s  10  -m  10 --verysensitive  /global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/assembly.fasta  /global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/SRR770018wriass/SRR770018-genome_sorted.bam  /global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/SRR770019wriass/SRR770019-genome_sorted.bam  /global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/SRR770574wriass/SRR770574-genome_sorted.bam ",
###Human microbiome from MetaBAT paper - assembly of all libs followed by bam file for each lib
"MetaHIT_human_microbiome" :
    """
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/kmernator-FR-342326.edique02.d-Ray-358992.edique02-k51-Ray-358992.edique02.d-Contigs.fasta 
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011094bam/ERR011094-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011144bam/ERR011144-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011149bam/ERR011149-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011129bam/ERR011129-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011095bam/ERR011095-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011162bam/ERR011162-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011137bam/ERR011137-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011143bam/ERR011143-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011099bam/ERR011099-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011098bam/ERR011098-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011103bam/ERR011103-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011102bam/ERR011102-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011111bam/ERR011111-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011151bam/ERR011151-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011124bam/ERR011124-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011161bam/ERR011161-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011125bam/ERR011125-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011131bam/ERR011131-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011126bam/ERR011126-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011120bam/ERR011120-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011140bam/ERR011140-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011159bam/ERR011159-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011154bam/ERR011154-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011105bam/ERR011105-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011110bam/ERR011110-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011121bam/ERR011121-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011104bam/ERR011104-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011118bam/ERR011118-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011116bam/ERR011116-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011096bam/ERR011096-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011106bam/ERR011106-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011157bam/ERR011157-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011109bam/ERR011109-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011113bam/ERR011113-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011123bam/ERR011123-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011122bam/ERR011122-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011087bam/ERR011087-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011138bam/ERR011138-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011093bam/ERR011093-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011119bam/ERR011119-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011164bam/ERR011164-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011091bam/ERR011091-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011088bam/ERR011088-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011090bam/ERR011090-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011097bam/ERR011097-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011145bam/ERR011145-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011150bam/ERR011150-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011163bam/ERR011163-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011133bam/ERR011133-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011146bam/ERR011146-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011089bam/ERR011089-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011139bam/ERR011139-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011101bam/ERR011101-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011135bam/ERR011135-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011100bam/ERR011100-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011127bam/ERR011127-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011134bam/ERR011134-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011132bam/ERR011132-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011117bam/ERR011117-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011128bam/ERR011128-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011165bam/ERR011165-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011108bam/ERR011108-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011136bam/ERR011136-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011148bam/ERR011148-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011147bam/ERR011147-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011115bam/ERR011115-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011112bam/ERR011112-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011141bam/ERR011141-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011152bam/ERR011152-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011158bam/ERR011158-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011156bam/ERR011156-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011160bam/ERR011160-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011130bam/ERR011130-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011142bam/ERR011142-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011153bam/ERR011153-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011092bam/ERR011092-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011114bam/ERR011114-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011107bam/ERR011107-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011166bam/ERR011166-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011155bam/ERR011155-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011168bam/ERR011168-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011169bam/ERR011169-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011167bam/ERR011167-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011170bam/ERR011170-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011175bam/ERR011175-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011171bam/ERR011171-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011172bam/ERR011172-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011173bam/ERR011173-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011174bam/ERR011174-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011176bam/ERR011176-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011177bam/ERR011177-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011178bam/ERR011178-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011180bam/ERR011180-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011179bam/ERR011179-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011181bam/ERR011181-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011182bam/ERR011182-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011184bam/ERR011184-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011183bam/ERR011183-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011190bam/ERR011190-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011189bam/ERR011189-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011186bam/ERR011186-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011193bam/ERR011193-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011187bam/ERR011187-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011194bam/ERR011194-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011191bam/ERR011191-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011192bam/ERR011192-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011195bam/ERR011195-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011185bam/ERR011185-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011188bam/ERR011188-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011196bam/ERR011196-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011198bam/ERR011198-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011201bam/ERR011201-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011199bam/ERR011199-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011200bam/ERR011200-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011197bam/ERR011197-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011205bam/ERR011205-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011203bam/ERR011203-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011206bam/ERR011206-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011204bam/ERR011204-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011202bam/ERR011202-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011209bam/ERR011209-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011210bam/ERR011210-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011211bam/ERR011211-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011212bam/ERR011212-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011216bam/ERR011216-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011213bam/ERR011213-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011207bam/ERR011207-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011214bam/ERR011214-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011219bam/ERR011219-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011224bam/ERR011224-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011221bam/ERR011221-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011226bam/ERR011226-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011227bam/ERR011227-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011220bam/ERR011220-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011223bam/ERR011223-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011228bam/ERR011228-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011222bam/ERR011222-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011218bam/ERR011218-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011230bam/ERR011230-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011225bam/ERR011225-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011229bam/ERR011229-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011217bam/ERR011217-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011215bam/ERR011215-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011208bam/ERR011208-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011233bam/ERR011233-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011237bam/ERR011237-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011234bam/ERR011234-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011232bam/ERR011232-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011231bam/ERR011231-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011236bam/ERR011236-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011235bam/ERR011235-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011238bam/ERR011238-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011239bam/ERR011239-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011240bam/ERR011240-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011242bam/ERR011242-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011246bam/ERR011246-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011243bam/ERR011243-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011248bam/ERR011248-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011241bam/ERR011241-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011247bam/ERR011247-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011245bam/ERR011245-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011249bam/ERR011249-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011250bam/ERR011250-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011251bam/ERR011251-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011244bam/ERR011244-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011261bam/ERR011261-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011258bam/ERR011258-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011253bam/ERR011253-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011256bam/ERR011256-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011262bam/ERR011262-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011260bam/ERR011260-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011257bam/ERR011257-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011255bam/ERR011255-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011259bam/ERR011259-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011254bam/ERR011254-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011252bam/ERR011252-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011267bam/ERR011267-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011263bam/ERR011263-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011266bam/ERR011266-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011268bam/ERR011268-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011272bam/ERR011272-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011264bam/ERR011264-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011265bam/ERR011265-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011270bam/ERR011270-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011269bam/ERR011269-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011271bam/ERR011271-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011280bam/ERR011280-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011281bam/ERR011281-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011274bam/ERR011274-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011282bam/ERR011282-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011278bam/ERR011278-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011277bam/ERR011277-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011279bam/ERR011279-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011285bam/ERR011285-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011284bam/ERR011284-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011283bam/ERR011283-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011273bam/ERR011273-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011275bam/ERR011275-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011276bam/ERR011276-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011286bam/ERR011286-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011287bam/ERR011287-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011292bam/ERR011292-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011291bam/ERR011291-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011288bam/ERR011288-genome_sorted.bam                                       
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011293bam/ERR011293-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011289bam/ERR011289-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011290bam/ERR011290-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011296bam/ERR011296-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011295bam/ERR011295-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011294bam/ERR011294-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011300bam/ERR011300-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011303bam/ERR011303-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011302bam/ERR011302-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011299bam/ERR011299-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011298bam/ERR011298-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011304bam/ERR011304-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011306bam/ERR011306-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011297bam/ERR011297-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011307bam/ERR011307-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011301bam/ERR011301-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011305bam/ERR011305-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011309bam/ERR011309-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011308bam/ERR011308-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011310bam/ERR011310-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011311bam/ERR011311-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011312bam/ERR011312-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011313bam/ERR011313-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011314bam/ERR011314-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011316bam/ERR011316-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011315bam/ERR011315-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011317bam/ERR011317-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011324bam/ERR011324-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011318bam/ERR011318-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011319bam/ERR011319-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011320bam/ERR011320-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011327bam/ERR011327-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011322bam/ERR011322-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011326bam/ERR011326-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011328bam/ERR011328-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011323bam/ERR011323-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011325bam/ERR011325-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011321bam/ERR011321-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011330bam/ERR011330-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011333bam/ERR011333-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011337bam/ERR011337-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011334bam/ERR011334-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011332bam/ERR011332-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011331bam/ERR011331-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011329bam/ERR011329-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011338bam/ERR011338-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011336bam/ERR011336-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011335bam/ERR011335-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011339bam/ERR011339-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011340bam/ERR011340-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011345bam/ERR011345-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011342bam/ERR011342-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011343bam/ERR011343-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011347bam/ERR011347-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011341bam/ERR011341-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011346bam/ERR011346-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011344bam/ERR011344-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011349bam/ERR011349-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011348bam/ERR011348-genome_sorted.bam
 /global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR011350bam/ERR011350-genome_sorted.bam
"""
#/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR033896_2/ERR033896_2-genome_sorted.bam
#/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR033897_1/ERR033897_1-genome_sorted.bam
#/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR033898_2/ERR033898_2-genome_sorted.bam
#/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR033896_1/ERR033896_1-genome_sorted.bam
#/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR033897_2/ERR033897_2-genome_sorted.bam
#/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/ERR033898_1/ERR033898_1-genome_sorted.bam "
},
binning_cmd_sup_nt_megan : {
"CAMI_high" : "cami.fa ",
"mockup_microbial" : "/global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/1030534.scaffolds.fasta ",
"Albertsen2013" : "/global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/assembly.fasta ",
"MetaHIT_human_microbiome" : "/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/kmernator-FR-342326.edique02.d-Ray-358992.edique02-k51-Ray-358992.edique02.d-Contigs.fasta "
}
}
for b in binning_cmd_sup:
    binning_cmd_params[b] = binning_cmd_params[binning_cmd_sup_nt_megan].copy()
for b in binning_cmd_unsup:
    binning_cmd_params[b] = binning_cmd_params[binning_cmd_unsup_metabat].copy()

binning_cmd_params[binning_cmd_sup_nt_megan]["mockup_microbial"] +=  "/global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/1030534.scaffolds.fasta.nt_bbdedupe_bbmasked_formatted.LOCALIZED.out "
binning_cmd_params[binning_cmd_sup_nt16s_megan]["mockup_microbial"] +=  "/global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/1030534.scaffolds.fasta.loose16S23S5SrRNAblastnreads.LOCALIZED.out "
binning_cmd_params[binning_cmd_sup_microb_megan]["mockup_microbial"] +=  "/global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/1030534.scaffolds.fasta.microbial.LOCALIZED.out "
binning_cmd_params[binning_cmd_sup_nt_megan]["Albertsen2013"] +=  "/global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/assembly.fasta.nt_bbdedupe_bbmasked_formatted.LOCALIZED.out "
binning_cmd_params[binning_cmd_sup_nt16s_megan]["Albertsen2013"] +=  "/global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/assembly.fasta.loose16S23S5SrRNAblastnreads.LOCALIZED.out "
binning_cmd_params[binning_cmd_sup_microb_megan]["Albertsen2013"] +=  "/global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/assembly.fasta.microbial.LOCALIZED.out "
binning_cmd_params[binning_cmd_sup_nt_megan]["MetaHIT_human_microbiome"] +=  "/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/kmernator-FR-342326.edique02.d-Ray-358992.edique02-k51-Ray-358992.edique02.d-Contigs.fasta.nt_bbdedupe_bbmasked_formatted.LOCALIZED.out "
binning_cmd_params[binning_cmd_sup_nt16s_megan]["MetaHIT_human_microbiome"] +=  "/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/kmernator-FR-342326.edique02.d-Ray-358992.edique02-k51-Ray-358992.edique02.d-Contigs.fasta.loose16S23S5SrRNAblastnreads.LOCALIZED.out "
binning_cmd_params[binning_cmd_sup_microb_megan]["MetaHIT_human_microbiome"] +=  "/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/kmernator-FR-342326.edique02.d-Ray-358992.edique02-k51-Ray-358992.edique02.d-Contigs.fasta.microbial.LOCALIZED.out "

'''
###The references of each project as fasta files are placed in an array.
gold_standard_bins = ["/global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/CONFIG.txt",
                "/global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/CONFIG.txt",
                "/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/CONFIG.txt"]
'''

###The output directory locations for each testrun. These dirs will be created under the outputdir specified by the user
testdata_outdir_goldbins = {
    "CAMI_high" : "/global/projectb/scratch/andreopo/binning/metabat/bin/cami_high/CONFIG.txt",
    "mockup_microbial" : "/global/projectb/scratch/andreopo/binning/metabat/bin/multiBAMS/CONFIG.txt",
    "Albertsen2013" : "/global/projectb/scratch/andreopo/binning/metabat/bin/Albertsen2013/CONFIG.txt",
    "MetaHIT_human_microbiome" : "/global/projectb/scratch/andreopo/binning/metabat/bin/MetaHIT_human_microbiome/CONFIG.txt"
}


#For all evals and all clusters + bins combos save all results
results = {}

#For all evals and all clusters + bins combos save all result_summaries
results_summaries = {}


'''
create a directory named after the supervised method used; if directory exists then rename it with a timestamp.
runs one supervised taxonomic classification binning command
returns: fastas, where the first fasta is the unclassfied contigs
'''
def run_sup_binning_test_datasets(binning_cmd, output_path, binning_cmd_params, log=None):
    #step 1: run binning_cmd_sup*
    #step 2: parse the output to get the list of contig names assigned to a genome/reference with high confidence
    #step 3: use findcontig script to print contigs to a separate file for each genome/reference and unclassified contigs also to a separate file.
    output_dir = output_path
    ###TODO cd output_dir, since binning_cmd will write output fasta files there
    if not os.path.isdir(output_dir):
        print "Cannot find %s, creating as new" % output_dir
        os.makedirs(output_dir)    
        os.chmod(output_dir, 0775)
    os.chdir(output_dir)
    clusters = []
    ###Run binning_cmd on test_data
    cmd = binning_cmd_sup.get(binning_cmd) + " " + binning_cmd_params ###.get(binning_cmd).get(test_data)
    if binning_cmd.find("_megan") > -1 or binning_cmd.find("_taxmapper") > -1:
        cmd = binning_cmd_sup.get(binning_cmd) + " " + binning_cmd_params + " " + output_dir
        
    std_out, std_err, exit_code = run_command(cmd, True, log)
    if exit_code != 0:
        print "CMD %s failed with OUTPUT %s ERROR %s" % (cmd, std_out, std_err)
        ###return
    ###Put clusters in BINS
    ###   make a new dir under the output_dir named BINS:
    BINS_dir = os.path.join(output_dir , "BINS" )
    ## Validate output dir. Create output_directory if it doesn't exist
    if not os.path.isdir(BINS_dir):
        print "Cannot find %s, creating as new" % BINS_dir
        os.makedirs(BINS_dir)    
        os.chmod(BINS_dir, 0775)
    ###   For all files under output directory suffixed with .fa
    ###   cp them under the BINS subdirectory
    BINS_src_files = os.path.join(output_dir , "*.keep.fasta" )
    ###cp BINS_src_files BINS_dir
    files_list  = open(os.path.join(output_dir , "files" ), 'w')
    i = 1
    for filename in glob.glob( BINS_src_files ):
        if filename.find("Not_assigned") < 0 and filename.find("No_hits") < 0:
            shutil.copy(filename, BINS_dir)
            shutil.move(os.path.join(BINS_dir, os.path.basename(filename)), os.path.join(BINS_dir, os.path.basename(filename + "." + str(i) + ".fa") ))
            files_list.write(os.path.basename(filename + "." + str(i) + ".fa") + "\n")
            i += 1
    files_list.close()
    '''
    BINS_file = os.path.join(output_dirs.get(binning_cmd1) , "BINS" )
    for r in result_sets1:
        BINS_fileh.write(r)
    for c in ccc_refs:
        clusters.append(c)
    clusters = BINS_dir
    i += 1
    '''
    result_bins =  glob.glob( os.path.join(BINS_dir, "*.fa") )
    unclass = glob.glob( os.path.join(output_dir, "unknown.fasta") )
    return (0, unclass, result_bins)
    
'''
runs all unsupervised and supervised commands. for each combination:
run supervised to get the classified contigs first; the first file is the unclassified contigs.
then run unsupervised on the unclassified contigs.
put the results in a special directory named after both unsupervised and classified method; if directory exists then rename it with a timestamp.
'''
def run_hybrid_binning_test_datasets(binning_sup, binning_unsup):
    pass


'''
create a directory named after the unsupervised method used; if directory exists then rename it with a timestamp.
'''
def run_unsup_binning_test_datasets(binning_cmd, output_path, binning_cmd_params, log=None):
    '''
    result_sets = []
    i = 0
    for test_data in test_datasets.get(binning_cmd):
        output_dir = output_dirs.get(binning_cmd)[i]
    '''
    output_dir = output_path
    ###TODO cd output_dir, since binning_cmd will write output fasta files there
    if not os.path.isdir(output_dir):
        print "Cannot find %s, creating as new" % output_dir
        os.makedirs(output_dir)    
        os.chmod(output_dir, 0775)
    os.chdir(output_dir)
    clusters = []
    ###Run binning_cmd on test_data
    cmd = binning_cmd_unsup.get(binning_cmd) + " " + binning_cmd_params ###.get(binning_cmd).get(test_data)
    std_out, std_err, exit_code = run_command(cmd, True, log)
    if exit_code != 0:
        print "CMD %s failed with OUT %s ERR %s" % (cmd, std_out, std_err)
        return
    ###Put clusters in BINS
    ###   make a new dir under the output_dir named BINS:
    BINS_dir = os.path.join(output_dir , "BINS" )
    ## Validate output dir. Create output_directory if it doesn't exist
    if not os.path.isdir(BINS_dir):
        print "Cannot find %s, creating as new" % BINS_dir
        os.makedirs(BINS_dir)    
        os.chmod(BINS_dir, 0775)
    ###   For all files under output directory suffixed with .fa
    ###   cp them under the BINS subdirectory
    BINS_src_files = os.path.join(output_dir , "*bins/*.fa" )
    ###cp BINS_src_files BINS_dir
    files_list  = open(os.path.join(output_dir , "files" ), 'w')
    for filename in glob.glob( BINS_src_files ):
        shutil.copy(filename, BINS_dir)
        files_list.write(os.path.basename(filename) + "\n")
    files_list.close()

    '''
    ###Run taxa finder
    for filename in glob.glob( BINS_src_files ):
        FASTA = filename
        cmd = "/usr/common/jgi/annotators/prodigal/2.50/bin/prodigal -a  %s.gene.faa -d  %s.gene.fasta  -i  %s  -o  %s.prodigal.out -p meta" % (FASTA, FASTA, FASTA, FASTA)
        std_out, std_err, exit_code = run_command(cmd, True, log)
        if exit_code != 0:
            print "CMD %s failed with OUT %s ERR %s" % (cmd, std_out, std_err)
            return
        cmd = "module load parallel; /projectb/sandbox/rqc/prod/pipelines/external_tools/Last/last-418/bin/parallel-fasta \"/projectb/sandbox/rqc/prod/pipelines/external_tools/Last/last-418/bin/lastal  -F 15 /projectb/sandbox/rqc/prod/pipelines/external_tools/Last/nr/nr\" <  %s.gene.fasta >  %s.gene.fasta.last.maf" % (FASTA, FASTA)
        std_out, std_err, exit_code = run_command(cmd, True, log)
        if exit_code != 0:
            print "CMD %s failed with OUT %s ERR %s" % (cmd, std_out, std_err)
            return
        cmd = "/projectb/sandbox/rqc/prod/pipelines/external_tools/Last/last2tophits.pl  %s.gene.fasta.last.maf  >   %s.gene.fasta.last.maf.parsed.tophit" % (FASTA, FASTA)
        std_out, std_err, exit_code = run_command(cmd, True, log)
        if exit_code != 0:
            print "CMD %s failed with OUT %s ERR %s" % (cmd, std_out, std_err)
            return
        cmd = "/projectb/sandbox/rqc/prod/pipelines/external_tools/taxMapper/DEFAULT/taxMapper.pl  %s.gene.fasta.last.maf.parsed.tophit -tax_only >  %s.gene.fasta.last.maf.parsed.tophit.taxlist" % (FASTA, FASTA)
        std_out, std_err, exit_code = run_command(cmd, True, log)
        if exit_code != 0:
            print "CMD %s failed with OUT %s ERR %s" % (cmd, std_out, std_err)
            return
    '''

    '''
    BINS_file = os.path.join(output_dirs.get(binning_cmd1) , "BINS" )
    for r in result_sets1:
        BINS_fileh.write(r)
    for c in ccc_refs:
        clusters.append(c)
    clusters = BINS_dir
    i += 1
    '''
    result_bins = glob.glob( os.path.join(BINS_dir, "*.fa") )  ### = glob.glob( BINS_src_files )  ###
    print("in run_unsup_binning_test_datasets result_bins: %s" % (result_bins))
    
    return (0, result_bins)




'''
read the number of files in the directory
'''
def num_files(directory):
    if not os.path.isdir(directory):
        return -1
    return len([name for name in os.listdir(directory) if os.path.isfile(name)])


###Jeff evaluation script
#eval_cmd1 = "qsub -cwd -b yes -j yes -w e -r y -N qc -l ram.c=5G -t 1:31 /global/projectb/scratch/andreopo/binning/metabat/github_jfroula/Metabat/binQC_array.sh "
eval_cmd1 = "export SGE_TASK_ID=%s ;  /global/projectb/scratch/andreopo/binning/metabat/github_jfroula/Metabat/binQC_array.sh  "
###TODO add more evaluation methods
eval_cmd2 = "nmi"
eval_cmd3 = "jeff eval normalized by contig length"

def run_eval(binning_cmd, eval_cmd, bins, output_path2, test_data):
    '''
    i = 0
    for gold_bins in gold_standard_bins.get(binning_cmd):
    '''
    gold_bins = testdata_outdir_goldbins.get(test_data)
    output_dir = output_path2 ###output_dirs.get(binning_cmd)[i]
    ###TODO cd output_dir, since eval_cmd is dependent on it to find the files and BINS
    os.chdir(output_dir)
    ###gold_bins = gold_standard_bins.get(i)
    num_bins = len(bins) ###num_files(bins[i])
    if num_bins <= 0:
        print "num_bins is " + str(num_bins)
        return
    for binn in range(num_bins):
        binnn = binn + 1
        print str(binnn)
        ###setenv to analyze the bin number
        eval_cmd2 = eval_cmd % ( binnn )
        cmd = eval_cmd2 + gold_bins
        print cmd
        std_out, std_err, exit_code = run_command(cmd, True)
        if exit_code != 0:
            print "CMD %s failed with OUT %s ERR %s" % (cmd, std_out, std_err)
            return
    ###Extract all evaluation results for this clusters + bins and save in results
    cmd = "find "+ output_dir +" -type f -name results.log  |xargs tail -n 100 | egrep 'Best reference|completeness:|precision:|==>' > " + os.path.join(output_dir, "BestRefs.Compl.Prec.txt")
    std_out, std_err, exit_code = run_command(cmd, True)
    if exit_code != 0:
        print "CMD %s failed with OUT %s ERR %s" % (cmd, std_out, std_err)
        exit(1)
    
    return
    ###Compute an average quality metric and save in results_summaries
    ###i += 1 
    ###TODO Compute an average quality metric over all clusters + bins combinations and print it out


'''
def compare_metabat_clustering_to_gold_standard_bins():
    result_sets1 = run_unsup_binning_test_datasets( binning_cmd1 )
    ###Write result_sets1 to BINS files and add abs path to BINS array.
    
    ###BINS = []
    ###BINS_file = os.path.join(output_dirs.get(binning_cmd1) , "BINS" )
    ###BINS.append( BINS_file )
    
    ###Read each dir in result_sets1 and list all *.fa files then put them in a files file under output_dir
    
    run_eval(binning_cmd1, eval_cmd1, result_sets1)
    ###run_eval(eval_cmd2, result_sets1)
    ###run_eval(eval_cmd3, result_sets1)
'''

'''
compare_older_clusterings_to_gold_standard_bins():
    result_sets2 = run_clustering_test_datasets( binning_cmd2 )
    result_sets3 = run_clustering_test_datasets( binning_cmd3 )
    result_sets4 = run_clustering_test_datasets( binning_cmd4 )
    
    run_eval(eval_cmd1, result_sets2)
    run_eval(eval_cmd2, result_sets2)
    run_eval(eval_cmd3, result_sets2)

    run_eval(eval_cmd1, result_sets3)
    run_eval(eval_cmd2, result_sets3)
    run_eval(eval_cmd3, result_sets3)

    run_eval(eval_cmd1, result_sets4)
    run_eval(eval_cmd2, result_sets4)
    run_eval(eval_cmd3, result_sets4)


#TODO
#compute_avgaccuracy_clustering_to_gold_standard_bins():
'''
if __name__ == "__main__":
    ###A run can be in test mode or live mode
    ###In case of test mode we run the tests defined above for metabat and possibly other software packages as well.
    
    desc = 'rqc_metagenome_binning'
    parser = argparse.ArgumentParser(description=desc)
    '''
    parser.add_argument("-f", "--fasta", dest="fasta", help="Fasta file name (absolute path) with contigs or scaffolds", required=False)
    parser.add_argument("-b", "--bam_list", dest="bam_list", help="List of bam file names (absolute paths) comma separated", required=False)
    '''
    parser.add_argument("-o", "--output-path", dest="output_path", help = "Output path to write to", required=True)
    parser.add_argument("-t", "--test-data", dest="test_data", help = "The test dataset to use as a name", required=True)
    parser.add_argument('-l', '--live', dest="live", action='store_true', default=False, help = "Run against a real case as opposed to evaluation on a test case")
    
    options = parser.parse_args()
    
    print "The following tools will be tested in order of unsupervised, supervised and combinations: " + str(holder_of_runs)
    
    '''
    fasta = None
    bam_list = []
    '''
    
    test_data = None
    output_path = None
    live = False
    
    '''
    if options.fasta:
        fasta = options.fasta
    
    if options.bam_list:
        bam_list = options.bam_list.split(',')
    '''
    
    if options.test_data:
        test_data = options.test_data
        
    if not test_data in testdata_outdir_goldbins.keys():
        print "Your test data must be one of: " + str(testdata_outdir_goldbins.keys())
        exit(1)

    ###Output directory
    ###TODO do not allow the output path to have "homes" in it.
    if options.output_path:
        output_path = options.output_path
    else:
        ###base = os.path.basename(fastq)
        output_path = os.getcwd()  ###os.path.join( os.getcwd() , base )
        
    # create output_directory if it doesn't exist
    if not os.path.isdir(output_path):
        os.makedirs(output_path)	



    timestamp = strftime("%m%d%Y-%H%M%S")



    if  options.live:
            live = True
        
    if not live:
        ###This is a test run
        ###iterate over all binning_cmds in sup, then in unsup, and do all hybrid combos in parallel.
        ###keep track in holder_of_runs of what has been run already.
        for bin_cmd in holder_of_runs:
            if bin_cmd in binning_cmd_unsup:
                output_path2 = os.path.join( output_path , os.path.join(bin_cmd, test_data) )
                
                if os.path.isdir(output_path2) and os.path.exists(output_path2):
                    shutil.move(output_path2, output_path2 + "_" + timestamp)
                    
                # create output_directory if it doesn't exist
                if not os.path.isdir(output_path2):
                    os.makedirs(output_path2)
        
                (exit_code, bins) = run_unsup_binning_test_datasets( bin_cmd, output_path2, binning_cmd_params.get(bin_cmd).get(test_data) )
                run_eval(bin_cmd, eval_cmd1, bins, output_path2, test_data )
            
        for bin_cmd in holder_of_runs:
            if bin_cmd in binning_cmd_sup:
                ###holder_of_runs.append(bin_cmd)
                output_path3 = os.path.join( output_path , os.path.join(bin_cmd, test_data) )
                
                if os.path.isdir(output_path3) and os.path.exists(output_path3):
                    shutil.move(output_path3, output_path3 + "_" + timestamp)
                    
                # create output_directory if it doesn't exist
                if not os.path.isdir(output_path3):
                    os.makedirs(output_path3)
        
                (exit_code, unclass, bins) = run_sup_binning_test_datasets( bin_cmd, output_path3, binning_cmd_params.get(bin_cmd).get(test_data) )
                run_eval(bin_cmd, eval_cmd1, bins, output_path3, test_data )
                
                ###Repeat a loop for all of the unsupervised binning commands taking as input just the unclassified contigs
                for bin_cmd_u in holder_of_runs:
                    if bin_cmd_u in binning_cmd_unsup:
                        ###holder_of_runs.append(bin_cmd + "U" + bin_cmd_u)
                        output_path4 = os.path.join( os.path.join( output_path , bin_cmd + "U" + bin_cmd_u ) , test_data )
                        
                        if os.path.isdir(output_path4) and os.path.exists(output_path4):
                            shutil.move(output_path4, output_path4 + "_" + timestamp)
                            
                        # create output_directory if it doesn't exist
                        if not os.path.isdir(output_path4):
                            os.makedirs(output_path4)
                
                        ind = binning_cmd_params.get(bin_cmd_u).get(test_data).find(".fasta ")
                        (exit_code, bins2) = run_unsup_binning_test_datasets( bin_cmd_u, output_path4, unclass[0] + "  " + binning_cmd_params.get(bin_cmd_u).get(test_data)[ind+6:])
                        run_eval(bin_cmd_u, eval_cmd1, bins2, output_path4, test_data)

    else:
        ###In case of live mode
        ###Input is set of bam files (from rnaseq-prok code run), as well as contigs or scaffolds for alignment
        exit(1)

    exit(0)
