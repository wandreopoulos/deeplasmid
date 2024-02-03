#!/usr/bin/python3
'''
# Read a fasta file.
# Each sequence consists of 2 rows: header line and the sequence itself.
#
# Objective: compute the statistcial features in each sequence.
#
# Output: a text (csv or space separated) array where rows are sequences/objects.
# Each column in the output is a signal/feature.
# The 1st column is the class (1 or 0) and the nth column is the name/header of a sequence.
# Columns 2 to n-1 are main features.
# Column n is the filename of the sequence, batch and whether it passed or failed assembly.
# A header line with column names is OPTIONAL.
# 
# The script also outputs a yaml file with all features.
#
#TODOs include in training model:
#- Pentamers - in progress
#- Genes or COGs - need plasmid-specific gene models
#- Repeat sequences? Code exists that do this, but training will take more time.
#Besides simple repeated sequences on a single strand, repeats should eventually include (the main question is computing this fast):
#Palindromes, inverted repeats that create hairpins or cruciforms, and mirror repeats in DNA.
#A palindromic sequence in DNA, such as the recognition site for the EcoRI restriction enzyme GAATTCGAATTC that reads GAATTCGAATTC from the opposite strand, is the same when reading from 5' to 3' on either strand of the DNA.
#Inverted repeats within a single strand of DNA can be converted into a hairpin, such as TGCGAT.....ATCGCA, whichh
 is double-stranded in the stem region. Inverted repeats within a double-stranded DNA sequence can also form a cruciform (double-hairpin) structure, e.g. consider the opposite strand: ACGCTA.....TAGCGT.
#DNA can also contain mirror repeats, such as TTAGCACCACGATT, which are not palindromic.
'''

__author__ = "Bill Andreopoulos"
__email__ = "wandreopoulos@lbl.gov"


import re
import os
import sys
import fileinput
import logging
import time
# from optparse import OptionParser
import argparse
from collections import defaultdict
from subprocess import Popen, call, PIPE, CalledProcessError
import subprocess
from multiprocessing import Process, Queue
import multiprocessing
### import datetime
from time import time, strftime
import yaml
import Constants
from threading import Timer
import csv
import operator
import collections
#https://riccomini.name/kill-subprocesses-linux-bash
#https://www.blog.pythonlibrary.org/2016/05/17/python-101-how-to-timeout-a-subprocess/

srcdir = os.path.dirname(__file__)

'''
sys.path.append(os.path.join(dir, './mpld3/mpld3'))

import mpld3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
'''
import numpy as np

###import khmer
import shutil

from JGI_Pipeline import JGI_Pipeline



DEBUG = 0




'''
Print all of the fatures to files for machine interpretation
'''
def print_read_features(output_path, id_run, gc_content, mingc, maxgc, longestHomopol, totalLongHomopol, longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, longest_revcompl_repeat_s2, occur_longest_revcompl_repeat_s2, len_sequence, longest_id_ecoli, longest_id_vector, energy, max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, \
                                    max_occur_pentamer, max_occur_pentamer_1hot, max_occur_same_pentamer, percent_occur_same_pentamer, \
            max_occur_hexamer, max_occur_same_hexamer, percent_occur_same_hexamer, \
            max_occur_heptamer, max_occur_same_heptamer, percent_occur_same_heptamer, \
            max_occur_octamer, max_occur_same_octamer, percent_occur_same_octamer, \
            max_occur_ninemer, max_occur_same_ninemer, percent_occur_same_ninemer, \
            max_occur_dekamer, max_occur_same_dekamer, percent_occur_same_dekamer, coghits, genesperMB, genecount, aalenavg, pfam_vector, plassketch, plasORIsketch, chromsketch, header):
    
    features_file = open(os.path.join(output_path, 'features.svn'), 'a')
    ###Space separated file (for possible input to excel :)
    features_file_excel = open(os.path.join(output_path, 'features.txt'), 'a')
    ###Yaml file
    if not os.path.exists(os.path.join(output_path, 'yml')):
        os.makedirs(os.path.join(output_path, 'yml'))
    header_start = header.split()[0].replace("/", "_")
    ###suff = header_start[-3:]
    ###if not os.path.exists(os.path.join(output_path, 'yml/' + suff)):
    ###    os.makedirs(os.path.join(output_path, 'yml/' + suff))
    features_file_yml = open(os.path.join(output_path, 'yml', header_start + '.yml'), 'a')
    yml_dict = {}
    
    line = ""
    line_excel = ""
    #print "filename " + filename
    ###TODO Change for each run
    ###Note the plasmid detection from a header isn't critical for this script.
    if header.find("plasm") > -1 or header.find("Plasm") > -1 or header.find("PLASM") > -1:
        line += "1"
        line_excel += "1"
    else:
        line += "0"
        line_excel += "0"
        
    line += " 1:" + str(gc_content)
    line_excel += " " + str(gc_content)
    yml_dict["gc_content"] = gc_content
    line += " 2:" + str(mingc)
    line_excel += " " + str(mingc)
    yml_dict["mingc"] = mingc
    line += " 3:" + str(maxgc)
    line_excel += " " + str(maxgc)
    yml_dict["maxgc"] = maxgc

    ###longestHomopolSum,totalLongHomopolSum
    feature_index = 4
    longestHomopolSum = 0
    totalLongHomopolSum = 0
    bases = ['A', 'C', 'G', 'T', 'N']
    for b in bases:
        if longestHomopol.has_key(b):
            line += " " + str(feature_index) + ":" + str(float(longestHomopol[b]) ) ###/ float(len_sequence) )
        feature_index += 1
        line_excel += " " + str(float(longestHomopol.get(b, 0)) ) ###/ float(len_sequence) )
        yml_dict[b+"_longestHomopol"] = float(longestHomopol.get(b, 0))
        longestHomopolSum += longestHomopol.get(b, 0)
        if totalLongHomopol.has_key(b):
            line += " " + str(feature_index) + ":" + str(float(totalLongHomopol[b]) ) ###/ float(len_sequence) )
        feature_index += 1
        line_excel += " " + str(float(totalLongHomopol.get(b, 0)) ) ###/ float(len_sequence) )
        yml_dict[b+"_totalLongHomopol"] = float(totalLongHomopol.get(b, 0))
        totalLongHomopolSum += totalLongHomopol.get(b, 0)
        
    if longestHomopolSum > 0:
        line += " " + str(feature_index) + ":" + str(float(longestHomopolSum) ) ###/ float(len_sequence) )
    if DEBUG == 1:
        print( "longestHomopolSum " + str(longestHomopolSum))
    feature_index += 1
    line_excel += " " + str(float(longestHomopolSum) ) ###/ float(len_sequence) )
    yml_dict["longestHomopolSum"] = float(longestHomopolSum)

    if totalLongHomopolSum > 0:
        line += " " + str(feature_index) + ":" + str(float(totalLongHomopolSum) ) ###/ float(len_sequence) )
    if DEBUG == 1:
        print( "totalLongHomopolSum " + str(totalLongHomopolSum))
    feature_index += 1
    line_excel += " " + str(float(totalLongHomopolSum) ) ###/ float(len_sequence) )
    yml_dict["totalLongHomopolSum"] = float(totalLongHomopolSum)

    ###The repeat featurs are not used atm (too time-consuming to compute)
    if len(longest_repeat) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_repeat))
    if DEBUG == 1:
        print( "longest_repeat " + longest_repeat)
    feature_index += 1
    line_excel += " " + str(len(longest_repeat))
    yml_dict["longest_repeat"] = len(longest_repeat)
    
    if occur_longest_repeat > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_repeat)
    if DEBUG == 1:
        print( "occur_longest_repeat " + str(occur_longest_repeat))
    feature_index += 1
    line_excel += " " + str(occur_longest_repeat)
    yml_dict["occur_longest_repeat"] = occur_longest_repeat

    if len(longest_rev_repeat) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_rev_repeat))
    if DEBUG == 1:
        print( "longest_rev_repeat " + longest_rev_repeat)
    feature_index += 1
    line_excel += " " + str(len(longest_rev_repeat))
    yml_dict["longest_rev_repeat"] = len(longest_rev_repeat)
    
    if occur_longest_rev_repeat > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_rev_repeat)
    if DEBUG == 1:
        print( "occur_longest_rev_repeat " + str(occur_longest_rev_repeat))
    feature_index += 1
    line_excel += " " + str(occur_longest_rev_repeat)
    yml_dict["occur_longest_rev_repeat"] = occur_longest_rev_repeat


    if len(longest_revcompl_repeat) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_revcompl_repeat))
    if DEBUG == 1:
        print( "longest_revcompl_repeat " + longest_revcompl_repeat)
    feature_index += 1
    line_excel += " " + str(len(longest_revcompl_repeat))
    yml_dict["longest_revcompl_repeat"] = len(longest_revcompl_repeat)

    
    if occur_longest_revcompl_repeat > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_revcompl_repeat)
    if DEBUG == 1:
        print( "occur_longest_revcompl_repeat " + str(occur_longest_revcompl_repeat))
    feature_index += 1
    line_excel += " " + str(occur_longest_revcompl_repeat)
    yml_dict["occur_longest_revcompl_repeat"] = occur_longest_revcompl_repeat



    if len(longest_revcompl_repeat_s2) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_revcompl_repeat_s2))
    if DEBUG == 1:
        print( "longest_revcompl_repeat_s2 " + longest_revcompl_repeat_s2)
    feature_index += 1
    line_excel += " " + str(len(longest_revcompl_repeat_s2))
    yml_dict["longest_revcompl_repeat_s2"] = len(longest_revcompl_repeat_s2)
    
    
    if occur_longest_revcompl_repeat_s2 > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_revcompl_repeat_s2)
    if DEBUG == 1:
        print( "occur_longest_revcompl_repeat_s2 " + str(occur_longest_revcompl_repeat_s2))
    feature_index += 1
    line_excel += " " + str(occur_longest_revcompl_repeat_s2)
    yml_dict["occur_longest_revcompl_repeat_s2"] = occur_longest_revcompl_repeat_s2


    line += " " + str(feature_index) + ":" + str(len_sequence)
    feature_index += 1
    line_excel += " " + str(len_sequence)
    yml_dict["len_sequence"] = len_sequence

    
    if longest_id_ecoli > 0:
        line += " " + str(feature_index) + ":" + str(longest_id_ecoli)
    if DEBUG == 1:
        print( "longest_id_ecoli " + str(longest_id_ecoli))
    feature_index += 1
    line_excel += " " + str(longest_id_ecoli)

    if longest_id_vector > 0:
        line += " " + str(feature_index) + ":" + str(longest_id_vector)
    if DEBUG == 1:
        print( "longest_id_vector " + str(longest_id_vector))
    feature_index += 1
    line_excel += " " + str(longest_id_vector)

    if energy < 0:
        line += " " + str(feature_index) + ":" + str(energy)
    if DEBUG == 1:
        print( "energy " + str(energy))
    feature_index += 1
    line_excel += " " + str(energy)

    ###max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, same for trimers and tetramers.
    if not max_occur_dimer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_dimer)
    if DEBUG == 1:
        print( "max_occur_dimer " + str(max_occur_dimer))
    feature_index += 1
    line_excel += " " + str(max_occur_dimer)
    yml_dict["max_occur_dimer"] = str(max_occur_dimer)

    if max_occur_same_dimer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_dimer)
    if DEBUG == 1:
        print( "max_occur_same_dimer " + str(max_occur_same_dimer))
    feature_index += 1
    line_excel += " " + str(max_occur_same_dimer)
    yml_dict["max_occur_same_dimer"] = str(max_occur_same_dimer)

    if percent_occur_same_dimer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_dimer)
    if DEBUG == 1:
        print( "percent_occur_same_dimer " + str(percent_occur_same_dimer))
    feature_index += 1
    line_excel += " " + str(percent_occur_same_dimer)
    yml_dict["percent_occur_same_dimer"] = str(percent_occur_same_dimer)

    if not max_occur_trimer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_trimer)
    if DEBUG == 1:
        print( "max_occur_trimer " + str(max_occur_trimer))
    feature_index += 1
    line_excel += " " + str(max_occur_trimer)
    yml_dict["max_occur_trimer"] = str(max_occur_trimer)

    if max_occur_same_trimer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_trimer)
    if DEBUG == 1:
        print( "max_occur_same_trimer " + str(max_occur_same_trimer))
    feature_index += 1
    line_excel += " " + str(max_occur_same_trimer)
    yml_dict["max_occur_same_trimer"] = str(max_occur_same_trimer)

    if percent_occur_same_trimer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_trimer)
    if DEBUG == 1:
        print( "percent_occur_same_trimer " + str(percent_occur_same_trimer))
    feature_index += 1
    line_excel += " " + str(percent_occur_same_trimer)
    yml_dict["percent_occur_same_trimer"] = str(percent_occur_same_trimer)

    if not max_occur_tetramer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_tetramer)
    if DEBUG == 1:
        print( "max_occur_tetramer " + str(max_occur_tetramer))
    feature_index += 1
    line_excel += " " + str(max_occur_tetramer)
    yml_dict["max_occur_tetramer"] = str(max_occur_tetramer)

    if max_occur_same_tetramer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_tetramer)
    if DEBUG == 1:
        print( "max_occur_same_tetramer " + str(max_occur_same_tetramer))
    feature_index += 1
    line_excel += " " + str(max_occur_same_tetramer)
    yml_dict["max_occur_same_tetramer"] = str(max_occur_same_tetramer)

    if percent_occur_same_tetramer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_tetramer)
    if DEBUG == 1:
        print( "percent_occur_same_tetramer " + str(percent_occur_same_tetramer))
    feature_index += 1
    line_excel += " " + str(percent_occur_same_tetramer)
    yml_dict["percent_occur_same_tetramer"] = str(percent_occur_same_tetramer)


    if not max_occur_pentamer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_pentamer)
    if DEBUG == 1:
        print( "max_occur_pentamer " + str(max_occur_pentamer))
    feature_index += 1
    line_excel += " " + str(max_occur_pentamer)
    yml_dict["max_occur_pentamer"] = str(max_occur_pentamer)

    if not max_occur_pentamer_1hot is None:
        line += " " + str(feature_index) + ":" + str(max_occur_pentamer_1hot)
    if DEBUG == 1:
        print( "max_occur_pentamer_1hot " + str(max_occur_pentamer_1hot))
    feature_index += 1
    line_excel += " " + str(max_occur_pentamer_1hot)
    yml_dict["max_occur_pentamer_1hot"] = str(max_occur_pentamer_1hot)


    if max_occur_same_pentamer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_pentamer)
    if DEBUG == 1:
        print( "max_occur_same_pentamer " + str(max_occur_same_pentamer))
    feature_index += 1
    line_excel += " " + str(max_occur_same_pentamer)
    yml_dict["max_occur_same_pentamer"] = str(max_occur_same_pentamer)

    if percent_occur_same_pentamer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_pentamer)
    if DEBUG == 1:
        print( "percent_occur_same_pentamer " + str(percent_occur_same_pentamer))
    feature_index += 1
    line_excel += " " + str(percent_occur_same_pentamer)
    yml_dict["percent_occur_same_pentamer"] = str(percent_occur_same_pentamer)


    if not max_occur_hexamer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_hexamer)
    if DEBUG == 1:
        print( "max_occur_hexamer " + str(max_occur_hexamer))
    feature_index += 1
    line_excel += " " + str(max_occur_hexamer)
    yml_dict["max_occur_hexamer"] = str(max_occur_hexamer)

    if max_occur_same_hexamer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_hexamer)
    if DEBUG == 1:
        print( "max_occur_same_hexamer " + str(max_occur_same_hexamer))
    feature_index += 1
    line_excel += " " + str(max_occur_same_hexamer)
    yml_dict["max_occur_same_hexamer"] = str(max_occur_same_hexamer)

    if percent_occur_same_hexamer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_hexamer)
    if DEBUG == 1:
        print( "percent_occur_same_hexamer " + str(percent_occur_same_hexamer))
    feature_index += 1
    line_excel += " " + str(percent_occur_same_hexamer)
    yml_dict["percent_occur_same_hexamer"] = str(percent_occur_same_hexamer)


    if not max_occur_heptamer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_heptamer)
    if DEBUG == 1:
        print( "max_occur_heptamer " + str(max_occur_heptamer))
    feature_index += 1
    line_excel += " " + str(max_occur_heptamer)
    yml_dict["max_occur_heptamer"] = str(max_occur_heptamer)

    if max_occur_same_heptamer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_heptamer)
    if DEBUG == 1:
        print( "max_occur_same_heptamer " + str(max_occur_same_heptamer))
    feature_index += 1
    line_excel += " " + str(max_occur_same_heptamer)
    yml_dict["max_occur_same_heptamer"] = str(max_occur_same_heptamer)

    if percent_occur_same_heptamer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_heptamer)
    if DEBUG == 1:
        print( "percent_occur_same_heptamer " + str(percent_occur_same_heptamer))
    feature_index += 1
    line_excel += " " + str(percent_occur_same_heptamer)
    yml_dict["percent_occur_same_heptamer"] = str(percent_occur_same_heptamer)

    #if USE_PRODIGAL: # It was deciced to print 0s for these fields if PRODIGAL is not installed or not used.
    '''
    for i in COGs:
        line += " " + str(feature_index) + ":" + str(coghits.get(i))
        if DEBUG == 1:
            print "COG " + i + " " + str(coghits.get(i))
        feature_index += 1
        line_excel += " " + str(coghits.get(i))
        yml_dict["COG_"+i] = str(coghits.get(i))
    '''
    line += " " + str(feature_index) + ":" + str(genesperMB)
    if DEBUG == 1:
        print( "genesperMB " + str(genesperMB))
    feature_index += 1
    line_excel += " " + str(genesperMB)
    yml_dict["genesperMB"] = float(genesperMB)
    
    line += " " + str(feature_index) + ":" + str(genecount)
    if DEBUG == 1:
        print( "genecount " + str(genecount))
    feature_index += 1
    line_excel += " " + str(genecount)
    yml_dict["genecount"] = float(genecount)

    line += " " + str(feature_index) + ":" + str(aalenavg)
    if DEBUG == 1:
        print( "aalenavg " + str(aalenavg))
    feature_index += 1
    line_excel += " " + str(aalenavg)
    yml_dict["aalenavg"] = float(aalenavg)

    line += " " + str(feature_index) + ":" + str(pfam_vector)
    if DEBUG == 1:
        print( "pfam_vector " + str(pfam_vector))
    feature_index += 1
    line_excel += " " + str(pfam_vector)
    yml_dict["pfam_vector"] = str(pfam_vector)

    #if USE_PROT_SKETCH: # It was deciced to print 0s for these fields if sketch is not installed or not used.
    #if plassketch is not None:
    line += " " + str(feature_index) + ":" + str(plassketch)
    if DEBUG == 1:
        print( "plassketch " + str(plassketch))
    feature_index += 1
    line_excel += " " + str(plassketch)
    yml_dict["plassketch"] = str(plassketch)

    #if plasORIsketch is not None:
    line += " " + str(feature_index) + ":" + str(plasORIsketch)
    if DEBUG == 1:
        print( "plasORIsketch " + str(plasORIsketch))
    feature_index += 1
    line_excel += " " + str(plasORIsketch)
    yml_dict["plasORIsketch"] = str(plasORIsketch)

    #if chromsketch is not None:
    line += " " + str(feature_index) + ":" + str(chromsketch)
    if DEBUG == 1:
        print( "chromsketch " + str(chromsketch))
    feature_index += 1
    line_excel += " " + str(chromsketch)
    yml_dict["chromsketch"] = str(chromsketch)


    '''
    if not max_occur_octamer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_octamer)
    if DEBUG == 1:
        print "max_occur_octamer " + str(max_occur_octamer)
    feature_index += 1
    line_excel += " " + str(max_occur_octamer)

    if max_occur_same_octamer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_octamer)
    if DEBUG == 1:
        print "max_occur_same_octamer " + str(max_occur_same_octamer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_octamer)

    if percent_occur_same_octamer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_octamer)
    if DEBUG == 1:
        print "percent_occur_same_octamer " + str(percent_occur_same_octamer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_octamer)


    if not max_occur_ninemer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_ninemer)
    if DEBUG == 1:
        print "max_occur_ninemer " + str(max_occur_ninemer)
    feature_index += 1
    line_excel += " " + str(max_occur_ninemer)

    if max_occur_same_ninemer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_ninemer)
    if DEBUG == 1:
        print "max_occur_same_ninemer " + str(max_occur_same_ninemer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_ninemer)

    if percent_occur_same_ninemer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_ninemer)
    if DEBUG == 1:
        print "percent_occur_same_ninemer " + str(percent_occur_same_ninemer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_ninemer)


    if not max_occur_dekamer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_dekamer)
    if DEBUG == 1:
        print "max_occur_dekamer " + str(max_occur_dekamer)
    feature_index += 1
    line_excel += " " + str(max_occur_dekamer)

    if max_occur_same_dekamer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_dekamer)
    if DEBUG == 1:
        print "max_occur_same_dekamer " + str(max_occur_same_dekamer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_dekamer)

    if percent_occur_same_dekamer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_dekamer)
    if DEBUG == 1:
        print "percent_occur_same_dekamer " + str(percent_occur_same_dekamer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_dekamer)
    '''

    line += " # sequence " + str(header) + " from file " + str(id_run)
    if DEBUG == 1:
        print( "line " + line)
    line_excel += " " + str(id_run) + " >" + str(header)
    yml_dict["id_run"] = id_run
    yml_dict["header"] = header
    yml_dict2 = {}
    yml_dict2["sequence"] = yml_dict
    
    print(line_excel)
    
    ###features_file is not used atm.
    features_file.write(line + "\n")
    features_file.flush()
    features_file.close()
    ###features_file_excel is the space separated file with all features!
    features_file_excel.write(line_excel + "\n")
    features_file_excel.flush()
    os.fsync(features_file_excel.fileno())
    features_file_excel.close()
    ###Write all features to the yml file!
    yaml.dump(yml_dict2, features_file_yml, default_flow_style=False)
    features_file_yml.flush()
    os.fsync(features_file_yml.fileno())
    features_file_yml.close()



#............................
# convert sequences to 1-hot 2-D arrays
#............................
def encode1hotBase(seqPad):
        hot2D=[]
        for y in seqPad:
            hot2D.append(Constants.oneHotBase[y])
        #''.join(map(str, a.flatten().tolist()))
        return ''.join(map(str, map(int, np.asarray(hot2D).flatten().tolist())))  #np.array(hot2D).astype(np.float32)


###############################
# RQC's runcommand function
###############################
def runCommand(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.strip(), stderr.strip(), process.returncode



'''
Compute total GC in seq
'''
def oneComputeGCtotal(seq):
    gc_content = 0
    total_bases = len(seq)
    idx = 0
    while idx < len(seq):
        seq_idx = seq[idx]
        if seq_idx == 'C' or seq_idx == 'G':
            gc_content += 1
        idx += 1
    gc_cont_perc = float(gc_content) / float(total_bases)
    return gc_cont_perc


'''
Compute min and max GC in any window of size windowSizeL in seq
'''
def twothreeComputeMinMaxGCinWindow(seq, windowSizeL, terminal=False):
    ###return 0,0,0,0
    #print str(seq)
    mingc = 10000.0
    maxgc = -10000.0
    mingcstart = -1
    maxgcstart = -1
    ###GCs = []
    idx = 0
    if terminal == True:
        subseq = seq[ : windowSizeL]
        gc_content = oneComputeGCtotal(subseq)
        if gc_content < mingc:
            mingc = gc_content
            mingcstart = 0
        if gc_content > maxgc:
            maxgc = gc_content
            maxgcstart = 0
            
        subseq = seq[ -windowSizeL : ]
        gc_content = oneComputeGCtotal(subseq)
        if gc_content < mingc:
            mingc = gc_content
            mingcstart = max( 0 , len(seq)-windowSizeL )
        if gc_content > maxgc:
            maxgc = gc_content
            maxgcstart = max( 0 , len(seq)-windowSizeL )
        return mingc, mingcstart, maxgc, maxgcstart

    while idx + windowSizeL <= len(seq):
        #print "idx " + str(idx) + " windowSizeL " + str(windowSizeL)
        subseq = seq[idx : idx + windowSizeL]
        gc_content = oneComputeGCtotal(subseq)
        if gc_content < mingc:
            mingc = gc_content
            mingcstart = idx
        if gc_content > maxgc:
            maxgc = gc_content
            maxgcstart = idx
        idx += 1

    #print str(GCs)
    return mingc, mingcstart, maxgc, maxgcstart


'''
Find the most freq homopolymers
'''
def fourFindLongestHomopolymer(seq):
    longestHomopol = {}
    totalLongHomopol = {}
    startLongestHomopol = {}
    min_length = 6
    seq = re.sub('[^ACGT]', 'N', seq)
    hits = re.findall(r'(([A-Z])\2\2+)', seq)
    for hit in hits:
        hit_len = len(hit[0])
        if longestHomopol.has_key(hit[1]):
            if hit_len > longestHomopol[hit[1]]:
                longestHomopol[hit[1]] = hit_len
                startLongestHomopol[hit[1]] = seq.find(hit[0])
            if hit_len >= min_length:
                totalLongHomopol[hit[1]] += hit_len
        else:
            longestHomopol[hit[1]] = hit_len
            startLongestHomopol[hit[1]] = seq.find(hit[0])
            if hit_len >= min_length:
                totalLongHomopol[hit[1]] = hit_len
            else:
                totalLongHomopol[hit[1]] = 0
        
    return longestHomopol, totalLongHomopol, startLongestHomopol




'''
  This function finds the most frequent 
  pentamer in a sequence. It only find the pentamer and converts it to 1hot encoding. No frequencies returned.
  Note computing the pentamer isn't too efficient compared to other features. Example:
RUNTIME oneComputeGCtotal: 0.00158214569092
RUNTIME twothreeComputeMinMaxGCinWindow100: 0.158836126328
RUNTIME fourFindLongestHomopolymer: 0.00112986564636
RUNTIME fivesixFindMers: 46.0283839703
  
RUNTIME oneComputeGCtotal: 0.00306391716003
RUNTIME twothreeComputeMinMaxGCinWindow100: 2.72661781311
RUNTIME fourFindLongestHomopolymer: 0.0174589157104
RUNTIME fivesixFindMers: 37.1501250267

This used to be computed with khmer (see an older commit)

  Parameters
  ----------
  seq : string
  seqin : string
  penalty_value : unused int

  Returns
  -----------
  max_occur_pentamer : most freq 5mer
  max_occur_pentamer_1hot : 1hot encoding of most freq 5mer
'''
def fivesixFindPentamer(seq, seqin, penalty_value):
    max_occur_pentamer = ""
    max_occur_same_pentamer = 0
    pos_occur_same_pentamer = 0

    cmd = [ os.path.join(srcdir, "run_pentamer.sh"), str(seqin)] #"shifter", "--image=bryce911/bbtools", "commonkmers.sh", "in=" + str(seqin), "out=stdout", "k=5", "display=3", "count"]
    proc1 = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = proc1.communicate()
    #std_out = subprocess.check_output(cmd)
    #subprocess.check_call(cmd)
    print("PENTAMER stdout: " + std_out)
    vals = ["AAAAA", 1]
    vals = std_out.split("\t")[1].split("=")
    print("PENTAMER vals: " + str(vals))
    max_occur_pentamer = vals[0]
    max_occur_same_pentamer = vals[1]

    percent_occur_same_pentamer = float(pos_occur_same_pentamer) / float(len(seq))

    if DEBUG == 1:
        print("max_occur_pentamer " + max_occur_pentamer + " max_occur_same_pentamer " + str(max_occur_same_pentamer) + " percent_occur_same_pentamer " + str(percent_occur_same_pentamer))

    #Convert max_occur_pentamer to one-hot encoding
    max_occur_pentamer_1hot = encode1hotBase(max_occur_pentamer)

    #proc1.kill()
    #proc1.terminate()
    #proc1.wait()

    return len(seq), max_occur_pentamer, max_occur_pentamer_1hot, max_occur_same_pentamer, percent_occur_same_pentamer, penalty_value




COGs = ["COG0009"]

'''
COGs not used atm because I only had chrom-specific COGs, not plasmid-specific COGs.
Also the COG computation wasn't that fast. Alevy said he'd send plasmid-specific COGs once they're ready.
Note use of COGs should improve accuracy of ML a lot.
'''
def prodigal(sequence, seqin, penalty_value):
        ###Run taxa finder
        ###for filename in glob.glob( BINS_src_files ):
        FASTA = seqin
        '''
        cmd = "/usr/common/jgi/annotators/prodigal/2.50/bin/prodigal -a  %s.gene.faa -d  %s.gene.fasta  -i  %s  -o  %s.prodigal.out -p meta" % (FASTA, FASTA, FASTA, FASTA)
        std_out, std_err, exit_code = run_command(cmd, True, log)
        if exit_code != 0:
            print "CMD %s failed with OUT %s ERR %s" % (cmd, std_out, std_err)
            return
        '''
        ###/usr/common/jgi/annotators/prodigal/2.50/bin/prodigal
        ###cmd = ["shifter", "--image=registry.services.nersc.gov/jgi/prodigal:latest", "prodigal", "-a", FASTA + ".gene.faa", "-d", FASTA + ".gene.fasta", "-i", FASTA, "-o", FASTA + ".prodigal.out", "-p", "meta" ]
        cmd = [ os.path.join(srcdir, "run_prodigal.sh"), FASTA]
        subprocess.check_call(cmd)

        prodigalFile = open(FASTA + ".prodigal.out", "r")
        x = {}
        genecount = 0
        contiglen = len(sequence)
        for i in xrange(0,contiglen):
           x[i] = 0
        for line in prodigalFile:
                if line.find("CDS") > -1:
                      [ start , end ] = [ int(i) for i in re.findall("\d+", line) ]
                      for i in xrange(min(start, end)-1, max(start, end)):
                           x[i] = 1
                      genecount += 1

        ones = x.values().count(1)
        genesperMB=0
        if contiglen > 0:
            genesperMB=(float(ones)/float(contiglen))
        prodigalFile.close()

        countaa = 0
        countprot = 0
        for line in fileinput.input( FASTA + ".gene.faa" ): ###READS####READS: ###:
            line = line.strip()
            if line.startswith(">"):
                 countprot += 1
            else:
                 countaa += len(line)
        aalenavg = 0
        if countprot > 0:
            aalenavg = countaa / float(countprot)

        #now run hmmsearch
        print ("HMMSearch Parsing...")
        cmd = [ os.path.join(srcdir, "run_hmmsearch.sh"), FASTA]
        subprocess.check_call(cmd)


        tblout_pfam = FASTA + ".domtblout"


        feature_table = get_table_from_tblout(tblout_pfam)
        feature_table = [i.strip().split(' ', 1) for i in feature_table]

        with open(FASTA + '.feature_table.txt', 'w') as output:
            writer = csv.writer(output, lineterminator='\n')
            writer.writerows(feature_table)


        feature_table_names=[]
        feature_table_genes=[]
        for i in feature_table:
            feature_table_names.append(i[0])
            feature_table_genes.append(i[1])


        print ("build_genehit_vector...")
        k = build_genehit_vector(feature_table_genes)



        if os.path.exists(FASTA + ".gene.faa"):  os.remove(FASTA + ".gene.faa")
        if os.path.exists(FASTA + ".prodigal.out"):  os.remove(FASTA + ".gene.fasta")
        if os.path.exists(FASTA + ".prodigal.out"):  os.remove(FASTA + ".prodigal.out")
        if os.path.exists(FASTA + ".domtblout"):  os.remove(FASTA + ".domtblout")
        if os.path.exists(FASTA + ".feature_table.txt"):  os.remove(FASTA + ".feature_table.txt")
        if os.path.exists(FASTA + ".out_pfam"):  os.remove(FASTA + ".out_pfam")

        return genesperMB, genecount, aalenavg, k



'''
Evaluate if the sequence hits one of the chromosome-specific aa sequences.
To make a separate sketch database from an aa fasta file, run this:
sketch.sh in=x.faa out=x.sketch amino persequence
Then to compare, run:
comparesketch.sh in=contigs.fa translate ref=x.sketch persequence
'''
def run_chromsketch(sequence, seqin, penalty_value):
        FASTA = seqin
        cmd = [os.path.join(srcdir, "run_chromsketch.sh"), FASTA]
        ####"shifter", "--image=bryce911/bbtools",  "comparesketch.sh", "-Xmx100M",  "in="+FASTA , "translate",  "ref=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_clean/classifier/dl/asafl_plasmidPred/protein_all.faa.sketch", "persequence"]
        print("CHROMFINDER cmd: " + str(cmd))
        #note because of the necessity to ensure there are no zombie processes left behind, I did not use RQC's runCommand.
        #subprocess.check_call(cmd)
        proc1 = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        kill = lambda process: process.kill()
        my_timer = Timer(50, kill, [proc1])
        try:
            my_timer.start()
            std_out, std_err = proc1.communicate()
        finally:
            my_timer.cancel()

        #std_out, std_err = proc1.communicate()
        print("CHROMFINDER std_out: " + str(std_out))
        #std_out = subprocess.check_output(cmd)

        results=1
        if std_out.find("No hits") >= 0:
            results=0

        #proc1.kill()
        #proc1.terminate()
        #proc1.wait()

        return results



'''
Evaluate if the sequence hits one of the plasmid-specific aa sequences.
To make a separate sketch database from an aa fasta file, run this:
sketch.sh in=x.faa out=x.sketch amino persequence
Then to compare, run:
comparesketch.sh in=contigs.fa translate ref=x.sketch persequence
'''
def run_plassketch(sequence, seqin, penalty_value):
        FASTA = seqin
        cmd = [ os.path.join(srcdir, "run_plassketch.sh"), FASTA]
        ####"shifter", "--image=bryce911/bbtools",  "comparesketch.sh", "-Xmx100M",  "in="+FASTA , "translate",  "ref=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_clean/classifier/dl/asafl_plasmidPred/protein_all.faa.sketch", "persequence"]
        print("PLASMIDFINDER cmd: " + str(cmd))
        #note because of the necessity to ensure there are no zombie processes left behind, I did not use RQC's runCommand.
        #subprocess.check_call(cmd)
        proc1 = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        kill = lambda process: process.kill()
        my_timer = Timer(50, kill, [proc1])
        try:
            my_timer.start()
            std_out, std_err = proc1.communicate()
        finally:
            my_timer.cancel()

        #std_out, std_err = proc1.communicate()
        print("PLASMIDFINDER std_out: " + str(std_out))
        #std_out = subprocess.check_output(cmd)

        results=1
        if std_out.find("No hits") >= 0:
            results=0

        #proc1.kill()
        #proc1.terminate()
        #proc1.wait()

        return results


'''
Evaluate if the sequence hits one of the plasmid-specific aa sequences.
In order to also sketch a nucleotide (not amino acid) fasta file, do these command parameters look ok? I removed the "amino" part.
sketch.sh in=x.fasta out=x.sketch persequence
Then to compare, run:
comparesketch.sh in=contigs.fa ref=x.sketch persequence
'''
def run_plasORIsketch(sequence, seqin, penalty_value):
        FASTA = seqin
        cmd = [ os.path.join(srcdir, "run_plasORIsketch.sh"), FASTA]
        ####"shifter", "--image=bryce911/bbtools",  "comparesketch.sh", "-Xmx100M",  "in="+FASTA , "translate",  "ref=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml_clean/classifier/dl/asafl_plasmidPred/protein_all.faa.sketch", "persequence"]
        print("PLASMIDORIFINDER cmd: " + str(cmd))
        #note because of the necessity to ensure there are no zombie processes left behind, I did not use RQC's runCommand.
        #subprocess.check_call(cmd)
        proc1 = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        kill = lambda process: process.kill()
        my_timer = Timer(50, kill, [proc1])
        try:
            my_timer.start()
            std_out, std_err = proc1.communicate()
        finally:
            my_timer.cancel()

        #std_out, std_err = proc1.communicate()
        print("PLASMIDORIFINDER std_out: " + str(std_out))
        #std_out = subprocess.check_output(cmd)

        results=1
        if std_out.find("No hits") >= 0:
            results=0

        #proc1.kill()
        #proc1.terminate()
        #proc1.wait()

        return results





def get_table_from_tblout(tblout_pfam):
    with open(tblout_pfam, "r") as infile:
        tblout_pfam=infile.readlines()
    tblout_pfam = [i.split() for i in tblout_pfam[3:-10]]
    for i in tblout_pfam:
        i[13] = float(i[13])
    tblout_pfam.sort(key = operator.itemgetter(0, 13,17), reverse = True)
    top_genes={}
    for i in tblout_pfam:
        if i[0] not in top_genes:
            top_genes[i[0]] = [[i[3],float(i[13]),float(i[17]),float(i[18])]]
        else:
            for j in top_genes[i[0]]:
                start_i, end_i, start_j, end_j = float(i[17]), float(i[18]), float(j[2]), float(j[3])
                if not ((end_i <= start_j) or (start_i >= end_j)):
                    break
                else:
                    top_genes[i[0]].append([i[3],float(i[13]),start_i,end_i])
                    break
    contigs = collections.OrderedDict()
    for i in top_genes:
        name = i.rsplit("_", 1)[0]
        if name not in contigs:
            contigs[name] = []
            for j in top_genes[i]:
                contigs[name].append(j[0])
        else:
            for j in top_genes[i]:
                contigs[name].append(j[0])
    out = []
    for key, value in contigs.items():
        out+=[str(key) + " "  +  " ".join(value)]
    return out


def build_genehit_vector(input_list):
    tr=os.path.dirname(os.path.abspath(__file__)) + "/pfams_discr.txt"
    hmm_dict = []
    with open(tr, 'r') as infile:
        table=infile.readlines()
        hmm_dict = map(hash, [i.strip() for i in table])

    # Calculate probabilities for each element of input list
    out_list = [0]*1538
    gene_list = []
    count = 0
    if len(input_list) > 0:
        gene_list = map(hash, input_list[0].split())
    for i in hmm_dict:
        if i in gene_list:
            out_list[count] = 1
        count += 1

    return out_list





def window_size(arr):
    mn = min(arr)
    mx = max(arr)
    return mx - mn + 1

def get_window(arr, seq):
    mn = min(arr)
    mx = max(arr)
    return seq[mn:mx + 1]

#
# Create a timestamp with the format YYYYMMDDHHMMSS.
#
def create_timestamp():
    return strftime("%m%d%Y-%H%M%S")


def process_seq(seqin, sequence, header, penalty_value, output_path, id_run, run_blast = None):
                    #if DEBUG == 1:
                    #    print "\n\n_____________________________________\n  -------> Computing features for fasta header : " + header + "\n"
                    
                    #TRIMREADS.write(line + "\n" + sequence + "\n")
                    #if DEBUG == 1:
                    #    print "sequence %s" % ( sequence )
                    #    print "length of sequence %s" % ( len(sequence) )

                    holder_bad_gc = []
                    startpointingtime = time()
                    gc_content = oneComputeGCtotal(sequence);
                    runtime = str(time() - startpointingtime)
                    #print "RUNTIME oneComputeGCtotal: " + runtime
                    #if DEBUG == 1:
                    #    print "oneComputeGCtotal for entire sequence"

                    window = min(100, len(sequence));
                    startpointingtime = time()
                    mingc, mingcstart, maxgc, maxgcstart = twothreeComputeMinMaxGCinWindow(sequence, window);
                    runtime = str(time() - startpointingtime)
                    #print "RUNTIME twothreeComputeMinMaxGCinWindow100: " + runtime
                    #if DEBUG == 1:
                    #    print "twothreeComputeMinMaxGCinWindow for window " + str(window)

                    ##TODO report where does the homopolymer start (coordinates for every homopolymer) and length, same for GC rules
                    ##TODO some reporting print statements show pointers to vars instead of the values, check all show values.
                    ##TODO check parallelism of writing output and errors, can you write results from parallel processes in order.
                    startpointingtime = time()
                    longestHomopol = {}
                    totalLongHomopol = {}
                    startLongestHomopol = {}
                    longestHomopol, totalLongHomopol, startLongestHomopol = fourFindLongestHomopolymer(sequence);
                    runtime = str(time() - startpointingtime)
                    #print "RUNTIME fourFindLongestHomopolymer: " + runtime
                    bases = ['A', 'C', 'G', 'T', 'N']
                    #if DEBUG == 1:
                    #    print "fourFindLongestHomopolymer"

                    startpointingtime = time()
                    longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, longest_revcompl_repeat_s2, occur_longest_revcompl_repeat_s2, len_sequence, max_occur_same_dimer, percent_occur_same_dimer, max_occur_same_trimer, percent_occur_same_trimer, penalty_value = ('', 0, '', 0, '', 0, '', 0, len(sequence), 0, 0, 0, 0, 0)
                    ###fivesixFindRepeatsInverted(sequence, penalty_value);
                    ###longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, longest_revcompl_repeat_s2, occur_longest_revcompl_repeat_s2, len_sequence, penalty_value = fivesixFindRepeats(sequence, penalty_value);
                    runtime = str(time() - startpointingtime)
                    ###print "RUNTIME fivesixFindRepeats: " + runtime
                    #if DEBUG == 1:
                    #    print "fivesixFindRepeats"
                        
                    startpointingtime = time()
                    len_sequence, max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, \
                    max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, \
                    max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, \
                    max_occur_pentamer, max_occur_pentamer_1hot, max_occur_same_pentamer, percent_occur_same_pentamer, \
                    max_occur_hexamer, max_occur_same_hexamer, percent_occur_same_hexamer, \
                    max_occur_heptamer, max_occur_same_heptamer, percent_occur_same_heptamer, \
                    max_occur_octamer, max_occur_same_octamer, percent_occur_same_octamer, \
                    max_occur_ninemer, max_occur_same_ninemer, percent_occur_same_ninemer, \
                    max_occur_dekamer, max_occur_same_dekamer, percent_occur_same_dekamer = (len(sequence), '-', 0, 0, '-', 0, 0, '-', 0, 0, '-', '-', 0, 0, '-', 0, 0, '-', 0, 0, '-', 0, 0, '-', 0, 0, '-', 0, 0)
                    len_sequence, max_occur_pentamer, max_occur_pentamer_1hot, max_occur_same_pentamer, percent_occur_same_pentamer, penalty_value = fivesixFindPentamer(sequence, seqin, penalty_value)
                    runtime = str(time() - startpointingtime)
                    #print "RUNTIME fivesixFindMers: " + runtime
                    #if DEBUG == 1:
                    #    print "fivesixFindMers"

                    startpointingtime = time()
                    coghits = None
                    genesperMB = 0
                    genecount = 0
                    aalenavg = 0
                    pfam_vector = []
                    if Constants.USE_PRODIGAL:
                       genesperMB, genecount, aalenavg, pfam_vector = prodigal(sequence, seqin, penalty_value);
                    runtime = str(time() - startpointingtime)
                    ###print "RUNTIME cogs: " + runtime
                    #if DEBUG == 1:
                    #    print "cogs"


                    startpointingtime = time()
                    plassketch = 0
                    plasORIsketch = 0
                    chromsketch = 0
                    if Constants.USE_PROT_SKETCH:
                       plassketch = run_plassketch(sequence, seqin, penalty_value);
                       plasORIsketch = run_plasORIsketch(sequence, seqin, penalty_value);
                       chromsketch = run_chromsketch(sequence, seqin, penalty_value);
                    runtime = str(time() - startpointingtime)
                    ###print "RUNTIME cogs: " + runtime
                    #if DEBUG == 1:
                    #    print "cogs"
                        
                    longest_id_ecoli = 0
                    longest_id_vector = 0
                    energy = 0
                    '''
                    seq_in_filename = os.path.join(output_path, 'seq.fasta')
                    seq_in = open(seq_in_filename, 'w')
                    seq_in.write('>' + header + '\n')
                    seq_in.write(sequence)
                    seq_in.close()
                    longest_id_ecoli = 0
                    longest_id_vector = 0
                    if run_blast != None and run_blast == 'p':
                        longest_id_ecoli, longest_id_vector, penalty_value = runBlastPlasmid(seq_in_filename, penalty_value); ###(0, 0) ###seveneightLongestAlignmentsEcoliVector(filename + '_' + header.replace(">",":"), sequence, ecoli_target_fasta, vector_target_fasta);
                    if DEBUG == 1:
                        print "seveneightLongestAlignmentsEcoliVector"
                    ###TODO remove
                    energy = 0 ###tenFoldingDeltaG(sequence)
                    if DEBUG == 1:
                        print "tenFoldingDeltaG"
                    '''
                        
                        
                    print_read_features(output_path, id_run, gc_content, mingc, maxgc, longestHomopol, totalLongHomopol, longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, longest_revcompl_repeat_s2, occur_longest_revcompl_repeat_s2, len_sequence, longest_id_ecoli, longest_id_vector, energy, max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, \
                                                    max_occur_pentamer, max_occur_pentamer_1hot, max_occur_same_pentamer, percent_occur_same_pentamer, \
            max_occur_hexamer, max_occur_same_hexamer, percent_occur_same_hexamer, \
            max_occur_heptamer, max_occur_same_heptamer, percent_occur_same_heptamer, \
            max_occur_octamer, max_occur_same_octamer, percent_occur_same_octamer, \
            max_occur_ninemer, max_occur_same_ninemer, percent_occur_same_ninemer, \
            max_occur_dekamer, max_occur_same_dekamer, percent_occur_same_dekamer, coghits, genesperMB, genecount, aalenavg, pfam_vector, plassketch, plasORIsketch, chromsketch, header);
                    
                    return penalty_value;




def multiproc_pool(pipeline):
        #print "remoteCommand: " + pipeline.getCommand();
        try:
                #print("Starting: " + pipeline.getCommand() + " ...with args: " + pipeline.params[0] + "  " +  pipeline.paramvalues[0] + "  " +  pipeline.params[1] + "  " +  pipeline.paramvalues[1] + "  " +  pipeline.params[2] + "  " +  pipeline.paramvalues[2] + "  " +  pipeline.params[3] + "  " +  pipeline.paramvalues[3] )
                #TODO append all params to cmd
                cmd = [ pipeline.getCommand() ] ###, pipeline.params[0] , pipeline.paramvalues[0] , pipeline.params[1] , pipeline.paramvalues[1] , pipeline.params[2] , pipeline.paramvalues[2] , pipeline.params[3] , pipeline.paramvalues[3] ]
                for i in range(0, len(pipeline.params)):
                        cmd.append(pipeline.params[i])
                        if i < len(pipeline.paramvalues):
                            cmd.append(pipeline.paramvalues[i])
                        
                proc1 = subprocess.Popen(cmd, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                kill = lambda process: process.kill()
                my_timer = Timer(1000, kill, [proc1])
                try:
                    my_timer.start()
                    #read sequence from stdin. https://stackoverflow.com/questions/8475290/how-do-i-write-to-a-python-subprocess-stdin
                    std_out, std_err = proc1.communicate(input=pipeline.sequence.encode('utf-8'))
                finally:
                    my_timer.cancel()

                #subprocess.check_call(cmd)

        except subprocess.CalledProcessError as e:
                result = "Failed to run job (%s)." % (e.output) ###.message   ###"+str(pipeline.paramvalues[2])+"
                print( result )
                ###p.kill()
                # store the result
                result_queue.put(result)

                return False

        result = 'Job finished with success.' ###% exit_code  ###'+str(pipeline.paramvalues[2])+'
        print( result )
        result_queue.put(result)

        #s.deleteJobTemplate(jt)
        #print 'Finished one Job '+str(pipeline.paramvalues[2])+' , Cleaning up.....'

        #s.exit()
        return True;

        #s = JGI_SGE_Submit(self.work_queue, self.result_queue)
        #s.run_job_set()
        #time.sleep(2)
        #return x*x


'''
Note: The same file is used for processing both the input fasta (passed with -i) 
and the individual sequences (passed with -s for the sequence file and -he for the header).
'''


if __name__ == "__main__":
    ###penalty_values = {}
    all_seqs = 0
    penalty_value = 0
    datetime = create_timestamp()
    ###id_run = None
    id_run = "DATETIME" + datetime

    #input-path and output-path are the only parameters used in the top call
    #in the child call, output-path header and seqin (a non-existent filename with the seqnumber) are specified and the sequence itself is passed through stdin
    desc = 'fasta_features'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-i", "--input-path", dest="input_path", help = "Input path to write to", required=False)
    parser.add_argument("-o", "--output-path", dest="output_path", help = "Output path to write to", required=False)
    parser.add_argument("-he", "--header", dest="header", help = "The header to use", required=False)
    parser.add_argument("-s", "--seqin", dest="seqin", help = "The sequence to use", required=False)
    
    options = parser.parse_args()

    ###if (options.input_path == None) != (options.output_path == None):
    ###    print "Either specify both input and output path or neither"
    ###    exit(1)

    output_path = None
    input_path = []
    run_blast = None
    header = "none"
    sequence = ""


    ###Output directory
    if options.output_path:
        output_path = options.output_path
    else:
        ###base = os.path.basename(fastq)
        output_path = os.getcwd()  ###os.path.join( os.getcwd() , base )

    if options.header:
        header = options.header

    #we are processing either the fasta with all sequences (-i) or a single sequence (-s and -he)
    if options.seqin and options.input_path:
        print("Please specify either seqin(-s) or input-path(-i) but not both")
        exit(2)
    
    if options.seqin and not options.header:
        print("If using seqin(-s) you may want to specify a header(-he)")
        
    ###Assume we are not passing a valid seqfilein to get rid of filesystem. Instead read sequence from stdin via a pipe.
    if options.seqin:
        '''
        if not os.path.exists(options.seqin):
            print("The dataset path does not exist: " + options.seqin)
            exit(2)
        if not os.path.isfile(options.seqin):
            print("The dataset file does not exist, is missing or is not readable: " + options.seqin)
            exit(2)
        seqfilein = open(options.seqin, 'r')
        '''
        for line in stdin: ###seqfilein: ###fileinput.input( [ options.seqin ] ): ###READS####READS: ###:
            ###line = CONTIGS.readline()
            ####sequence = options.seqin
            if not line.startswith(">"):
                sequence += line.strip()
        ###Remove file options.seqin
        #seqfilein.close()

    ###Do not allow the output path to have "homes" in it.
    if output_path.find("homes") > -1:
        print("The output path can not be under homes. The output_path given: " + output_path)
        exit(2)

    # create output_directory if it doesn't exist
    if not os.path.isdir(output_path):
        print("Cannot find %s, creating as new" % output_path)
        os.makedirs(output_path)    
        os.chmod(output_path, 777)

    ###Input directory
    if options.input_path:
        input_path = options.input_path
        if not os.path.isabs(input_path):
            ###TODO Convert a relative path to abs path
            input_path = os.path.join(os.getcwd(), input_path)
        if not os.path.exists(input_path):
            print("The dataset path does not exist: " + input_path)
            exit(2)
        if not os.path.isfile(input_path):
            print("The dataset file does not exist, is missing or is not readable: " + input_path)
            exit(2)
        #TODO verify number of contigs in input_path is <10k. This is file is only for contigs, not reads!
        input_path = [ input_path ]
        
    
    colheaders = "id_run,gc_content,mingc,maxgc,AlongestHomopol,AtotalLongHomopol,ClongestHomopol,CtotalLongHomopol,GlongestHomopol,GtotalLongHomopol,TlongestHomopol,TtotalLongHomopol,longestHomopolSum,totalLongHomopolSum,longest_repeat,occur_longest_repeat,longest_rev_repeat,occur_longest_rev_repeat,longest_revcompl_repeat,occur_longest_revcompl_repeat,longest_revcompl_repeat_s2,occur_longest_revcompl_repeat_s2,len_sequence,longest_id_ecoli,longest_id_vector,energy,max_occur_dimer,max_occur_same_dimer,percent_occur_same_dimer,max_occur_trimer,max_occur_same_trimer,percent_occur_same_trimer,max_occur_tetramer,max_occur_same_tetramer,percent_occur_same_tetramer,max_occur_pentamer,max_occur_pentamer_1hot,max_occur_same_pentamer,percent_occur_same_pentamer,max_occur_hexamer,max_occur_same_hexamer,percent_occur_same_hexamer,max_occur_heptamer,max_occur_same_heptamer,percent_occur_same_heptamer,COGS,genesperMB,genecount,aalenavg,plassketch,plasORIsketch,chromsketch"
    if len(input_path) >0:
        features_file = open(os.path.join(output_path, 'features.svn'), 'w')
        features_file.write(colheaders + "\n")
        features_file.flush()
        features_file.close()
        features_file_excel = open(os.path.join(output_path, 'features.txt'), 'w')
        features_file_excel.write(colheaders + "\n")
        features_file_excel.flush()
        os.fsync(features_file_excel.fileno())
        features_file_excel.close()

    ###for filename in fasta_files:

    ###READS = open(os.path.join(out_subdir , filename), 'r')
    ###TRIMREADS = open(os.path.join(out_subdir , 'TRIMREADS_' + datetime + '_' + filename), 'w')

    startpointingtimetotal = time()
    jobs_to_run = []
    num_jobs = 0
    seq_count = 0

    # Initialize the main queue of pending jobs.
    # Grab a job from the queue and submit it.
    # run
    # load up work queue
    # result_queue : a queue to pass to workers to store the results
    result_queue = multiprocessing.Queue()
    #work_queue is not used, the jobs_to_run list is used instead.
    work_queue = multiprocessing.Queue()

    #while True:
    if len(input_path) > 0:
        for line in fileinput.input( files = ( input_path ) , encoding="utf-8"): ###READS####READS: ###:
            ###line = CONTIGS.readline()
            line = line.strip()

            ###if line == '':
            ###    break
                
            if line.startswith(">"):
                #if (sequence == ""):
                #    header = line
                ### =~ m/(>\S+)\_length\_\d+\_cov\_(\S+)/o) {
                #elif (sequence != ""):
                if (sequence != ""):
                    #if len(input_path) <1:
                    #    print "penalty_value = process_seq(sequence %s, header %s, penalty_value %s, output_path %s, id_run %s, run_blast %s)" % (sequence, header, penalty_value, output_path, id_run, run_blast)
                    #    penalty_value = process_seq(sequence, header, penalty_value, output_path, id_run, run_blast)
                    #else:
                    ###rqc_command_dir = os.getcwd() ###os.path.join(os.path.abspath(os.pardir) , "categorical" ) ###current_dir
                    pars = ['-o', '-he', '-s']
                    seq_filename = os.path.join(output_path, str(seq_count) + '.fa')
                    '''
                    seq_file = open(seq_filename, 'w')
                    seq_file.write(header + "\n")
                    seq_file.write(sequence)
                    seq_file.flush()
                    os.fsync(seq_file.fileno())
                    seq_file.close()
                    '''
                    seq_count += 1
                    paramvals = [ str(output_path), str(header).replace(">",""), str(seq_filename) ]
                    rqc_command = os.path.realpath(__file__)  ###os.path.join(os.path.dirname(__file__), "read_fastaB.py") ###os.path.realpath(__file__)  ###os.path.join(rqc_command_dir , 'read_fasta.py')
                    p = JGI_Pipeline(rqc_command, pars, paramvals, sequence, header)
                    work_queue.put(p)
                    jobs_to_run.append(p)
                    num_jobs += 1
                    if num_jobs % 100 == 0:
                        print("job number processed: %s " % (num_jobs))
                    #num_processes += 1

                sequence = "";
                header = line   ###.replace(" ", "_").replace(",", ".");
                if DEBUG == 1: print("header %s" % (header))
                all_seqs += 1
                penalty_value = 0
            elif len(line) > 2:
                line = line.upper()
                ###line = re.sub('[^ATCG]', '', line)  ###Removed because we want Ns to be counted in homopolymer function
                sequence += line
        ###rqc_command_dir = os.getcwd() ###os.path.join(os.path.abspath(os.pardir) , "categorical" ) ###current_dir
        pars = ['-o', '-he', '-s']
        seq_filename = os.path.join(output_path, str(seq_count) + '.fa')
        '''
        seq_file = open(seq_filename, 'w')
        seq_file.write(header + "\n")
        seq_file.write(sequence)
        seq_file.flush()
        os.fsync(seq_file.fileno())
        seq_file.close()
        '''
        seq_count += 1
        paramvals = [ str(output_path), str(header).replace(">",""), str(seq_filename) ]
        rqc_command = os.path.realpath(__file__)  ###os.path.join(os.path.dirname(__file__), "read_fastaB.py") ###os.path.realpath(__file__)  ###os.path.join(rqc_command_dir , 'read_fasta.py')
        p = JGI_Pipeline(rqc_command, pars, paramvals, sequence, header)
        work_queue.put(p)
        jobs_to_run.append(p)
        num_jobs += 1
        if num_jobs % 100 == 0:
            print("job number processed: %s " % (num_jobs))
        #num_processes += 1
    else:
        ###if DEBUG == 1: print "\n\nSEQUENCE %s\n\n" % sequence
        ###if len(input_path) <1:
        ###print "penalty_value = process_seq(options.seqin %s, sequence %s, header %s, penalty_value %s, output_path %s, id_run %s, run_blast %s)" % (options.seqin, sequence, header, penalty_value, output_path, id_run, run_blast)
        penalty_value = process_seq(options.seqin, sequence, header, penalty_value, output_path, id_run, run_blast)
        if os.path.exists(options.seqin):  os.remove(options.seqin)

    #
    #
    ###MULTI-processing part - this is run only from the input fasta process.
    #
    #
    if len(input_path) >0:
        print('cpu_count() = %d\n' % multiprocessing.cpu_count())
        #
        # Create pool
        #
        print('Creating pool with %d processes\n' % Constants.PROCESSES)
        pool = multiprocessing.Pool(Constants.PROCESSES)
        #print 'pool = %s' % pool
        #print
        #pool = mp.Pool()
        #for i in range(10):
        #    pool.apply_async(multiproc_pool, args = (i, ), callback = log_result)
        #pool.map_async(multiproc_pool, range(100000), 10, callback = log_result)
        #pool.close()
        #pool.join()
        N = len(jobs_to_run)
        #print 'def pow3(x): return x**3'
        t = time()
        print('\tlist1(pool.imap(multiproc_pool, jobs_to_run(%d), chunksize=%d)):\n\t\t%s' \
            ' seconds' % (N, N//(Constants.PROCESSES), time() - t))
        # This is the multiproc pool and jobs_to_run is the list with all jobs that need to be run.
        C = list(pool.imap(multiproc_pool, jobs_to_run, chunksize=max(2, N//(Constants.PROCESSES))))
        print('\tlist2(pool.imap(multiproc_pool, jobs_to_run(%d), chunksize=%d)):\n\t\t%s' \
            ' seconds' % (N, N//(Constants.PROCESSES), time() - t))
        #print(result_list)
        while not result_queue.empty():
            print("Status of a done job: " + result_queue.get())
            
            
    runtime = str(time() - startpointingtimetotal)
    print("TOTAL_RUNTIME: " + runtime)
    #runtimes.write(runtime + "\n")
    print("exiting........")
    ###sys.exit(0)
    os._exit(0)

