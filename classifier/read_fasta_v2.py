#!/usr/bin/env python

# Read a fasta file.
# Each sequence consists of 2 rows: header line and the sequence itself.
#
# Objective: correlate signals and find associations and do a clustering.
#
# Output: a text (csv or space separated) array where rows are sequences/objects.
# Each column in the output is a signal/feature.
# The 1st column is the name of each sequence.
# Columns 2 to n-1 are main features.
# Column n is the filename of the sequence, batch and whether it passed or failed assembly.
# A header line with column names is OPTIONAL.
#
# Need to replace 4 variables:
# out_subdir : the directory with your fasta file(s)
# fasta_files : the array of you fasta file names to process.
# cmd under seveneightLongestAlignmentsEcoliVector: the ssw command for Smith-Waterman alignment.
#     This step did not reveal much info in disambiguating pass from fail so it could also be skipped.
# cmd under tenFoldingDeltaG : the location to the ViennaRNA package.
#
#

import re
import os
import sys
import fileinput
import logging
import time
# from optparse import OptionParser
import argparse
from collections import defaultdict
from subprocess import Popen, call, PIPE
from multiprocessing import Process, Queue
### import datetime
from time import time, strftime

dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir, './mpld3/mpld3'))

import mpld3
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

import khmer
import shutil



def print_read_features(output_path, id_run, gc_content, mingc, maxgc, longestHomopol, totalLongHomopol, longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, len_sequence, longest_id_ecoli, longest_id_vector, energy, max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, header):
    
    features_file = open(os.path.join(output_path, 'features.svn'), 'a')
    features_file_excel = open(os.path.join(output_path, 'features.txt'), 'a')
    
    line = ""
    line_excel = ""
    #print "filename " + filename
    ###TODO Change for each run
    if id_run.find("perfect") > -1 or id_run.find("pass") > -1 or id_run.find("single") > -1:
        line += "1"
        line_excel += "1"
    ###TODO Change for each run
    elif id_run.find("fail") > -1:
        line += "-1"
        line_excel += "-1"
    else:
        line += "-1"
        line_excel += "-1"
        
    line += " 1:" + str(gc_content)
    line_excel += " " + str(gc_content)
    line += " 2:" + str(mingc)
    line_excel += " " + str(mingc)
    line += " 3:" + str(maxgc)
    line_excel += " " + str(maxgc)

    ###longestHomopolSum,totalLongHomopolSum
    feature_index = 4
    longestHomopolSum = 0
    totalLongHomopolSum = 0
    bases = ['A', 'C', 'G', 'T']
    for b in bases:
        if longestHomopol.has_key(b):
            line += " " + str(feature_index) + ":" + str(float(longestHomopol[b]) ) ###/ float(len_sequence) )
        feature_index += 1
        line_excel += " " + str(float(longestHomopol.get(b, 0)) ) ###/ float(len_sequence) )
        longestHomopolSum += longestHomopol.get(b, 0)
        if totalLongHomopol.has_key(b):
            line += " " + str(feature_index) + ":" + str(float(totalLongHomopol[b]) ) ###/ float(len_sequence) )
        feature_index += 1
        line_excel += " " + str(float(totalLongHomopol.get(b, 0)) ) ###/ float(len_sequence) )
        totalLongHomopolSum += totalLongHomopol.get(b, 0)
        
    if longestHomopolSum > 0:
        line += " " + str(feature_index) + ":" + str(float(longestHomopolSum) ) ###/ float(len_sequence) )
    if DEBUG == 1:
        print "longestHomopolSum " + str(longestHomopolSum)
    feature_index += 1
    line_excel += " " + str(float(longestHomopolSum) ) ###/ float(len_sequence) )

    if totalLongHomopolSum > 0:
        line += " " + str(feature_index) + ":" + str(float(totalLongHomopolSum) ) ###/ float(len_sequence) )
    if DEBUG == 1:
        print "totalLongHomopolSum " + str(totalLongHomopolSum)
    feature_index += 1
    line_excel += " " + str(float(totalLongHomopolSum) ) ###/ float(len_sequence) )


    if len(longest_repeat) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_repeat))
    if DEBUG == 1:
        print "longest_repeat " + longest_repeat
    feature_index += 1
    line_excel += " " + str(len(longest_repeat))
    
    if occur_longest_repeat > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_repeat)
    if DEBUG == 1:
        print "occur_longest_repeat " + str(occur_longest_repeat)
    feature_index += 1
    line_excel += " " + str(occur_longest_repeat)

    if len(longest_rev_repeat) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_rev_repeat))
    if DEBUG == 1:
        print "longest_rev_repeat " + longest_rev_repeat
    feature_index += 1
    line_excel += " " + str(len(longest_rev_repeat))
    
    if occur_longest_rev_repeat > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_rev_repeat)
    if DEBUG == 1:
        print "occur_longest_rev_repeat " + str(occur_longest_rev_repeat)
    feature_index += 1
    line_excel += " " + str(occur_longest_rev_repeat)


    if len(longest_revcompl_repeat) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_revcompl_repeat))
    if DEBUG == 1:
        print "longest_revcompl_repeat " + longest_revcompl_repeat
    feature_index += 1
    line_excel += " " + str(len(longest_revcompl_repeat))
    
    
    if occur_longest_revcompl_repeat > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_revcompl_repeat)
    if DEBUG == 1:
        print "occur_longest_revcompl_repeat " + str(occur_longest_revcompl_repeat)
    feature_index += 1
    line_excel += " " + str(occur_longest_revcompl_repeat)


    line += " " + str(feature_index) + ":" + str(len_sequence)
    feature_index += 1
    line_excel += " " + str(len_sequence)
    
    
    if longest_id_ecoli > 0:
        line += " " + str(feature_index) + ":" + str(longest_id_ecoli)
    if DEBUG == 1:
        print "longest_id_ecoli " + str(longest_id_ecoli)
    feature_index += 1
    line_excel += " " + str(longest_id_ecoli)

    if longest_id_vector > 0:
        line += " " + str(feature_index) + ":" + str(longest_id_vector)
    if DEBUG == 1:
        print "longest_id_vector " + str(longest_id_vector)
    feature_index += 1
    line_excel += " " + str(longest_id_vector)

    if energy < 0:
        line += " " + str(feature_index) + ":" + str(energy)
    if DEBUG == 1:
        print "energy " + str(energy)
    feature_index += 1
    line_excel += " " + str(energy)

    ###max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, same for trimers and tetramers.
    if not max_occur_dimer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_dimer)
    if DEBUG == 1:
        print "max_occur_dimer " + str(max_occur_dimer)
    feature_index += 1
    line_excel += " " + str(max_occur_dimer)

    if max_occur_same_dimer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_dimer)
    if DEBUG == 1:
        print "max_occur_same_dimer " + str(max_occur_same_dimer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_dimer)

    if percent_occur_same_dimer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_dimer)
    if DEBUG == 1:
        print "percent_occur_same_dimer " + str(percent_occur_same_dimer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_dimer)

    if not max_occur_trimer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_trimer)
    if DEBUG == 1:
        print "max_occur_trimer " + str(max_occur_trimer)
    feature_index += 1
    line_excel += " " + str(max_occur_trimer)

    if max_occur_same_trimer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_trimer)
    if DEBUG == 1:
        print "max_occur_same_trimer " + str(max_occur_same_trimer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_trimer)

    if percent_occur_same_trimer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_trimer)
    if DEBUG == 1:
        print "percent_occur_same_trimer " + str(percent_occur_same_trimer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_trimer)

    if not max_occur_tetramer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_tetramer)
    if DEBUG == 1:
        print "max_occur_tetramer " + str(max_occur_tetramer)
    feature_index += 1
    line_excel += " " + str(max_occur_tetramer)

    if max_occur_same_tetramer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_tetramer)
    if DEBUG == 1:
        print "max_occur_same_tetramer " + str(max_occur_same_tetramer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_tetramer)

    if percent_occur_same_tetramer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_tetramer)
    if DEBUG == 1:
        print "percent_occur_same_tetramer " + str(percent_occur_same_tetramer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_tetramer)

    line += " # sequence " + str(header) + " from file " + str(id_run)
    if DEBUG == 1:
        print "line " + line
    line_excel += " " + str(id_run) + " " + str(header)
    
    print line_excel
    
    features_file.write(line + "\n")
    features_file.flush()
    features_file.close()
    features_file_excel.write(line_excel + "\n")
    features_file_excel.flush()
    features_file_excel.close()
    




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

def twothreeComputeMinMaxGCinWindow(seq, windowSizeL, terminal=False):
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

def fourFindLongestHomopolymer(seq):
    longestHomopol = {}
    totalLongHomopol = {}
    startLongestHomopol = {}
    min_length = 6
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

def fivesixFindMers(seq, penalty_value):

    lmers_dimer = khmer.new_ktable(2)
    lmers_indexes_dimer = {}
    lmers_trimer = khmer.new_ktable(3)
    lmers_indexes_trimer = {}
    lmers_tetramer = khmer.new_ktable(4)
    lmers_indexes_tetramer = {}
    
    lmers_dimer.consume(seq)
    lmers_trimer.consume(seq)
    lmers_tetramer.consume(seq)

    '''
    for windowSizeL in [4]:
        for idx_start in range(0, windowSizeL):
            idx = idx_start
            #consec_identical_mers = ''
            count_consec_identical_mers = 0
            prev_mer = ''
            while idx + windowSizeL <= len(seq):
                subseq = seq[idx : idx + windowSizeL]
                if subseq == len(subseq) * subseq[0]: ###is a homopolymer:
                    if DEBUG == 1:
                        print "subseq is a homopolymer " + subseq
                    idx += windowSizeL
                    prev_mer = subseq
                    continue
                lmers  = lmers_dimer
                lmers_indexes = lmers_indexes_dimer
                if windowSizeL == 3:
                    lmers  = lmers_trimer
                    lmers_indexes = lmers_indexes_trimer
                if windowSizeL == 4:
                    lmers  = lmers_tetramer
                    lmers_indexes = lmers_indexes_tetramer
                lmers[subseq] = lmers.get(subseq, 0) + 1
                lmers_indexes[subseq] = lmers_indexes.get(subseq, []) + range(idx , idx + windowSizeL)
                ###print "subseq " + subseq + " prev_mer " + prev_mer
                
                if not subseq == prev_mer:
                    key = str(windowSizeL) + "repeat"
                    if count_consec_identical_mers > rules.get(key).get("X"):
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Sequence=prev_mer) + rules.get(key).get("explanation") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, "    (...) explanation: consec_identical_mers " + prev_mer + " count_consec_identical_mers " + str(count_consec_identical_mers)
                        if windowSizeL == 2:
                            if count_consec_identical_mers > max_occur_same_dimer:
                                max_occur_same_dimer = count_consec_identical_mers
                                pos_occur_same_dimer = count_consec_identical_mers * 2
                                max_occur_dimer = prev_mer
                        elif windowSizeL == 3:
                            if count_consec_identical_mers > max_occur_same_trimer:
                                max_occur_same_trimer = count_consec_identical_mers
                                pos_occur_same_trimer = count_consec_identical_mers * 3
                                max_occur_trimer = prev_mer
                        elif windowSizeL == 4:
                            if count_consec_identical_mers > max_occur_same_tetramer:
                                max_occur_same_tetramer = count_consec_identical_mers
                                pos_occur_same_tetramer = count_consec_identical_mers * 4
                                max_occur_tetramer = prev_mer
                    count_consec_identical_mers = 1
                    #consec_identical_mers = subseq
                else:
                    count_consec_identical_mers += 1
                    #consec_identical_mers = subseq
                    
                idx += windowSizeL
                prev_mer = subseq
    '''

    max_occur_dimer = ""
    max_occur_same_dimer = 0
    pos_occur_same_dimer = 0
    '''
    for l in lmers_dimer.keys():
        if lmers_dimer.get(l) > max_occur_same_dimer and lmers_dimer.get(l) > 1:
            max_occur_dimer = l
            max_occur_same_dimer = lmers_dimer.get(l)
            pos_occur_same_dimer = len(set(lmers_indexes_dimer.get(l)))
    '''
    # run through all entries. if they have nonzero presence, print.
    for i in range(0, lmers_dimer.n_entries()):
       n = lmers_dimer.get(i)
       if n > max_occur_same_dimer and n > 1:
          if DEBUG == 1:  print lmers_dimer.reverse_hash(i), "is present", n, "times."
          max_occur_dimer = lmers_dimer.reverse_hash(i)
          max_occur_same_dimer = n
          
    max_occur_trimer = ""
    max_occur_same_trimer = 0
    pos_occur_same_trimer = 0
    '''
    for l in lmers_trimer.keys():
        if lmers_trimer.get(l) > max_occur_same_trimer and lmers_trimer.get(l) > 1:
            max_occur_trimer = l
            max_occur_same_trimer = lmers_trimer.get(l)
            pos_occur_same_trimer = len(set(lmers_indexes_trimer.get(l)))
    '''
    # run through all entries. if they have nonzero presence, print.
    for i in range(0, lmers_trimer.n_entries()):
       n = lmers_trimer.get(i)
       if n > max_occur_same_trimer and n > 1:
          if DEBUG == 1:  print lmers_trimer.reverse_hash(i), "is present", n, "times."
          max_occur_trimer = lmers_trimer.reverse_hash(i)
          max_occur_same_trimer = n

    max_occur_tetramer = ""
    max_occur_same_tetramer = 0
    pos_occur_same_tetramer = 0
    '''
    for l in lmers_tetramer.keys():
        if lmers_tetramer.get(l) > max_occur_same_tetramer and lmers_tetramer.get(l) > 1:
            max_occur_tetramer = l
            max_occur_same_tetramer = lmers_tetramer.get(l)
            pos_occur_same_tetramer = len(set(lmers_indexes_tetramer.get(l)))
    '''
    # run through all entries. if they have nonzero presence, print.
    for i in range(0, lmers_tetramer.n_entries()):
       n = lmers_tetramer.get(i)
       if n > max_occur_same_tetramer and n > 1:
          if DEBUG == 1:  print lmers_tetramer.reverse_hash(i), "is present", n, "times."
          max_occur_tetramer = lmers_tetramer.reverse_hash(i)
          max_occur_same_tetramer = n

    percent_occur_same_dimer = float(pos_occur_same_dimer) / float(len(seq))
    percent_occur_same_trimer = float(pos_occur_same_trimer) / float(len(seq))
    percent_occur_same_tetramer = float(pos_occur_same_tetramer) / float(len(seq))

    if DEBUG == 1:
        print "max_occur_dimer " + max_occur_dimer + " max_occur_same_dimer " + str(max_occur_same_dimer) + " percent_occur_same_dimer " + str(percent_occur_same_dimer)
        print "max_occur_trimer " + max_occur_trimer + " max_occur_same_trimer " + str(max_occur_same_trimer) + " percent_occur_same_trimer " + str(percent_occur_same_trimer)
        print "max_occur_tetramer " + max_occur_tetramer + " max_occur_same_tetramer " + str(max_occur_same_tetramer) + " percent_occur_same_tetramer " + str(percent_occur_same_tetramer)

    return len(seq), max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, penalty_value


def fivesixFindRepeatsInverted(seq, penalty_value):
    ###Use http://tandem.bu.edu/trf/trf.html
    ###No, things that look for repeats in human genome are prob looking for long sequences
    ###whereas our entire sequence is shorter, in comparision to human
    ###if you hash all lenght l mers and count the number of entries with count>1,
    ###that's an easy way to measure repeats (exact repeats)
    lmers = {}
    lmers_indexes = {}
    min_l = 4
    max_l = 31
    for windowSizeL in range(max_l, min_l, -1):
        idx = 0
        while idx + windowSizeL <= len(seq):
            subseq = seq[idx : idx + windowSizeL]
            # TODO Check it does not intersect with a previous repeat pattern found? It does not matter because we only
            # care about the indexes of a mer in the seq and those are computed based on sets, so if there is overlap with previous lmers it is ok.
            #if len( list(set(range(idx , idx + windowSizeL)) & set(lmers_indexes.get(subseq, [])) ) ) == 0:
            lmers[subseq] = lmers.get(subseq, 0) + 1
            lmers_indexes[subseq] = lmers_indexes.get(subseq, []) + range(idx , idx + windowSizeL)
            idx += 1

    longest_lmer = ""
    for l in lmers.keys():
        if len(l) > len(longest_lmer) and lmers[l] > 1:
            longest_lmer = l

    max_occur_same_dimer = 0
    pos_occur_same_dimer = 0
    max_occur_dimer = ""
    max_occur_same_trimer = 0
    pos_occur_same_trimer = 0
    max_occur_trimer = ""
    max_occur_same_tetramer = 0
    pos_occur_same_tetramer = 0
    max_occur_tetramer = ""
    
    '''
    ###This logic for max dimers and trimers was redone below.
        if len(l) == 2 and lmers[l] > max_occur_same_dimer:
            max_occur_same_dimer = lmers[l]
            max_occur_dimer = l
            pos_occur_same_dimer = len(set(lmers_indexes[l]))
        if len(l) == 3 and lmers[l] > max_occur_same_trimer:
            max_occur_same_trimer = lmers[l]
            max_occur_trimer = l
            pos_occur_same_trimer = len(set(lmers_indexes[l]))
    '''
    ##TODO dimers and trimers should not include homopolymers.
    ##TODO generalize the code for dimers trimers and homopolymers with regexs.
    ###TODO regex python for dimer and trimer repeats, but not homopolymers. Exclude homopolymers. OK check below and continue if homopolymer.
    for windowSizeL in [2, 3, 4]:
        for idx_start in range(0, windowSizeL):
            idx = idx_start
            #consec_identical_mers = ''
            count_consec_identical_mers = 0
            prev_mer = ''
            while idx + windowSizeL <= len(seq):
                subseq = seq[idx : idx + windowSizeL]
                if subseq == len(subseq) * subseq[0]: ###is a homopolymer:
                    if DEBUG == 1:
                        print "subseq is a homopolymer " + subseq
                    idx += windowSizeL
                    prev_mer = subseq
                    continue
                ###print "subseq " + subseq + " prev_mer " + prev_mer
                if not subseq == prev_mer:
                    key = str(windowSizeL) + "repeat"
                    if count_consec_identical_mers > rules.get(key).get("X"):
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Sequence=prev_mer) + rules.get(key).get("explanation") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, "    (...) explanation: consec_identical_mers " + prev_mer + " count_consec_identical_mers " + str(count_consec_identical_mers)
                        if windowSizeL == 2:
                            if count_consec_identical_mers > max_occur_same_dimer:
                                max_occur_same_dimer = count_consec_identical_mers
                                pos_occur_same_dimer = count_consec_identical_mers * 2
                                max_occur_dimer = prev_mer
                        elif windowSizeL == 3:
                            if count_consec_identical_mers > max_occur_same_trimer:
                                max_occur_same_trimer = count_consec_identical_mers
                                pos_occur_same_trimer = count_consec_identical_mers * 3
                                max_occur_trimer = prev_mer
                        elif windowSizeL == 4:
                            if count_consec_identical_mers > max_occur_same_tetramer:
                                max_occur_same_tetramer = count_consec_identical_mers
                                pos_occur_same_tetramer = count_consec_identical_mers * 4
                                max_occur_tetramer = prev_mer
                    count_consec_identical_mers = 1
                    #consec_identical_mers = subseq
                else:
                    count_consec_identical_mers += 1
                    #consec_identical_mers = subseq
                idx += windowSizeL
                prev_mer = subseq


    percent_occur_same_dimer = float(pos_occur_same_dimer) / float(len(seq))
    percent_occur_same_trimer = float(pos_occur_same_trimer) / float(len(seq))
    percent_occur_same_tetramer = float(pos_occur_same_tetramer) / float(len(seq))

    if DEBUG == 1:
        print "max_occur_dimer " + max_occur_dimer + " max_occur_same_dimer " + str(max_occur_same_dimer) + " percent_occur_same_dimer " + str(percent_occur_same_dimer)
        print "max_occur_trimer " + max_occur_trimer + " max_occur_same_trimer " + str(max_occur_same_trimer) + " percent_occur_same_trimer " + str(percent_occur_same_trimer)
        print "max_occur_tetramer " + max_occur_tetramer + " max_occur_same_tetramer " + str(max_occur_same_tetramer) + " percent_occur_same_tetramer " + str(percent_occur_same_tetramer)
    
    #for i in lmers.keys():
    #    lmers[i] = 1
    
    rev_lmers = {}
    revcompl_lmers = {}
    union_bases_in_repeats = []
    #for windowSizeL in range(5, 61):
        #idx = 0
        #while idx + windowSizeL <= len(seq):
            #subseq = seq[idx : idx + windowSizeL]
    for subseq in lmers.keys():
        if len(subseq) > min_l:
            lmers_indexes_subseq = lmers_indexes.get(subseq)
            #TODO rules 1 only need 5mers on 30base terminals , compare first 30 bases to last 30 bases and see if they have repeats. Find overlapping kmers(5mers) in beginning and end.
            #3: check if any terminal direct repeats exist for subseq at both ends of the sequence.
            # (len(set( range(0, 25) ) & set(lmers_indexes_subseq)) > 0) and (len(set( range(len(seq)-25, len(seq)) ) & set(lmers_indexes_subseq)) > 0)
            key = "Terminal Repeats longer than 5 bp"
            if lmers.get(subseq) > 1 and len(subseq) >= rules.get(key).get("X") and (29 < len(seq)-30) and (len(set( range(0, 30) ) & set(lmers_indexes_subseq)) > min_l) and (len(set( range(len(seq)-30, len(seq)) ) & set(lmers_indexes_subseq)) > min_l):
                ###OLD check (0 in lmers_indexes_subseq or 1 in lmers_indexes_subseq) and ( (len(seq)-2) in lmers_indexes_subseq or (len(seq)-1) in lmers_indexes_subseq):
                penalty_value += float( rules.get(key).get("value") )
                if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(ViolatingSequence=subseq) + rules.get(key).get("explanation") )
                        
            rev_subseq = subseq[::-1]
            #TODO CHECK
            rev_subseq_indexes = set([m.start() for m in re.finditer(rev_subseq, seq)]) ###[seq.find(rev_subseq) , seq.rfind(rev_subseq) ]) ###[m.start() for m in re.finditer(rev_subseq, seq) ]
            #Remove -1 for not found case
            rev_subseq_indexes.discard(-1)
            ### print str(rev_subseq_indexes)
            if len(rev_subseq_indexes) > 1 and rev_subseq == subseq:
                rev_lmers[rev_subseq] = lmers.get(subseq)
            if len(rev_subseq_indexes) > 0 and rev_subseq != subseq:  ### seq.find(rev_subseq) > -1:
                rev_lmers[rev_subseq] = len(rev_subseq_indexes) + lmers.get(subseq) ### rev_lmers.get(rev_subseq, lmers.get(subseq) ) + 1
                #print subseq + " " + rev_subseq + " rev_lmers[rev_subseq] " + str(rev_lmers[rev_subseq]) + " lmers.get(subseq, 0) " + str(lmers.get(subseq, 0))
            ###Update rev_subseq_indexes to include all indexes
            rev_subseq_indexes_tmp = []
            for i in rev_subseq_indexes:
                rev_subseq_indexes_tmp = rev_subseq_indexes_tmp + range(i, i+len(rev_subseq))
            rev_subseq_indexes = rev_subseq_indexes_tmp
            
            revcompl_subseq = ""
            for i in range(len(rev_subseq)):
                if rev_subseq[i] == 'A':
                    revcompl_subseq += 'T'
                elif rev_subseq[i] == 'T':
                    revcompl_subseq += 'A'
                elif rev_subseq[i] == 'C':
                    revcompl_subseq += 'G'
                elif rev_subseq[i] == 'G':
                    revcompl_subseq += 'C'
            #TODO CHECK
            revcompl_subseq_indexes = set([m.start() for m in re.finditer(revcompl_subseq, seq)]) ###[seq.find(revcompl_subseq) , seq.rfind(revcompl_subseq) ]) ###[m.start() for m in re.finditer(revcompl_subseq, seq) ]
            revcompl_subseq_indexes.discard(-1)
            if len(revcompl_subseq_indexes) > 1 and revcompl_subseq == subseq:
                revcompl_lmers[revcompl_subseq] = lmers.get(subseq)
            if len(revcompl_subseq_indexes) > 0 and revcompl_subseq != subseq:  ### seq.find(revcompl_subseq) > -1:
                revcompl_lmers[revcompl_subseq] = len(revcompl_subseq_indexes) + lmers.get(subseq) ### revcompl_lmers.get(revcompl_subseq, lmers.get(subseq) ) + 1
                #print subseq + " " + revcompl_subseq + " revcompl_lmers[revcompl_subseq] " + str(revcompl_lmers[revcompl_subseq]) + " lmers.get(subseq, 0) " + str(lmers.get(subseq, 0))
            ###Update revcompl_subseq_indexes to include all indexes
            revcompl_subseq_indexes_tmp = []
            for i in revcompl_subseq_indexes:
                revcompl_subseq_indexes_tmp = revcompl_subseq_indexes_tmp + range(i, i+len(revcompl_subseq))
            revcompl_subseq_indexes = revcompl_subseq_indexes_tmp

            
            if revcompl_lmers.has_key(revcompl_subseq):
                ###These all involve hairpins, so we only consider the revcompl.
                #6: check if len(subseq) >= 10 and subseq and revcompl_subseq are in window of 40 bases and GC content >= 80%
                key = "Windowed Hairpin > 80% GC of Length > 10 bp"
                if len(subseq) >= rules.get(key).get("Y") and len(revcompl_subseq) >= rules.get(key).get("Y"):
                    if window_size(lmers_indexes.get(subseq) + revcompl_subseq_indexes) <= rules.get(key).get("X"):
                        if oneComputeGCtotal( get_window(lmers_indexes.get(subseq) + revcompl_subseq_indexes, seq) ) >= rules.get(key).get("Z"):
                            penalty_value += float( rules.get(key).get("value") )
                            if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Sequence=subseq, RevComplSequence=revcompl_subseq) + rules.get(key).get("explanation") )

                #5: check if len(subseq) > 19 and len(revcompl_subseq) > 19
                key = "Hairpin anywhere stems"
                if len(subseq) > rules.get(key).get("X") and len(revcompl_subseq) > rules.get(key).get("X"):
                    penalty_value += float( rules.get(key).get("value") )
                    if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Sequence=subseq, RevComplSequence=revcompl_subseq) + rules.get(key).get("explanation") )

                #4: check if len(subseq) > 15 and len(revcompl_subseq) > 15 and subseq and revcompl_subseq are in window of 100 bases
                key = "Hairpin close stems"
                if len(subseq) > rules.get(key).get("X") and len(revcompl_subseq) > rules.get(key).get("X"):
                    win_size = window_size(lmers_indexes.get(subseq) + revcompl_subseq_indexes)
                    if win_size <= rules.get(key).get("Y"):
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Sequence=subseq, ForwardLocations=win_size) + rules.get(key).get("explanation") )
    
            #idx += 1
            ###Check the last failure cases here
            #21: window of 70b with repeats involving at least 89% of the bases
            #TODO CHECK correctness
            key = "Repeat Percentage (total windowed repeats)"
            if len(subseq) >= rules.get(key).get("Z") and lmers.get(subseq) > 1:
                win_size = window_size(lmers_indexes.get(subseq) + revcompl_subseq_indexes + rev_subseq_indexes)
                if win_size >= rules.get(key).get("X"):
                    threshold_bases = int(rules.get(key).get("Y") * win_size )
                    bases = len( set( lmers_indexes.get(subseq) + revcompl_subseq_indexes + rev_subseq_indexes ) )
                    if bases >= threshold_bases:
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Segment=subseq, MinRepeatLength=len(subseq), Percentage=bases, WindowSize=win_size) + rules.get(key).get("explanation") )
                    
            #22: total bases in repeats are at least 40% of the sequence
            key = "Repeat Percentage > 35% (Individual repeats)"
            #TODO CHECK > or >= 8
            if len(subseq) > rules.get(key).get("X"):
                threshold_bases = int( rules.get(key).get("Y") * len(seq) )
                bases = len( set( lmers_indexes.get(subseq) + revcompl_subseq_indexes + rev_subseq_indexes ) )
                if bases >= threshold_bases:
                    penalty_value += float( rules.get(key).get("value") )
                    if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Segment=subseq, Percentage=bases, ThresholdPercentage=threshold_bases) + rules.get(key).get("explanation") )

            #23: sum up bases involved in repeats. Afterwards check if they account for at least 69% of sequence
            key = "Repeat Percentage (total repeats)"
            if len(subseq) >= rules.get(key).get("Y"):
                union_bases_in_repeats = union_bases_in_repeats + lmers_indexes.get(subseq) + revcompl_subseq_indexes + rev_subseq_indexes
                
    #23: sum up bases involved in repeats. Afterwards check if they account for at least 69% of sequence
    key = "Repeat Percentage (total repeats)"
    threshold_bases = int( rules.get(key).get("X") * len(seq) )
    bases = len( set( union_bases_in_repeats ) )
    if bases >= threshold_bases:
        penalty_value += float( rules.get(key).get("value") )
        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Percentage=bases, ThresholdPercentage=threshold_bases) + rules.get(key).get("explanation") )
        
    longest_rev_lmer = ""
    for l in rev_lmers.keys():
        if len(l) > len(longest_rev_lmer) and rev_lmers[l] > 1:
            longest_rev_lmer = l
        
    if not longest_rev_lmer == "":
        if DEBUG == 1:
            print "longest_rev_lmer " + longest_rev_lmer + " rev_lmers " + str(rev_lmers[longest_rev_lmer])
    
    longest_revcompl_lmer = ""
    for l in revcompl_lmers.keys():
        if len(l) > len(longest_revcompl_lmer) and revcompl_lmers[l] > 1:
            longest_revcompl_lmer = l

    if not longest_revcompl_lmer == "":
        if DEBUG == 1:
            print "longest_revcompl_lmer " + longest_revcompl_lmer + " revcompl_lmers " + str(revcompl_lmers[longest_revcompl_lmer])
    
    if not longest_lmer == "":
        if DEBUG == 1:
            print "longest_lmer " + longest_lmer + " lmers " + str(lmers[longest_lmer])
            print "penalty_value " + str(penalty_value)

    return longest_lmer, lmers.get(longest_lmer, 0), longest_rev_lmer, rev_lmers.get(longest_rev_lmer, 0), longest_revcompl_lmer, revcompl_lmers.get(longest_revcompl_lmer, 0), len(seq), max_occur_same_dimer, percent_occur_same_dimer, max_occur_same_trimer, percent_occur_same_trimer, penalty_value






def fivesixFindRepeats(seq, penalty_value):
    ###Use http://tandem.bu.edu/trf/trf.html
    ###No, things that look for repeats in human genome are prob looking for long sequences
    ###whereas our entire sequence is shorter, in comparision to human
    ###if you hash all lenght l mers and count the number of entries with count>1,
    ###that's an easy way to measure repeats (exact repeats)
    lmers = {}
    ###lmers_indexes = {}
    revcompl_lmers = {}
    min_l = 4
    max_l = 31
    for windowSizeL in range(max_l, min_l, -1):
        idx = 0
        while idx + windowSizeL <= len(seq):
            subseq = seq[idx : idx + windowSizeL]
            # TODO Check it does not intersect with a previous repeat pattern found? It does not matter because we only
            # care about the indexes of a mer in the seq and those are computed based on sets, so if there is overlap with previous lmers it is ok.
            #if len( list(set(range(idx , idx + windowSizeL)) & set(lmers_indexes.get(subseq, [])) ) ) == 0:
            lmers[subseq] = lmers.get(subseq, 0) + 1
            ###Find the reverse complement
            rev_subseq = subseq[::-1]
            revcompl_subseq = ""
            for i in range(len(rev_subseq)):
                if rev_subseq[i] == 'A':
                    revcompl_subseq += 'T'
                elif rev_subseq[i] == 'T':
                    revcompl_subseq += 'A'
                elif rev_subseq[i] == 'C':
                    revcompl_subseq += 'G'
                elif rev_subseq[i] == 'G':
                    revcompl_subseq += 'C'
            if lmers.has_key(revcompl_subseq):
                ###Store the subseq so we do not have to recompute the revcompl when searching for the longest revcompl
                revcompl_lmers[subseq] = revcompl_lmers.get(subseq, 0) + 1

            ###lmers_indexes[subseq] = lmers_indexes.get(subseq, []) + range(idx , idx + windowSizeL)
            idx += 1

    longest_lmer = ""
    longest_revcompl_lmer = ""
    most_freq_lmer = ""
    most_freq_lmer_count = 1
    for l in lmers.keys():
        if len(l) > len(longest_lmer) and lmers[l] > 1:
            longest_lmer = l
        if len(l) > len(longest_revcompl_lmer) and revcompl_lmers.get(l, 0) > 0:
            longest_revcompl_lmer = l
        if lmers[l] > most_freq_lmer_count:
            most_freq_lmer = l
            most_freq_lmer_count = lmers[l]

    if not longest_lmer == "":
        if DEBUG == 1:
            print "longest_lmer " + longest_lmer + " lmers " + str(lmers[longest_lmer])
            print "penalty_value " + str(penalty_value)

    return longest_lmer, lmers.get(longest_lmer, 0), most_freq_lmer, most_freq_lmer_count, longest_revcompl_lmer, revcompl_lmers.get(longest_revcompl_lmer, 0), len(seq), penalty_value




def runBlastMito(seq_in_filename, penalty_value):
    #print seq to seq.fasta
    #f = open('seq.fasta', 'w')
    #f.write(">SEQUENCE\n")
    #f.write(seq)
    #f.close()
    if os.path.exists(seq_in_filename + '.main.vs.refseq.mitochondrion.parsed.tophit'):
        shutil.move(seq_in_filename + '.main.vs.refseq.mitochondrion.parsed.tophit', seq_in_filename + '.main.vs.refseq.mitochondrion.parsed.tophit.BAK')

    ##############################
    #run run_mito.sh
    outfile = seq_in_filename + ".out" ###os.path.join("tmp", timestamp+'.out')
    cmd = ["/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/run_blast_mito.sh", seq_in_filename]
    if DEBUG == 1:
        print "cmd " + str(cmd)
    with open(outfile, 'w') as output:
                p = Popen(cmd, stdout = output, stderr = PIPE)
                std_out, std_err = p.communicate()
                #print "std_out " + std_out
                #print "std_err " + std_err
                #print "retcode " + str(p.returncode)

    #parse tax_list file
    if not os.path.exists(seq_in_filename + '.main.vs.refseq.mitochondrion.parsed.tophit'):
        return 0, 0, penalty_value
    
    RESULTS = open(seq_in_filename + '.main.vs.refseq.mitochondrion.parsed.tophit', 'r')
    seq_id = 0
    num_hits = 0
    while True:
                line = RESULTS.readline()
                if DEBUG == 1:
                    print "line " + line
                if line == '':
                    break
                if not line.startswith("#"):
                    vals = line.split("\t")
                    if float(vals[4]) > seq_id:
                        seq_id = float(vals[4])
                        num_hits = float(vals[3]) / min(float(vals[5]), float(vals[6]))
                    
    return seq_id, num_hits, penalty_value

def runBlastPlasmid(seq_in_filename, penalty_value):
    #print seq to seq.fasta
    #f = open('seq.fasta', 'w')
    #f.write(">SEQUENCE\n")
    #f.write(seq)
    #f.close()
    if os.path.exists(seq_in_filename + '.main.vs.refseq.plasmid.parsed.tophit'):
        shutil.move(seq_in_filename + '.main.vs.refseq.plasmid.parsed.tophit', seq_in_filename + '.main.vs.refseq.plasmid.parsed.tophit.BAK')

    ##############################
    #run run_mito.sh
    outfile = seq_in_filename + ".out" ###os.path.join("tmp", timestamp+'.out')
    cmd = ["/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/run_blast_plasmid.sh", seq_in_filename]
    if DEBUG == 1:
        print "cmd " + str(cmd)
    with open(outfile, 'w') as output:
                p = Popen(cmd, stdout = output, stderr = PIPE)
                std_out, std_err = p.communicate()
                #print "std_out " + std_out
                #print "std_err " + std_err
                #print "retcode " + str(p.returncode)

    #parse tax_list file
    if not os.path.exists(seq_in_filename + '.main.vs.refseq.plasmid.parsed.tophit'):
        return 0, 0, penalty_value
    
    RESULTS = open(seq_in_filename + '.main.vs.refseq.plasmid.parsed.tophit', 'r')
    seq_id = 0
    num_hits = 0
    while True:
                line = RESULTS.readline()
                if DEBUG == 1:
                    print "line " + line
                if line == '':
                    break
                if not line.startswith("#"):
                    vals = line.split("\t")
                    if float(vals[4]) > seq_id:
                        seq_id = float(vals[4])
                        num_hits = float(vals[3]) / min(float(vals[5]), float(vals[6]))
                    
    return seq_id, num_hits, penalty_value



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
    # 14 digits YYYYMMDDHHMMSS
    #year = datetime.datetime.now().year
    #month = datetime.datetime.now().month
    #day = datetime.datetime.now().day
    #hour = datetime.datetime.now().hour
    #minute = datetime.datetime.now().minute
    #second = datetime.datetime.now().second
    #if (month < 10):
    #    month = "0" + str(month);
    #if (day < 10):
    #    day = "0" + str(day);
    #if (hour < 10):
    #    hour = "0" + str(hour);
    #if (minute < 10):
    #    minute = "0" + str(minute);
    #if (second < 10):
    #    second = "0" + str(second);
    #res = str(year) + str(month) + str(day) + str(hour) + str(minute) + str(second);
    #return res;


def seveneightLongestAlignmentsEcoliVector(query_sequence_name, query_sequence, ecoli_target_fasta, vector_target_fasta):
    timestamp = create_timestamp()
    query_fasta = os.path.join("tmp", "query_seq.fasta")
    QUER_FA = open(query_fasta, 'w')
    QUER_FA.write(">" + query_sequence_name + "\n")
    QUER_FA.write(query_sequence + "\n")
    QUER_FA.close()

    longest_id_ecoli = 0
    longest_id_vector = 0

    ######################ECOLI BLAST######################
    #defaults = "-D 2 -a 8 -p 60 -e 1e-3 -J F -f T -I T -W 11";
    
    #ext = defaults.replace("-e 1e-", "E").replace("\"", "").replace(" ", "").replace("-", "")

    #output the .log and megablast raw file, which will be deleted at the end of the pipeline.
    #outfile = os.path.join("tmp", query_sequence_name + '_ecoli_MEGABLASTresults_' + timestamp + '.out') ###+ '_' + ext
    #cmd_str = "/jgi/tools/bin/megablast " + defaults + " -i " + query_fasta + " -d " + ecoli_target_fasta  ### + " > " + outfile + " 2> " + outfile + ".log"
    #print cmd_str
    #cmd = cmd_str.split(" ")
    #with open(outfile, 'w') as output:
    #    p = Popen(cmd, stdout = output, stderr = PIPE)
    #    std_out, std_err = p.communicate()


    ### This is a fast Smith-Waterman aligner, which also finds reverse complement read alignments: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
    #TODO Run Smith-Waterman on query (key_query_timestamp.fasta) and target (key_target_timestamp.fasta) and output results to file key_results_timestamp.out
    #############ECOLI Smith-Waterman#################
    outfile = os.path.join("tmp", query_sequence_name + '_ecoli_SSWresults_'+timestamp+'.out')
    cmd = ["/global/homes/a/andreopo/test_smith_waterman/Complete-Striped-Smith-Waterman-Library-master/src/ssw_test", "-cr", ecoli_target_fasta, query_fasta] ###, ">" , outfile]
    if DEBUG == 1:
        print "cmd " + str(cmd)
    with open(outfile, 'w') as output:
                p = Popen(cmd, stdout = output, stderr = PIPE)
                std_out, std_err = p.communicate()
                #print "std_out " + std_out
                #print "std_err " + std_err
                #print "retcode " + str(p.returncode)

    RESULTS = open(outfile, 'r')
    longest_id_ecoli = 0
    count_id = 0
    while True:
                line = RESULTS.readline()
                if DEBUG == 1:
                    print "line " + line
                if line == '':
                    break
                if line.startswith("target_name"):
                        # We reached a new target, so we save the longest hit.
                        if count_id and count_id > 0:
                            perc_id = (count_id + count_id) / ((target_end - target_begin + 1) + (query_end - query_begin + 1))
                            if perc_id >= 0.98 and count_id >= 18:
                                if count_id > longest_id_ecoli:
                                    longest_id_ecoli = count_id
                            if DEBUG == 1:
                                print "ECOLI PERC_ID " + str(perc_id) + " COUNT_ID " + str(count_id)
                            count_id = 0

                        key2 = line.strip().replace("target_name: ", "")
                        ###print "line " + line + " key2 " + key2
                        line = RESULTS.readline()
                        line = RESULTS.readline()
                        n = 0
                        if line.startswith("optimal_alignment_score"):
                            myarray = re.findall("(\d+)", line)
                            alignment = float(myarray[0])
                            target_begin = float(myarray[-4])
                            target_end = float(myarray[-3])
                            query_begin = float(myarray[-2])
                            query_end = float(myarray[-1])
                            #n = alignment / ((target_end - target_begin + 1) + (query_end - query_begin + 1))
                            #n = int(n*100)
                            #if n > longest_id_ecoli:
                            #    longest_id_ecoli = n
                            ###n = myarray[0]
                else:
                    count_id += line.count("|")

    if count_id and count_id > 0:
        perc_id = (count_id + count_id) / ((target_end - target_begin + 1) + (query_end - query_begin + 1))
        if perc_id >= 0.98 and count_id >= 18:
            if count_id > longest_id_ecoli:
                longest_id_ecoli = count_id
        if DEBUG == 1:
            print "ECOLI PERC_ID " + str(perc_id) + " COUNT_ID " + str(count_id)
        count_id = 0
    #############VECTOR#################
    outfile = os.path.join("tmp", query_sequence_name + '_vector_SSWresults_'+timestamp+'.out')
    cmd = ["/global/homes/a/andreopo/test_smith_waterman/Complete-Striped-Smith-Waterman-Library-master/src/ssw_test", "-cr", vector_target_fasta, query_fasta] ###, ">" , outfile]
    if DEBUG == 1:
        print "cmd " + str(cmd)
    with open(outfile, 'w') as output:
                p = Popen(cmd, stdout = output, stderr = PIPE)
                std_out, std_err = p.communicate()
                #print "std_out " + std_out
                #print "std_err " + std_err
                #print "retcode " + str(p.returncode)

    RESULTS = open(outfile, 'r')
    longest_id_vector = 0
    count_id = 0
    while True:
                line = RESULTS.readline()
                if DEBUG == 1:
                    print "line " + line
                if line == '':
                    break
                if line.startswith("target_name"):
                        # We reached a new target, so we save the longest hit.
                        if count_id and count_id > 0:
                            perc_id = (count_id + count_id) / ((target_end - target_begin + 1) + (query_end - query_begin + 1))
                            if perc_id >= 0.98 and count_id >= 18:
                                if count_id > longest_id_vector:
                                    longest_id_vector = count_id
                            if DEBUG == 1:
                                print "VECTOR PERC_ID " + str(perc_id) + " COUNT_ID " + str(count_id)
                            count_id = 0

                        key2 = line.strip().replace("target_name: ", "")
                        ###print "line " + line + " key2 " + key2
                        line = RESULTS.readline()
                        line = RESULTS.readline()
                        n = 0
                        if line.startswith("optimal_alignment_score"):
                            myarray = re.findall("(\d+)", line)
                            alignment = float(myarray[0])
                            target_begin = float(myarray[-4])
                            target_end = float(myarray[-3])
                            query_begin = float(myarray[-2])
                            query_end = float(myarray[-1])
                            #n = alignment / ((target_end - target_begin + 1) + (query_end - query_begin + 1))
                            #n = int(n*100)
                            #if n > longest_id_vector:
                            #    longest_id_vector = n
                            ###n = myarray[0]
                else:
                    count_id += line.count("|")

    if count_id and count_id > 0:
        perc_id = (count_id + count_id) / ((target_end - target_begin + 1) + (query_end - query_begin + 1))
        if perc_id >= 0.98 and count_id >= 18:
            if count_id > longest_id_vector:
                longest_id_vector = count_id
        if DEBUG == 1:
            print "VECTOR PERC_ID " + str(perc_id) + " COUNT_ID " + str(count_id)
        count_id = 0

    return longest_id_ecoli, longest_id_vector


def tenFoldingDeltaG(seq):
    ###Use ViennaRNA or http://mfold.rna.albany.edu/?q=mfold/DNA-Folding-Form
    cmd = ["/global/homes/a/andreopo/ViennaRNA/ViennaRNA-2.0.7h/Progs/RNAfold", "-P", "/global/homes/a/andreopo/ViennaRNA/ViennaRNA-2.0.7h/dna_mathews1999.par"] ###, ">" , outfile]
    if DEBUG == 1:
        print "cmd " + str(cmd)
    p = Popen(cmd, stdout = PIPE, stderr = PIPE, stdin = PIPE)
    p.stdin.write(seq)
    std_out, std_err = p.communicate()
    if DEBUG == 1:
        print "ENERGY_std_out " + std_out
    #print "ENERGY_std_err " + std_err
    #print "ENERGY_retcode " + str(p.returncode)

    energy = 0    
    #Parse the free energy number and return it
    m = re.match(".*\(-(.+)\).*", std_out.replace("\n", " "))
    if not m is None and len(m.groups()) > 0:
        energy = m.groups()[0]

    if DEBUG == 1:
        print "ENERGY " + str(energy)
    return -1 * float(energy)


#def classifySVM:
    ###  ./svm_learn seq_features seq_features.model
    ###  ./svm_classify seq_features seq_features.model seq_features.predict

rules = {
    "Terminal Repeats longer than 5 bp" : {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: This sequence cannot be made as gBlocks Gene Fragments because of the direct repeat {ViolatingSequence} that is found at each end of the gBlock, located at the ends of the sequence.  Solution: Change the sequence of one of the repeats or add supplementary bases to one or both ends of the sequence (Note: Restriction enzymes tend to perform better when 4-7 additional bases are added to the ends of the sequences). ",
        "explanation": "Explanation: Measures repeats that are 5 bases or longer on both terminal 30 base ends of the sequence. ",
        "X": 5,
        "value": 11
    },
    "Hairpin close stems" : {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: This gBlock cannot be made because a hairpin with the stem sequence of {Sequence} within a window size of {ForwardLocations} bases into the sequence exists. ",
        "explanation": "Explanation: Hairpin with a stem length greater than 15 bases where both stems are within 100 bases of each other. ",
        "X": 15,
        "Y": 100,
        "value": 11
    },
    "Hairpin anywhere stems": {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: This gBlock cannot be made because a hairpin with the stem sequence of {Sequence} and {RevComplSequence} exists in the sequence. ",
        "explanation": "Explanation: Hairpin with a stem length greater than 19 bases. ",
        "X": 19,
        "value": 11
    },
    "Windowed Hairpin > 80% GC of Length > 10 bp": {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: This sequence features a hairpin with the stem sequence {Sequence} and {RevComplSequence} that has a GC% >80%. ",
        "explanation": "Explanation: Hairpin within a window of 40 bases length that contains a stem with a length of at least 10. That stem has a GC content => 80%. ",
        "X": 40,
        "Y": 10,
        "Z": 0.80,
        "value": 9
    },
    "Homopolymeric A": {
        "name": fourFindLongestHomopolymer,
        "text": "Risk: This sequence cannot be synthesized as gBlocks Gene Fragments because a homopolymeric run of {SequenceLength} A bases is present starting at position {SequenceStart}. Solution: Redesign to reduce the homopolymeric run or contact our Gene specialists to assist you. Note that this sequence may be made as an IDT Gene or Minigene as the homopolymeric run of As is less than 21 bases. ",
        "explanation": "Explanation: 12 As in a row. ",
        "X": 12,
        "value": 11
    },
    "Homopolymeric C": {
        "name": fourFindLongestHomopolymer,
        "text": "Risk: This sequence cannot be synthesized as gBlocks Gene Fragments because a homopolymeric run of {SequenceLength} C bases is present starting at position {SequenceStart}. Solution: Redesign to reduce the homopolymeric run or contact our Gene specialists to assist you. Note that this sequence may be made as an IDT Gene or Minigene as the homopolymeric run of Cs is less than 21 bases. ",
        "explanation": "Explanation: 8 Cs in a row. ",
        "X": 8,
        "value": 11
    },
    "Homopolymeric G": {
        "name": fourFindLongestHomopolymer,
        "text": "Risk: This sequence cannot be synthesized as gBlocks Gene Fragments because a homopolymeric run of {SequenceLength} G bases is present starting at position {SequenceStart}. Solution: Redesign to reduce the homopolymeric run or contact our Gene specialists to assist you. Note that this sequence may be made as an IDT Gene or Minigene as the homopolymeric run of Gs is less than 21 bases. ",
        "explanation": "Explanation: 8 Gs in a row. ",
        "X": 8,
        "value": 11
    },
    "Homopolymeric T": {
        "name": fourFindLongestHomopolymer,
        "text": "Risk: This sequence cannot be synthesized as gBlocks Gene Fragments because a homopolymeric run of {SequenceLength} T bases is present starting at position {SequenceStart}. Solution: Redesign to reduce the homopolymeric run or contact our Gene specialists to assist you. Note that this sequence may be made as an IDT Gene or Minigene as the homopolymeric run of Ts is less than 21 bases. ",
        "explanation": "Explanation: 12 Ts in a row. ",
        "X": 12,
        "value": 11
    },
    "High GC % of whole gBlock": {
        "name": oneComputeGCtotal,
        "text": "Risk: This sequence cannot be made as gBlocks Gene Fragments because the total GC% is {Percentage} starting at position {SequenceStart}. Solution:  Redesign to reduce GC content below 69% or contact our Gene specialists to assist you. ",
        "explanation": "Explanation: Overall GC% cannot exceed 69%. ",
        "X": 0.69,
        "value": 11
    },
    "Low GC % of whole gBlock": {
        "name": oneComputeGCtotal,
        "text": "Risk: This sequence cannot be made as gBlocks Gene Fragments because the total GC% is {Percentage} starting at position {SequenceStart}. Solution:  Redesign to increase GC content above 25% or contact our Gene specialists to assist you. ",
        "explanation": "Explanation: Overall GC% must be at least 25%. ",
        "X": 0.25,
        "value": 11
    },
    "Windowed GC 100 bp high": {
        "name": twothreeComputeMinMaxGCinWindow,
        "text": "Risk: This sequence has a window of 100 bases with a GC% of {Percentage} starting at position {SequenceStart} . ",
        "explanation": "Explanation: GC% in any window of 100 bases cannot exceed 77 % . ",
        "X": 100,
        "Y": 0.77,
        "value": 11
    },
    "Windowed GC 100 bp low": {
        "name": twothreeComputeMinMaxGCinWindow,
        "text": "Risk: This sequence has a window of 100 bases with a GC% of {Percentage} starting at position {SequenceStart} . ",
        "explanation": "Explanation: GC% in any window of 100 bases cannot be below 28 % . ",
        "X": 100,
        "Y": 0.28,
        "value": 11
    },
    "Windowed GC 20 bp high": {
        "name": twothreeComputeMinMaxGCinWindow,
        "text": "Risk: This sequence has a window of 20 bases with a GC% of {Percentage} starting at position {SequenceStart} . ",
        "explanation": "Explanation: GC% in any window of 20 bases cannot exceed 90 % . ",
        "X": 20,
        "Y": 0.90,
        "value": 3
    },
    "Windowed GC 20 bp low": {
        "name": twothreeComputeMinMaxGCinWindow,
        "text": "Risk: This sequence has a window of 20 bases with a GC% of {Percentage} starting at position {SequenceStart} . ",
        "explanation": "Explanation: GC% in any window of 20 bases cannot be below 15 % . ",
        "X": 20,
        "Y": 0.15,
        "value": 3
    },
    "Terminal GC percentage high": {
        "name": twothreeComputeMinMaxGCinWindow,
        "text": "Risk: This sequence has a GC content of {Percentage} starting at position {SequenceStart} in one or both of the terminal 30 bases. Please redesign the terminals to have a GC content less than 76% . ",
        "explanation": "Explanation: The terminal 30 bases of either end cannont exceed 76 % GC . ",
        "X": 0.76,
        "value": 3
    },
    "Terminal GC percentage low": {
        "name": twothreeComputeMinMaxGCinWindow,
        "text": "Risk: This sequence has a GC content of {Percentage} starting at position {SequenceStart} in one or both of the terminal 30 bases. Please redesign the terminals to have a GC content greater than 24% . ",
        "explanation": "Explanation: The terminal 30 bases of either end cannot be below 24 % GC . ",
        "X": 0.24,
        "value": 3
    },
    "4repeat": {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: This gBlock cannot be made because more than 2 consecutive copies with the 4 nucleotide repeat {Sequence} exist in the sequence. ",
        "explanation": "Explanation: More than 2 concesecutive copies of a tetranucleotide repeat are present. ",
        "X": 2,
        "value": 11
    },
    "3repeat": {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: This gBlock cannot be made because more than 5 consecutive copies with the 3 nucleotide repeat {Sequence} exist in the sequence. ",
        "explanation": "Explanation: More than 5 concesecutive copies of a trinucleotide repeat are present. ",
        "X": 2,
        "value": 11
    },
    "2repeat": {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: This gBlock cannot be made because more than 9 consecutive copies with the 2 nucleotide repeat {Sequence} exist in the sequence. ",
        "explanation": "Explanation: More than 9 concesecutive copies of a dinucleotide repeat are present. ",
        "X": 2,
        "value": 11
    },
    "Repeat Percentage (total windowed repeats)": {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: Repeated sequences {Segment} of {MinRepeatLength} bases (more than 8 bases) comprise a combined {Percentage} bases of the sequence within a window of {WindowSize} bases.  Please redesign to reduce or eliminate these repeats. ",
        "explanation": "Explanation: Any window of 70 bases length that contains more than 89 % of the bases involved in repeats of 8 bases or longer that are repeated at least one time within the entire sequence. ",
        "X": 70,
        "Y": 0.89,
        "Z": 8,
        "value": 11
    },
    "Repeat Percentage > 35% (Individual repeats)": {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: The repeated sequence {Segment} constitutes {Percentage} bases of the gBlock. Please redesign to reduce this amount to be below {ThresholdPercentage} bases. ",
        "explanation": "Explanation: Any single repeat sequence that is over 8 bases in length and makes up more than 40 % of the sequence. ",
        "X": 8,
        "Y": 0.40,
        "value": 11
    },
    "Repeat Percentage (total repeats)": {
        "name": fivesixFindRepeatsInverted,
        "text": "Risk: Repeated sequences 8 bp and longer comprise {Percentage} bases of the sequence.  Please redesign to reduce to below {ThresholdPercentage} bases. ",
        "explanation": "Explanation: Any sequence contains more that 69 % of bases that are part of repeats 8 bases or longer. ",
        "X": 0.69,
        "Y": 8,
        "value": 11
    }
}


def process_seq(sequence, header, penalty_value, output_path, id_run, run_blast = None):
                    if DEBUG == 1:
                        print "\n\n_____________________________________\n  -------> Computing features for fasta header : " + header + "\n"
                    #Remove the vector sequence from all reads
                    '''
                    sequence = sequence.replace(vectorsequence5endA, "")
                    sequence = sequence.replace(rev_vectorsequence5endA, "")
                    sequence = sequence.replace(vectorsequence5endB, "")
                    sequence = sequence.replace(rev_vectorsequence5endB, "")
                    sequence = sequence.replace(vectorsequence5endC, "")
                    sequence = sequence.replace(rev_vectorsequence5endC, "")
                    sequence = sequence.replace(vectorsequence3end, "")
                    sequence = sequence.replace(rev_vectorsequence3end, "")
                    '''
                    
                    #TRIMREADS.write(line + "\n" + sequence + "\n")
                    if DEBUG == 1:
                        print "sequence %s" % ( sequence )
                        print "length of sequence %s" % ( len(sequence) )

                    holder_bad_gc = []
                    gc_content = oneComputeGCtotal(sequence);
                    key = "High GC % of whole gBlock"
                    if float(gc_content) > rules.get(key).get("X"):
                        holder_bad_gc.append(0)
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Percentage=gc_content, SequenceStart=0) + rules.get(key).get("explanation") )
                    key = "Low GC % of whole gBlock"
                    if float(gc_content) < rules.get(key).get("X"):
                        holder_bad_gc.append(0)
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Percentage=gc_content, SequenceStart=0) + rules.get(key).get("explanation") )
                    if DEBUG == 1:
                        print "oneComputeGCtotal for entire sequence"

                    #TODO coalesce windows/consolidate so the same does not get reported many times, especially when moving to the larger windows.
                    window = min(20, len(sequence));
                    mingc, mingcstart, maxgc, maxgcstart = twothreeComputeMinMaxGCinWindow(sequence, window);
                    key = "Windowed GC 20 bp high"
                    if maxgc > rules.get(key).get("Y") and not maxgcstart in holder_bad_gc:
                        holder_bad_gc.append(maxgcstart)
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Percentage=maxgc, SequenceStart=maxgcstart) + rules.get(key).get("explanation") )
                    key = "Windowed GC 20 bp low"
                    if mingc < rules.get(key).get("Y") and not mingcstart in holder_bad_gc:
                        holder_bad_gc.append(mingcstart)
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Percentage=mingc, SequenceStart=mingcstart) + rules.get(key).get("explanation") )
                    if DEBUG == 1:
                        print "twothreeComputeMinMaxGCinWindow for window " + str(window)

                    window = min(30, len(sequence));
                    mingc, mingcstart, maxgc, maxgcstart = twothreeComputeMinMaxGCinWindow(sequence, window, True);
                    key = "Terminal GC percentage high"
                    if maxgc > rules.get(key).get("X") and not maxgcstart in holder_bad_gc:
                        holder_bad_gc.append(maxgcstart)
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Percentage=maxgc, SequenceStart=maxgcstart) + rules.get(key).get("explanation") )
                    key = "Terminal GC percentage low"
                    if mingc < rules.get(key).get("X") and not mingcstart in holder_bad_gc:
                        holder_bad_gc.append(mingcstart)
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Percentage=mingc, SequenceStart=mingcstart) + rules.get(key).get("explanation") )
                    if DEBUG == 1:
                        print "twothreeComputeMinMaxGCinWindow for terminalend window " + str(window)

                    window = min(100, len(sequence));
                    mingc, mingcstart, maxgc, maxgcstart = twothreeComputeMinMaxGCinWindow(sequence, window);
                    key = "Windowed GC 100 bp high"
                    if maxgc > rules.get(key).get("Y") and not maxgcstart in holder_bad_gc:
                        holder_bad_gc.append(maxgcstart)
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Percentage=maxgc, SequenceStart=maxgcstart) + rules.get(key).get("explanation") )
                    key = "Windowed GC 100 bp low"
                    if mingc < rules.get(key).get("Y") and not mingcstart in holder_bad_gc:
                        holder_bad_gc.append(mingcstart)
                        penalty_value += float( rules.get(key).get("value") )
                        if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(Percentage=mingc, SequenceStart=mingcstart) + rules.get(key).get("explanation") )
                    if DEBUG == 1:
                        print "twothreeComputeMinMaxGCinWindow for window " + str(window)

                    ##TODO report where does the homopolymer start (coordinates for every homopolymer) and length, same for GC rules
                    ##TODO some reporting print statements show pointers to vars instead of the values, check all show values.
                    ##TODO check parallelism of writing output and errors, can you write results from parallel processes in order.
                    longestHomopol, totalLongHomopol, startLongestHomopol = fourFindLongestHomopolymer(sequence);
                    bases = ['A', 'C', 'G', 'T']
                    for b in bases:
                        key = "Homopolymeric "+b
                        if float(longestHomopol.get(b, 0)) > rules.get(key).get("X"):
                            penalty_value += float( rules.get(key).get("value") )
                            if  PRINT_SUCC_FAIL: print >>sys.stderr, str( key + " . " + rules.get(key).get("text").format(SequenceLength=longestHomopol.get(b, 0), SequenceStart=startLongestHomopol.get(b, 0)) + rules.get(key).get("explanation") )
                    if DEBUG == 1:
                        print "fourFindLongestHomopolymer"

                    longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, len_sequence, max_occur_same_dimer, percent_occur_same_dimer, max_occur_same_trimer, percent_occur_same_trimer, penalty_value = ('', 0, '', 0, '', 0, len(sequence), 0, 0, 0, 0, 0)
                    ###fivesixFindRepeatsInverted(sequence, penalty_value);
                    longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, len_sequence, penalty_value = fivesixFindRepeats(sequence, penalty_value);
                    if DEBUG == 1:
                        print "fivesixFindRepeats"
                        
                    len_sequence, max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, penalty_value = fivesixFindMers(sequence, penalty_value);
                    if DEBUG == 1:
                        print "fivesixFindMers"
                        
                    seq_in_filename = os.path.join(output_path, 'seq.fasta')
                    seq_in = open(seq_in_filename, 'w')
                    seq_in.write('>' + header + '\n')
                    seq_in.write(sequence)
                    seq_in.close()
                    longest_id_ecoli = 0
                    longest_id_vector = 0
                    if run_blast != None and run_blast == 'm':
                        longest_id_ecoli, longest_id_vector, penalty_value = runBlastMito(seq_in_filename, penalty_value); ###(0, 0) ###seveneightLongestAlignmentsEcoliVector(filename + '_' + header.replace(">",":"), sequence, ecoli_target_fasta, vector_target_fasta);
                    if run_blast != None and run_blast == 'p':
                        longest_id_ecoli, longest_id_vector, penalty_value = runBlastPlasmid(seq_in_filename, penalty_value); ###(0, 0) ###seveneightLongestAlignmentsEcoliVector(filename + '_' + header.replace(">",":"), sequence, ecoli_target_fasta, vector_target_fasta);
                    if DEBUG == 1:
                        print "seveneightLongestAlignmentsEcoliVector"
                    ###TODO remove
                    energy = 0 ###tenFoldingDeltaG(sequence)
                    if DEBUG == 1:
                        print "tenFoldingDeltaG"
                        
                        
                    print_read_features(output_path, id_run, gc_content, mingc, maxgc, longestHomopol, totalLongHomopol, longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, 0, 0, 0, 0, max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, header);
                    
                    return penalty_value;




THRESHOLD_VALUE_FAIL = 10
DEBUG = 0
PRINT_SUCC_FAIL = 0


if __name__ == "__main__":
    ###penalty_values = {}
    all_seqs = 0
    penalty_value = 0
    datetime = create_timestamp()
    ###id_run = None
    id_run = "DATETIME" + datetime

    succ = []
    fail = []
    
    desc = 'fasta_features'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-i", "--input-path", dest="input_path", help = "Input path to write to", required=False)
    parser.add_argument("-o", "--output-path", dest="output_path", help = "Output path to write to", required=False)
    parser.add_argument("-b", "-run-blast", dest="run_blast", help = "Specify either p for refseq.plasmid or m for refseq.mito", required=False)
    
    options = parser.parse_args()

    if (options.input_path == None) != (options.output_path == None):
        print "Either specify both input and output path or neither"
        exit(1)

    output_path = None
    input_path = None
    run_blast = None

    if (options.run_blast != None) and ((options.run_blast == 'p') or (options.run_blast == 'm')):
        run_blast = options.run_blast
    elif (options.run_blast != None):
        print("Specify either p for refseq.plasmid or m for refseq.mito")
        exit(1)
    
    ###Output directory
    if options.output_path:
        output_path = options.output_path
    else:
        ###base = os.path.basename(fastq)
        output_path = os.getcwd()  ###os.path.join( os.getcwd() , base )

    ###Do not allow the output path to have "homes" in it.
    if output_path.find("homes") > -1:
        print("The output path can not be under homes. The output_path given: " + output_path)
        exit(2)

    # create output_directory if it doesn't exist
    if not os.path.isdir(output_path):
        print "Cannot find %s, creating as new" % output_path
        os.makedirs(output_path)    
        os.chmod(output_path, 0775)

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
        
    '''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-o", "--output-path", dest="output_path", help = "Output path to write to", required=False)
    parser.add_argument("-i", "--id-run", dest="id_run", help = "The id for this run. You may want to use the input filename as id to distinguish the feature output file", required=False)
    
    options = parser.parse_args()
    
    output_path = None
    output_path = os.getcwd()

    ###TODO do not allow the output path to have "homes" in it.
    if options.output_path:
            output_path = options.output_path
    
    if options.id_run:
            id_run = options.id_run

    print >>sys.stdout, "Output path: %s" % (output_path)
    print >>sys.stdout, "Id for the run: %s" % (id_run)
    
    ###out_subdir = "/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/metagenomics_shoudan"
    ###1st batch
    ###fasta_files = [ "OOTO_4.fasta" ]
    ###fasta_files = ['B30_perfect_frags.txt', 'B30_no_screening_frags.txt'] ###'B30_single_error_frags_fix.txt', 'B30_multiple_errors_frags.txt',
    ###2nd batch
    ###fasta_files = ['no_assembly_frags_9_12_2013.txt', 'perfect_frags_9_12_2013.txt'] ###'1+_error_frags_9_12_2013.txt',
    ###fasta_files = ['B30_no_screening_frags.txt', 'no_assembly_frags_9_12_2013.txt']
    ###fasta_files = ['11_19_2013_No_evaluation_reference_seqs_pass.txt']
    ###fasta_files += ['11_19_2013_No_evaluation_reference_seqs_failed_assembly.txt']
    ###TODO Change for each run
    ###fasta_files = ['Brady_set3_CSPFY14_tests.fasta']
    ###fff = '_'.join(fasta_files)
    svn_file = os.path.join(output_path , "out_seq_features_" + str(id_run) + ".svn")
    txt_file = os.path.join(output_path , "out_seq_features_" + str(id_run) + ".txt")

    print >>sys.stdout, "Output text file with features: %s" % (txt_file)
    print >>sys.stdout, "Output text file with features for svn training: %s" % (svn_file)
    '''

    '''
    vectorsequence5endA = "CCGCCTTGTTTAACTTTAAGAAGGAGCCCTTCCCCCC"
    rev_vectorsequence5endA = vectorsequence5endA[::-1]

    vectorsequence5endB = "CCGCCTTGTTTAACTTTAAGAAGGAGCCCTTCCCCC"
    rev_vectorsequence5endB = vectorsequence5endB[::-1]

    vectorsequence5endC = "CCGCCTTGTTTAACTTTAAGAAGGAGCCCTTCCCC"
    rev_vectorsequence5endC = vectorsequence5endC[::-1]

    vectorsequence3end = "GGGAAGGGTGGGCGCGCCGACCC"
    rev_vectorsequence3end = vectorsequence3end[::-1]
    
    ecoli_target_fasta = "ecoli_target_fasta.fasta"
    vector_target_fasta = "vector_target_fasta.fasta"
    '''
    
    ###Print the header line.
    '''
    features_file = open(svn_file, 'w')
    features_file_excel = open(txt_file, 'w')
    features_file.write("pf,gccontent,mingc,maxgc,longestAhomopol,totalLongAhomopol,longestChomopol,totalLongChomopol,longestGhomopol,totalLongGhomopol,longestThomopol,totalLongThomopol,longestHomopolSum,totalLongHomopolSum,longestdirectrepeat,occurdirectrepeat,longestreverserepeat,occurreverserepeat,longestreversecomplrepeat,occurreversecomplrepeat,sequencelengthbases,ecolialign,vectoralign,freenergy,max_occur_same_dimer,percent_occur_same_dimer,max_occur_same_trimer,percent_occur_same_trimer,class,faheader\n")
    features_file.close()
    features_file_excel.write("pf,gccontent,mingc,maxgc,longestAhomopol,totalLongAhomopol,longestChomopol,totalLongChomopol,longestGhomopol,totalLongGhomopol,longestThomopol,totalLongThomopol,longestHomopolSum,totalLongHomopolSum,longestdirectrepeat,occurdirectrepeat,longestreverserepeat,occurreverserepeat,longestreversecomplrepeat,occurreversecomplrepeat,sequencelengthbases,ecolialign,vectoralign,freenergy,max_occur_same_dimer,percent_occur_same_dimer,max_occur_same_trimer,percent_occur_same_trimer,class,faheader\n")
    features_file_excel.close()
    '''
    colheaders = "id_run,gc_content,mingc,maxgc,AlongestHomopol,AtotalLongHomopol,ClongestHomopol,CtotalLongHomopol,GlongestHomopol,GtotalLongHomopol,TlongestHomopol,TtotalLongHomopol,longestHomopolSum,totalLongHomopolSum,longest_repeat,occur_longest_repeat,longest_rev_repeat,occur_longest_rev_repeat,longest_revcompl_repeat,occur_longest_revcompl_repeat,len_sequence,longest_id_ecoli,longest_id_vector,energy,max_occur_dimer,max_occur_same_dimer,percent_occur_same_dimer,max_occur_trimer,max_occur_same_trimer,percent_occur_same_trimer,max_occur_tetramer,max_occur_same_tetramer,percent_occur_same_tetramer"
    features_file = open(os.path.join(output_path, 'features.svn'), 'w')
    features_file.write(colheaders + "\n")
    features_file.flush()
    features_file.close()
    features_file_excel = open(os.path.join(output_path, 'features.txt'), 'w')
    features_file_excel.write(colheaders + "\n")
    features_file_excel.flush()
    features_file_excel.close()

    ###for filename in fasta_files:

    ###READS = open(os.path.join(out_subdir , filename), 'r')
    ###TRIMREADS = open(os.path.join(out_subdir , 'TRIMREADS_' + datetime + '_' + filename), 'w')

    sequence = ""

    #while True:
    for line in fileinput.input(input_path): ###READS####READS: ###:
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
                    penalty_value = process_seq(sequence, header, penalty_value, output_path, id_run, run_blast)

                    ###We just went through a sequence so print the prediction
                    if  PRINT_SUCC_FAIL:
                        print >>sys.stderr, "\nPenalty_value : " + str(penalty_value) + "    Threshold_value_fail : " + str(THRESHOLD_VALUE_FAIL) + "\n"
                        if penalty_value > THRESHOLD_VALUE_FAIL:
                            print >>sys.stderr, "\nFAILURE synthesis prediction for fasta header : " + header
                        else:
                            print >>sys.stderr, "\nSUCCESS synthesis prediction for fasta header : " + header
                    '''
                    if penalty_values.has_key(penalty_value):
                        penalty_values.get(penalty_value).append(header)
                    else:
                        penalty_values[penalty_value] = [ header ]
                    '''
                    if header.find("FAIL") > -1:
                        fail.append(penalty_value)
                        #if  PRINT_SUCC_FAIL: print >>sys.stderr, "APPENDED FAIL " + str(penalty_value) + " to " + str(fail)
                    else:
                        succ.append(penalty_value)
                        #if  PRINT_SUCC_FAIL: print >>sys.stderr, "APPENDED SUCC " + str(penalty_value) + " to " + str(succ)

                sequence = "";
                header = line.replace(" ", "_").replace(",", ".");
                if DEBUG == 1: print "header %s" % (header)
                if  PRINT_SUCC_FAIL: print >>sys.stderr, "\n\n_____________________________\nPredicting synthesis for fasta header : " + header + "\n"
                all_seqs += 1
                penalty_value = 0
            elif len(line) > 2:
                sequence += line.upper();

    penalty_value = process_seq(sequence, header, penalty_value, output_path, id_run, run_blast)
    
    ###We just went through a sequence so print the prediction
    if  PRINT_SUCC_FAIL:
        print >>sys.stderr, "\nPenalty_value : " + str(penalty_value) + "    Threshold_value_fail : " + str(THRESHOLD_VALUE_FAIL) + "\n"
        if penalty_value > THRESHOLD_VALUE_FAIL:
            print >>sys.stderr, "\nFAILURE synthesis prediction for fasta header : " + header
        else:
            print >>sys.stderr, "\nSUCCESS synthesis prediction for fasta header : " + header
    '''
    if penalty_values.has_key(penalty_value):
        penalty_values.get(penalty_value).append(header)
    else:
        penalty_values[penalty_value] = [ header ]
    '''
    if header.find("FAIL") > -1:
        fail.append(penalty_value)
        #if  PRINT_SUCC_FAIL: print >>sys.stderr, "APPENDED FAIL " + str(penalty_value) + " to " + str(fail)
    else:
        succ.append(penalty_value)
        #if  PRINT_SUCC_FAIL: print >>sys.stderr, "APPENDED SUCC " + str(penalty_value) + " to " + str(succ)

    '''    
    succ_file = open('LIST_success_class_penalties.txt', 'a')
    succ_file.write("Processed %s sequences: " % (all_seqs) )
    succ_file.write( str(succ) + " " )
    succ_file.flush()
    succ_file.close()
    
    fail_file = open('LIST_failure_class_penalties.txt', 'a')
    fail_file.write("Processed %s sequences: " % (all_seqs) )
    fail_file.write( str(fail) + " " )
    fail_file.flush()
    fail_file.close()
    '''
    
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, axisbg='#EEEEEE')
    ax.grid(color='white', linestyle='solid')
    ax.hist(np.array(succ), 10, histtype='stepfilled', fc='lightblue', alpha=0.5);
    ax.hist(np.array(fail), 10, histtype='stepfilled', fc='red', alpha=0.5);
    mpld3.save_html(fig, 'hist.html')
    '''
    
    '''
    segment = 0
    print "\n\nSequences ranked by likelihood to fail\n"
    for pen in sorted( penalty_values.iterkeys(), reverse=True ):
        segment += len ( penalty_values.get(pen) )
        for i in penalty_values.get(pen):
            print "%s is at the top %s of most likely to fail sequences" % ( i , float(segment) / float(all_seqs) )
    '''
    
    ###READS.close();
    ###TRIMREADS.close();
#andreopo@gpint108:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/mpld3$ git clone https://github.com/jakevdp/mpld3.git
#andreopo@gpint108:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/mpld3/mpld3$ python setup.py submodule
#andreopo@gpint108:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/mpld3/mpld3$ python setup.py build
#
#andreopo@gpint108:/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/mpld3/mpld3$ export PYTHONPATH=/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-rqc-synbio/io/mpld3/mpld3:$PYTHONPATH
#
#import mpld3
#import matplotlib
#matplotlib.use('Agg')
#import numpy as np
#import matplotlib.pyplot as plt
#fig, ax = plt.subplots()
#fig = plt.figure()
#
#ax = fig.add_subplot(111, axisbg='#EEEEEE')
#ax.grid(color='white', linestyle='solid')
#
#x = np.random.normal(size=1000)
#ax.hist(x, 30, histtype='stepfilled', fc='lightblue', alpha=0.5);
#
#np.array([[1,2,3,4,5,6,7,8,9], [5,6,7,12,15,20]])
#(array([ 2,  2,  2,  6,  8, 12, 23, 26, 31, 38, 47, 56, 64, 69, 79, 94, 80,
#       74, 74, 57, 44, 38, 19, 17, 13,  7,  9,  5,  2,  2]), array([-3.03001246, -2.83051633, -2.63102021, -2.43152408, -2.23202795,
#       -2.03253183, -1.8330357 , -1.63353957, -1.43404345, -1.23454732,
#       -1.03505119, -0.83555507, -0.63605894, -0.43656281, -0.23706668,
#       -0.03757056,  0.16192557,  0.3614217 ,  0.56091782,  0.76041395,
#        0.95991008,  1.1594062 ,  1.35890233,  1.55839846,  1.75789459,
#        1.95739071,  2.15688684,  2.35638297,  2.55587909,  2.75537522,
#        2.95487135]), <a list of 1 Patch objects>)
#>>> ax.hist(np.array([1,2,3,4,5,6,7,8,9]), 5, histtype='stepfilled', fc='lightblue', alpha=0.5);
#(array([1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1,
#       0, 0, 0, 1, 0, 0, 1]), array([ 1.        ,  1.26666667,  1.53333333,  1.8       ,  2.06666667,
#        2.33333333,  2.6       ,  2.86666667,  3.13333333,  3.4       ,
#        3.66666667,  3.93333333,  4.2       ,  4.46666667,  4.73333333,
#        5.        ,  5.26666667,  5.53333333,  5.8       ,  6.06666667,
#        6.33333333,  6.6       ,  6.86666667,  7.13333333,  7.4       ,
#        7.66666667,  7.93333333,  8.2       ,  8.46666667,  8.73333333,  9.        ]), <a list of 1 Patch objects>)
#>>>
#>>>
#>>> ax.hist(np.array([5,6,7,12,15,20]), 5, histtype='stepfilled', fc='red', alpha=0.5);
#(array([1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
#       0, 0, 0, 0, 0, 0, 1]), array([  5. ,   5.5,   6. ,   6.5,   7. ,   7.5,   8. ,   8.5,   9. ,
#         9.5,  10. ,  10.5,  11. ,  11.5,  12. ,  12.5,  13. ,  13.5,
#        14. ,  14.5,  15. ,  15.5,  16. ,  16.5,  17. ,  17.5,  18. ,
#        18.5,  19. ,  19.5,  20. ]), <a list of 1 Patch objects>)
#mpld3.fig_to_html(fig)
#mpld3.save_html(fig, '/tmp/test.html')
