#!/usr/bin/env python
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


'''
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir, './mpld3/mpld3'))

import mpld3
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
'''

###import khmer
import shutil

from JGI_Pipeline import JGI_Pipeline


#COGs = [ "gnl|CDD|30359", "gnl|CDD|30362", "gnl|CDD|30363", "gnl|CDD|30366", "gnl|CDD|30367", "gnl|CDD|30368", "gnl|CDD|30373", "gnl|CDD|30374", "gnl|CDD|30386", "gnl|CDD|30391", "gnl|CDD|30397", "gnl|CDD|30398", "gnl|CDD|30399", "gnl|CDD|30400", "gnl|CDD|30401", "gnl|CDD|30413", "gnl|CDD|30421", "gnl|CDD|30422", "gnl|CDD|30429", "gnl|CDD|30430", "gnl|CDD|30434", "gnl|CDD|30435", "gnl|CDD|30436", "gnl|CDD|30437", "gnl|CDD|30438", "gnl|CDD|30439", "gnl|CDD|30440", "gnl|CDD|30441", "gnl|CDD|30442", "gnl|CDD|30443", "gnl|CDD|30445", "gnl|CDD|30446", "gnl|CDD|30447", "gnl|CDD|30448", "gnl|CDD|30449", "gnl|CDD|30450", "gnl|CDD|30451", "gnl|CDD|30452", "gnl|CDD|30465", "gnl|CDD|30479", "gnl|CDD|30493", "gnl|CDD|30511", "gnl|CDD|30513", "gnl|CDD|30521", "gnl|CDD|30522", "gnl|CDD|30529", "gnl|CDD|30533", "gnl|CDD|30534", "gnl|CDD|30535", "gnl|CDD|30536", "gnl|CDD|30537", "gnl|CDD|30542", "gnl|CDD|30544", "gnl|CDD|30546", "gnl|CDD|30547", "gnl|CDD|30548", "gnl|CDD|30549", "gnl|CDD|30550", "gnl|CDD|30551", "gnl|CDD|30552", "gnl|CDD|30560", "gnl|CDD|30565", "gnl|CDD|30566", "gnl|CDD|30568", "gnl|CDD|30569", "gnl|CDD|30571", "gnl|CDD|30572", "gnl|CDD|30576", "gnl|CDD|30577", "gnl|CDD|30579", "gnl|CDD|30580", "gnl|CDD|30582", "gnl|CDD|30587", "gnl|CDD|30593", "gnl|CDD|30599", "gnl|CDD|30601", "gnl|CDD|30603", "gnl|CDD|30604", "gnl|CDD|30605", "gnl|CDD|30606", "gnl|CDD|30610", "gnl|CDD|30613", "gnl|CDD|30616", "gnl|CDD|30617", "gnl|CDD|30621", "gnl|CDD|30623", "gnl|CDD|30638", "gnl|CDD|30639", "gnl|CDD|30640", "gnl|CDD|30641", "gnl|CDD|30661", "gnl|CDD|30667", "gnl|CDD|30672", "gnl|CDD|30676", "gnl|CDD|30681", "gnl|CDD|30683", "gnl|CDD|30684", "gnl|CDD|30689", "gnl|CDD|30690", "gnl|CDD|30691", "gnl|CDD|30697", "gnl|CDD|30706", "gnl|CDD|30708", "gnl|CDD|30709", "gnl|CDD|30710", "gnl|CDD|30766", "gnl|CDD|30772", "gnl|CDD|30790", "gnl|CDD|30791", "gnl|CDD|30792", "gnl|CDD|30794", "gnl|CDD|30829", "gnl|CDD|30830", "gnl|CDD|30833", "gnl|CDD|30868", "gnl|CDD|30871", "gnl|CDD|30878", "gnl|CDD|30879", "gnl|CDD|30885", "gnl|CDD|30887", "gnl|CDD|30898", "gnl|CDD|30910", "gnl|CDD|30911", "gnl|CDD|30912", "gnl|CDD|30921", "gnl|CDD|30930", "gnl|CDD|30935", "gnl|CDD|30937", "gnl|CDD|30938", "gnl|CDD|30939", "gnl|CDD|30940", "gnl|CDD|30948", "gnl|CDD|30962", "gnl|CDD|30966", "gnl|CDD|30983", "gnl|CDD|30998", "gnl|CDD|31033", "gnl|CDD|31034", "gnl|CDD|31065", "gnl|CDD|31068", "gnl|CDD|31075", "gnl|CDD|31085", "gnl|CDD|31092", "gnl|CDD|31093", "gnl|CDD|31094", "gnl|CDD|31095", "gnl|CDD|31122", "gnl|CDD|31123", "gnl|CDD|31124", "gnl|CDD|31125", "gnl|CDD|31142", "gnl|CDD|31145", "gnl|CDD|31148", "gnl|CDD|31149", "gnl|CDD|31151", "gnl|CDD|31170", "gnl|CDD|31199", "gnl|CDD|31233", "gnl|CDD|31241", "gnl|CDD|31243", "gnl|CDD|31289", "gnl|CDD|31290", "gnl|CDD|31291", "gnl|CDD|31292", "gnl|CDD|31293", "gnl|CDD|31294", "gnl|CDD|31300", "gnl|CDD|31352", "gnl|CDD|31353", "gnl|CDD|31355", "gnl|CDD|31356", "gnl|CDD|31357", "gnl|CDD|31372", "gnl|CDD|31378", "gnl|CDD|31379", "gnl|CDD|31380", "gnl|CDD|31381", "gnl|CDD|31382", "gnl|CDD|31383", "gnl|CDD|31391", "gnl|CDD|31399", "gnl|CDD|31407", "gnl|CDD|31412", "gnl|CDD|31413", "gnl|CDD|31415", "gnl|CDD|31427", "gnl|CDD|31429", "gnl|CDD|31438", "gnl|CDD|31450", "gnl|CDD|31484", "gnl|CDD|31499", "gnl|CDD|31502", "gnl|CDD|31505", "gnl|CDD|31516", "gnl|CDD|31532", "gnl|CDD|31549", "gnl|CDD|31560", "gnl|CDD|31561", "gnl|CDD|31565", "gnl|CDD|31572", "gnl|CDD|31573", "gnl|CDD|31574", "gnl|CDD|31575", "gnl|CDD|31590", "gnl|CDD|31595", "gnl|CDD|31602", "gnl|CDD|31628", "gnl|CDD|31633", "gnl|CDD|31656", "gnl|CDD|31660", "gnl|CDD|31679", "gnl|CDD|31687", "gnl|CDD|31688", "gnl|CDD|31689", "gnl|CDD|31692", "gnl|CDD|31703", "gnl|CDD|31719", "gnl|CDD|31723", "gnl|CDD|31733", "gnl|CDD|31741", "gnl|CDD|31759", "gnl|CDD|31764", "gnl|CDD|31776", "gnl|CDD|31778", "gnl|CDD|31782", "gnl|CDD|31787", "gnl|CDD|31788", "gnl|CDD|31789", "gnl|CDD|31791", "gnl|CDD|31804", "gnl|CDD|31818", "gnl|CDD|31819", "gnl|CDD|31829", "gnl|CDD|31830", "gnl|CDD|31836", "gnl|CDD|31844", "gnl|CDD|31861", "gnl|CDD|31862", "gnl|CDD|31876", "gnl|CDD|31903", "gnl|CDD|31913", "gnl|CDD|31916", "gnl|CDD|31922", "gnl|CDD|31932", "gnl|CDD|31942", "gnl|CDD|31944", "gnl|CDD|31947", "gnl|CDD|31983", "gnl|CDD|31998", "gnl|CDD|32003", "gnl|CDD|32010", "gnl|CDD|32011", "gnl|CDD|32026", "gnl|CDD|32044", "gnl|CDD|32047", "gnl|CDD|32052", "gnl|CDD|32073", "gnl|CDD|32074", "gnl|CDD|32083", "gnl|CDD|32085", "gnl|CDD|32095", "gnl|CDD|32104", "gnl|CDD|32116", "gnl|CDD|32122", "gnl|CDD|32127", "gnl|CDD|32135", "gnl|CDD|32159", "gnl|CDD|32173", "gnl|CDD|32179", "gnl|CDD|32180", "gnl|CDD|32181", "gnl|CDD|32184", "gnl|CDD|32187", "gnl|CDD|32190", "gnl|CDD|32195", "gnl|CDD|32199", "gnl|CDD|32206", "gnl|CDD|32207", "gnl|CDD|32225", "gnl|CDD|32234", "gnl|CDD|32236", "gnl|CDD|32241", "gnl|CDD|32258", "gnl|CDD|32275", "gnl|CDD|32276", "gnl|CDD|32280", "gnl|CDD|32284", "gnl|CDD|32285", "gnl|CDD|32300", "gnl|CDD|32306", "gnl|CDD|32308", "gnl|CDD|32309", "gnl|CDD|32310", "gnl|CDD|32319", "gnl|CDD|32322", "gnl|CDD|32330", "gnl|CDD|32340", "gnl|CDD|32346", "gnl|CDD|32350", "gnl|CDD|32357", "gnl|CDD|32359", "gnl|CDD|32401", "gnl|CDD|32410", "gnl|CDD|32419", "gnl|CDD|32441", "gnl|CDD|32443", "gnl|CDD|32445", "gnl|CDD|32446", "gnl|CDD|32450", "gnl|CDD|32530", "gnl|CDD|32540", "gnl|CDD|32550", "gnl|CDD|32556", "gnl|CDD|32564", "gnl|CDD|32567", "gnl|CDD|32575", "gnl|CDD|32581", "gnl|CDD|32589", "gnl|CDD|32590", "gnl|CDD|32595", "gnl|CDD|32598", "gnl|CDD|32642", "gnl|CDD|32653", "gnl|CDD|32699", "gnl|CDD|32715", "gnl|CDD|32717", "gnl|CDD|32728", "gnl|CDD|32757", "gnl|CDD|32781", "gnl|CDD|32844", "gnl|CDD|32850", "gnl|CDD|32890", "gnl|CDD|32892", "gnl|CDD|32915", "gnl|CDD|32943", "gnl|CDD|32944", "gnl|CDD|33081", "gnl|CDD|33087", "gnl|CDD|33088", "gnl|CDD|33152", "gnl|CDD|33440", "gnl|CDD|33475", "gnl|CDD|33514", "gnl|CDD|33580", "gnl|CDD|33780", "gnl|CDD|33782", "gnl|CDD|33838", "gnl|CDD|33845", "gnl|CDD|33865", "gnl|CDD|33878", "gnl|CDD|33880", "gnl|CDD|34070", "gnl|CDD|34125", "gnl|CDD|34213", "gnl|CDD|34310", "gnl|CDD|34439", "gnl|CDD|34510", "gnl|CDD|34528", "gnl|CDD|34545", "gnl|CDD|34853", "gnl|CDD|34854", "gnl|CDD|34867", "gnl|CDD|34868", "gnl|CDD|34906", "gnl|CDD|34946", "gnl|CDD|34964", "gnl|CDD|34966", "gnl|CDD|35018", "gnl|CDD|35139", "gnl|CDD|58281", \
#      "gnl|CDD|30350", "gnl|CDD|30352", "gnl|CDD|30356", "gnl|CDD|30358", "gnl|CDD|30359", "gnl|CDD|30360", "gnl|CDD|30362", "gnl|CDD|30363", "gnl|CDD|30365", "gnl|CDD|30366", "gnl|CDD|30368", "gnl|CDD|30369", "gnl|CDD|30370", "gnl|CDD|30372", "gnl|CDD|30374", "gnl|CDD|30378", "gnl|CDD|30380", "gnl|CDD|30381", "gnl|CDD|30383", "gnl|CDD|30385", "gnl|CDD|30386", "gnl|CDD|30387", "gnl|CDD|30388", "gnl|CDD|30390", "gnl|CDD|30391", "gnl|CDD|30393", "gnl|CDD|30397", "gnl|CDD|30398", "gnl|CDD|30399", "gnl|CDD|30400", "gnl|CDD|30401", "gnl|CDD|30402", "gnl|CDD|30406", "gnl|CDD|30409", "gnl|CDD|30410", "gnl|CDD|30412", "gnl|CDD|30413", "gnl|CDD|30414", "gnl|CDD|30415", "gnl|CDD|30420", "gnl|CDD|30421", "gnl|CDD|30424", "gnl|CDD|30425", "gnl|CDD|30427", "gnl|CDD|30428", "gnl|CDD|30429", "gnl|CDD|30430", "gnl|CDD|30431", "gnl|CDD|30433", "gnl|CDD|30435", "gnl|CDD|30436", "gnl|CDD|30437", "gnl|CDD|30438", "gnl|CDD|30439", "gnl|CDD|30440", "gnl|CDD|30441", "gnl|CDD|30442", "gnl|CDD|30443", "gnl|CDD|30445", "gnl|CDD|30446", "gnl|CDD|30447", "gnl|CDD|30448", "gnl|CDD|30449", "gnl|CDD|30450", "gnl|CDD|30451", "gnl|CDD|30452", "gnl|CDD|30453", "gnl|CDD|30457", "gnl|CDD|30459", "gnl|CDD|30460", "gnl|CDD|30461", "gnl|CDD|30464", "gnl|CDD|30468", "gnl|CDD|30473", "gnl|CDD|30474", "gnl|CDD|30475", "gnl|CDD|30476", "gnl|CDD|30477", "gnl|CDD|30478", "gnl|CDD|30485", "gnl|CDD|30491", "gnl|CDD|30492", "gnl|CDD|30493", "gnl|CDD|30494", "gnl|CDD|30495", "gnl|CDD|30496", "gnl|CDD|30497", "gnl|CDD|30498", "gnl|CDD|30500", "gnl|CDD|30501", "gnl|CDD|30503", "gnl|CDD|30505", "gnl|CDD|30509", "gnl|CDD|30510", "gnl|CDD|30511", "gnl|CDD|30513", "gnl|CDD|30514", "gnl|CDD|30515", "gnl|CDD|30516", "gnl|CDD|30517", "gnl|CDD|30518", "gnl|CDD|30520", "gnl|CDD|30521", "gnl|CDD|30523", "gnl|CDD|30524", "gnl|CDD|30526", "gnl|CDD|30527", "gnl|CDD|30528", "gnl|CDD|30529", "gnl|CDD|30532", "gnl|CDD|30533", "gnl|CDD|30534", "gnl|CDD|30535", "gnl|CDD|30536", "gnl|CDD|30537", "gnl|CDD|30538", "gnl|CDD|30539", "gnl|CDD|30540", "gnl|CDD|30544", "gnl|CDD|30546", "gnl|CDD|30547", "gnl|CDD|30548", "gnl|CDD|30549", "gnl|CDD|30550", "gnl|CDD|30551", "gnl|CDD|30553", "gnl|CDD|30555", "gnl|CDD|30557", "gnl|CDD|30559", "gnl|CDD|30561", "gnl|CDD|30564", "gnl|CDD|30566", "gnl|CDD|30572", "gnl|CDD|30574", "gnl|CDD|30575", "gnl|CDD|30580", "gnl|CDD|30583", "gnl|CDD|30584", "gnl|CDD|30585", "gnl|CDD|30586", "gnl|CDD|30588", "gnl|CDD|30591", "gnl|CDD|30592", "gnl|CDD|30593", "gnl|CDD|30596", "gnl|CDD|30599", "gnl|CDD|30600", "gnl|CDD|30604", "gnl|CDD|30605", "gnl|CDD|30611", "gnl|CDD|30614", "gnl|CDD|30619", "gnl|CDD|30625", "gnl|CDD|30629", "gnl|CDD|30632", "gnl|CDD|30633", "gnl|CDD|30634", "gnl|CDD|30636", "gnl|CDD|30642", "gnl|CDD|30648", "gnl|CDD|30651", "gnl|CDD|30652", "gnl|CDD|30653", "gnl|CDD|30660", "gnl|CDD|30666", "gnl|CDD|30677", "gnl|CDD|30678", "gnl|CDD|30680", "gnl|CDD|30682", "gnl|CDD|30686", "gnl|CDD|30688", "gnl|CDD|30690", "gnl|CDD|30691", "gnl|CDD|30694", "gnl|CDD|30695", "gnl|CDD|30700", "gnl|CDD|30701", "gnl|CDD|30707", "gnl|CDD|30710", "gnl|CDD|30714", "gnl|CDD|30715", "gnl|CDD|30716", "gnl|CDD|30731", "gnl|CDD|30737", "gnl|CDD|30738", "gnl|CDD|30743", "gnl|CDD|30744", "gnl|CDD|30745", "gnl|CDD|30748", "gnl|CDD|30751", "gnl|CDD|30753", "gnl|CDD|30755", "gnl|CDD|30756", "gnl|CDD|30759", "gnl|CDD|30760", "gnl|CDD|30761", "gnl|CDD|30769", "gnl|CDD|30774", "gnl|CDD|30780", "gnl|CDD|30782", "gnl|CDD|30785", "gnl|CDD|30786", "gnl|CDD|30787", "gnl|CDD|30790", "gnl|CDD|30791", "gnl|CDD|30792", "gnl|CDD|30793", "gnl|CDD|30795", "gnl|CDD|30798", "gnl|CDD|30799", "gnl|CDD|30800", "gnl|CDD|30801", "gnl|CDD|30803", "gnl|CDD|30804", "gnl|CDD|30806", "gnl|CDD|30807", "gnl|CDD|30808", "gnl|CDD|30809", "gnl|CDD|30810", "gnl|CDD|30811", "gnl|CDD|30812", "gnl|CDD|30813", "gnl|CDD|30815", "gnl|CDD|30816", "gnl|CDD|30817", "gnl|CDD|30819", "gnl|CDD|30820", "gnl|CDD|30821", "gnl|CDD|30822", "gnl|CDD|30823", "gnl|CDD|30824", "gnl|CDD|30828", "gnl|CDD|30830", "gnl|CDD|30831", "gnl|CDD|30832", "gnl|CDD|30834", "gnl|CDD|30835", "gnl|CDD|30837", "gnl|CDD|30838", "gnl|CDD|30839", "gnl|CDD|30840", "gnl|CDD|30841", "gnl|CDD|30844", "gnl|CDD|30847", "gnl|CDD|30849", "gnl|CDD|30850", "gnl|CDD|30851", "gnl|CDD|30853", "gnl|CDD|30854", "gnl|CDD|30858", "gnl|CDD|30859", "gnl|CDD|30861", "gnl|CDD|30862", "gnl|CDD|30865", "gnl|CDD|30866", "gnl|CDD|30868", "gnl|CDD|30869", "gnl|CDD|30870", "gnl|CDD|30871", "gnl|CDD|30872", "gnl|CDD|30873", "gnl|CDD|30874", "gnl|CDD|30877", "gnl|CDD|30878", "gnl|CDD|30879", "gnl|CDD|30880", "gnl|CDD|30881", "gnl|CDD|30883", "gnl|CDD|30886", "gnl|CDD|30887", "gnl|CDD|30888", "gnl|CDD|30889", "gnl|CDD|30892", "gnl|CDD|30896", "gnl|CDD|30898", "gnl|CDD|30904", "gnl|CDD|30905", "gnl|CDD|30906", "gnl|CDD|30907", "gnl|CDD|30909", "gnl|CDD|30910", "gnl|CDD|30912", "gnl|CDD|30914", "gnl|CDD|30915", "gnl|CDD|30919", "gnl|CDD|30928", "gnl|CDD|30929", "gnl|CDD|30931", "gnl|CDD|30932", "gnl|CDD|30934", "gnl|CDD|30935", "gnl|CDD|30936", "gnl|CDD|30937", "gnl|CDD|30940", "gnl|CDD|30941", "gnl|CDD|30942", "gnl|CDD|30943", "gnl|CDD|30944", "gnl|CDD|30945", "gnl|CDD|30946", "gnl|CDD|30949", "gnl|CDD|30952", "gnl|CDD|30953", "gnl|CDD|30954", "gnl|CDD|30955", "gnl|CDD|30957", "gnl|CDD|30958", "gnl|CDD|30959", "gnl|CDD|30961", "gnl|CDD|30962", "gnl|CDD|30964", "gnl|CDD|30965", "gnl|CDD|30966", "gnl|CDD|30967", "gnl|CDD|30969", "gnl|CDD|30970", "gnl|CDD|30971", "gnl|CDD|30972", "gnl|CDD|30973", "gnl|CDD|30974", "gnl|CDD|30975", "gnl|CDD|30978", "gnl|CDD|30980", "gnl|CDD|30982", "gnl|CDD|30985", "gnl|CDD|30987", "gnl|CDD|30989", "gnl|CDD|30996", "gnl|CDD|30997", "gnl|CDD|30999", "gnl|CDD|31000", "gnl|CDD|31001", "gnl|CDD|31002", "gnl|CDD|31004", "gnl|CDD|31006", "gnl|CDD|31008", "gnl|CDD|31009", "gnl|CDD|31011", "gnl|CDD|31012", "gnl|CDD|31017", "gnl|CDD|31018", "gnl|CDD|31019", "gnl|CDD|31023", "gnl|CDD|31025", "gnl|CDD|31027", "gnl|CDD|31028", "gnl|CDD|31031", "gnl|CDD|31037", "gnl|CDD|31039", "gnl|CDD|31041", "gnl|CDD|31046", "gnl|CDD|31048", "gnl|CDD|31049", "gnl|CDD|31051", "gnl|CDD|31052", "gnl|CDD|31055", "gnl|CDD|31058", "gnl|CDD|31059", "gnl|CDD|31061", "gnl|CDD|31063", "gnl|CDD|31069", "gnl|CDD|31070", "gnl|CDD|31074", "gnl|CDD|31076", "gnl|CDD|31078", "gnl|CDD|31082", "gnl|CDD|31083", "gnl|CDD|31088", "gnl|CDD|31090", "gnl|CDD|31107", "gnl|CDD|31108", "gnl|CDD|31109", "gnl|CDD|31111", "gnl|CDD|31112", "gnl|CDD|31115", "gnl|CDD|31116", "gnl|CDD|31119", "gnl|CDD|31121", "gnl|CDD|31125", "gnl|CDD|31127", "gnl|CDD|31128", "gnl|CDD|31132", "gnl|CDD|31134", "gnl|CDD|31138", "gnl|CDD|31141", "gnl|CDD|31146", "gnl|CDD|31153", "gnl|CDD|31164", "gnl|CDD|31166", "gnl|CDD|31168", "gnl|CDD|31176", "gnl|CDD|31177", "gnl|CDD|31182", "gnl|CDD|31183", "gnl|CDD|31184", "gnl|CDD|31185", "gnl|CDD|31186", "gnl|CDD|31188", "gnl|CDD|31189", "gnl|CDD|31200", "gnl|CDD|31203", "gnl|CDD|31206", "gnl|CDD|31213", "gnl|CDD|31215", "gnl|CDD|31216", "gnl|CDD|31217", "gnl|CDD|31218", "gnl|CDD|31221", "gnl|CDD|31223", "gnl|CDD|31227", "gnl|CDD|31231", "gnl|CDD|31235", "gnl|CDD|31251", "gnl|CDD|31252", "gnl|CDD|31253", "gnl|CDD|31255", "gnl|CDD|31260", "gnl|CDD|31261", "gnl|CDD|31262", "gnl|CDD|31263", "gnl|CDD|31264", "gnl|CDD|31268", "gnl|CDD|31277", "gnl|CDD|31279", "gnl|CDD|31284", "gnl|CDD|31285", "gnl|CDD|31288", "gnl|CDD|31301", "gnl|CDD|31305", "gnl|CDD|31306", "gnl|CDD|31313", "gnl|CDD|31314", "gnl|CDD|31317", "gnl|CDD|31318", "gnl|CDD|31319", "gnl|CDD|31320", "gnl|CDD|31323", "gnl|CDD|31325", "gnl|CDD|31326", "gnl|CDD|31327", "gnl|CDD|31331", "gnl|CDD|31351", "gnl|CDD|31354", "gnl|CDD|31361", "gnl|CDD|31364", "gnl|CDD|31365", "gnl|CDD|31366", "gnl|CDD|31368", "gnl|CDD|31369", "gnl|CDD|31370", "gnl|CDD|31371", "gnl|CDD|31373", "gnl|CDD|31374", "gnl|CDD|31380", "gnl|CDD|31385", "gnl|CDD|31401", "gnl|CDD|31402", "gnl|CDD|31408", "gnl|CDD|31409", "gnl|CDD|31418", "gnl|CDD|31421", "gnl|CDD|31423", "gnl|CDD|31428", "gnl|CDD|31440", "gnl|CDD|31441", "gnl|CDD|31442", "gnl|CDD|31444", "gnl|CDD|31445", "gnl|CDD|31454", "gnl|CDD|31464", "gnl|CDD|31467", "gnl|CDD|31469", "gnl|CDD|31471", "gnl|CDD|31483", "gnl|CDD|31492", "gnl|CDD|31495", "gnl|CDD|31500", "gnl|CDD|31510", "gnl|CDD|31526", "gnl|CDD|31534", "gnl|CDD|31540", "gnl|CDD|31548", "gnl|CDD|31550", "gnl|CDD|31564", "gnl|CDD|31566", "gnl|CDD|31583", "gnl|CDD|31586", "gnl|CDD|31604", "gnl|CDD|31639", "gnl|CDD|31643", "gnl|CDD|31648", "gnl|CDD|31661", "gnl|CDD|31662", "gnl|CDD|31663", "gnl|CDD|31664", "gnl|CDD|31665", "gnl|CDD|31673", "gnl|CDD|31675", "gnl|CDD|31676", "gnl|CDD|31677", "gnl|CDD|31691", "gnl|CDD|31707", "gnl|CDD|31711", "gnl|CDD|31714", "gnl|CDD|31718", "gnl|CDD|31727", "gnl|CDD|31730", "gnl|CDD|31754", "gnl|CDD|31757", "gnl|CDD|31761", "gnl|CDD|31781", "gnl|CDD|31783", "gnl|CDD|31786", "gnl|CDD|31797", "gnl|CDD|31816", "gnl|CDD|31825", "gnl|CDD|31837", "gnl|CDD|31839", "gnl|CDD|31848", "gnl|CDD|31855", "gnl|CDD|31856", "gnl|CDD|31858", "gnl|CDD|31866", "gnl|CDD|31868", "gnl|CDD|31881", "gnl|CDD|31919", "gnl|CDD|31923", "gnl|CDD|31930", "gnl|CDD|31948", "gnl|CDD|31974", "gnl|CDD|31978", "gnl|CDD|31987", "gnl|CDD|31989", "gnl|CDD|32004", "gnl|CDD|32020", "gnl|CDD|32025", "gnl|CDD|32026", "gnl|CDD|32030", "gnl|CDD|32031", "gnl|CDD|32038", "gnl|CDD|32049", "gnl|CDD|32063", "gnl|CDD|32078", "gnl|CDD|32082", "gnl|CDD|32086", "gnl|CDD|32098", "gnl|CDD|32101", "gnl|CDD|32106", "gnl|CDD|32120", "gnl|CDD|32123", "gnl|CDD|32125", "gnl|CDD|32126", "gnl|CDD|32142", "gnl|CDD|32143", "gnl|CDD|32144", "gnl|CDD|32158", "gnl|CDD|32171", "gnl|CDD|32185", "gnl|CDD|32193", "gnl|CDD|32203", "gnl|CDD|32208", "gnl|CDD|32209", "gnl|CDD|32213", "gnl|CDD|32233", "gnl|CDD|32238", "gnl|CDD|32240", "gnl|CDD|32242", "gnl|CDD|32243", "gnl|CDD|32255", "gnl|CDD|32259", "gnl|CDD|32262", "gnl|CDD|32263", "gnl|CDD|32267", "gnl|CDD|32269", "gnl|CDD|32296", "gnl|CDD|32297", "gnl|CDD|32307", "gnl|CDD|32311", "gnl|CDD|32315", "gnl|CDD|32318", "gnl|CDD|32324", "gnl|CDD|32329", "gnl|CDD|32331", "gnl|CDD|32339", "gnl|CDD|32342", "gnl|CDD|32344", "gnl|CDD|32367", "gnl|CDD|32369", "gnl|CDD|32371", "gnl|CDD|32379", "gnl|CDD|32381", "gnl|CDD|32382", "gnl|CDD|32386", "gnl|CDD|32389", "gnl|CDD|32393", "gnl|CDD|32398", "gnl|CDD|32399", "gnl|CDD|32402", "gnl|CDD|32405", "gnl|CDD|32408", "gnl|CDD|32425", "gnl|CDD|32440", "gnl|CDD|32446", "gnl|CDD|32452", "gnl|CDD|32455", "gnl|CDD|32456", "gnl|CDD|32458", "gnl|CDD|32464", "gnl|CDD|32473", "gnl|CDD|32487", "gnl|CDD|32488", "gnl|CDD|32514", "gnl|CDD|32519", "gnl|CDD|32525", "gnl|CDD|32535", "gnl|CDD|32552", "gnl|CDD|32594", "gnl|CDD|32600", "gnl|CDD|32607", "gnl|CDD|32615", "gnl|CDD|32632", "gnl|CDD|32633", "gnl|CDD|32635", "gnl|CDD|32636", "gnl|CDD|32643", "gnl|CDD|32654", "gnl|CDD|32670", "gnl|CDD|32684", "gnl|CDD|32711", "gnl|CDD|32715", "gnl|CDD|32740", "gnl|CDD|32753", "gnl|CDD|32755", "gnl|CDD|32767", "gnl|CDD|32771", "gnl|CDD|32783", "gnl|CDD|32836", "gnl|CDD|32853", "gnl|CDD|32855", "gnl|CDD|32891", "gnl|CDD|32904", "gnl|CDD|32907", "gnl|CDD|32933", "gnl|CDD|32935", "gnl|CDD|32986", "gnl|CDD|32990", "gnl|CDD|32994", "gnl|CDD|33001", "gnl|CDD|33006", "gnl|CDD|33016", "gnl|CDD|33019", "gnl|CDD|33022", "gnl|CDD|33034", "gnl|CDD|33043", "gnl|CDD|33078", "gnl|CDD|33095", "gnl|CDD|33102", "gnl|CDD|33125", "gnl|CDD|33130", "gnl|CDD|33137", "gnl|CDD|33144", "gnl|CDD|33145", "gnl|CDD|33153", "gnl|CDD|33185", "gnl|CDD|33192", "gnl|CDD|33193", "gnl|CDD|33198", "gnl|CDD|33202", "gnl|CDD|33219", "gnl|CDD|33221", "gnl|CDD|33242", "gnl|CDD|33254", "gnl|CDD|33267", "gnl|CDD|33288", "gnl|CDD|33289", "gnl|CDD|33304", "gnl|CDD|33307", "gnl|CDD|33308", "gnl|CDD|33341", "gnl|CDD|33346", "gnl|CDD|33349", "gnl|CDD|33351", "gnl|CDD|33352", "gnl|CDD|33372", "gnl|CDD|33393", "gnl|CDD|33408", "gnl|CDD|33417", "gnl|CDD|33428", "gnl|CDD|33434", "gnl|CDD|33452", "gnl|CDD|33453", "gnl|CDD|33455", "gnl|CDD|33464", "gnl|CDD|33466", "gnl|CDD|33474", "gnl|CDD|33479", "gnl|CDD|33482", "gnl|CDD|33492", "gnl|CDD|33496", "gnl|CDD|33497", "gnl|CDD|33499", "gnl|CDD|33506", "gnl|CDD|33531", "gnl|CDD|33537", "gnl|CDD|33539", "gnl|CDD|33586", "gnl|CDD|33624", "gnl|CDD|33628", "gnl|CDD|33630", "gnl|CDD|33631", "gnl|CDD|33633", "gnl|CDD|33636", "gnl|CDD|33637", "gnl|CDD|33692", "gnl|CDD|33702", "gnl|CDD|33706", "gnl|CDD|33752", "gnl|CDD|33872", "gnl|CDD|33875", "gnl|CDD|33906", "gnl|CDD|33908", "gnl|CDD|33916", "gnl|CDD|33917", "gnl|CDD|33942", "gnl|CDD|33946", "gnl|CDD|33952", "gnl|CDD|33968", "gnl|CDD|33987", "gnl|CDD|34022", "gnl|CDD|34132", "gnl|CDD|34135", "gnl|CDD|34204", "gnl|CDD|34207", "gnl|CDD|34209", "gnl|CDD|34222", "gnl|CDD|34223", "gnl|CDD|34229", "gnl|CDD|34234", "gnl|CDD|34257", "gnl|CDD|34262", "gnl|CDD|34263", "gnl|CDD|34265", "gnl|CDD|34287", "gnl|CDD|34296", "gnl|CDD|34300", "gnl|CDD|34349", "gnl|CDD|34384", "gnl|CDD|34399", "gnl|CDD|34400", "gnl|CDD|34401", "gnl|CDD|34454", "gnl|CDD|34495", "gnl|CDD|34556", "gnl|CDD|34565", "gnl|CDD|34566", "gnl|CDD|34569", "gnl|CDD|34579", "gnl|CDD|34580", "gnl|CDD|34583", "gnl|CDD|34611", "gnl|CDD|34865", "gnl|CDD|34891", "gnl|CDD|34899", "gnl|CDD|35023", "gnl|CDD|35048", "gnl|CDD|35051", "gnl|CDD|35058", "gnl|CDD|35076", "gnl|CDD|35086", "gnl|CDD|35093", "gnl|CDD|35128", "gnl|CDD|35165", "gnl|CDD|35175", "gnl|CDD|35201", "gnl|CDD|35213", "gnl|CDD|35214", "gnl|CDD|35218", "gnl|CDD|58280", "gnl|CDD|58285"]
###COGs = [ "gnl|CDD|9885", "gnl|CDD|9888", "gnl|CDD|9889", "gnl|CDD|9892", "gnl|CDD|9893", "gnl|CDD|9894", "gnl|CDD|9899", "gnl|CDD|9900", "gnl|CDD|9912", "gnl|CDD|9917", "gnl|CDD|9923", "gnl|CDD|9924", "gnl|CDD|9925", "gnl|CDD|9926", "gnl|CDD|9927", "gnl|CDD|9939", "gnl|CDD|9947", "gnl|CDD|9948", "gnl|CDD|9955", "gnl|CDD|9956", "gnl|CDD|9960", "gnl|CDD|9961", "gnl|CDD|9962", "gnl|CDD|9963", "gnl|CDD|9964", "gnl|CDD|9965", "gnl|CDD|9966", "gnl|CDD|9967", "gnl|CDD|9968", "gnl|CDD|9969", "gnl|CDD|9971", "gnl|CDD|9972", "gnl|CDD|9973", "gnl|CDD|9974", "gnl|CDD|9975", "gnl|CDD|9976", "gnl|CDD|9977", "gnl|CDD|9978", "gnl|CDD|9991", "gnl|CDD|10005", "gnl|CDD|10019", "gnl|CDD|10037", "gnl|CDD|10039", "gnl|CDD|10047", "gnl|CDD|10048", "gnl|CDD|10055", "gnl|CDD|10059", "gnl|CDD|10060", "gnl|CDD|10061", "gnl|CDD|10062", "gnl|CDD|10063", "gnl|CDD|10068", "gnl|CDD|10070", "gnl|CDD|10072", "gnl|CDD|10073", "gnl|CDD|10074", "gnl|CDD|10075", "gnl|CDD|10076", "gnl|CDD|10077", "gnl|CDD|10078", "gnl|CDD|10086", "gnl|CDD|10091", "gnl|CDD|10092", "gnl|CDD|10094", "gnl|CDD|10095", "gnl|CDD|10097", "gnl|CDD|10098", "gnl|CDD|10102", "gnl|CDD|10103", "gnl|CDD|10105", "gnl|CDD|10106", "gnl|CDD|10108", "gnl|CDD|10113", "gnl|CDD|10119", "gnl|CDD|10125", "gnl|CDD|10127", "gnl|CDD|10129", "gnl|CDD|10130", "gnl|CDD|10131", "gnl|CDD|10132", "gnl|CDD|10136", "gnl|CDD|10139", "gnl|CDD|10142", "gnl|CDD|10143", "gnl|CDD|10147", "gnl|CDD|10149", "gnl|CDD|10164", "gnl|CDD|10165", "gnl|CDD|10166", "gnl|CDD|10167", "gnl|CDD|10187", "gnl|CDD|10193", "gnl|CDD|10198", "gnl|CDD|10202", "gnl|CDD|10207", "gnl|CDD|10209", "gnl|CDD|10210", "gnl|CDD|10215", "gnl|CDD|10216", "gnl|CDD|10217", "gnl|CDD|10223", "gnl|CDD|10231", "gnl|CDD|10233", "gnl|CDD|10234", "gnl|CDD|10235", "gnl|CDD|10291", "gnl|CDD|10297", "gnl|CDD|10315", "gnl|CDD|10316", "gnl|CDD|10317", "gnl|CDD|10319", "gnl|CDD|10343", "gnl|CDD|10354", "gnl|CDD|10355", "gnl|CDD|10358", "gnl|CDD|10393", "gnl|CDD|10396", "gnl|CDD|10403", "gnl|CDD|10404", "gnl|CDD|10410", "gnl|CDD|10412", "gnl|CDD|10423", "gnl|CDD|10435", "gnl|CDD|10436", "gnl|CDD|10437", "gnl|CDD|10446", "gnl|CDD|10455", "gnl|CDD|10460", "gnl|CDD|10462", "gnl|CDD|10463", "gnl|CDD|10464", "gnl|CDD|10465", "gnl|CDD|10473", "gnl|CDD|10487", "gnl|CDD|10491", "gnl|CDD|10508", "gnl|CDD|10523", "gnl|CDD|10558", "gnl|CDD|10559", "gnl|CDD|10590", "gnl|CDD|10593", "gnl|CDD|10600", "gnl|CDD|10610", "gnl|CDD|10617", "gnl|CDD|10618", "gnl|CDD|10619", "gnl|CDD|10620", "gnl|CDD|10647", "gnl|CDD|10648", "gnl|CDD|10649", "gnl|CDD|10650", "gnl|CDD|10667", "gnl|CDD|10670", "gnl|CDD|10673", "gnl|CDD|10674", "gnl|CDD|10676", "gnl|CDD|10695", "gnl|CDD|10724", "gnl|CDD|10758", "gnl|CDD|10766", "gnl|CDD|10768", "gnl|CDD|10814", "gnl|CDD|10815", "gnl|CDD|10816", "gnl|CDD|10817", "gnl|CDD|10818", "gnl|CDD|10819", "gnl|CDD|10825", "gnl|CDD|10877", "gnl|CDD|10878", "gnl|CDD|10880", "gnl|CDD|10881", "gnl|CDD|10882", "gnl|CDD|10897", "gnl|CDD|10903", "gnl|CDD|10904", "gnl|CDD|10905", "gnl|CDD|10906", "gnl|CDD|10907", "gnl|CDD|10908", "gnl|CDD|10916", "gnl|CDD|10924", "gnl|CDD|10932", "gnl|CDD|10937", "gnl|CDD|10938", "gnl|CDD|10940", "gnl|CDD|10952", "gnl|CDD|10954", "gnl|CDD|10963", "gnl|CDD|10975", "gnl|CDD|11009", "gnl|CDD|11024", "gnl|CDD|11027", "gnl|CDD|11030", "gnl|CDD|11041", "gnl|CDD|11057", "gnl|CDD|11074", "gnl|CDD|11085", "gnl|CDD|11086", "gnl|CDD|11090", "gnl|CDD|11097", "gnl|CDD|11098", "gnl|CDD|11099", "gnl|CDD|11100", "gnl|CDD|11115", "gnl|CDD|11120", "gnl|CDD|11127", "gnl|CDD|11153", "gnl|CDD|11158", "gnl|CDD|11181", "gnl|CDD|11185", "gnl|CDD|11204", "gnl|CDD|11212", "gnl|CDD|11213", "gnl|CDD|11214", "gnl|CDD|11217", "gnl|CDD|11228", "gnl|CDD|11244", "gnl|CDD|11248", "gnl|CDD|11258", "gnl|CDD|11266", "gnl|CDD|11284", "gnl|CDD|11289", "gnl|CDD|11301", "gnl|CDD|11303", "gnl|CDD|11307", "gnl|CDD|11312", "gnl|CDD|11313", "gnl|CDD|11314", "gnl|CDD|11316", "gnl|CDD|11329", "gnl|CDD|11343", "gnl|CDD|11344", "gnl|CDD|11354", "gnl|CDD|11355", "gnl|CDD|11361", "gnl|CDD|11369", "gnl|CDD|11386", "gnl|CDD|11387", "gnl|CDD|11401", "gnl|CDD|11428", "gnl|CDD|11438", "gnl|CDD|11441", "gnl|CDD|11447", "gnl|CDD|11457", "gnl|CDD|11467", "gnl|CDD|11469", "gnl|CDD|11472", "gnl|CDD|11508", "gnl|CDD|11523", "gnl|CDD|11528", "gnl|CDD|11535", "gnl|CDD|11536", "gnl|CDD|11551", "gnl|CDD|11569", "gnl|CDD|11572", "gnl|CDD|11577", "gnl|CDD|11598", "gnl|CDD|11599", "gnl|CDD|11608", "gnl|CDD|11610", "gnl|CDD|11620", "gnl|CDD|11629", "gnl|CDD|11641", "gnl|CDD|11647", "gnl|CDD|11652", "gnl|CDD|11660", "gnl|CDD|11684", "gnl|CDD|11698", "gnl|CDD|11704", "gnl|CDD|11705", "gnl|CDD|11706", "gnl|CDD|11709", "gnl|CDD|11712", "gnl|CDD|11715", "gnl|CDD|11720", "gnl|CDD|11724", "gnl|CDD|11731", "gnl|CDD|11732", "gnl|CDD|11750", "gnl|CDD|11759", "gnl|CDD|11761", "gnl|CDD|11766", "gnl|CDD|11783", "gnl|CDD|11800", "gnl|CDD|11801", "gnl|CDD|11805", "gnl|CDD|11809", "gnl|CDD|11810", "gnl|CDD|11825", "gnl|CDD|11831", "gnl|CDD|11833", "gnl|CDD|11834", "gnl|CDD|11835", "gnl|CDD|11844", "gnl|CDD|11847", "gnl|CDD|11855", "gnl|CDD|11865", "gnl|CDD|11871", "gnl|CDD|11875", "gnl|CDD|11882", "gnl|CDD|11884", "gnl|CDD|11926", "gnl|CDD|11935", "gnl|CDD|11944", "gnl|CDD|11966", "gnl|CDD|11968", "gnl|CDD|11970", "gnl|CDD|11971", "gnl|CDD|11975", "gnl|CDD|12055", "gnl|CDD|12065", "gnl|CDD|12075", "gnl|CDD|12081", "gnl|CDD|12089", "gnl|CDD|12092", "gnl|CDD|12100", "gnl|CDD|12106", "gnl|CDD|12114", "gnl|CDD|12115", "gnl|CDD|12120", "gnl|CDD|12123", "gnl|CDD|12167", "gnl|CDD|12178", "gnl|CDD|12224", "gnl|CDD|12240", "gnl|CDD|12242", "gnl|CDD|12253", "gnl|CDD|12282", "gnl|CDD|12306", "gnl|CDD|12369", "gnl|CDD|12375", "gnl|CDD|12415", "gnl|CDD|12417", "gnl|CDD|12440", "gnl|CDD|12468", "gnl|CDD|12469", "gnl|CDD|12606", "gnl|CDD|12612", "gnl|CDD|12613", "gnl|CDD|12677", "gnl|CDD|12965", "gnl|CDD|13000", "gnl|CDD|13039", "gnl|CDD|13105", "gnl|CDD|13305", "gnl|CDD|13307", "gnl|CDD|13363", "gnl|CDD|13370", "gnl|CDD|13390", "gnl|CDD|13403", "gnl|CDD|13405", "gnl|CDD|13595", "gnl|CDD|13650", "gnl|CDD|13738", "gnl|CDD|13835", "gnl|CDD|13964", "gnl|CDD|14035", "gnl|CDD|14053", "gnl|CDD|14070", "gnl|CDD|14378", "gnl|CDD|14379", "gnl|CDD|14392", "gnl|CDD|14393", "gnl|CDD|14431", "gnl|CDD|14471", "gnl|CDD|14489", "gnl|CDD|14491", "gnl|CDD|14543", "gnl|CDD|14664" ]

COGs = ["COG0009", "COG0012", "COG0013", "COG0016", "COG0017", "COG0018", "COG0023", "COG0024", "COG0037", "COG0042", "COG0048", "COG0049", "COG0050", "COG0051", "COG0052", "COG0064", "COG0072", "COG0073", "COG0080", "COG0081", "COG0085", "COG0086", "COG0087", "COG0088", "COG0089", "COG0090", "COG0091", "COG0092", "COG0093", "COG0094", "COG0096", "COG0097", "COG0098", "COG0099", "COG0100", "COG0101", "COG0102", "COG0103", "COG0116", "COG0130", "COG0144", "COG0162", "COG0164", "COG0172", "COG0173", "COG0180", "COG0184", "COG0185", "COG0186", "COG0187", "COG0188", "COG0193", "COG0195", "COG0197", "COG0198", "COG0199", "COG0200", "COG0201", "COG0202", "COG0203", "COG0211", "COG0216", "COG0217", "COG0219", "COG0220", "COG0222", "COG0223", "COG0227", "COG0228", "COG0230", "COG0231", "COG0233", "COG0238", "COG0244", "COG0250", "COG0252", "COG0254", "COG0255", "COG0256", "COG0257", "COG0261", "COG0264", "COG0267", "COG0268", "COG0272", "COG0275", "COG0290", "COG0291", "COG0292", "COG0293", "COG0313", "COG0319", "COG0324", "COG0328", "COG0333", "COG0335", "COG0336", "COG0341", "COG0342", "COG0343", "COG0349", "COG0357", "COG0359", "COG0360", "COG0361", "COG0417", "COG0423", "COG0441", "COG0442", "COG0443", "COG0445", "COG0470", "COG0481", "COG0482", "COG0486", "COG0522", "COG0525", "COG0532", "COG0533", "COG0539", "COG0541", "COG0552", "COG0564", "COG0565", "COG0566", "COG0576", "COG0585", "COG0590", "COG0592", "COG0593", "COG0594", "COG0595", "COG0603", "COG0617", "COG0621", "COG0638", "COG0653", "COG0689", "COG0690", "COG0721", "COG0724", "COG0731", "COG0742", "COG0749", "COG0750", "COG0751", "COG0752", "COG0779", "COG0780", "COG0781", "COG0782", "COG0799", "COG0802", "COG0805", "COG0806", "COG0809", "COG0828", "COG0858", "COG1030", "COG1039", "COG1041", "COG1092", "COG1093", "COG1094", "COG1095", "COG1096", "COG1097", "COG1103", "COG1158", "COG1159", "COG1161", "COG1162", "COG1163", "COG1179", "COG1185", "COG1186", "COG1187", "COG1188", "COG1189", "COG1190", "COG1198", "COG1206", "COG1214", "COG1219", "COG1220", "COG1222", "COG1234", "COG1236", "COG1245", "COG1258", "COG1293", "COG1308", "COG1311", "COG1314", "COG1325", "COG1341", "COG1358", "COG1369", "COG1370", "COG1374", "COG1382", "COG1383", "COG1384", "COG1385", "COG1400", "COG1405", "COG1412", "COG1439", "COG1444", "COG1467", "COG1471", "COG1490", "COG1498", "COG1499", "COG1500", "COG1503", "COG1514", "COG1530", "COG1534", "COG1544", "COG1552", "COG1571", "COG1576", "COG1588", "COG1590", "COG1594", "COG1599", "COG1600", "COG1601", "COG1603", "COG1617", "COG1631", "COG1632", "COG1643", "COG1644", "COG1650", "COG1658", "COG1675", "COG1676", "COG1690", "COG1717", "COG1727", "COG1730", "COG1736", "COG1746", "COG1756", "COG1758", "COG1761", "COG1798", "COG1813", "COG1818", "COG1825", "COG1826", "COG1841", "COG1859", "COG1862", "COG1867", "COG1889", "COG1890", "COG1899", "COG1901", "COG1911", "COG1921", "COG1933", "COG1939", "COG1944", "COG1952", "COG1976", "COG1990", "COG1996", "COG1997", "COG1998", "COG2001", "COG2004", "COG2007", "COG2012", "COG2016", "COG2023", "COG2024", "COG2042", "COG2051", "COG2053", "COG2058", "COG2075", "COG2092", "COG2093", "COG2097", "COG2101", "COG2102", "COG2117", "COG2123", "COG2125", "COG2126", "COG2127", "COG2136", "COG2139", "COG2147", "COG2157", "COG2163", "COG2167", "COG2174", "COG2176", "COG2219", "COG2229", "COG2238", "COG2260", "COG2262", "COG2264", "COG2265", "COG2269", "COG2384", "COG2404", "COG2419", "COG2428", "COG2443", "COG2451", "COG2501", "COG2511", "COG2519", "COG2520", "COG2603", "COG2606", "COG2813", "COG2825", "COG2872", "COG2890", "COG2892", "COG2904", "COG2933", "COG2961", "COG3028", "COG3036", "COG3076", "COG3078", "COG3101", "COG3129", "COG3130", "COG3270", "COG3276", "COG3277", "COG3343", "COG3642", "COG3678", "COG3719", "COG3785", "COG4021", "COG4023", "COG4080", "COG4088", "COG4108", "COG4121", "COG4123", "COG4352", "COG4445", "COG4575", "COG4694", "COG4830", "COG4901", "COG4919", "COG4937", "COG5256", "COG5257", "COG5270", "COG5271", "COG5324", "COG5384", "COG5405", "COG5407", "COG5459", "COG5580"]

'''
COGs not used atm because I only had chrom-specific COGs, not plasmid-specific COGs.
Also the COG computation wasn't that fast. Alevy said he'd send plasmid-specific COGs once they're ready.
Note use of COGs should improve accuracy of ML a lot.
'''
def cogs(sequence, seqin, penalty_value):
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
        cmd = ["shifter", "--image=registry.services.nersc.gov/jgi/prodigal:latest", "prodigal", "-a", FASTA + ".gene.faa", "-d", FASTA + ".gene.fasta", "-i", FASTA, "-o", FASTA + ".prodigal.out", "-p", "meta" ]
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
        genesperMB=(float(ones)/float(contiglen))
        prodigalFile.close()

        ###prodigal input contig file output genes as aas
        ###rpsblast on genes as aas
        if os.path.getsize(FASTA + ".gene.faa") > 0:
            Cog = "/global/dna/projectdirs/microbial/img/databases/cog/Cog"  ###"/global/projectb/sandbox/IMG/img/databases/cog"
            ###cmd2 = ["/usr/bin/rpsblast", "-m", "8", "-d", Cog, "-i", FASTA + ".gene.faa", "-b", "1", "-a", "1", "-e", "1e-2", "-l", FASTA + ".gene.faa.rpsblast.log" ]    ####
            cmd2 = ["shifter", "--image=registry.services.nersc.gov/jgi/hmmer:latest", "hmmsearch", "--cut_ga", "--domtblout", FASTA + ".COG.hmm.hmmsearch.domtblout.txt", "/global/projectb/scratch/andreopo/GAA-1290_plasmids/Dingl/microbial/NEW3_ALL_FEATURES/tmp170.refseq.bacteria.NOplasmids_NOmito.fasta.1/COG_sub2.hmm",  FASTA + ".gene.faa"] ###COG_sub2  ###"-E", "100"
            std_out2 = subprocess.check_output(cmd2)
            print "std_out2_________ " + std_out2
        '''
        cmd = "/usr/bin/rpsblast -m 8 -d %s -i %s.gene.faa -o %s.gene.faa.out -b 1 -a 1 -e 1e-2 -l %s.gene.faa.log" % (Cog, FASTA, FASTA, FASTA)
        std_out, std_err, exit_code = run_command(cmd, True, log)
        if exit_code != 0:
            print "CMD %s failed with OUT %s ERR %s" % (cmd, std_out, std_err)
            return
        '''
            
        '''
        std_out = subprocess.check_output(cmd4)
        cmd2 = ["/usr/common/jgi/utilities/bbtools/prod-v35.50/bin/commonkmers.sh", "in=" + str(seqin), "out=stdout", "k=2", "display=3", "count"]
            ###print "CMD " + str(cmd2)
            ###run
            ###p = Popen(cmd2)
            ###std_out, std_err = p.communicate()
            std_out = subprocess.check_output(cmd2)
            ###read first line from stdout
            ###print "STDOUT " + std_out
            ###print std_err
            ###for line in std_out:
            ###    print "LINE " + line
            vals = std_out.split()[1].split("=")
            ###print "VALS " + str(vals)
            max_occur_dimer = vals[0]
            max_occur_same_dimer = vals[1]
        '''
        
        ###save the second_column objects only in hits
        ###For each gene 'contig' there should be just one COG hit with the highest bit score
        ###This assume the same COGs for a specific gene are consecutive in order.
        hits = {}
        if os.path.getsize(FASTA + ".gene.faa") > 0:
            tbloutFile = open(FASTA + ".COG.hmm.hmmsearch.domtblout.txt", "r")
            for s in tbloutFile:
                if len(s) > 1 and not s.startswith("#"):
                    ###hits.append(s.split()[2])
                    s_array = s.split()
                    contig = s_array[0]
                    cog_hit = s_array[3]
                    cog_len = int(s_array[5])
                    bit_score = float(s_array[7])
                    hit_start = int(s_array[19])
                    hit_end = int(s_array[20])
                    hit_len = float(hit_end - hit_start) / cog_len
                    if not contig in hits:
                        hits[contig] = [cog_hit, bit_score, hit_len]
                    else:
                        cog_bit_len = hits.get(contig)
                        if cog_bit_len[0] == cog_hit:
                                cog_bit_len[2] += hit_len
                                hits[contig] = cog_bit_len
                        elif bit_score > cog_bit_len[1]:
                                hits[contig] = [cog_hit, bit_score, hit_len]
        
        print "HITS___________ " + str(hits)
        
        ###check if hits contain the COGs unique to chromosomes - if yes set the corresponding bit field to 1.
        ###The same COGs may be encounterd for multiple genes in the fille, so keep one with highest hit_len.
        results = {}
        hits_summary = {}
        for k, v in hits.items():
                cog_hit = v[0]
                bit_score = v[1]
                hit_len = v[2]
                if cog_hit in hits_summary:
                        if hits_summary[cog_hit] < hit_len:
                             hits_summary[cog_hit] = hit_len
                else:
                        hits_summary[cog_hit] = hit_len
                        
        for i in COGs:
            if i in hits_summary and hits_summary.get(i) >= 0.80:
                results[i] = "1"
            else:
                results[i] = "0"

        os.remove(FASTA + ".prodigal.out")
        ###os.remove(FASTA + ".gene.faa.rpsblast.log")
        if os.path.getsize(FASTA + ".gene.faa") > 0:
            os.remove(FASTA + ".COG.hmm.hmmsearch.domtblout.txt")
        os.remove(FASTA + ".gene.faa")
        os.remove(FASTA + ".gene.fasta")

        
        return results, genesperMB, genecount


'''
Print all of the fatures to files for machine interpretation
'''
def print_read_features(output_path, id_run, gc_content, mingc, maxgc, longestHomopol, totalLongHomopol, longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, longest_revcompl_repeat_s2, occur_longest_revcompl_repeat_s2, len_sequence, longest_id_ecoli, longest_id_vector, energy, max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, \
                                    max_occur_pentamer, max_occur_pentamer_1hot, max_occur_same_pentamer, percent_occur_same_pentamer, \
            max_occur_hexamer, max_occur_same_hexamer, percent_occur_same_hexamer, \
            max_occur_heptamer, max_occur_same_heptamer, percent_occur_same_heptamer, \
            max_occur_octamer, max_occur_same_octamer, percent_occur_same_octamer, \
            max_occur_ninemer, max_occur_same_ninemer, percent_occur_same_ninemer, \
            max_occur_dekamer, max_occur_same_dekamer, percent_occur_same_dekamer, coghits, genesperMB, genecount, header):
    
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
        print "longestHomopolSum " + str(longestHomopolSum)
    feature_index += 1
    line_excel += " " + str(float(longestHomopolSum) ) ###/ float(len_sequence) )
    yml_dict["longestHomopolSum"] = float(longestHomopolSum)

    if totalLongHomopolSum > 0:
        line += " " + str(feature_index) + ":" + str(float(totalLongHomopolSum) ) ###/ float(len_sequence) )
    if DEBUG == 1:
        print "totalLongHomopolSum " + str(totalLongHomopolSum)
    feature_index += 1
    line_excel += " " + str(float(totalLongHomopolSum) ) ###/ float(len_sequence) )
    yml_dict["totalLongHomopolSum"] = float(totalLongHomopolSum)

    ###The repeat featurs are not used atm (too time-consuming to compute)
    if len(longest_repeat) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_repeat))
    if DEBUG == 1:
        print "longest_repeat " + longest_repeat
    feature_index += 1
    line_excel += " " + str(len(longest_repeat))
    yml_dict["longest_repeat"] = len(longest_repeat)
    
    if occur_longest_repeat > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_repeat)
    if DEBUG == 1:
        print "occur_longest_repeat " + str(occur_longest_repeat)
    feature_index += 1
    line_excel += " " + str(occur_longest_repeat)
    yml_dict["occur_longest_repeat"] = occur_longest_repeat

    if len(longest_rev_repeat) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_rev_repeat))
    if DEBUG == 1:
        print "longest_rev_repeat " + longest_rev_repeat
    feature_index += 1
    line_excel += " " + str(len(longest_rev_repeat))
    yml_dict["longest_rev_repeat"] = len(longest_rev_repeat)
    
    if occur_longest_rev_repeat > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_rev_repeat)
    if DEBUG == 1:
        print "occur_longest_rev_repeat " + str(occur_longest_rev_repeat)
    feature_index += 1
    line_excel += " " + str(occur_longest_rev_repeat)
    yml_dict["occur_longest_rev_repeat"] = occur_longest_rev_repeat


    if len(longest_revcompl_repeat) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_revcompl_repeat))
    if DEBUG == 1:
        print "longest_revcompl_repeat " + longest_revcompl_repeat
    feature_index += 1
    line_excel += " " + str(len(longest_revcompl_repeat))
    yml_dict["longest_revcompl_repeat"] = len(longest_revcompl_repeat)

    
    if occur_longest_revcompl_repeat > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_revcompl_repeat)
    if DEBUG == 1:
        print "occur_longest_revcompl_repeat " + str(occur_longest_revcompl_repeat)
    feature_index += 1
    line_excel += " " + str(occur_longest_revcompl_repeat)
    yml_dict["occur_longest_revcompl_repeat"] = occur_longest_revcompl_repeat



    if len(longest_revcompl_repeat_s2) > 0:
        line += " " + str(feature_index) + ":" + str(len(longest_revcompl_repeat_s2))
    if DEBUG == 1:
        print "longest_revcompl_repeat_s2 " + longest_revcompl_repeat_s2
    feature_index += 1
    line_excel += " " + str(len(longest_revcompl_repeat_s2))
    yml_dict["longest_revcompl_repeat_s2"] = len(longest_revcompl_repeat_s2)
    
    
    if occur_longest_revcompl_repeat_s2 > 0:
        line += " " + str(feature_index) + ":" + str(occur_longest_revcompl_repeat_s2)
    if DEBUG == 1:
        print "occur_longest_revcompl_repeat_s2 " + str(occur_longest_revcompl_repeat_s2)
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
    yml_dict["max_occur_dimer"] = str(max_occur_dimer)

    if max_occur_same_dimer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_dimer)
    if DEBUG == 1:
        print "max_occur_same_dimer " + str(max_occur_same_dimer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_dimer)
    yml_dict["max_occur_same_dimer"] = str(max_occur_same_dimer)

    if percent_occur_same_dimer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_dimer)
    if DEBUG == 1:
        print "percent_occur_same_dimer " + str(percent_occur_same_dimer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_dimer)
    yml_dict["percent_occur_same_dimer"] = str(percent_occur_same_dimer)

    if not max_occur_trimer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_trimer)
    if DEBUG == 1:
        print "max_occur_trimer " + str(max_occur_trimer)
    feature_index += 1
    line_excel += " " + str(max_occur_trimer)
    yml_dict["max_occur_trimer"] = str(max_occur_trimer)

    if max_occur_same_trimer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_trimer)
    if DEBUG == 1:
        print "max_occur_same_trimer " + str(max_occur_same_trimer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_trimer)
    yml_dict["max_occur_same_trimer"] = str(max_occur_same_trimer)

    if percent_occur_same_trimer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_trimer)
    if DEBUG == 1:
        print "percent_occur_same_trimer " + str(percent_occur_same_trimer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_trimer)
    yml_dict["percent_occur_same_trimer"] = str(percent_occur_same_trimer)

    if not max_occur_tetramer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_tetramer)
    if DEBUG == 1:
        print "max_occur_tetramer " + str(max_occur_tetramer)
    feature_index += 1
    line_excel += " " + str(max_occur_tetramer)
    yml_dict["max_occur_tetramer"] = str(max_occur_tetramer)

    if max_occur_same_tetramer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_tetramer)
    if DEBUG == 1:
        print "max_occur_same_tetramer " + str(max_occur_same_tetramer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_tetramer)
    yml_dict["max_occur_same_tetramer"] = str(max_occur_same_tetramer)

    if percent_occur_same_tetramer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_tetramer)
    if DEBUG == 1:
        print "percent_occur_same_tetramer " + str(percent_occur_same_tetramer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_tetramer)
    yml_dict["percent_occur_same_tetramer"] = str(percent_occur_same_tetramer)


    if not max_occur_pentamer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_pentamer)
    if DEBUG == 1:
        print "max_occur_pentamer " + str(max_occur_pentamer)
    feature_index += 1
    line_excel += " " + str(max_occur_pentamer)
    yml_dict["max_occur_pentamer"] = str(max_occur_pentamer)

    if not max_occur_pentamer_1hot is None:
        line += " " + str(feature_index) + ":" + str(max_occur_pentamer_1hot)
    if DEBUG == 1:
        print "max_occur_pentamer_1hot " + str(max_occur_pentamer_1hot)
    feature_index += 1
    line_excel += " " + str(max_occur_pentamer_1hot)
    yml_dict["max_occur_pentamer_1hot"] = str(max_occur_pentamer_1hot)


    if max_occur_same_pentamer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_pentamer)
    if DEBUG == 1:
        print "max_occur_same_pentamer " + str(max_occur_same_pentamer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_pentamer)
    yml_dict["max_occur_same_pentamer"] = str(max_occur_same_pentamer)

    if percent_occur_same_pentamer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_pentamer)
    if DEBUG == 1:
        print "percent_occur_same_pentamer " + str(percent_occur_same_pentamer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_pentamer)
    yml_dict["percent_occur_same_pentamer"] = str(percent_occur_same_pentamer)


    if not max_occur_hexamer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_hexamer)
    if DEBUG == 1:
        print "max_occur_hexamer " + str(max_occur_hexamer)
    feature_index += 1
    line_excel += " " + str(max_occur_hexamer)
    yml_dict["max_occur_hexamer"] = str(max_occur_hexamer)

    if max_occur_same_hexamer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_hexamer)
    if DEBUG == 1:
        print "max_occur_same_hexamer " + str(max_occur_same_hexamer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_hexamer)
    yml_dict["max_occur_same_hexamer"] = str(max_occur_same_hexamer)

    if percent_occur_same_hexamer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_hexamer)
    if DEBUG == 1:
        print "percent_occur_same_hexamer " + str(percent_occur_same_hexamer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_hexamer)
    yml_dict["percent_occur_same_hexamer"] = str(percent_occur_same_hexamer)


    if not max_occur_heptamer is None:
        line += " " + str(feature_index) + ":" + str(max_occur_heptamer)
    if DEBUG == 1:
        print "max_occur_heptamer " + str(max_occur_heptamer)
    feature_index += 1
    line_excel += " " + str(max_occur_heptamer)
    yml_dict["max_occur_heptamer"] = str(max_occur_heptamer)

    if max_occur_same_heptamer > 0:
        line += " " + str(feature_index) + ":" + str(max_occur_same_heptamer)
    if DEBUG == 1:
        print "max_occur_same_heptamer " + str(max_occur_same_heptamer)
    feature_index += 1
    line_excel += " " + str(max_occur_same_heptamer)
    yml_dict["max_occur_same_heptamer"] = str(max_occur_same_heptamer)

    if percent_occur_same_heptamer > 0:
        line += " " + str(feature_index) + ":" + str(percent_occur_same_heptamer)
    if DEBUG == 1:
        print "percent_occur_same_heptamer " + str(percent_occur_same_heptamer)
    feature_index += 1
    line_excel += " " + str(percent_occur_same_heptamer)
    yml_dict["percent_occur_same_heptamer"] = str(percent_occur_same_heptamer)

    if coghits:
        for i in COGs:
            line += " " + str(feature_index) + ":" + str(coghits.get(i))
            if DEBUG == 1:
                print "COG " + i + " " + str(coghits.get(i))
            feature_index += 1
            line_excel += " " + str(coghits.get(i))
            yml_dict["COG_"+i] = str(coghits.get(i))
            
        line += " " + str(feature_index) + ":" + str(genesperMB)
        if DEBUG == 1:
            print "genesperMB " + str(genesperMB)
        feature_index += 1
        line_excel += " " + str(genesperMB)
        yml_dict["genesperMB"] = str(genesperMB)
        
        line += " " + str(feature_index) + ":" + str(genecount)
        if DEBUG == 1:
            print "genecount " + str(genecount)
        feature_index += 1
        line_excel += " " + str(genecount)
        yml_dict["genecount"] = str(genecount)

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
        print "line " + line
    line_excel += " " + str(id_run) + " >" + str(header)
    yml_dict["id_run"] = id_run
    yml_dict["header"] = header
    yml_dict2 = {}
    yml_dict2["sequence"] = yml_dict
    
    print line_excel
    
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




#............................
# convert sequences to 1-hot 2-D arrays
#............................
def encode1hotBase(seqPad):
        hot2D=[]
        for y in seqPad:
            hot2D.append(Constants.oneHotBase[y])
        return hot2D #np.array(hot2D).astype(np.float32)


###############################
# RQC's runcommand function
###############################
def runCommand(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.strip(), stderr.strip(), process.returncode


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

    cmd = ["shifter", "--image=bryce911/bbtools", "commonkmers.sh", "in=" + str(seqin), "out=stdout", "k=5", "display=3", "count"]
    std_out = subprocess.check_output(cmd)
    vals = std_out.split("\t")[1].split("=")
    max_occur_pentamer = vals[0]
    max_occur_same_pentamer = vals[1]

    percent_occur_same_pentamer = float(pos_occur_same_pentamer) / float(len(seq))

    if DEBUG == 1:
        print "max_occur_pentamer " + max_occur_pentamer + " max_occur_same_pentamer " + str(max_occur_same_pentamer) + " percent_occur_same_pentamer " + str(percent_occur_same_pentamer)

    #Convert max_occur_pentamer to one-hot encoding
    max_occur_pentamer_1hot = encode1hotBase(max_occur_pentamer)

    return len(seq), max_occur_pentamer, max_occur_pentamer_1hot, max_occur_same_pentamer, percent_occur_same_pentamer, penalty_value



'''
Unused at the moment
'''
def fivesixFindMers_BBtoolsCommonKmers(seq, seqin, penalty_value):
    '''
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
    '''
    for i in range(0, lmers_dimer.n_entries()):
       n = lmers_dimer.get(i)
       if n > max_occur_same_dimer and n > 1:
          if DEBUG == 1:  print lmers_dimer.reverse_hash(i), "is present", n, "times."
          max_occur_dimer = lmers_dimer.reverse_hash(i)
          max_occur_same_dimer = n
    '''
    cmd2 = ["/usr/common/jgi/utilities/bbtools/prod-v35.50/bin/commonkmers.sh", "in=" + str(seqin), "out=stdout", "k=2", "display=3", "count"]
    ###print "CMD " + str(cmd2)
    ###run
    ###p = Popen(cmd2)
    ###std_out, std_err = p.communicate()
    std_out = subprocess.check_output(cmd2)
    ###read first line from stdout
    ###print "STDOUT " + std_out
    ###print std_err
    ###for line in std_out:
    ###    print "LINE " + line
    vals = std_out.split()[1].split("=")
    ###print "VALS " + str(vals)
    max_occur_dimer = vals[0]
    max_occur_same_dimer = vals[1]

          
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
    '''
    for i in range(0, lmers_trimer.n_entries()):
       n = lmers_trimer.get(i)
       if n > max_occur_same_trimer and n > 1:
          if DEBUG == 1:  print lmers_trimer.reverse_hash(i), "is present", n, "times."
          max_occur_trimer = lmers_trimer.reverse_hash(i)
          max_occur_same_trimer = n
    '''
    cmd3 = ["/usr/common/jgi/utilities/bbtools/prod-v35.50/bin/commonkmers.sh", "in=" + str(seqin), "out=stdout", "k=3", "display=3", "count"]
    ###print "CMD " + str(cmd3)
    ###run
    ###p = Popen(cmd3)
    ###std_out, std_err = p.communicate()
    std_out = subprocess.check_output(cmd3)
    ###read first line from stdout
    ###print "STDOUT " + std_out
    ###print std_err
    ###for line in std_out:
    ###    print "LINE " + line
    vals = std_out.split()[1].split("=")
    ###print "VALS " + str(vals)
    max_occur_trimer = vals[0]
    max_occur_same_trimer = vals[1]

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
    '''
    for i in range(0, lmers_tetramer.n_entries()):
       n = lmers_tetramer.get(i)
       if n > max_occur_same_tetramer and n > 1:
          if DEBUG == 1:  print lmers_tetramer.reverse_hash(i), "is present", n, "times."
          max_occur_tetramer = lmers_tetramer.reverse_hash(i)
          max_occur_same_tetramer = n
    '''
    cmd4 = ["/usr/common/jgi/utilities/bbtools/prod-v35.50/bin/commonkmers.sh", "in=" + str(seqin), "out=stdout", "k=4", "display=3", "count"]
    ###print "CMD " + str(cmd4)
    ###run
    ###p = Popen(cmd4)  ###stdout = PIPE, stderr = PIPE
    ###std_out, std_err = p.communicate()
    std_out = subprocess.check_output(cmd4)
    ###read first line from stdout
    ###print "STDOUT " + std_out
    ###print std_err
    ###for line in std_out:
    ###    print "LINE " + line
    vals = std_out.split()[1].split("=")
    ###print "VALS " + str(vals)
    max_occur_tetramer = vals[0]
    max_occur_same_tetramer = vals[1]

    percent_occur_same_dimer = float(pos_occur_same_dimer) / float(len(seq))
    percent_occur_same_trimer = float(pos_occur_same_trimer) / float(len(seq))
    percent_occur_same_tetramer = float(pos_occur_same_tetramer) / float(len(seq))

    if DEBUG == 1:
        print "max_occur_dimer " + max_occur_dimer + " max_occur_same_dimer " + str(max_occur_same_dimer) + " percent_occur_same_dimer " + str(percent_occur_same_dimer)
        print "max_occur_trimer " + max_occur_trimer + " max_occur_same_trimer " + str(max_occur_same_trimer) + " percent_occur_same_trimer " + str(percent_occur_same_trimer)
        print "max_occur_tetramer " + max_occur_tetramer + " max_occur_same_tetramer " + str(max_occur_same_tetramer) + " percent_occur_same_tetramer " + str(percent_occur_same_tetramer)

    return len(seq), max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, penalty_value


'''
Unused at the moment
'''
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
    '''
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
    '''




'''
Unused at the moment
'''
def fivesixFindRepeats(seq, penalty_value):
    return "", 0, "", 0, "", 0, "", 0, len(seq), penalty_value
    '''
    ###Use http://tandem.bu.edu/trf/trf.html
    ###No, things that look for repeats in human genome are prob looking for long sequences
    ###whereas our entire sequence is shorter, in comparision to human
    ###if you hash all lenght l mers and count the number of entries with count>1,
    ###that's an easy way to measure repeats (exact repeats)
    
    ###These are the repeats on the same strand. We return the longest lmer (sorted alphbetically to break ties) and the most frequent lmer with the frequencies.
    ###For the longest lmer and the most frequent lmer we return the sequences and frequencies but not the length. Print out lengths and frequencies for both seqs.
    lmers = {}
    ###lmers_indexes = {}
    ####revcompl_lmers = {}
    revcompls = []
    ###Palindrome on both strands matching anywhere. E.g. AAGCTT-TTCGAA but at non-oppposite positions of the two strands. We return the longest one (sorted alphbetically to break ties).
    ###For the longest palindrome we return the sequence and the length but not the frequencies. Print out length only.
    palindrome2strands_anywhere = {}
    ###Palindrome on both strands matching at same positions. E.g. AAGCTT-TTCGAA at the same positions (opposite) of the two strands. We return the longest one (sorted alphbetically to break ties).
    ###For the longest palindrome we return the sequence and the length but not the frequencies. Print out length only.
    palindrome2strands = {}
    min_l = 15
    max_l = 31
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    startpointingtime = time()
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
            #if subseq == rev_subseq:
            #    palindrome1strand[subseq] = palindrome1strand.get(subseq, 0) + 1
            bases = list(rev_subseq) 
            bases = [complement[base] for base in bases]
            revcompl_subseq = ''.join(bases)
            #for i in range(len(rev_subseq)):
            #    if rev_subseq[i] == 'A':
            #        revcompl_subseq += 'T'
            #    elif rev_subseq[i] == 'T':
            #        revcompl_subseq += 'A'
            #    elif rev_subseq[i] == 'C':
            #        revcompl_subseq += 'G'
            #    elif rev_subseq[i] == 'G':
            #        revcompl_subseq += 'C'
            if subseq == revcompl_subseq:
                palindrome2strands[subseq] = palindrome2strands.get(subseq, 0) + 1
            ###if lmers.has_key(revcompl_subseq):
            revcompls.append(revcompl_subseq)
            idx += 1
    runtime = str(time() - startpointingtime)
    print "RUNTIME fivesixFindRepeatsA: " + runtime

    ###revcomplset = set(revcompls)
    ###lmerset = set(lmers.keys())
    ###lmerrevcomplint = lmerset & revcomplset
    startpointingtime = time()
    for subseq in revcompls:
        ###if seq.find(revcompl_subseq) > -1:
        ###Store the subseq so we do not have to recompute the revcompl when searching for the longest revcompl
        if subseq in lmers:
            palindrome2strands_anywhere[subseq] = palindrome2strands_anywhere.get(subseq, 0) + 1
        ###revcompl_lmers[subseq] = 1 ###revcompl_lmers.get(subseq, 0) + 1
        ###lmers_indexes[subseq] = lmers_indexes.get(subseq, []) + range(idx , idx + windowSizeL)
    runtime = str(time() - startpointingtime)
    print "RUNTIME fivesixFindRepeatsB: " + runtime

    longest_lmer = ""
    most_freq_lmer = ""
    most_freq_lmer_count = 1
    startpointingtime = time()
    for l in lmers.keys():
        if len(l) > len(longest_lmer) and lmers[l] > 1:
            longest_lmer = l
        elif len(l) == len(longest_lmer) and lmers[l] > 1:
            longest_lmer = sorted([l, longest_lmer])[0]
        if lmers[l] > most_freq_lmer_count:
            most_freq_lmer = l
            most_freq_lmer_count = lmers[l]
    runtime = str(time() - startpointingtime)
    print "RUNTIME fivesixFindRepeatsC: " + runtime
    
    longest_revcompl_lmer_1s = ""
    startpointingtime = time()
    for l in palindrome2strands_anywhere.keys():
        if len(l) > len(longest_revcompl_lmer_1s): ### and revcompl_lmers.get(l, 0) > 0:
            longest_revcompl_lmer_1s = l
        elif len(l) == len(longest_revcompl_lmer_1s):
            longest_revcompl_lmer_1s = sorted([l, longest_revcompl_lmer_1s])[0]
    runtime = str(time() - startpointingtime)
    print "RUNTIME fivesixFindRepeatsD: " + runtime

    longest_revcompl_lmer_2s = ""
    startpointingtime = time()
    for l in palindrome2strands.keys():
        if len(l) > len(longest_revcompl_lmer_2s): ### and revcompl_lmers.get(l, 0) > 0:
            longest_revcompl_lmer_2s = l
        elif len(l) == len(longest_revcompl_lmer_2s):
            longest_revcompl_lmer_2s = sorted([l, longest_revcompl_lmer_2s])[0]
    runtime = str(time() - startpointingtime)
    print "RUNTIME fivesixFindRepeatsE: " + runtime

    if not longest_lmer == "":
        if DEBUG == 1:
            print "longest_lmer " + longest_lmer + " lmers " + str(lmers[longest_lmer])
            print "penalty_value " + str(penalty_value)

    return longest_lmer, lmers.get(longest_lmer, 0), most_freq_lmer, most_freq_lmer_count, longest_revcompl_lmer_1s, len(longest_revcompl_lmer_1s), longest_revcompl_lmer_2s, len(longest_revcompl_lmer_2s), len(seq), penalty_value
    '''


'''
Unused at the moment, run blast against a plasmid database
'''
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



#def classifySVM:
    ###  ./svm_learn seq_features seq_features.model
    ###  ./svm_classify seq_features seq_features.model seq_features.predict


def process_seq(seqin, sequence, header, penalty_value, output_path, id_run, run_blast = None):
                    if DEBUG == 1:
                        print "\n\n_____________________________________\n  -------> Computing features for fasta header : " + header + "\n"
                    
                    #TRIMREADS.write(line + "\n" + sequence + "\n")
                    if DEBUG == 1:
                        print "sequence %s" % ( sequence )
                        print "length of sequence %s" % ( len(sequence) )

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
                    #len_sequence, max_occur_pentamer, max_occur_pentamer_1hot, max_occur_same_pentamer, percent_occur_same_pentamer, penalty_value = fivesixFindPentamer(sequence, seqin, penalty_value)
                    runtime = str(time() - startpointingtime)
                    #print "RUNTIME fivesixFindMers: " + runtime
                    #if DEBUG == 1:
                    #    print "fivesixFindMers"

                    startpointingtime = time()
                    coghits = {}
                    genesperMB = 0
                    genecount = 0
                    if USE_COGs:
                       coghits, genesperMB, genecount = cogs(sequence, seqin, penalty_value);
                    runtime = str(time() - startpointingtime)
                    ###print "RUNTIME cogs: " + runtime
                    #if DEBUG == 1:
                    #    print "cogs"

                        
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
                        
                        
                    print_read_features(output_path, id_run, gc_content, mingc, maxgc, longestHomopol, totalLongHomopol, longest_repeat, occur_longest_repeat, longest_rev_repeat, occur_longest_rev_repeat, longest_revcompl_repeat, occur_longest_revcompl_repeat, longest_revcompl_repeat_s2, occur_longest_revcompl_repeat_s2, len_sequence, longest_id_ecoli, longest_id_vector, energy, max_occur_dimer, max_occur_same_dimer, percent_occur_same_dimer, max_occur_trimer, max_occur_same_trimer, percent_occur_same_trimer, max_occur_tetramer, max_occur_same_tetramer, percent_occur_same_tetramer, \
                                                    max_occur_pentamer, max_occur_pentamer_1hot, max_occur_same_pentamer, percent_occur_same_pentamer, \
            max_occur_hexamer, max_occur_same_hexamer, percent_occur_same_hexamer, \
            max_occur_heptamer, max_occur_same_heptamer, percent_occur_same_heptamer, \
            max_occur_octamer, max_occur_same_octamer, percent_occur_same_octamer, \
            max_occur_ninemer, max_occur_same_ninemer, percent_occur_same_ninemer, \
            max_occur_dekamer, max_occur_same_dekamer, percent_occur_same_dekamer, coghits, genesperMB, genecount, header);
                    
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
                        
                subprocess.check_call(cmd)

        except subprocess.CalledProcessError as e:
                result = "Failed to run job (%s)." % (e.output) ###.message   ###"+str(pipeline.paramvalues[2])+"
                print result
                ###p.kill()
                # store the result
                result_queue.put(result)

                return False

        result = 'Job finished with success.' ###% exit_code  ###'+str(pipeline.paramvalues[2])+'
        print result
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

DEBUG = 0
USE_COGs = 0

if __name__ == "__main__":
    ###penalty_values = {}
    all_seqs = 0
    penalty_value = 0
    datetime = create_timestamp()
    ###id_run = None
    id_run = "DATETIME" + datetime

    desc = 'fasta_features'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-i", "--input-path", dest="input_path", help = "Input path to write to", required=False)
    parser.add_argument("-o", "--output-path", dest="output_path", help = "Output path to write to", required=False)
    parser.add_argument("-he", "--header", dest="header", help = "The header to use", required=False)
    parser.add_argument("-s", "--seqin", dest="seqin", help = "The sequence to use", required=False)
    parser.add_argument("-b", "-run-blast", dest="run_blast", help = "Specify either p for refseq.plasmid or m for refseq.mito", required=False)
    
    options = parser.parse_args()

    ###if (options.input_path == None) != (options.output_path == None):
    ###    print "Either specify both input and output path or neither"
    ###    exit(1)

    output_path = None
    input_path = []
    run_blast = None
    header = "none"
    sequence = ""

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

    if options.header:
        header = options.header

    #we are processing either the fasta with all sequences (-i) or a single sequence (-s and -he)
    if options.seqin and options.input_path:
        print("Please specify either seqin(-s) or input-path(-i) but not both")
        exit(2)
    
    if options.seqin and not options.header:
        print("If using seqin(-s) you may want to specify a header(-he)")
        
    if options.seqin:
        if not os.path.exists(options.seqin):
            print("The dataset path does not exist: " + options.seqin)
            exit(2)
        if not os.path.isfile(options.seqin):
            print("The dataset file does not exist, is missing or is not readable: " + options.seqin)
            exit(2)
        seqfilein = open(options.seqin, 'r')
        for line in seqfilein: ###fileinput.input( [ options.seqin ] ): ###READS####READS: ###:
            ###line = CONTIGS.readline()
            ####sequence = options.seqin
            if not line.startswith(">"):
                sequence += line.strip()
        ###Remove file options.seqin
        seqfilein.close()

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
        #TODO verify number of contigs in input_path is <10k. This is file is only for contigs, not reads!
        input_path = [ input_path ]
        
    
    colheaders = "id_run,gc_content,mingc,maxgc,AlongestHomopol,AtotalLongHomopol,ClongestHomopol,CtotalLongHomopol,GlongestHomopol,GtotalLongHomopol,TlongestHomopol,TtotalLongHomopol,longestHomopolSum,totalLongHomopolSum,longest_repeat,occur_longest_repeat,longest_rev_repeat,occur_longest_rev_repeat,longest_revcompl_repeat,occur_longest_revcompl_repeat,longest_revcompl_repeat_s2,occur_longest_revcompl_repeat_s2,len_sequence,longest_id_ecoli,longest_id_vector,energy,max_occur_dimer,max_occur_same_dimer,percent_occur_same_dimer,max_occur_trimer,max_occur_same_trimer,percent_occur_same_trimer,max_occur_tetramer,max_occur_same_tetramer,percent_occur_same_tetramer,max_occur_pentamer,max_occur_pentamer_1hot,max_occur_same_pentamer,percent_occur_same_pentamer,max_occur_hexamer,max_occur_same_hexamer,percent_occur_same_hexamer,max_occur_heptamer,max_occur_same_heptamer,percent_occur_same_heptamer,COGS,genesperMB,genecount"
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
        for line in fileinput.input( input_path ): ###READS####READS: ###:
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
                    seq_file = open(seq_filename, 'w')
                    seq_file.write(header + "\n")
                    seq_file.write(sequence)
                    seq_file.flush()
                    os.fsync(seq_file.fileno())
                    seq_file.close()
                    seq_count += 1
                    paramvals = [ str(output_path), str(header).replace(">",""), str(seq_filename) ]
                    rqc_command = os.path.realpath(__file__)  ###os.path.join(os.path.dirname(__file__), "read_fastaB.py") ###os.path.realpath(__file__)  ###os.path.join(rqc_command_dir , 'read_fasta.py')
                    p = JGI_Pipeline(rqc_command, pars, paramvals, sequence, header)
                    work_queue.put(p)
                    jobs_to_run.append(p)
                    num_jobs += 1
                    if num_jobs % 100 == 0:
                        print "job number processed: %s " % (num_jobs)
                    #num_processes += 1

                sequence = "";
                header = line   ###.replace(" ", "_").replace(",", ".");
                if DEBUG == 1: print "header %s" % (header)
                all_seqs += 1
                penalty_value = 0
            elif len(line) > 2:
                line = line.upper()
                ###line = re.sub('[^ATCG]', '', line)  ###Removed because we want Ns to be counted in homopolymer function
                sequence += line
        ###rqc_command_dir = os.getcwd() ###os.path.join(os.path.abspath(os.pardir) , "categorical" ) ###current_dir
        pars = ['-o', '-he', '-s']
        seq_filename = os.path.join(output_path, str(seq_count) + '.fa')
        seq_file = open(seq_filename, 'w')
        seq_file.write(header + "\n")
        seq_file.write(sequence)
        seq_file.flush()
        os.fsync(seq_file.fileno())
        seq_file.close()
        seq_count += 1
        paramvals = [ str(output_path), str(header).replace(">",""), str(seq_filename) ]
        rqc_command = os.path.realpath(__file__)  ###os.path.join(os.path.dirname(__file__), "read_fastaB.py") ###os.path.realpath(__file__)  ###os.path.join(rqc_command_dir , 'read_fasta.py')
        p = JGI_Pipeline(rqc_command, pars, paramvals, sequence, header)
        work_queue.put(p)
        jobs_to_run.append(p)
        num_jobs += 1
        if num_jobs % 100 == 0:
            print "job number processed: %s " % (num_jobs)
        #num_processes += 1
    else:
        ###print "\n\nSEQUENCE %s\n\n" % sequence
        ###if len(input_path) <1:
        ###print "penalty_value = process_seq(options.seqin %s, sequence %s, header %s, penalty_value %s, output_path %s, id_run %s, run_blast %s)" % (options.seqin, sequence, header, penalty_value, output_path, id_run, run_blast)
        penalty_value = process_seq(options.seqin, sequence, header, penalty_value, output_path, id_run, run_blast)
        os.remove(options.seqin)

    #
    #
    ###MULTI-processing part - this is run only from the input fasta process.
    #
    #
    if len(input_path) >0:
        print 'cpu_count() = %d\n' % multiprocessing.cpu_count()
        #
        # Create pool
        #
        print 'Creating pool with %d processes\n' % Constants.PROCESSES
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
        print '\tlist1(pool.imap(multiproc_pool, jobs_to_run(%d), chunksize=%d)):\n\t\t%s' \
            ' seconds' % (N, N//(Constants.PROCESSES), time() - t)
        # This is the multiproc pool and jobs_to_run is the list with all jobs that need to be run.
        C = list(pool.imap(multiproc_pool, jobs_to_run, chunksize=max(2, N//(Constants.PROCESSES))))
        print '\tlist2(pool.imap(multiproc_pool, jobs_to_run(%d), chunksize=%d)):\n\t\t%s' \
            ' seconds' % (N, N//(Constants.PROCESSES), time() - t)
        #print(result_list)
        while not result_queue.empty():
            print("Status of a done job: " + result_queue.get())
            
            
    runtime = str(time() - startpointingtimetotal)
    print "TOTAL_RUNTIME: " + runtime
    #runtimes.write(runtime + "\n")
    print "exiting........"
    ###sys.exit(0)
    os._exit(0)

