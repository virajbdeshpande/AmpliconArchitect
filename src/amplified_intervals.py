# This software is Copyright 2014 The Regents of the University of
# California. All Rights Reserved.
#
# Permission to copy, modify, and distribute this software and its
# documentation for educational, research and non-profit purposes, without fee,
# and without a written agreement is hereby granted, provided that the above
# copyright notice, this paragraph and the following three paragraphs appear
# in all copies.
#
# Permission to make commercial use of this software may be obtained by
# contacting:
# Technology Transfer Office
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# invent@ucsd.edu
#
# This software program and documentation are copyrighted by The Regents of the
# University of California. The software program and documentation are supplied
# "as is", without any accompanying services from The Regents. The Regents does
# not warrant that the operation of the program will be uninterrupted or
# error-free. The end-user understands that the program was developed for
# research purposes and is advised not to rely exclusively on the program for
# any reason.
#
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO
# ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
# CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING
# OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
# EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
# CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESSFOR A PARTICULAR PURPOSE.
# THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
# CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
# ENHANCEMENTS, OR MODIFICATIONS.

#Author: Viraj Deshpande
#Contact: virajbdeshpande@gmail.com


import copy
from collections import defaultdict
import sys
from sets import Set
import numpy as np
import matplotlib.pyplot as plt
import re
sys.setrecursionlimit(10000)
import argparse

import hg19util as hg

GAIN = 5
CNSIZE_MIN = 100000

parser = argparse.\
ArgumentParser(description="Filter and merge amplified intervals")
parser.add_argument('--bed', dest='bed',
                    help="Bed file with list of amplified intervals", metavar='FILE',
                    action='store', type=str, nargs=1)
args = parser.parse_args()
rdAlts = args.bed[0]
rdList0 = hg.interval_list(rdAlts + '/output/alts.dat', 'bed')
rdList = hg.interval_list([r for r in rdList0 if float(r.info[1]) > GAIN ])

genome_features = hg.oncogene_list

amplicon_listl = rdList
amplicon_listl = hg.interval_list([a for a in amplicon_listl if a.size() > CNSIZE_MIN]) 
amplicon_listl.sort()


cr = hg.conserved_regions
uc_list = hg.interval_list([])
for a in amplicon_listl:
    if (len(hg.interval_list([a]).intersection(cr)) == 0 or
        a.size() > max(1000000, 10 * sum([a.intersection(ci[1]).size() for ci in hg.interval_list([a]).intersection(cr)])) or
       a.size() - sum([a.intersection(ci[1]).size() for ci in hg.interval_list([a]).intersection(cr)]) > 2000000):
        if (len(hg.interval_list([a]).intersection(cr))) == 0:
            uc_list.append(a)
        else:
            cra = hg.interval_list([a]).intersection(cr)
            cpos = a.start
            for crai in cra:
                if cpos < crai[1].start - 1000000:
                    uc_list.append(hg.interval(a.chrom, cpos, crai[1].start - 1000000, info=a.info))
                cpos = crai[1].end + 1000000
            if a.end > cpos:
                uc_list.append(hg.interval(a.chrom, cpos, a.end, info = a.info))

uc_list = hg.interval_list([a for a in uc_list if float(a.info[1]) * a.segdup_uniqueness() > 5.0 and a.rep_content() < 2.5])

uc_merge = uc_list.merge_clusters(extend=300000)

all_uc = hg.interval_list([a[0] for a in uc_merge if sum([ai.size() for ai in a[1]]) > 100000] )


for a in uc_merge:
    if sum([ai.size() for ai in a[1]]) > 100000:
        print str(a[0]), sum([ai.size() * float(ai.info[1]) for ai in a[1]]) / sum([ai.size() for ai in a[1]]), rdAlts
exit()


