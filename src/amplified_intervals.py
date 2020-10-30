#!/usr/bin/env python2

# This software is Copyright 2017 The Regents of the University of California. All Rights Reserved. Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies. Permission to make commercial use of this software may be obtained by contacting:
#
# Office of Innovation and Commercialization
#
# University of California
#
# La Jolla, CA 92093-0910
#
# (858) 534-5815
#
# invent@ucsd.edu
#
# This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
#
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

#Author: Viraj Deshpande
#Contact: virajbdeshpande@gmail.com


import copy
from collections import defaultdict
import sys
from sets import Set
import numpy as np
import re
sys.setrecursionlimit(10000)
import argparse
import os
import pysam
import global_names

if sys.version_info >= (3,0):
    sys.stderr.write("AA must be run with python2. Python3 support is under development.\n")
    sys.exit(1)

GAIN = 5.0
CNSIZE_MIN = 100000

parser = argparse.\
ArgumentParser(description="Filter and merge amplified intervals")
parser.add_argument('--bed', dest='bed',
                    help="Input bed file with list of amplified intervals", metavar='FILE',
                    action='store', type=str, required=True)
parser.add_argument('--out', dest='out',
                    help="OPTIONAL: Prefix filename for output bed file. Default: <INPUT_BED_BASENAME>_amplified.bed", metavar='FILE',
                    action='store', type=str, default='')
parser.add_argument('--bam', dest='bam',
                    help="OPTIONAL: Bamfile, used to avoid large aneuploidies", metavar='FILE',
                    action='store', type=str, default='')
parser.add_argument('--gain', dest='gain',
                    help="OPTIONAL: CN gain threshold for interval to be considered as a seed. Default: 5",
                    action='store', type=float, default=GAIN)
parser.add_argument('--cnsize_min', dest='cnsize_min',
                    help="OPTIONAL: Minimum size (in bp) for interval to be considered as a seed. Default: 100000",
                    action='store', type=int, default=CNSIZE_MIN)
parser.add_argument('--ref', dest='ref',
                    help="Values: [hg19, GRCh37, GRCh38, None]. \"hg19\"(default) & \"GRCh38\" : chr1, .. chrM etc / \"GRCh37\" : '1', '2', .. 'MT' etc/ \"None\" : Do not use any annotations. AA can tolerate additional chromosomes not stated but accuracy and annotations may be affected. Default: hg19", metavar='STR',
                    action='store', type=str, default='hg19')
args = parser.parse_args()

global_names.REF = args.ref
import hg19util as hg


if args.bed != '':
    rdAlts = args.bed

if args.out != '':
    outname= args.out + ".bed"
else:
    outname = os.path.splitext(rdAlts)[0] + "_amplified.bed"

GAIN,CNSIZE_MIN = args.gain,args.cnsize_min

rdList0 = hg.interval_list(rdAlts, 'bed')
if rdList0:
    try:
        if len(rdList0[0].info) == 0:
            sys.stderr.write("ERROR: CNV estimate bed file had too few columns.\n"
                             "Must contain: chr  pos1  pos2  cnv_estimate\n")
            sys.exit(1)
        _ = float(rdList0[0].info[-1])

    except ValueError:
        sys.stderr.write("ERROR: CNV estimates must be in last column of bed file.\n")
        sys.exit(1)

rdList = hg.interval_list([r for r in rdList0 if float(r.info[-1]) > GAIN ])

if args.bam != "":
    import bam_to_breakpoint as b2b
    if os.path.splitext(args.bam)[-1] == '.cram':
        bamFile = pysam.Samfile(args.bam, 'rc')
    else:
        bamFile = pysam.Samfile(args.bam, 'rb')
    cstats = None
    cb = bamFile
    if os.path.exists(os.path.join(hg.DATA_REPO, "coverage.stats")):
        coverage_stats_file = open(os.path.join(hg.DATA_REPO, "coverage.stats"))
        for l in coverage_stats_file:
            ll = l.strip().split()
            if ll[0] == os.path.abspath(cb.filename):
                cstats = tuple(map(float, ll[1:]))
        coverage_stats_file.close()
    bamFileb2b = b2b.bam_to_breakpoint(bamFile, coverage_stats=cstats)
    rdList = hg.interval_list([r for r in rdList if float(r.info[-1]) >
                               GAIN + 2 * max(1.0, bamFileb2b.median_coverage(refi=r)[0] / bamFileb2b.median_coverage()[0]) - 2
                               and bamFileb2b.median_coverage(refi=r)[0] / bamFileb2b.median_coverage()[0] > 0])

genome_features = hg.oncogene_list
amplicon_listl = rdList

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

uc_list = hg.interval_list([a for a in uc_list if float(a.info[-1]) * a.segdup_uniqueness() > GAIN and a.rep_content() < 2.5])
uc_merge = uc_list.merge_clusters(extend=300000)

with open(outname,"w") as outfile:
    for a in uc_merge:
        if sum([ai.size() for ai in a[1]]) > CNSIZE_MIN:
            outfile.write('\t'.join([str(a[0]), str(sum([ai.size() * float(ai.info[-1]) for ai in a[1]]) / sum([ai.size() for ai in a[1]])), rdAlts]) + '\n')

