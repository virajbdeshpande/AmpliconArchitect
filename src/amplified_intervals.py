#This file is part of AmpliconArchitect.
#AmpliconArchitect is software which can use whole genome sequencing data to reconstruct the structure of focal amplifications.
#Copyright (C) 2018 Viraj Deshpande
#
#AmpliconArchitect is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#AmpliconArchitect is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with AmpliconArchitect.  If not, see <http://www.gnu.org/licenses/>


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

import hg19util as hg

GAIN = 5
CNSIZE_MIN = 100000

parser = argparse.\
ArgumentParser(description="Filter and merge amplified intervals")
parser.add_argument('--bed', dest='bed',
                    help="Input bed file with list of amplified intervals", metavar='FILE',
                    action='store', type=str, nargs=1, default='')
parser.add_argument('--bam', dest='bam',
                    help="OPTIONAL: Bamfile, used to avoid large aneuploidies", metavar='FILE',
                    action='store', type=str, nargs=1, default='')
parser.add_argument('--bedlist', dest='bedlist',
                    help="List of bed files with amplified intervals in each sample", metavar='FILE',
                    action='store', type=str, nargs=1, default=[])
args = parser.parse_args()
rdAltsl = []
if args.bed[0] != '':
    rdAltsl.append(args.bed[0])
elif len(args.bedlist) != 0 and args.bedlist[0] != '':
    for l in open(args.bedlist[0]):
        rdAltsl.append(l.strip())


for rdAlts in rdAltsl:

    rdList0 = hg.interval_list(rdAlts, 'bed')
    rdList = hg.interval_list([r for r in rdList0 if float(r.info[1]) > GAIN ])

    if args.bam != "":
        import bam_to_breakpoint as b2b
        if os.path.splitext(args.bam[0])[-1] == '.cram':
            bamFile = pysam.Samfile(args.bam[0], 'rc')
        else:
            bamFile = pysam.Samfile(args.bam[0], 'rb')
        coverage_stats_file = open(hg.DATA_REPO + "/coverage.stats")
        cstats = None
        cb = bamFile
        for l in coverage_stats_file:
            ll = l.strip().split()
            if ll[0] == os.path.abspath(cb.filename):
                cstats = tuple(map(float, ll[1:]))
        coverage_stats_file.close()
        bamFileb2b = b2b.bam_to_breakpoint(bamFile, coverage_stats=cstats)
        rdList = hg.interval_list([r for r in rdList if float(r.info[1]) > GAIN + 2 * max(1.0, bamFileb2b.median_coverage(refi=r)[0] / bamFileb2b.median_coverage()[0]) - 2])

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

    uc_list = hg.interval_list([a for a in uc_list if float(a.info[1]) * a.segdup_uniqueness() > 5.0 and a.rep_content() < 2.5])
    uc_merge = uc_list.merge_clusters(extend=300000)
    all_uc = hg.interval_list([a[0] for a in uc_merge if sum([ai.size() for ai in a[1]]) > CNSIZE_MIN] )

    if len(args.bedlist) != 0 and args.bedlist[0] != '':
        outfile = open(rdAlts[:-4] + "_amplified.bed")
    for a in uc_merge:
        if sum([ai.size() for ai in a[1]]) > 100000:
            if len(args.bedlist) != 0 and args.bedlist[0] != '':
                outfile.write('\t'.join([str(a[0]), sum([ai.size() * float(ai.info[1]) for ai in a[1]]) / sum([ai.size() for ai in a[1]]), rdAlts]) + '\n')
            else:
                print str(a[0]), sum([ai.size() * float(ai.info[1]) for ai in a[1]]) / sum([ai.size() for ai in a[1]]), rdAlts
exit()


