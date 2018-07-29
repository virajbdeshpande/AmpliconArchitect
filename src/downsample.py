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


import pysam
import argparse
import math
from time import clock
from collections import defaultdict
from sets import Set
from cStringIO import StringIO
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import logging
import JobNotifier
import random
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

import global_names

parser = argparse.\
ArgumentParser(description="Reconstruct Amplicons connected to listed intervals.")
parser.add_argument('--bam', dest='bam',
                    help="Coordinate sorted BAM file with index mapped to hg19 reference sequence", metavar='FILE',
                    action='store', type=str, nargs=1)
parser.add_argument('--final', dest='final',
                    help="Optional Final coverage. Default is 10. If initial coverage is less than final, do nothing.", metavar='FLOAT',
                    action='store', type=float, default=10.0)
parser.add_argument('--downsample_dir', dest='downsample_dir',
                    help="Optional directory to output. Default is same as original bamfile", metavar='DIR',
                    action='store', type=str, default='')
parser.add_argument('--cbam', dest='cbam',
                    help="Optional bamfile to use for coverage calculation. Also generates new coverage bam file in downsample_dir.", metavar='FILE',
                    action='store', type=str, default=None)
parser.add_argument('--cbed', dest='cbed',
                    help="Optional bedfile defining 1000 10kbp genomic windows for coverage calcualtion", metavar='FILE',
                    action='store', type=str, default=None)
parser.add_argument('--ref', dest='ref',
                    help="Values: [hg19, GRCh37, None]. \"hg19\"(default) : chr1, .. chrM etc / \"GRCh37\" : '1', '2', .. 'MT' etc/ \"None\" : Do not use any annotations. AA can tolerate additional chromosomes not stated but accuracy and annotations may be affected. Default: hg19", metavar='STR',
                    action='store', type=str, default='hg19')

args = parser.parse_args()

global_names.REF = args.ref



import hg19util as hg
import bam_to_breakpoint as b2b
from breakpoint_graph import *


if os.path.splitext(args.bam[0])[-1] == '.cram':
    bamFile = pysam.Samfile(args.bam[0], 'rc')
else:
    bamFile = pysam.Samfile(args.bam[0], 'rb')
cbam = None
if args.cbam is not None:
    if os.path.splitext(args.cbam[0])[-1] == '.cram':
        cbam = pysam.Samfile(args.cbam, 'rc')
    else:
        cbam = pysam.Samfile(args.cbam, 'rb')
cbed = args.cbed




coverage_stats_file = open(hg.DATA_REPO + "/coverage.stats")
cstats = None
cb = bamFile
if cbam is not None:
    cb = cbam
for l in coverage_stats_file:
    ll = l.strip().split()
    if ll[0] == os.path.abspath(cb.filename):
        cstats = tuple(map(float, ll[1:]))
coverage_stats_file.close()
coverage_windows=None
if cbed is not None:
    coverage_windows=hg.interval_list(cbed, 'bed')
    coverage_windows.sort()
if cstats is None and cbam is not None:
    cbam2b = b2b.bam_to_breakpoint(cbam, coverage_stats=cstats, coverage_windows=coverage_windows)
    cstats = cbam2b.basic_stats
elif cstats is None:
    bamFileb2b = b2b.bam_to_breakpoint(bamFile, coverage_stats=cstats, coverage_windows=coverage_windows)
    cstats = bamFileb2b.basic_stats


final = args.final

if cstats[0] <= final:
    exit()    
ratio = float(final) / cstats[0]


downsample_dir = os.path.dirname(os.path.abspath(args.bam[0]))
if args.downsample_dir != '':
    downsample_dir = args.downsample_dir

b2 = pysam.Samfile(downsample_dir + '/' + os.path.basename(args.bam[0])[:-4] + '.DS.bam', 'wb', template = bamFile)
for a in bamFile.fetch():
    random.seed(a.qname)
    if random.uniform(0, 1) < ratio:
        b2.write(a)
b2.close()
pysam.index(downsample_dir + '/' + os.path.basename(args.bam[0])[:-4] + '.DS.bam')

# if args.cbam is not None and not os.path.exists(downsample_dir + '/' + os.path.basename(args.cbam)[:-4] + '.DS.bam'):
#     c2 = pysam.Samfile(downsample_dir + '/' + os.path.basename(args.cbam)[:-4] + '.DS.bam', 'wb', template = cbam)
#     for a in cbam.fetch():
#         random.seed(a.qname)
#         if random.uniform(0, 1) < ratio:
#             c2.write(a)
#     c2.close()
#     pysam.index(downsample_dir + '/' + os.path.basename(args.cbam)[:-4] + '.DS.bam')


