#!/usr/bin/env python

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

from time import time
TSTART = time()
import pysam
import argparse
from time import time
import os
import matplotlib
matplotlib.use('Agg')
import random

import global_names

parser = argparse.\
ArgumentParser(description="Reconstruct Amplicons connected to listed intervals.")
parser.add_argument('--bam', dest='bam',
                    help="Coordinate sorted BAM file with index", metavar='FILE',
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
                    help="Values: [hg19, GRCh37, GRCh38, GRCh38_viral, mm10, None]. \"hg19\", \"mm10\", \"GRCh38\" : chr1, .. chrM etc / \"GRCh37\" : '1', '2', .. 'MT' etc/ \"None\" : Do not use any annotations. AA can tolerate additional chromosomes not stated but accuracy and annotations may be affected.", metavar='STR',
                    action='store', type=str, required=True)
parser.add_argument('--cstats_only', help="Compute the coverage statistics for the BAM file and exit. Do not perform any downsampling.", action='store_true')
parser.add_argument('--random_seed', dest="random_seed",
                    help="Set flag to use the numpy default random seed (sets np.random.seed(seed=None)), otherwise will use seed=0",
                    action='store_true', default=False)

args = parser.parse_args()

global_names.REF = args.ref
global_names.TSTART = TSTART
if args.random_seed:
    global_names.SEED = None

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
    bamfile_pathname = str(cb.filename.decode())
    if ll[0] == os.path.abspath(bamfile_pathname):
        bamfile_filesize = os.path.getsize(bamfile_pathname)

        cstats = tuple(map(float, ll[1:]))
        if len(cstats) < 15 or cstats[13] != 3 or bamfile_filesize != int(cstats[14]):  # 3 is default sdevs
            cstats = None

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

print("Estimated bamfile coverage is ", str(cstats[0]))
if args.cstats_only:
    sys.exit(0)

final = args.final

if cstats[0] <= final:
    exit()    
ratio = float(final) / float(cstats[0])

print("Downsampling:", args.bam[0], "Estimated original coverage:", float(cstats[0]), "Desired final coverage:", final, "DS ratio:", ratio)

downsample_dir = os.path.dirname(os.path.abspath(args.bam[0]))
if args.downsample_dir != '':
    downsample_dir = args.downsample_dir

i=0
rulist = []
t0 = time()
b2 = pysam.Samfile(downsample_dir + '/' + os.path.basename(args.bam[0])[:-4] + '.DS.bam', 'wb', template=bamFile)

seed_shift = str(t0)
if global_names.SEED is not None:
    seed_shift = str(global_names.SEED)

for a in bamFile.fetch():
    random.seed(a.query_name + seed_shift)

    ru = random.uniform(0, 1)
    if ru < ratio:
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


