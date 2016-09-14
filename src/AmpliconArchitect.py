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

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

import hg19util as hg
import bam_to_breakpoint as b2b


parser = argparse.\
ArgumentParser(description="Reconstruct Amplicons connected to listed intervals.")
parser.add_argument('--bed', dest='rdAlts',
                    help="Bed file with putative list of amplified intervals", metavar='FILE',
                    action='store', type=str, nargs=1)
parser.add_argument('--bam', dest='bam',
                    help="Coordinate sorted BAM file with index mapped to hg19 reference sequence", metavar='FILE',
                    action='store', type=str, nargs=1)
parser.add_argument('--out', dest='outName',
                    help="Prefix for output files", metavar='FILE',
                    action='store', type=str, nargs=1)
args = parser.parse_args()
rdAlts = args.rdAlts[0]
bamFile = pysam.Samfile(args.bam[0], 'rb')
outName = args.outName[0]
logging.basicConfig(filename=outName+'.log',level=logging.DEBUG)
logging.info("#TIME " + str(clock()) + " import done")
summary_logger = logging.getLogger('summary')
summary_logger.addHandler(logging.FileHandler(outName + '_summary.txt', 'w'))
graph_logger = logging.getLogger('graph')
cycle_logger = logging.getLogger('cycle')

rdList0 = hg.interval_list(rdAlts, 'bed')
rdList = hg.interval_list([r for r in rdList0])
coverage_stats_file = open(hg.DATA_REPO + "/coverage.stats")
cstats = None
for l in coverage_stats_file:
    ll = l.strip().split()
    if ll[0] == os.path.abspath(bamFile.filename):
        cstats = tuple(map(float, ll[1:]))
coverage_stats_file.close()
coverage_windows=None
#coverage_windows=hg.interval_list('coverage_estimation_seq_coords_hg19.txt', 'bed')
#coverage_windows.sort()
bamFileb2b = b2b.bam_to_breakpoint(bamFile, coverage_stats=cstats, coverage_windows=coverage_windows)
ecolor = {'interchromosomal' : 'blue',
          'concordant' : 'black',
          'everted' : 'yellow',
          'forward' : 'magenta',
          'reverse' : 'cyan',
          'discordant' : 'red'}
irdhops = []
irddict = {}
irdSets = Set([Set([ird]) for ird in rdList])
irdgroupdict = {ird:Set([ird]) for ird in rdList}

logging.info("#TIME " + str(clock()) + " initiated bamFileb2b")

for ird in rdList:
#    print bamFileb2b.gc_scaling()
#    de = bamFileb2b.interval_discordant_edges(ird)
#    logging.info("#TIME " + str(clock()) + " discordant edges done")
#    bamFileb2b.meanshift_refined(ird)
#    exit()
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()
#    ilist = bamFileb2b.interval_hops(ird)
    ilist = hg.interval_list([ird])
#    ilist = hg.interval_list([bamFileb2b.interval_extend(ird)])
    irdhops.append((ird,ilist))
    for i in ilist:
        irddict[i] = ird
#    iout = open(outName + '.' + ird.chrom + ":" + str(ird.start) + '-' + str(ird.end) + '.out', 'w')
#    iout.write(mystdout.getvalue())
#    iout.close()
    sys.stdout = old_stdout
#exit()
allhops = hg.interval_list(reduce(lambda x, y: x + y, [irdh[1] for irdh in irdhops], []))
allhops.sort()
allmerge = allhops.merge_clusters()
for am in allmerge:
    nset = Set([])
    for ami in am[1]:
        nset.update(irdgroupdict[irddict[ami]])
        if irdgroupdict[irddict[ami]] in irdSets:
            irdSets.remove(irdgroupdict[irddict[ami]])
    for ird in nset:
        irdgroupdict[ird] = nset
    irdSets.add(nset)
irdgroups = []
for nset in irdSets:
    ngroup = hg.interval_list([])
    for am in allmerge:
        if irddict[am[1][0]] in nset:
            ngroup.append(am[0])
    ngroup.sort()
    irdgroups.append(ngroup)







summary_logger.info('#Amplicons = ' + str(len(irdgroups)))
summary_logger.info('-----------------------------------------------------------------------------------------')

amplicon_id = 1

for ig in irdgroups:
    ilist = ig
    ird = ig[0]
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()
    summary_logger.info("AmpliconID = " + str(amplicon_id))
    graph_handler = logging.FileHandler(outName + '_amplicon' + str(amplicon_id) + '_graph.txt', 'w')
    cycle_handler = logging.FileHandler(outName + '_amplicon' + str(amplicon_id) + '_cycles.txt', 'w')
    graph_logger.addHandler(graph_handler)
    cycle_logger.addHandler(cycle_handler)
    bamFileb2b.interval_filter_vertices(ilist)
    summary_logger.info('-----------------------------------------------------------------------------------------')
    graph_logger.removeHandler(graph_handler)
    cycle_logger.removeHandler(cycle_handler)
#    iout = open(outName + '.' + ird.chrom + ":" + str(ird.start) + '-' + str(ird.end) + '.out2', 'w')
#    iout.write(mystdout.getvalue())
#    iout.close()
    sys.stdout = old_stdout
    amplicon_id += 1
    continue

