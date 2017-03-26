#!/usr/bin/python
# This software is Copyright 2017 The Regents of the University of
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
import JobNotifier
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')



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
parser.add_argument('--cbam', dest='cbam',
                    help="Optional bamfile to use for coverage calculation", metavar='FILE',
                    action='store', type=str, default=None)
parser.add_argument('--cbed', dest='cbed',
                    help="Optional bedfile defining 1000 10kbp genomic windows for coverage calcualtion", metavar='FILE',
                    action='store', type=str, default=None)
parser.add_argument('--extendmode', dest='extendmode',
                    help="EXPLORE : Search for all connected intervals /CLUSTERED : Analyse input intervals as single amplicon only. /UNCLUSTERED : Analyse input intervals as independent amplicons only.", metavar='STR',
                    action='store', type=str, default='EXPLORE')
parser.add_argument('--sensitivems', dest='sensitivems',
                    help="Set \"True\" only if expected copy counts to vary by orders of magnitude, .e.g viral integration. Default: False", metavar='STR',
                    action='store', type=str, default='False')
parser.add_argument('--ref', dest='ref',
                    help="\"hg19\"(default) : chr1, .. chrM etc / \"GRCh37\" : '1', '2', .. 'MT' etc/ \"None\" : Do not use any annotations. AA can tolerate additional chromosomes not stated but accuracy and annotations may be affected. Default: hg19", metavar='STR',
                    action='store', type=str, default='hg19')
parser.add_argument('--downsample', dest='downsample',
                    help="Downsample the bam file during analysis (Alternatively pre-process $AA_SRC/downsample.py). Values: -1 : Do not downsample / 0 (default): Downsample to 10X coverage if larger / >0 : Downsample to stated float if larger", metavar='FLOAT',
                    action='store', type=float, default=0)
args = parser.parse_args()


logging.basicConfig(filename=args.outName[0] + '.log',level=logging.DEBUG)
# # output logs to stdout
root = logging.getLogger()
# root.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter('[%(name)s:%(levelname)s]\t%(message)s')
ch.setFormatter(formatter)
root.addHandler(ch)
summary_logger = logging.getLogger('summary')
summary_logger.propagate = False
summary_logger.addHandler(logging.FileHandler(args.outName[0] + '_summary.txt', 'w'))
graph_logger = logging.getLogger('graph')
graph_logger.propagate = False
cycle_logger = logging.getLogger('cycle')
cycle_logger.propagate = False
class PrefixAdapter(logging.LoggerAdapter):
    def process(self, msg, kwargs):
        return '[%s] %s' % (self.extra['prefix'], msg), kwargs


rdAlts = args.rdAlts[0]
if os.path.splitext(args.bam[0])[-1] == '.cram':
    bamFile = pysam.Samfile(args.bam[0], 'rc')
else:
    bamFile = pysam.Samfile(args.bam[0], 'rb')
outName = args.outName[0]
cbam = None
if args.cbam is not None:
    if os.path.splitext(args.cbam[0])[-1] == '.cram':
        cbam = pysam.Samfile(args.cbam, 'rc')
    else:
        cbam = pysam.Samfile(args.cbam, 'rb')
cbed = args.cbed
try:
    DATA_REPO = os.environ["AA_DATA_REPO"]
except:
    logging.warning("#TIME " + '%.3f\t'%clock() + "unable to set AA_DATA_REPO variable. Setting to working directory")
    DATA_REPO = '.'
if DATA_REPO == '.' or DATA_REPO == '':
    logging.warning("#TIME " + '%.3f\t'%clock() + "AA_DATA_REPO not set or empy. Setting to working directory")
    DATA_REPO = '.'
try:
    reffile = open(DATA_REPO + "/reference.txt", 'w')
    reffile.write(args.ref)
    reffile.close()
except:
    logging.warning("#TIME " + '%.3f\t'%clock() + "unable to set reference in $AA_DATA_REPO/reference.txt. Setting in working directory.")



logging.info("#TIME " + '%.3f\t'%clock() + " Loading libraries and reference annotations for: " + args.ref)
import hg19util as hg
import bam_to_breakpoint as b2b
from breakpoint_graph import *


logging.info("#TIME " + '%.3f\t'%clock() + " Initiating bam_to_breakpoint object for: " + args.bam[0])
rdList0 = hg.interval_list(rdAlts, 'bed')
rdList = hg.interval_list([r for r in rdList0])
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
bamFileb2b = b2b.bam_to_breakpoint(bamFile, coverage_stats=cstats, coverage_windows=coverage_windows, downsample=args.downsample, sensitivems=(args.sensitivems=='True'))
irdhops = []
irddict = {}
irdSets = Set([Set([ird]) for ird in rdList])
irdgroupdict = {ird:Set([ird]) for ird in rdList}


# bamFileb2b.plot_segmentation(rdList, outName + '_amplicon' + str(1))
# e1 = breakpoint_edge('chr18:27112747+->chr18:24718814-')
# c1 = bamFileb2b.refine_discordant_edge(e1)
# print c1[0], c1[1], len(c1[2])
# exit()


if args.extendmode == 'EXPLORE':
    for ird in rdList:
        logging.info("#TIME " + '%.3f\t'%clock() + " Exploring interval: " + str(ird))
        # print bamFileb2b.gc_scaling()
        # de = bamFileb2b.interval_discordant_edges(ird)
        # logging.info("#TIME " + '%.3f\t'%clock() + " discordant edges done")
        # bamFileb2b.meanshift_refined(ird)
        # exit()    
        old_stdout = sys.stdout
        sys.stdout = mystdout = StringIO()
        ilist = bamFileb2b.interval_hops(ird)
        # ilist = hg.interval_list([ird])
        # ilist = hg.interval_list([bamFileb2b.interval_extend(ird)])
        irdhops.append((ird,ilist))
        for i in ilist:
            irddict[i] = ird
        iout = open(outName + '.' + ird.chrom + ":" + str(ird.start) + '-' + str(ird.end) + '.out', 'w')
        iout.write(mystdout.getvalue())
        iout.close()
        sys.stdout = old_stdout
    # exit()

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
elif args.extendmode == 'CLUSTERED':
    irdgroups = [rdList]
else:
    irdgroups = [hg.interval_list([r]) for r in rdList]


logging.info("#TIME " + '%.3f\t'%clock() + " Interval sets for amplicons determined: ")
for il in enumerate(irdgroups):
    logging.info("[amplicon" + str(il[0] + 1) + ']\t' + ','.join([i.chrom + ':' + str(i.start) + '-' + str(i.end) for i in il[1]]))



summary_logger.info('#Amplicons = ' + str(len(irdgroups)))
summary_logger.info('-----------------------------------------------------------------------------------------')

amplicon_id = 1

for ig in irdgroups:
    logging.info("#TIME " + '%.3f\t'%clock() + " Reconstructing amplicon" + str(amplicon_id))
    ilist = ig
    ird = ig[0]
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()
    adapter = PrefixAdapter(summary_logger, {'prefix': str(amplicon_id)})
    summaryFormatter = logging.Formatter('[amplicon' + str(amplicon_id) + '] %(message)s')
    for handler in summary_logger.handlers:
        handler.setFormatter(summaryFormatter)
        # summary_logger.addHandler(handler)
    summary_logger.info("AmpliconID = " + str(amplicon_id))
    graph_handler = logging.FileHandler(outName + '_amplicon' + str(amplicon_id) + '_graph.txt', 'w')
    cycle_handler = logging.FileHandler(outName + '_amplicon' + str(amplicon_id) + '_cycles.txt', 'w')
    graph_logger.addHandler(graph_handler)
    cycle_logger.addHandler(cycle_handler)
    bamFileb2b.interval_filter_vertices(ilist)
    summary_logger.info('-----------------------------------------------------------------------------------------')
    bamFileb2b.plot_segmentation(ilist, outName + '_amplicon' + str(amplicon_id))
    graph_logger.removeHandler(graph_handler)
    cycle_logger.removeHandler(cycle_handler)
    iout = open(outName + '.' + ird.chrom + ":" + str(ird.start) + '-' + str(ird.end) + '.out2', 'w')
    iout.write(mystdout.getvalue())
    iout.close()
    sys.stdout = old_stdout
    amplicon_id += 1
    continue
# JobNotifier.sendMessage("subject","message","virajbdeshpande@gmail.com")
logging.info("#TIME " + '%.3f\t'%clock() + " Total Runtime")  