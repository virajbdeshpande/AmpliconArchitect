#!/usr/bin/python

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


from time import clock, time
TSTART = time()
import pysam
import argparse
import math
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
# import JobNotifier
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

import global_names


parser = argparse.\
ArgumentParser(description="Reconstruct Amplicons connected to listed intervals.")
parser.add_argument('--bed', dest='rdAlts',
                    help="Bed file with putative list of amplified intervals", metavar='FILE',
                    action='store', type=str)
parser.add_argument('--bam', dest='bam',
                    help="Coordinate sorted BAM file with index mapped to hg19 reference sequence", metavar='FILE',
                    action='store', type=str)
parser.add_argument('--out', dest='outName',
                    help="Prefix for output files", metavar='FILE',
                    action='store', type=str, nargs=1)
parser.add_argument('--runmode', dest='runmode',
                    help="Values: [FULL/BPGRAPH/CYCLES/SVVIEW]. This option determines which stages of AA will be run. FULL: Run the full reconstruction including breakpoint graph, cycles as well as SV visualization. BPGRAPH: Only reconstruct the breakpoint graph and estimate copy counts, but do not reconstruct the amplicon cycles. CYCLES: Only reconstruct the breakpoint graph and cycles, but do not create the output for SV visualization. SVVIEW: Only create the SV visualization, but do not reconstruct the breakpoint graph or cycles", metavar='STR',
                    action='store', type=str, default='FULL')
parser.add_argument('--extendmode', dest='extendmode',
                    help="Values: [EXPLORE/CLUSTERED/UNCLUSTERED/VIRAL]. This determines how the input intervals in bed file are treated. EXPLORE : Search for all connected intervals in genome that may be connected to input intervals. CLUSTERED : Input intervals are treated as part of a single connected amplicon and no new connected intervals are added. UNCLUSTERED : Input intervals are treated as part of a single connected amplicon and no new connected intervals are added.", metavar='STR',
                    action='store', type=str, default='EXPLORE')
parser.add_argument('--sensitivems', dest='sensitivems',
                    help="Values: [True, False]. Set \"True\" only if expected copy counts to vary by orders of magnitude, .e.g viral integration. Default: False", metavar='STR',
                    action='store', type=str, default='False')
parser.add_argument('--ref', dest='ref',
                    help="Values: [hg19, GRCh37, None]. \"hg19\"(default) : chr1, .. chrM etc / \"GRCh37\" : '1', '2', .. 'MT' etc/ \"None\" : Do not use any annotations. AA can tolerate additional chromosomes not stated but accuracy and annotations may be affected. Default: hg19", metavar='STR',
                    action='store', type=str, default='hg19')
parser.add_argument('--downsample', dest='downsample',
                    help="Values: [-1, 0, C(>0)]. Decide how to downsample the bamfile during reconstruction. Reads are automatically downsampled in real time for speedup. Alternatively pre-process bam file using $AA_SRC/downsample.py. -1 : Do not downsample bam file, use full coverage. 0 (default): Downsample bamfile to 10X coverage if original coverage larger then 10. C (>0) : Downsample bam file to coverage C if original coverage larger than C", metavar='FLOAT',
                    action='store', type=float, default=0)
parser.add_argument('--cbam', dest='cbam',
                    help="Optional bamfile to use for coverage calculation", metavar='FILE',
                    action='store', type=str, default=None)
parser.add_argument('--cbed', dest='cbed',
                    help="Optional bedfile defining 1000 10kbp genomic windows for coverage calcualtion", metavar='FILE',
                    action='store', type=str, default=None)
args = parser.parse_args()

global_names.REF = args.ref

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


rdAlts = args.rdAlts
if os.path.splitext(args.bam)[-1] == '.cram':
    bamFile = pysam.Samfile(args.bam, 'rc')
else:
    bamFile = pysam.Samfile(args.bam, 'rb')
outName = args.outName[0]
cbam = None
if args.cbam is not None:
    if os.path.splitext(args.cbam)[-1] == '.cram':
        cbam = pysam.Samfile(args.cbam, 'rc')
    else:
        cbam = pysam.Samfile(args.cbam, 'rb')
cbed = args.cbed
try:
    DATA_REPO = os.environ["AA_DATA_REPO"]
except:
    logging.warning("#TIME " + '%.3f\t'%(time() - TSTART) + "unable to set AA_DATA_REPO variable. Setting to working directory")
    DATA_REPO = '.'
if DATA_REPO == '.' or DATA_REPO == '':
    logging.warning("#TIME " + '%.3f\t'%(time() - TSTART) + "AA_DATA_REPO not set or empy. Setting to working directory")
    DATA_REPO = '.'


logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Loading libraries and reference annotations for: " + args.ref)
import hg19util as hg
import bam_to_breakpoint as b2b


logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Initiating bam_to_breakpoint object for: " + args.bam)
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
    cbam2b = b2b.bam_to_breakpoint(cbam, sample_name=outName, coverage_stats=cstats, coverage_windows=coverage_windows)
    cstats = cbam2b.basic_stats
bamFileb2b = b2b.bam_to_breakpoint(bamFile, sample_name=outName, coverage_stats=cstats, coverage_windows=coverage_windows, downsample=args.downsample, sensitivems=(args.sensitivems == 'True'), span_coverage=(args.cbam is None))



segments = []
# segments=hg.interval_list(rdAlts.replace('.bed', '_segments.bed'), 'bed')

# bandsfile="karyotype.HK359.EGFR.txt"
# segments = [(l[2], hg.interval(l[1], int(l[4]), int(l[5])).intersection(i), l[6]) for l in [ll.strip().split() for ll in open(bandsfile) if 'band' in ll and ll.strip().split()[1][:3] == 'chr'] if hg.interval(l[1], int(l[4]), int(l[5])).intersects(i)]
# segments = [('', hg.interval(l[1], int(l[4]), int(l[5])), l[6]) for l in [ll.strip().split() for ll in open(bandsfile) if 'band' in ll and ll.strip().split()[1][:3] == 'chr']]


if args.extendmode == 'VIRAL':
    logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Finding integration sites: " + str(rdList[0]))
    de = bamFileb2b.interval_discordant_edges(rdList)
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()
    amplist = bamFileb2b.interval_hops(rdList, explore=False)
    alist = hg.interval_list([hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos) for e in de] + [hg.interval(e[0].v2.chrom, e[0].v2.pos, e[0].v2.pos) for e in de] + rdList)
    alist.sort()
    rdList = hg.interval_list([i[0] for i in alist.merge_clusters(extend=5000000) if len(hg.interval_list([i[0]]).intersection(amplist) + hg.interval_list([i[0]]).intersection(rdList)) > 0 ])
    rdList = hg.interval_list([hg.interval(i.chrom, max(0, i.start - 10000), min(i.end + 10000, hg.chrLen[hg.chrNum(i.chrom)])) for i in rdList])
    iout = open(outName + '.integration_search.out', 'w')
    iout.write(mystdout.getvalue())
    iout.close()
    sys.stdout = old_stdout



irdhops = []
irddict = {}
irdSets = Set([Set([ird]) for ird in rdList])
irdgroupdict = {ird:Set([ird]) for ird in rdList}
if args.extendmode == 'EXPLORE' or args.extendmode == 'VIRAL':
    for ird in rdList:
        logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Exploring interval: " + str(ird))
        old_stdout = sys.stdout
        sys.stdout = mystdout = StringIO()
        ilist = bamFileb2b.interval_hops(ird)
        irdhops.append((ird, ilist))
        for i in ilist:
            irddict[i] = ird
        iout = open(outName + '.' + ird.chrom + ":" + str(ird.start) + '-' + str(ird.end) + '.out', 'w')
        iout.write(mystdout.getvalue())
        iout.close()
        sys.stdout = old_stdout

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
elif args.extendmode == 'CLUSTERED' or args.extendmode == 'VIRAL_CLUSTERED':
    irdgroups = [rdList]
else:
    irdgroups = [hg.interval_list([r]) for r in rdList]


logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Interval sets for amplicons determined: ")
for il in enumerate(irdgroups):
    logging.info("[amplicon" + str(il[0] + 1) + ']\t' + ','.join([i.chrom + ':' + str(i.start) + '-' + str(i.end) for i in il[1]]))

summary_logger.info('#Amplicons = ' + str(len(irdgroups)))
summary_logger.info('-----------------------------------------------------------------------------------------')

if args.extendmode == 'VIRAL':
    amplicon_id = 0
else:
    amplicon_id = 1

for ig in irdgroups:
    ilist = ig
    ird = ig[0]
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()
    adapter = PrefixAdapter(summary_logger, {'prefix': str(amplicon_id)})
    summaryFormatter = logging.Formatter('[amplicon' + str(amplicon_id) + '] %(message)s')
    for handler in summary_logger.handlers:
        handler.setFormatter(summaryFormatter)
    summary_logger.info("AmpliconID = " + str(amplicon_id))
    summary_logger.info("#Intervals = " + str(len(ilist)))
    ilist1 = hg.interval_list([a[0] for a in ilist.merge_clusters()])
    istr = ','.join([i.chrom + ':' + str(i.start) + '-' + str(i.end) for i in ilist1])
    summary_logger.info("Intervals = " + str(istr))
    oncolist = ','.join(Set([a[1].info['Name'] for a in ilist1.intersection(hg.oncogene_list)]))+','
    summary_logger.info("OncogenesAmplified = " + str(oncolist))
    amplicon_name = outName + '_amplicon' + str(amplicon_id)
    if args.runmode in ['FULL', 'CYCLES', 'BPGRAPH']:
        logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Reconstructing amplicon" + str(amplicon_id))
        graph_handler = logging.FileHandler(amplicon_name + '_graph.txt', 'w')
        cycle_handler = logging.FileHandler(amplicon_name + '_cycles.txt', 'w')
        graph_logger.addHandler(graph_handler)
        cycle_logger.addHandler(cycle_handler)
        bamFileb2b.interval_filter_vertices(ilist, runmode=args.runmode)
        graph_logger.removeHandler(graph_handler)
        cycle_logger.removeHandler(cycle_handler)
    if args.runmode in ['FULL', 'SVVIEW']:
        logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Plotting SV View for amplicon" + str(amplicon_id))
        bamFileb2b.plot_segmentation(
            ilist, amplicon_name, segments=segments)
    summary_logger.info(
            '-----------------------------------------------------------------------------------------')
    iout = open(amplicon_name + '_logs.txt', 'w')
    iout.write(mystdout.getvalue())
    iout.close()
    sys.stdout = old_stdout
    amplicon_id += 1
    continue


if (args.extendmode in ['VIRAL', 'VIRAL_CLUSTERED']) and (args.runmode in ['FULL', 'SVVIEW', 'VIRALVIEW']):
    for i in irdgroups[0]:
        logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Plotting viral view for interval " + str(i))
        if i.intersects(rdList0[-1]) or len(hg.interval_list([i]).intersection(rdList)) == 0:
            continue
        bamFileb2b.plot_segmentation(hg.interval_list([i, rdList0[-1]]), outName + '_amplicon' + str(amplicon_id), scale_list=hg.interval_list([i]), font='large')
        amplicon_id += 1


logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Total Runtime")  
