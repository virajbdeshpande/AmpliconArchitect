#!/usr/bin/env python


# This software is Copyright 2017 The Regents of the University of California. All Rights Reserved. Permission to copy,
# modify, and distribute this software and its documentation for educational, research and non-profit purposes, without
# fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and
# the following three paragraphs appear in all copies. Permission to make commercial use of this software may be obtained
# by contacting:
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
# This software program and documentation are copyrighted by The Regents of the University of California. The software
# program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents
# does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that
# the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
#
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN
# IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY
# DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO
# OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.



# Author: Viraj Deshpande, virajbdeshpande@gmail.com
# Maintained by Jens Luebeck, jluebeck@ucsd.edu


from time import time
TSTART = time()
from datetime import datetime
launch_datetime = datetime.now().strftime("%B %d, %Y %H:%M:%S")

import argparse
import copy
from functools import reduce
import logging
import os
import sys

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pysam


if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO

import global_names

__version__ = "1.4.r1"


parser = argparse.\
ArgumentParser(description="Reconstruct Amplicons connected to listed intervals.")
parser.add_argument('--bed', dest='rdAlts',
                    help="Bed file with putative list of amplified intervals", metavar='FILE',
                    action='store', type=str, required=True)
parser.add_argument('--bam', dest='bam',
                    help="Coordinate sorted BAM file with index.", metavar='FILE',
                    action='store', type=str, required=True)
parser.add_argument('-o', '--out', dest='outName',
                    help="Prefix for output files", metavar='FILE',
                    action='store', type=str, nargs=1, required=True)
parser.add_argument('--runmode', dest='runmode', choices=['FULL', 'BPGRAPH', 'CYCLES', 'SVVIEW'],
                    help="Values: [FULL/BPGRAPH/CYCLES/SVVIEW]. This option determines which stages of AA will be run. "
                         "FULL: Run the full reconstruction including breakpoint graph, cycles as well as SV visualization."
                         " BPGRAPH: Only reconstruct the breakpoint graph and estimate copy counts, but do not reconstruct"
                         " the amplicon cycles. CYCLES: Only reconstruct the breakpoint graph and cycles, but do not create"
                         " the output for SV visualization. SVVIEW: Only create the SV visualization, but do not reconstruct"
                         " the breakpoint graph or cycles", metavar='STR',
                    action='store', type=str, default='FULL')
parser.add_argument('--extendmode', dest='extendmode', choices=['EXPLORE', 'CLUSTERED', 'UNCLUSTERED', 'VIRAL'],
                    help="Values: [EXPLORE/CLUSTERED/UNCLUSTERED/VIRAL]. This determines how the input intervals in bed"
                         " file are treated. EXPLORE : Search for all connected intervals in genome that may be connected"
                         " to input intervals. CLUSTERED : Input intervals are treated as part of a single connected amplicon"
                         " and no new connected intervals are added. UNCLUSTERED : Each input interval is treated as a distinct"
                         " single interval amplicon and no new intervals are added.", metavar='STR',
                    action='store', type=str, default='EXPLORE')
parser.add_argument('--sensitivems', dest='sensitivems',
                    help="Values: [True, False]. Set \"True\" only if expected copy counts to vary by orders of magnitude,"
                         " .e.g viral integration. Default: False", metavar='STR',
                    action='store', type=str, default='False')
parser.add_argument('--plotstyle', dest='plotstyle',
                    help="Values: [small large, all_amplicons]. \"small\": small font, \"all_amplicons\": display a large"
                         " number of intervals in a single plot, recommeded for visualizing multiple amplicons in CLUSTERED"
                         " mode. Default: \"large\"", metavar='STR',
                    action='store', type=str, default="small")
parser.add_argument('--ref', dest='ref',
                    help="Values: [hg19, GRCh37, GRCh38, GRCh38_viral, mm10, GRCm38]. \"hg19\", \"GRCh38\", \"mm10\" : "
                         "chr1, .. chrM etc / \"GRCh37\", \"GRCm38\" : '1', '2', .. 'MT' etc/ \"None\" : Do not use any"
                         " annotations. AA can tolerate additional chromosomes not stated but accuracy and annotations"
                         " may be affected.", metavar='STR',
                    action='store', type=str, choices=["hg19", "GRCh37", "GRCh38", "GRCh38_viral", "mm10", "GRCm38"], required=True)
parser.add_argument('--downsample', dest='downsample',
                    help="Values: [-1, 0, C(>0)]. Decide how to downsample the bamfile during reconstruction. Reads are"
                         " automatically downsampled in real time for speedup. Alternatively pre-process bam file using"
                         " $AA_SRC/downsample.py. -1 : Do not downsample bam file, use full coverage. 0 (default): Downsample"
                         " bamfile to 10X coverage if original coverage larger then 10. C (>0) : Downsample bam file to"
                         " coverage C if original coverage larger than C", metavar='FLOAT',
                    action='store', type=float, default=0)
parser.add_argument('--cbam', dest='cbam',
                    help="Optional bamfile to use for coverage calculation", metavar='FILE',
                    action='store', type=str, default=None)
parser.add_argument('--cbed', dest='cbed',
                    help="Optional bedfile defining 1000 10kbp genomic windows for coverage calcualtion", metavar='FILE',
                    action='store', type=str, default=None)
parser.add_argument('--sv_vcf', help='Provide a VCF file of externally-called SVs to augment SVs identified by AA internally.',
                    metavar='FILE', action='store', type=str)
parser.add_argument('--sv_vcf_no_filter', help='Use all external SV calls from the --sv_vcf arg, even those without "PASS"'
                    ' in the FILTER column.', action='store_true', default=False)
parser.add_argument('--insert_sdevs', dest='insert_sdevs',
                    help="Number of standard deviations around the insert size. May need to increase for sequencing runs"
                    " with high variance after insert size selection step. (default 3.0)", metavar='FLOAT',
                    action='store', type=float, default=3)
parser.add_argument('--pair_support_min', dest='pair_support_min',
                    help="Number of read pairs for minimum breakpoint support (default 2 but typically becomes higher due"
                    " to coverage-scaled cutoffs)", metavar='INT', action='store', type=int, default=2)
parser.add_argument('--foldback_pair_support_min', help="Number of read pairs for minimum foldback SV support "
                    "(default 2 but typically becomes higher due to coverage-scaled cutoffs). Used value will be the maximum"
                    " of pair_support and this argument. Raising to 3 will help dramatically in heavily artifacted samples.",
                    metavar='INT', action='store', type=int, default=2)

parser.add_argument('--no_cstats', dest='no_cstats', help="Do not re-use coverage statistics from coverage.stats.",
                    action='store_true', default=False)
parser.add_argument('--random_seed', dest="random_seed",
                    help="Set flag to use the numpy default random seed (sets np.random.seed(seed=None)), otherwise will"
                    " use seed=0", action='store_true', default=False)

parser.add_argument('--max_seed_len', dest="max_seed_len", help="Maximum length (in bp) of seed regions to allow AA to"
                    " take as input. Will print an error if in excess. This argument helps ensure users are giving seed"
                    " regions to AA, not whole-genome CNV calls to AA.", type=int, default=500000000)
parser.add_argument('--reuse_intermediate_files', dest="reuse_intermediate_files", help="Reuse and do not delete intermediate"
                        " files of same name from previous runs in this directory of the same name."
                    " Reusing intermediate files can cause crashes if the bam or bed inputs have changed. Set this flag"
                    " if providing own _edges.txt or _cnseg.txt files.", action='store_true')
parser.add_argument("-v", "--version", action='version', version='AmpliconArchitect version {version} \n'.format(version=__version__))

args = parser.parse_args()
global_names.REF = args.ref
global_names.TSTART = TSTART
if args.random_seed:
    global_names.SEED = None

outName = args.outName[0]
logging.basicConfig(filename=outName + '.log',level=logging.DEBUG)
logging.getLogger('fontTools.subset').level = logging.WARN
logging.getLogger('fontTools.ttLib').level = logging.WARN
logging.getLogger('matplotlib.backends').level = logging.WARN
logging.getLogger('matplotlib.font_manager').level = logging.WARN

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
summary_logger.addHandler(logging.FileHandler(outName + '_summary.txt', 'w'))
graph_logger = logging.getLogger('graph')
graph_logger.propagate = False
graph_logger.info("# AmpliconArchitect " + __version__ + " " + launch_datetime)
cycle_logger = logging.getLogger('cycle')
cycle_logger.propagate = False
cycle_logger.info("# AmpliconArchitect " + __version__ + " " + launch_datetime)
class PrefixAdapter(logging.LoggerAdapter):
    def process(self, msg, kwargs):
        return '[%s] %s' % (self.extra['prefix'], msg), kwargs


commandstring = 'Commandline: '

for arg in sys.argv:
    if ' ' in arg:
        commandstring += '"{}" '.format(arg)
    else:
        commandstring+="{} ".format(arg)

logging.info(commandstring)

logging.info("AmpliconArchitect version " + __version__)
logging.info("Python version " + sys.version + "\n")
rdAlts = args.rdAlts
if not rdAlts.endswith("_AA_CNV_SEEDS.bed"):
    logging.warning("\nWARNING! The --bed file of seed regions does not appear to have been generated by AmpliconSuite-"
                    "pipeline. Failure to properly filter seed regions for AmpliconArchitect may result in errors,"
                    " uninterpretable results, long runtimes, and other serious issues!\n For more information on using"
                    " AmpliconSuite-pipeline, please visit https://github.com/AmpliconSuite/AmpliconSuite-pipeline\n")

if os.path.splitext(args.bam)[-1] == '.cram':
    bamFile = pysam.Samfile(args.bam, 'rc')
else:
    bamFile = pysam.Samfile(args.bam, 'rb')
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

# check if there are files from previous runs of the same name in the output directory. If so, print a warning.
abs_out_dir = os.path.abspath(os.path.dirname(outName))
curr_files = os.listdir(abs_out_dir)
bname = os.path.basename(outName)
holdovers = [abs_out_dir + "/" + f for f in curr_files if f.startswith(bname + "_") and (f.endswith("_edges.txt") or f.endswith("_cnseg.txt"))]
if holdovers and not args.reuse_intermediate_files and not args.runmode == "SVVIEW":
    logging.warning("#TIME " + '%.3f\t'%(time() - TSTART) + "Intermediate files from a previous run of same sample name"
                    " in same output directory were detected. These will be removed unless '--reuse_intermediate_files'"
                    " was set.\n")
    for x in holdovers:
        os.remove(x)

logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Loading libraries and reference annotations for: " + args.ref)
import ref_util as hg
import bam_to_breakpoint as b2b
import external_data_handler as exdh

logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Initiating bam_to_breakpoint object for: " + args.bam)
rdList0 = hg.interval_list(rdAlts, 'bed', exclude_info_string=True)

# check the size of the seed intervals to ensure it is within what AA can analyze.
sumlen = sum([abs(r.end - r.start) for r in rdList0])
logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Total length of seed regions: " + str(sumlen) + "bp")
if sumlen > args.max_seed_len:
    logging.error("The total length of focal amp seed regions was in excess of maximum default length allowed by AA.\n"
                  "Please ensure that the file specified by --bed is a valid AA_CNV_SEEDS.bed file, NOT whole genome CNV"
                  " calls.\n For more information on producing seed regions for AA, please see AmpliconSuite-pipeline.\n"
                  "https://github.com/jluebeck/AmpliconSuite-pipeline\n\n"
                  "To bypass this error message (not recommended!), set a different value for --max_seed_len (default "
                  "500000000).\n")
    sys.exit(1)

rdList = hg.interval_list([r for r in rdList0])
cb = bamFile
if cbam is not None:
    cb = cbam

cstats = None
if args.no_cstats:
    logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "--no_cstats was set. Will not attempt to re-use coverage.stats info")

if os.path.exists(os.path.join(hg.DATA_REPO, "coverage.stats")) and not args.no_cstats:
    coverage_stats_file = open(os.path.join(hg.DATA_REPO, "coverage.stats"))
    for l in coverage_stats_file:
        ll = l.strip().split()
        if not ll:
            continue
        bamfile_pathname = str(cb.filename.decode())
        if ll[0] == os.path.abspath(bamfile_pathname):
            bamfile_filesize = os.path.getsize(bamfile_pathname)
            cstats = tuple(map(float, ll[1:]))
            if len(cstats) < 15 or int(round(cstats[11])) < args.pair_support_min:
                cstats = None
            elif cstats[13] != args.insert_sdevs or bamfile_filesize != int(cstats[14]) or any(np.isnan(cstats)):
                cstats = None


    coverage_stats_file.close()

if cstats:
    logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Reusing cstats from " + str(os.path.join(hg.DATA_REPO, "coverage.stats")))
else:
    logging.debug("#TIME " + '%.3f\t'%(time() - TSTART) + "cstats not found, generating coverage statistics... ")


coverage_windows = None
if cbed is not None:
    coverage_windows = hg.interval_list(cbed, 'bed')
    coverage_windows.sort()
if cstats is None and cbam is not None:
    cbam2b = b2b.bam_to_breakpoint(cbam, sample_name=outName, num_sdevs=args.insert_sdevs, pair_support_min=args.pair_support_min,
                                   coverage_stats=cstats, coverage_windows=coverage_windows, foldback_pair_support_min=args.foldback_pair_support_min)
    cstats = cbam2b.basic_stats

if args.sv_vcf:
    ext_dnlist = exdh.sv_vcf_to_bplist(args.sv_vcf, filter_by_pass=not args.sv_vcf_no_filter)
else:
    ext_dnlist = []


bamFileb2b = b2b.bam_to_breakpoint(bamFile, sample_name=outName, num_sdevs=args.insert_sdevs, pair_support_min=args.pair_support_min,
                                   coverage_stats=cstats, coverage_windows=coverage_windows, downsample=args.downsample,
                                   sensitivems=(args.sensitivems == 'True'), span_coverage=(args.cbam is None), tstart=TSTART,
                                   ext_dnlist=ext_dnlist, foldback_pair_support_min=args.foldback_pair_support_min)
logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Pair support requirement: " + str(bamFileb2b.pair_support))
segments = []

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

all_ilist = copy.copy(rdList)
irdhops = []
irddict = {}
irdSets = set([frozenset([ird]) for ird in rdList])
irdgroupdict = {ird: frozenset([ird]) for ird in rdList}
if args.extendmode == 'EXPLORE' or args.extendmode == 'VIRAL':
    for ird in rdList:
        logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Exploring interval: " + str(ird))
        old_stdout = sys.stdout
        sys.stdout = mystdout = StringIO()
        ilist = bamFileb2b.interval_hops(ird, rdlist=all_ilist)
        irdhops.append((ird, ilist))
        for i in ilist:
            irddict[i] = ird

        # These files are for debugging purposes and can be used to store information in further development.
        iout = open(outName + '.' + ird.chrom + "-" + str(ird.start) + '-' + str(ird.end) + '.out', 'w')
        iout.write(mystdout.getvalue())
        iout.close()

        sys.stdout = old_stdout
        all_ilist += ilist
        all_ilist.sort()

    allhops = hg.interval_list(reduce(lambda x, y: x + y, [irdh[1] for irdh in irdhops], []))
    allhops.sort()
    allmerge = allhops.merge_clusters()
    for am in allmerge:
        nset = set()
        for ami in am[1]:
            nset.update(irdgroupdict[irddict[ami]])
            if irdgroupdict[irddict[ami]] in irdSets:
                irdSets.remove(irdgroupdict[irddict[ami]])
        for ird in nset:
            irdgroupdict[ird] = nset
        irdSets.add(frozenset(nset))
    irdgroups = []
    for nset in irdSets:
        ngroup = hg.interval_list([])
        for am in allmerge:
            if irddict[am[1][0]] in nset:
                ngroup.append(am[0])
        ngroup.sort()
        irdgroups.append(ngroup)

    irdgroups.sort()


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
    oncolist = ','.join(set([a[1].info['Name'] for a in ilist1.intersection(hg.oncogene_list)]))+','
    summary_logger.info("OncogenesAmplified = " + str(oncolist))
    amplicon_name = outName + '_amplicon' + str(amplicon_id)
    if args.runmode in ['FULL', 'CYCLES', 'BPGRAPH']:
        logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Reconstructing amplicon" + str(amplicon_id))
        graph_handler = logging.FileHandler(amplicon_name + '_graph.txt', 'w')
        cycle_handler = logging.FileHandler(amplicon_name + '_cycles.txt', 'w')
        graph_logger.addHandler(graph_handler)
        cycle_logger.addHandler(cycle_handler)
        bamFileb2b.interval_filter_vertices(ilist, amplicon_name=amplicon_name, runmode=args.runmode)
        graph_logger.removeHandler(graph_handler)
        cycle_logger.removeHandler(cycle_handler)

    if args.runmode in ['FULL', 'SVVIEW']:
        logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Plotting SV View for amplicon" + str(amplicon_id))
        bamFileb2b.plot_segmentation(
            ilist, amplicon_name, segments=segments, font=args.plotstyle)
    summary_logger.info(
            '-----------------------------------------------------------------------------------------')
    iout = open(amplicon_name + '_logs.txt', 'w')
    iout.write(mystdout.getvalue())
    iout.close()
    sys.stdout = old_stdout
    amplicon_id += 1
    continue


if (args.extendmode in ['VIRAL', 'VIRAL_CLUSTERED']) and (args.runmode in ['FULL', 'SVVIEW', 'VIRALVIEW']):
    amplicon_id = 1
    for i in irdgroups[0]:
        if i.intersects(rdList0[-1]) or len(hg.interval_list([i]).intersection(rdList)) == 0:
            continue
        logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Plotting viral view for interval " + str(i))
        bamFileb2b.plot_segmentation(hg.interval_list([i, rdList0[-1]]), outName + '_amplicon' + str(amplicon_id), scale_list=hg.interval_list([i]), font='large')
        amplicon_id += 1


logging.info("#TIME " + '%.3f\t'%(time() - TSTART) + "Total Runtime")
