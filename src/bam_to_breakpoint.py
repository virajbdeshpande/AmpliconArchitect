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

import itertools
from time import clock
import pysam
import math
import copy
from collections import defaultdict
import mosek
import sys
import numpy as np
from scipy import stats
import heapq
import random
import os
import logging
import bisect
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle, Arc
import matplotlib.ticker as ticker
from matplotlib import gridspec
from cStringIO import StringIO
import random
import re

from breakpoint_graph import *
import hg19util as hg

from mycolors import *

# matplotlib.rcParams.update({'font.size': 28})

# use Arial font if you have it. will fall back to default if not available.
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']

summary_logger = logging.getLogger('summary')
graph_logger = logging.getLogger('graph')
cycle_logger = logging.getLogger('cycle')

# suppress some specific harmless numpy warnings during AA
np.seterr(divide='ignore', invalid='ignore', )

# check mosek version
mosek_major_version = mosek.Env.getversion()[0]
if mosek_major_version > 8:
    logging.warning("Mosek version is " + '.'.join([str(x) for x in mosek.Env.getversion()]) +
                    " AA requires version 8\n")

class breakpoint_cluster:
    def __init__(self, edge, bamfile, max_insert):
        self.edge = edge


class bam_to_breakpoint():
    def __init__(self, bamfile, sample_name='', read_length=100, max_insert=400, insert_size=300,
        window_size=10000, min_coverage=30, pair_support=-1, downsample=-1,
        secondary_index=None, coverage_stats=None, coverage_windows=None,
        sensitivems=False, span_coverage=True, tstart=0):
        self.bamfile = bamfile
        self.sample_name = sample_name
        self.window_size = window_size
        self.min_coverage = min_coverage
        self.max_insert = max_insert
        self.insert_size = insert_size
        self.read_length = read_length
        self.secondary_index=secondary_index
        self.gc_scale = defaultdict(lambda: 1.0) #self.gc_scaling()
        self.gc_set = False
        self.ms_window_size = 10000
        self.downsample = downsample
        self.downsample_ratio = 1
        self.sensitivems = sensitivems
        self.span_coverage = span_coverage
        self.mapping_quality_cutoff = 5
        self.breakpoint_mapping_quality_cutoff = 20
        self.breakpoint_entropy_cutoff = 0.75
        hg.update_chrLen([(c['SN'], c['LN']) for c in self.bamfile.header['SQ']])
        self.discordant_edge_calls = {}
        self.interval_coverage_calls = {}
        self.tstart = tstart if tstart != 0 else clock()
        if coverage_stats is None:
            self.basic_stats_set = False
            self.median_coverage(window_list=coverage_windows)
        else:
            (wc_10000_median, wc_10000_avg, wc_10000_std, wc_300_median, wc_300_avg, wc_300_std, self.read_length, self.insert_size, self.insert_std, self.min_insert, self.max_insert, self.pair_support, self.percent_proper) = coverage_stats
            self.basic_stats = coverage_stats
            self.basic_stats_set = True
            r = coverage_stats

            if self.downsample < 0 or self.downsample > self.basic_stats[0]:
                self.downsample_ratio = 1
            elif self.downsample == 0:
                self.downsample_ratio = 10.0 / self.basic_stats[0] if self.basic_stats[0] > 10 else 1
            else:
                self.downsample_ratio = float(self.downsample) / self.basic_stats[0] if self.basic_stats[0] > float(self.downsample) else 1

            if self.downsample_ratio != 1:
                rr = self.downsample_ratio
                rsq = math.sqrt(rr)
                r = [i[0] * i[1] for i in zip([rr, rr, rsq, rr, rr, rsq, 1, 1, 1, 1, 1, 1, 1], r)]
                r[11] = max((r[4] / 10.0)  * ((r[7] - r[6]) / 2 / r[6])*r[12], 2)
                self.pair_support = r[11]
                self.downsample_stats = r
            else:
                self.downsample_stats = self.basic_stats
        self.coverage_logs = {}

        if pair_support != -1:
            self.pair_support = pair_support

    # Methods to find coverage and other statistics of bam file

    def fetch(self, c, s, e):
        if s > e:
            (s, e) = (e, s)
        if s < 0:
            s = 1
        if s > hg.chrLen[hg.chrNum(c)]:
            s = hg.chrLen[hg.chrNum(c)] - 1
            e = hg.chrLen[hg.chrNum(c)] - 1
        if e < 0:
            s = 1
            e = 1
        if e > hg.chrLen[hg.chrNum(c)]:
            e = hg.chrLen[hg.chrNum(c)] - 1
        if self.downsample_ratio == 1:
            for a in self.bamfile.fetch(c, s, e + 1):
                    yield a
        else:                    
            for a in self.bamfile.fetch(c, s, e + 1):
                random.seed(a.query_name.encode('hex'))
                if random.uniform(0, 1) < self.downsample_ratio:
                        yield a

    def interval_coverage(self, i, clip=False, gcc=False):
        call_args = (i.chrom, i.start, i.end, clip, gcc)
        if call_args in self.interval_coverage_calls:
            return self.interval_coverage_calls[call_args]
        if gcc:
            wc_raw = self.window_coverage(i)
            wc_corrected = 0
            j = 0
            for w in wc_raw:
                alist = [a for a in self.fetch(w[0].chrom, w[0].start, w[0].end)]
                wc_corrected += w[0].size() * w[1] / self.gc_scaling()[int(w[0].gc_content() * 10)  / 10.0]
                if w[0].size() * w[1] / self.gc_scaling()[int(w[0].gc_content() * 10)  / 10.0]  > 10 * len(alist) * self.read_length:
                    print str(i).strip(), str(w[0]).strip(), wc_corrected, len(alist), w[1], self.gc_scaling()[int(w[0].gc_content() * 10)  / 10.0], w[0].gc_content(), w[0].sequence()
                    j += 1
                    if j > 100:
                        raise ValueError("j>100")
                self.interval_coverage_calls[call_args] = wc_corrected / i.size()
            return self.interval_coverage_calls[call_args]
        s2 = i.start
        e2 = i.end
        if i.start < 0:
            s2 = 0
        if i.end > hg.chrLen[hg.chrNum(i.chrom)]:
            e2 = hg.chrLen[hg.chrNum(i.chrom)]
        if s2 >= e2:
            return 0
        alist = [a for a in self.fetch(i.chrom, s2, e2)
                 if not a.is_unmapped]
        # if e2 - s2 >= window_size and clip == False and gcc == False:
        #     return len(alist) * self.read_length / float(e2 - s2)
        sumb = 0
        # if clip == True or (clip == False and i.size() >= 100 * self.read_length):
        #     return len(alist) * self.read_length / float(i.size())
        if clip == True or (clip is None and e2-s2 < 1000):
            icc = sum([sum(a) for a in self.bamfile.count_coverage(i.chrom, s2, e2)]) / max(1.0, float(e2-s2 + 1))
            self.interval_coverage_calls[call_args] = icc
            return self.interval_coverage_calls[call_args]
        else:
            self.interval_coverage_calls[call_args] = len(
                [a for a in alist if a.reference_end - 1 <= e2]) * self.read_length / max(1.0, float(e2 - s2 + 1))
            return self.interval_coverage_calls[call_args]
        for a in alist:
            ai = hg.interval(a, bamfile=self.bamfile).intersection(i)
            if ai is not None:
                sumb += ai.size()
        if sumb / float(i.size()) > 10 * len(alist) * self.read_length / float(i.size()):
            print str(i), sumb, len(alist)
            raise ValueError("isize exception")
        self.interval_coverage_calls[call_args] = sumb / float(i.size())
        return self.interval_coverage_calls[call_args]

    def window_coverage_stats(self, i, window_size=-1, gcc=False):
        if window_size == -1:
            window_size = self.max_insert - self.read_length
        j = range(i.start, i.end, window_size)
        jj = [hg.interval(i.chrom, k, k + window_size) for k in j]
        cc = [self.interval_coverage(k, gcc=gcc) for k in jj]
        dd = [abs(cc[j + 1] - cc[j]) for j in range(len(jj) - 1)]
        return (sum(cc)/len(cc), sum(dd)/len(dd))

    def window_coverage(self, i, window_size=-1, gcc=False, clip=None, exact=True):
        # print str(i)
        if window_size == -1:
            window_size = self.max_insert - self.read_length
            
        def win_breakup(i, window_size):
            if (exact):
                (istart, iend) = (i.start, i.end)
            else:
                istart = window_size * int(round(float(i.start) / window_size))
                iend = window_size * int(round(float(i.end) / window_size))
            for k in xrange(istart, iend, window_size):
                yield hg.interval(i.chrom, k, k + window_size - 1)
        for k in win_breakup(i, window_size):
            yield (k, self.interval_coverage(k, gcc=gcc, clip=clip))
        # return [(k, self.interval_coverage(k, gcc=gcc)) for k in jj]

    def median_coverage(self, window_size=-1, gcc=False, refi=-1, window_list=None):
        if (window_size == 10000 or window_size == -1) and self.basic_stats_set and refi == -1:
            return self.downsample_stats
        if window_size == 300 and self.basic_stats_set and refi == -1:
            return self.downsample_stats[3:6]

        num_iter = 1000
        iteri = 0
        chroffset = 0
        sumchrLen = sum([l for l in hg.chrLen.values()])
        if refi != -1:
            if type(refi) == str:
                sumchrLen = hg.chrLen[hg.chrNum(refi)]
                chroffset = hg.absPos(refi, 1)
            elif type(refi) == hg.interval:
                if len([i for i in hg.centromere_list if i.chrom == refi.chrom]) == 0:
                    chr_cent = None
                else:
                    chr_cent = [i for i in hg.centromere_list if i.chrom == refi.chrom][0]
                if chr_cent is None:
                    sumchrLen = hg.chrLen[hg.chrNum(refi.chrom)]
                    chroffset = hg.absPos(refi.chrom, 1)
                elif chr_cent.end > refi.end and chr_cent.start > refi.start:
                    sumchrLen = chr_cent.start
                    chroffset = hg.absPos(refi.chrom, 1)
                elif chr_cent.start < refi.start and chr_cent.end < refi.end:
                    sumchrLen = hg.chrLen[hg.chrNum(refi.chrom)] - chr_cent.end
                    chroffset = hg.absPos(refi.chrom, 1) + chr_cent.end
                else:
                    sumchrLen = hg.chrLen[hg.chrNum(refi.chrom)]
                    chroffset = hg.absPos(refi.chrom, 1)
        # if hg.chrPos(chroffset) is None:
        if refi != -1:
            cp = hg.chrPos(chroffset)
            if cp is not None:
                ii = hg.interval(cp[0], cp[1], cp[1] + sumchrLen)
                unconserved_len = sumchrLen - sum([i[0].intersection(i[1]).size() for i in hg.interval_list([ii]).intersection(hg.conserved_regions)])
                if (sumchrLen < 1000000 or (refi != -1 and unconserved_len < 1000000)) and window_size == -1:
                    return self.downsample_stats

        elif (sumchrLen < 1000000) and window_size == -1:
            return self.downsample_stats

        if (refi != -1 or window_size != -1) and (chroffset, sumchrLen, window_size) in self.coverage_logs:
            return self.coverage_logs[(chroffset, sumchrLen, window_size)]
        # logging.info("Calculating median arm coverage " + str(refi) + " " + str(window_size))

        if not self.basic_stats_set:
            read_length = []
            insert_size = []
            window_list_index = 0
            non_mapping = 0
            while (window_list is not None and window_list_index < len(window_list)) or (window_list is None and iteri <= num_iter):
                if window_list is None:
                    newpos = int(random.random() * sumchrLen) + chroffset
                else:
                    cwindow = window_list[window_list_index]
                    window_list_index += 1
                    if cwindow.end - cwindow.start < 10000:
                        continue
                    newpos = hg.absPos(cwindow.chrom, ((cwindow.end + cwindow.start) / 2) - 5000)
                if hg.chrPos(newpos) is None:
                    logging.debug("Unable to locate reference position: " + refi.chrom + " " + str(refi.start) + " "
                                 + str(refi.end) + " " + str(newpos) + " " + str(sumchrLen))
                    iteri+=1
                    continue
                (c,p) = hg.chrPos(newpos)
                if c not in self.bamfile.references or p < 10000 or hg.chrLen[hg.chrNum(c)] < p + 10000 or len(hg.interval_list([hg.interval(c, p, p+10000)]).intersection(hg.conserved_regions, extend=10000)) > 0 or len(hg.interval_list([hg.interval(c, p, p+10000)]).intersection(hg.centromere_list, extend=10000)) > 0:
                    continue
                read_length += [a.infer_query_length(always=False) for a in self.fetch(c, p, p+10000) if not a.is_unmapped]
                insert_size += [a.template_length for a in self.fetch(c, p, p+10000) if a.is_proper_pair and not a.is_reverse and a.template_length < 10000 and a.template_length > 0]
                iteri += 1
            self.read_length = np.average(read_length)
            self.insert_size = np.average(insert_size)
            percent_proper = len(insert_size) * 2.0 / (len(read_length) + non_mapping)
            self.percent_proper = percent_proper
            self.insert_std = np.std(insert_size)
            self.max_insert = self.insert_size + 3*self.insert_std
            self.min_insert = max(0, self.insert_size - 3*self.insert_std)

        if window_size not in [-1, 300, 10000]:
            ws_list = [window_size]
        else:
            ws_list = [10000, 300]

        wc_median = []
        wc_avg = []
        wc_std = []
        for ws in ws_list:
            wc_ws = []
            iteri = 0
            window_list_index = 0
            while (window_list is not None and window_list_index < len(window_list)) or (window_list is None and iteri <= num_iter):
                if window_list is None:
                    newpos = int(random.random() * sumchrLen) + chroffset
                else:
                    cwindow = window_list[window_list_index]
                    window_list_index += 1
                    if cwindow.end - cwindow.start < 10000:
                        continue
                    newpos = hg.absPos(cwindow.chrom, ((cwindow.end + cwindow.start) / 2) - 5000)
                if hg.chrPos(newpos) is None:
                    logging.warning("Unable to locate reference position: " + refi.chrom + " " + str(refi.start) + " "
                                 + str(refi.end) + " " + str(newpos) + " " + str(sumchrLen))
                    iteri+=1
                    continue
                (c,p) = hg.chrPos(newpos)
                if c not in self.bamfile.references or p < ws or hg.chrLen[hg.chrNum(c)] < p + ws or len(hg.interval_list([hg.interval(c, p, p+ws)]).intersection(hg.conserved_regions, extend=ws)) > 0 or len(hg.interval_list([hg.interval(c, p, p+ws)]).intersection(hg.centromere_list, extend=ws)) > 0:
                    continue
                wc_ws.append(self.interval_coverage(hg.interval(c,p, p+ws), gcc=gcc))
                iteri += 1
            wc_ws.sort()
            wc_ws_median = np.median(wc_ws)
            wc_ws_filter = [c for c in wc_ws if c < 5 * wc_ws_median and c > 0]
            if len(wc_ws_filter) == 0:
                print len(wc_ws_filter), len(wc_ws), len([c for c in wc_ws if c > 0]), wc_ws_median
                wc_median.append(0)
                wc_avg.append(0)
                wc_std.append(0)
            else:
                wc_median.append(wc_ws_filter[len(wc_ws_filter) / 2])
                wc_avg.append(np.average(wc_ws_filter))
                wc_std.append(np.std(wc_ws_filter))

        if window_size not in [-1, 300, 10000] or refi != -1:
            self.coverage_logs[(chroffset, sumchrLen, window_size)] = (wc_median[0], wc_avg[0], wc_std[0])
            return (wc_median[0], wc_avg[0], wc_std[0])
        
        (wc_10000_median, wc_10000_avg, wc_10000_std) = (wc_median[0], wc_avg[0], wc_std[0])
        (wc_300_median, wc_300_avg, wc_300_std) = (wc_median[1], wc_avg[1], wc_std[1])
        self.pair_support = max((wc_300_avg / 10.0)  * ((self.insert_size - self.read_length) / 2 / self.read_length)*percent_proper, 2)
        rstats = (wc_10000_median, wc_10000_avg, wc_10000_std, wc_300_median, wc_300_avg, wc_300_std, self.read_length, self.insert_size, self.insert_std, self.min_insert, self.max_insert, self.pair_support, self.percent_proper)
        if refi == -1:
            self.basic_stats = rstats
            self.basic_stats_set = True
            print "read length", self.read_length, self.insert_size, self.insert_std, self.max_insert, percent_proper
            print "coverage stats", self.basic_stats, len(wc_ws_filter)
            print "pair support", self.pair_support
            coverage_stats_file = open(hg.DATA_REPO + "/coverage.stats", 'a')
            coverage_stats_file.write(os.path.abspath(self.bamfile.filename) + '\t' + '\t'.join(map(str, rstats)) + '\n')
            coverage_stats_file.close()

        r = rstats
        if self.downsample < 0 or self.downsample > self.basic_stats[0]:
            self.downsample_ratio = 1
        elif self.downsample == 0:
            self.downsample_ratio = 10.0 / self.basic_stats[0] if self.basic_stats[0] > 10 else 1
        else:
            self.downsample_ratio = float(self.downsample) / self.basic_stats[0] if self.basic_stats[0] > float(self.downsample) else 1
        if self.downsample_ratio != 1:
            rr = self.downsample_ratio
            rsq = math.sqrt(rr)         
            r = [i[0] * i[1] for i in zip([rr, rr, rsq, rr, rr, rsq, 1, 1, 1, 1, 1, 1, 1], r)]
            r[11] = max((r[4] / 10.0)  * ((r[7] - r[6]) / 2 / r[6])*r[12], 2)
            self.pair_support = r[11]
            self.downsample_stats = r
        else:
            self.downsample_stats = self.basic_stats

        return rstats

    def gc_scaling(self):
        if self.gc_set:
            return self.gc_scale
        gc_total_rd = {i / 10.0:0 for i in range(11)}
        gc_num_windows = {i / 10.0:0 for i in range(11)}
        for ri in range(len(self.bamfile.references)):
            # print self.bamfile.references[ri]
            # print hg.chrLen
            if hg.chrNum(self.bamfile.references[ri]) not in hg.chrLen:
                continue
            # print self.bamfile.references[ri]
            wc = self.window_coverage(hg.interval(self.bamfile.references[ri], 0, self.bamfile.lengths[ri]))
            lwc = 0
            for w in wc:
                lwc += 1
                gc_total_rd[int(w[0].gc_content() * 10.0) / 10.0] += w[1]
                gc_num_windows[int(w[0].gc_content() * 10.0) / 10.0] += 1
            # print gc_num_windows, gc_total_rd, lwc
            break
        sc_factor = sum(gc_total_rd.values()) / sum(gc_num_windows.values())
        scale = {}
        for i in gc_total_rd:
            if gc_num_windows[i] == 0:
                scale[i] = 1.0
            else:
                scale[i] = gc_total_rd[i] / gc_num_windows[i] / sc_factor
        self.gc_scale = scale
        self.gc_set = True
        print "GC scale:", scale
        return scale


    # Methods to find all coverage shifts in amplicon
    def meanshift(self, i, window_size=-1, hb=2, cov=None, rd_global=-1, h0=-1, gcc=False, n=-1):
        if window_size == -1:
           window_size = self.max_insert - self.read_length
        if rd_global == -1:
            rd_global = self.median_coverage(window_size, gcc)[0]
        if h0 == -1:
            h0 = self.median_coverage(window_size, gcc)[2]
        if rd_global == 0:
            rd_global = self.median_coverage()[0]
            h0 = self.median_coverage()[2]
        if n==-1:
            n = min(max(100, 10 * hb), 10000000 / window_size)
        j = range(len(cov))
        if cov is None:
            s2 = i.start - window_size * n
            e2 = i.end + window_size * n
            if s2 < 0:
                s2 = 0
            if e2 > hg.chrLen[hg.chrNum(i.chrom)]:
                e2 = hg.chrLen[hg.chrNum(i.chrom)]
            j = range(s2, e2, window_size)
            # j = range(i.start, i.end, window_size)
            jj = [hg.interval(i.chrom, k, k + window_size) for k in j]
            i2 = hg.interval(i.chrom, s2, e2)
            cov = [c for c in self.window_coverage(i2, window_size, gcc, exact=False)]
            # cov = [self.interval_coverage(k) for k in jj]
        # print window_size, len(cov), str(cov[0][0]).strip(), cov[0][1], str(cov[1][0]).strip(), cov[1][1]
        def hr(wi):
            if cov[wi][1] < rd_global / 4.0:
                return h0 / 2.0
            else:
                # return math.sqrt(cov[wi][1] * self.read_length / window_size)
                return math.sqrt(cov[wi][1] / rd_global) * h0
        dfi = [(cov[wi][0], sum([wj * math.exp(-0.5 * wj ** 2 / hb ** 2) * math.exp(-0.5 * (cov[wi + wj][1] - cov[wi][1])** 2 / hr(wi)**2) for wj in range(-1 * n, n + 1)])) for wi in range(n, len(j) - n) if cov[wi][0] is not None]
        # print 'meanshift', str(i), len(cov), len(dfi), len([c for c in cov if c[0] is not None]), [(str(c[0][0]), c[0][1], c[1][1]) for c in zip([cc for cc in cov if cc[0] is not None], dfi) ]
        #[(interval,ms)]
        return dfi

    
    def meanshift_pval(self, s1, s2):
        if len(s1) <= 1 and len(s2) <= 1:
           return 1.0
        if len(s1) > 1 and len(s2) > 1:
            return stats.ttest_ind(s1, s2, equal_var=False)[1]
        elif len(s1) == 1:
            zscore = abs(s1[0] - np.average(s1 + s2)) / np.std(s1 + s2)
            return stats.norm.sf(zscore)
        elif len(s2) == 1:
            zscore = abs(s2[0] - np.average(s1 + s2)) / np.std(s1 + s2)
            return stats.norm.sf(zscore)
        return 1.0
 

    def meanshift_segmentation(self, i, window_size=-1, gcc=False, pvalue=0.01):
        if window_size == -1:
            window_size = 10000
        i = hg.interval(i.chrom, window_size * int(round(float(i.start) / window_size)), window_size * int(round(float(i.end) / window_size)))
        mc = self.median_coverage(window_size, gcc)
        rd_global = mc[0]
        h0 = mc[2]
        hb_profile = [2, 5, 10, 50, 100]
        # hb_profile = [2]
        n = min(max(100, 10 * hb_profile[-1]), 10000000 / window_size) #number of windows used to calculate meanshift
        s2 = i.start - window_size * n
        e2 = i.end + window_size * n
        startskip = 0
        endskip = 0
        if s2 < 0:
            s2 = i.start % window_size
            startskip = n - (i.start - s2) / window_size
        if e2 > hg.chrLen[hg.chrNum(i.chrom)]:
            hgl = hg.chrLen[hg.chrNum(i.chrom)]
            e2 = hgl - (hgl - i.end) % window_size
            endskip = n - (hgl - i.end) / window_size
        i2 = hg.interval(i.chrom, s2, e2)
        cov = [c for c in self.window_coverage(i2, window_size, gcc, exact=False)]
        cov = [(None, 0) for ni in range(startskip)] + cov + [(None, 0) for ni in range(endskip)]
        frozen = []
        def hr(c, wlen):
            if c < rd_global / 4.0:
                return h0 / 2.0
            else:
                return math.sqrt(c / rd_global) * h0

        for hb in hb_profile:
            cov2 = copy.copy(cov)
            for ms_iterate in range(1):
                fi = -1
                if len(frozen) > 0:
                    fi = 0
                ms = [w for w in self.meanshift(i, window_size, hb, cov2, rd_global=rd_global, h0=h0, gcc=gcc, n=n)]
                segs = []
                new_seg = []
                msi = 0
                # print 'THIS0', len(frozen), fi #, frozen[0][0][0].start
                # for ff in range(len(frozen)):
                #     print "THIS", ff, frozen[ff][0][0][0].start
                while msi < len(ms):
                    if fi >= 0 and fi < len(frozen) and ms[msi][0].start == frozen[fi][0][0][0].start:
                        if len(new_seg) > 0 and (frozen[fi][1] % 2 == 1 or (ms[msi][1] > 0 and ms[msi -1] <= 0)):
                            segs.append(new_seg)
                            new_seg = ms[msi: msi + len(frozen[fi][0])]
                        else:
                            new_seg += ms[msi: msi + len(frozen[fi][0])]
                        # segs.append(ms[msi: msi + len(frozen[fi][0])])
                        msi += len(frozen[fi][0])
                        fi += 1
                        continue
                    elif ms[msi][1] > 0 and ms[msi - 1][1] <= 0 and len(new_seg) > 0:
                        segs.append(new_seg)
                        new_seg = []
                    new_seg.append(ms[msi])
                    msi += 1
                if len(new_seg) > 0:
                   segs.append(new_seg)
                cov2 = copy.copy(cov[:n])
                covi = n
                for msi in range(len(segs)):
                    s = segs[msi]
                    c = np.average([cc[1] for cc in cov[covi: covi + len(s)]])
                    cov2 += [(ss[0], c) for ss in segs[msi]]
                    covi += len(segs[msi])
                cov2 += cov[-n:]
            ci = n
            frozen = []
            cpi = n
            for si in range(len(segs)):
                c = cov2[ci][1]
                c0 = cov[ci][1]
                lseg = segs[si][-1][0].end - segs[si][0][0].start
                freeze = 0
                # if segs[si][0][0].start < 54857402 and segs[si][-1][0].end > 54857402:
                #     print (segs[si][0][0].start, segs[si][-1][0].end), (segs[si-1][0][0].start, segs[si-1][-1][0].end), (segs[si+1][0][0].start, segs[si+1][-1][0].end)
                #     print stats.ttest_ind([cc[1] for cc in cov[ci:ci + len(segs[si])]], [cs[1] for cs in cov[ci - len (segs[si - 1]):ci]], equal_var=False) 
                #     print abs(cp - c), 3 * math.sqrt(max(cp, c) / rd_global) * h0
                #     print abs(cn - c), 3 * math.sqrt(max(cn, c) / rd_global) * h0
                #     print  [cs[1] for cs in cov[ci - len (segs[si - 1]):ci]]
                #     print [cc[1] for cc in cov[ci:ci + len(segs[si])]]
                #     print  [cs[1] for cs in cov[ci + len (segs[si]):ci + len(segs[si]) + len(segs[si + 1])]]
                if si > 0:
                    if (len(segs[si]) < 15 or len(segs[si - 1]) < 15):
                        cp = cov2[ci - 1][1]
                        if abs(cp - c) > 3 * math.sqrt(max(cp, c) / rd_global) * h0:
                        # if abs(cp - c) > 2 * hr(c, window_size * len(segs[si])):
                            freeze |= 1
                    if len(segs[si]) > 1 and len(segs[si - 1]) > 1:
                        if self.meanshift_pval([cc[1] for cc in cov[ci:ci + len(segs[si])]], [cs[1] for cs in cov[ci - len (segs[si - 1]):ci]]) < pvalue:
                            freeze |= 1
                if si < len(segs) - 1:
                    if (len(segs[si]) < 15 or len(segs[si + 1]) < 15):
                        cn = cov2[ci + len(segs[si])][1]
                        if abs(cn - c) > 3 * math.sqrt(max(cn, c) / rd_global) * h0:
                        # if abs(cn - c) > 2 * hr(c, window_size * len(segs[si])):
                            freeze |= 2
                    if self.meanshift_pval([cc[1] for cc in cov[ci: ci + len(segs[si])]], [cs[1] for cs in cov[ci + len(segs[si]):ci + len(segs[si]) + len(segs[si + 1])]]) < pvalue:
                        freeze |= 2
                # if freeze > 0:
                frozen.append((segs[si], freeze, c, cov2[cpi:ci+len(segs[si])], cov[cpi:ci+len(segs[si])]))
                ci += len(segs[si])
                if freeze > 0:
                    cpi = ci
        # for f in frozen:
        #     print str(hg.interval(f[0][0][0].chrom, f[0][0][0].start, f[0][-1][0].end)), f[1], f[2], str(i)
        # print '----...-------------...------------...-------------...-----------------...----------------------'
        # (list of windows[(windowinterval,ms)], left%2/right%4freeze, avg_coverage)

        plist = []
        ms1list = []
        ms2list = []
        cms = []
        c1list = []
        c2list = []
        for msi in range(len(frozen)):
            cms.append(frozen[msi])
            if frozen[msi][1] % 4 >= 2:
                plist.append(frozen[msi][0][-1][0].end)
                avgc = np.average(reduce(lambda x,y:x+y, [[c[1] for c in cc[4]] for cc in cms], []))
                ms1list.append(avgc * 2 / self.median_coverage(window_size, gcc)[0])
                c1list.append(reduce(lambda x,y:x+y, [[c[1] for c in cc[4]] for cc in cms], []))
                if len(ms1list) > 1:
                    ms2list.append(avgc * 2 / self.median_coverage(window_size, gcc)[0])
                    c2list.append((reduce(lambda x,y:x+y, [[c[1] for c in cc[4]] for cc in cms], [])))
                cms = []
        if len(cms) > 0:
            avgc = np.average(reduce(lambda x,y:x+y, [[c[1] for c in cc[4]] for cc in cms], []))
            ms2list.append(avgc * 2 / self.median_coverage(window_size, gcc)[0])
            c2list.append((reduce(lambda x,y:x+y, [[c[1] for c in cc[4]] for cc in cms], [])))
        shifts = zip(plist, ms1list, ms2list, c1list, c2list)
        # for a in shifts:
        #    print a[0], a[1], a[2], len(a[3]), len(a[4]), str(i)
        # print '---------------------------------------------------------'

        if len(shifts) > 0:
            merge = True
        else:
            merge = False
        while merge:
            merge = False
            mergelist = []
            for shiftsi in range(len(shifts)):
                s3 = [shifts[shiftsi][3], shifts[shiftsi][3][1:], shifts[shiftsi][3][:-1], shifts[shiftsi][3][1:-1]]
                s4 = [shifts[shiftsi][4], shifts[shiftsi][4][1:], shifts[shiftsi][4][:-1], shifts[shiftsi][4][1:-1]]
                min_ttest_val = 1.0
                for s3i in s3:
                    for s4i in s4:
                        p = self.meanshift_pval(s3i, s4i)
                        min_ttest_val = min(min_ttest_val, p)
                if min_ttest_val > pvalue:
                    mergelist.append(shiftsi)
            if len(mergelist) > 0:
                merge = True
                plist = []
                ms1list = []
                ms2list = []
                c1list = []
                c2list = []
                c1 = []
                for shiftsi in range(len(shifts)):
                    c1.extend(shifts[shiftsi][3])
                    if shiftsi not in mergelist:
                        plist.append(shifts[shiftsi][0])
                        avgc = np.average(c1)
                        ms1list.append(avgc * 2 / self.median_coverage(window_size, gcc)[0])
                        c1list.append(c1)
                        if len(plist) > 1:
                            c2list.append(c1)
                            ms2list.append(avgc * 2 / self.median_coverage(window_size, gcc)[0])
                        c1 = []
                if len(plist) > 0:
                    c1.extend(shifts[-1][4])
                    avgc = np.average(c1)
                    ms2list.append(avgc * 2 / self.median_coverage(window_size, gcc)[0])
                    c2list.append(c1)
                shifts = zip(plist, ms1list, ms2list, c1list, c2list)
            # for a in shifts:
            #     print a[0], a[1], a[2], len(a[3]), len(a[4])
            # print '---------------------------------------------------------'
        if self.sensitivems:
            shifts_select = [s for s in shifts if abs(s[2] - s[1]) >= 1]
        else:
            shifts_select = [s for s in shifts if abs(s[2] - s[1]) >= max(1, min(max(s[2], s[1]) / 10.0, math.sqrt(max(s[2], s[1]))))]
        if len(shifts_select) == 0:
            return hg.interval_list([hg.interval(i.chrom, i.start, i.end, info={'cn': np.average([c[1] for c in cov[n:-n]]) * 2 / self.median_coverage(window_size, gcc)[0]})])
        else:
            shift_intervals = hg.interval_list([])
            start = i.start
            for si in shifts_select:
                shift_intervals.append(hg.interval(i.chrom, start, si[0], info={'cn':si[1]}))
                start = si[0] + 1
            shift_intervals.append(hg.interval(i.chrom, start, i.end, info={'cn':shifts_select[-1][2]}))
            return shift_intervals

    def meanshift_refined(self, i, window_size0=10000, window_size1=300, gcc=False, shifts_unrefined=None):
        if hg.chrLen[hg.chrNum(i.chrom)] < 3 * window_size0:
            ms_ws1 = self.meanshift_segmentation(i, window_size1, gcc)
            for ii in ms_ws1:
                ii.info['start_refined'] = True
                ii.info['end_refined'] = True
            return ms_ws1
        if shifts_unrefined is None:
            shifts0 = self.meanshift_segmentation(i, window_size0, gcc, pvalue=0.0027)
        else:
            shifts0 = shifts_unrefined

        shift1_intervals = hg.interval_list(hg.interval(msi.chrom, msi.end, msi.end) for msi in shifts0[:-1])
        shift1_intervals = [msi[0] for msi in shift1_intervals.merge_clusters(extend=3 * window_size0)]
        shifts1 = reduce(lambda x,y: x+y, [self.meanshift_segmentation(hg.interval(i.chrom, s.start - 3 * window_size0, s.start + 3 * window_size0), window_size1, gcc, pvalue=0.05) for s in shift1_intervals], [])

        matched_shifts = []
        prev_end = None
        for s0i in range(len(shifts0[:-1])):
            cndiff0 = shifts0[s0i + 1].info['cn'] - shifts0[s0i].info['cn']
            bests1i = None
            bestscore = 0
            for s1i in range(len(shifts1) - 1):
                if shifts1[s1i].end < i.start or shifts1[s1i].end >= i.end:
                    continue
                if abs(shifts0[s0i].end - shifts1[s1i].end) >= window_size0:
                    continue
                cndiff1 = shifts1[s1i + 1].info['cn']-shifts1[s1i].info['cn']
                if cndiff0 * cndiff1 < 0 or cndiff0 / cndiff1 <= 0.5 and cndiff0 / cndiff1 >= 2:
                    continue
                if bests1i is None:
                    bests1i = s1i
                    bestscore = abs(cndiff0 - cndiff1)
                elif abs(cndiff0 - cndiff1) < bestscore:
                    bestscore = abs(cndiff0 - cndiff1)
                    bests1i = s1i
            best_start = prev_end + \
                1 if prev_end is not None else shifts0[s0i].start
            best_end = shifts1[bests1i].end if bests1i is not None else shifts0[s0i].end
            matched_shifts.append(hg.interval(i.chrom, best_start, best_end, info={'cn':shifts0[s0i].info['cn'], 'start_refined':prev_end is not None, 'end_refined': bests1i is not None}))
            prev_end = shifts1[bests1i].end if bests1i is not None else None
        if len(shifts0) > 1:
            s0i = -1
            best_start = prev_end + \
                1 if prev_end is not None else shifts0[s0i].start
            best_end = shifts0[s0i].end
            matched_shifts.append(hg.interval(i.chrom, best_start, best_end, info={'cn':shifts0[s0i].info['cn'], 'start_refined':prev_end is not None, 'end_refined': False}))
        else:
            matched_shifts.append(hg.interval(i.chrom, i.start, i.end, info={'cn':shifts0[0].info['cn'], 'start_refined': False, 'end_refined': False}))
        return matched_shifts

    def get_meanshift(self, i, window_size0=10000, window_size1=300, gcc=False):
        file_name = "%s_%s_%s_%s_cnseg.txt" % (self.sample_name, i.chrom, i.start, i.end)
        if os.path.exists(file_name):
            msfile = open(file_name)
            msr = []
            for line in msfile:
                if len(line) == 0 or line[0] == '#':
                    continue
                ll = line.strip().split()
                msi = hg.interval(str(ll[0]), int(ll[1]), int(ll[2]), info={'cn': float(ll[3]), 'start_refined': bool(ll[4]), 'end_refined': bool(ll[5])})
                msr.append(msi)
        else:
            msr = self.meanshift_refined(i, window_size0=window_size0, window_size1=window_size1, gcc=gcc)
            msfile = open(file_name, 'w')
            msfile.write('#chrom\tstart\tend\tcn\tstart_refined\tend_refined\n')
            for ms in msr:
                msfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ms.chrom, ms.start, ms.end, ms.info['cn'], ms.info['start_refined'], ms.info['end_refined']))
            msfile.close()
        return msr

    def interval_crossing_arcs(self, chrom, start, end, strand, ilist):
            if strand == -1:
                return [a for a in self.fetch(chrom, max(0, start), min(end, hg.chrLen[hg.chrNum(chrom)]))
                        if not a.is_unmapped and a.is_reverse
                        and (a.mate_is_unmapped or a.next_reference_id == -1 or len(ilist.intersection([hg.interval(a.next_reference_name, a.next_reference_start, a.next_reference_start)])) == 0)]
            else:
                return [a for a in self.fetch(chrom, max(0, start), min(end, hg.chrLen[hg.chrNum(chrom)]))
                        if not a.is_unmapped and not a.is_reverse 
                        and (a.mate_is_unmapped or a.next_reference_id == -1 or len(ilist.intersection([hg.interval(a.next_reference_name, a.next_reference_start, a.next_reference_start)])) == 0)]


    # Methods to find breakpoint edges in amplicon
    def get_mates(self, a):
        gmt = clock()
        self.get_mates_num_calls += 1
        try:
            miter = self.secondary_index.find(a.query_name)
            retval = [m for m in miter if m.is_read1 != a.is_read1]
            self.get_mates_time += clock() - gmt
            return retval
        except:
            # print clock(), 'get_mates', str(a)
            retval = [a2 for a2 in self.fetch(a.next_reference_name, a.next_reference_start, a.next_reference_start + 1) if a2.query_name == a.query_name and a2.is_read1 != a.is_read1]
            # retval = [self.bamfile.mate(a)]
            # print clock(), 'got_mates'
            self.get_mates_time += clock() - gmt
            return retval

    def pair_support_count(self, chrom, position, strand, meanshift, foldup=False, sensitivems=True):
        #str(hg.interval(f[0][0][0].chrom, f[0][0][0].start, f[0][-1][0].end)), f[1], f[2]
        cd = 1
        for fi in range(len(meanshift)):
            f = meanshift[fi]
            if len(f) == 0:
                continue
            if not hg.interval(f[0].chrom, f[0].start, f[-1].end).intersects(hg.interval(chrom, position, position), extend=self.ms_window_size):
                continue
            for pi in range(len(f)):
                if f[pi].start + self.ms_window_size >= position:
                    break
            # pi = bisect.bisect_left(f[1], (position,))
            if pi > 0 and pi < len(f) and f[pi].start - self.ms_window_size <= position and (f[pi].info['cn'] - f[pi - 1].info['cn']) / strand > 0:
                cd = abs(f[pi].info['cn'] - f[pi - 1].info['cn'])
            elif pi > 0:
                cd = f[pi - 1].info['cn']
            else:
                cd = f[0].info['cn']
        mc = self.median_coverage()
        cd = max(1, cd)
        if self.sensitivems and sensitivems:
            cd = min(cd, 10)
        pcount = max(mc[4] * cd /20.0  * ((self.insert_size - self.read_length) / 2 / self.read_length)*mc[12], 2)
        pmincount = mc[11]
        if pcount < mc[11]:
            pcount = pmincount
        return pcount 

    def concordant_edge(self, v, bp_margin=0):
        if v.pos == 0:
            return (None, [])
        elif v.strand == 1:
            dlist = [a for a in self.fetch(v.chrom,
                                                   max(1, v.pos - self.max_insert),
                                                   v.pos)
                     if not a.is_unmapped and not a.is_reverse and a.is_proper_pair
                     and a.next_reference_name == v.chrom and a.next_reference_start >= v.pos and a.reference_start < v.pos - bp_margin
                     and a.next_reference_start < a.reference_start + self.max_insert - self.read_length]
            if len(dlist) > self.pair_support:
                v2= breakpoint_vertex(v.chrom,
                                         max(v.pos + 1, min([a.next_reference_start for a in dlist])), -1)
                logging.debug("#TIME " + '%.3f\t'%clock() + " concordant edges " + str(v) + " " + str(len(dlist)))
                return (breakpoint_edge(v, v2), dlist)
        else:
            dlist = [a for a in self.fetch(v.chrom,
                                                   max(1, v.pos - self.max_insert),
                                                   v.pos)
                     if not a.is_reverse and a.is_proper_pair and not a.is_unmapped
                     and a.next_reference_name == v.chrom and a.next_reference_start >= v.pos and a.reference_start < v.pos - bp_margin
                     and a.next_reference_start < a.reference_start + self.max_insert - self.read_length]
            if len(dlist) > self.pair_support:
                v2 = breakpoint_vertex(v.chrom,
                                         min(v.pos - 1, max([a.reference_end - 1 for a in dlist])), 1)
                logging.debug("#TIME " + '%.3f\t'%clock() + " concordant edges " + str(v) + " " + str(len(dlist)))
                return (breakpoint_edge(v, v2), dlist)
        logging.debug("#TIME " + '%.3f\t'%clock() + " concordant edges " + str(v) + " not found")
        return (None, dlist)
            
    def foldup_count(self, chrom, position, strand, cdiff=-1):
        interval = hg.interval(chrom, max(1, position - self.ms_window_size), min(hg.chrLen[hg.chrNum(chrom)], position + self.ms_window_size))
        if strand == 1:
            dlist = [a for a in self.fetch(interval.chrom,
                                                   interval.start,
                                                   interval.end)
                     if not a.is_unmapped and not a.is_reverse and a.is_paired
                     and not a.is_proper_pair and not a.mate_is_unmapped
                     and not a.mate_is_reverse and a.reference_name == a.next_reference_name
                     and abs(a.next_reference_start - a.reference_start) < 100000]#self.ms_window_size]
        else:
            dlist = [a for a in self.fetch(interval.chrom,
                                                   interval.start,
                                                   interval.end)
                     if not a.is_unmapped and a.is_reverse and a.is_paired
                     and not a.is_proper_pair and not a.mate_is_unmapped
                     and a.mate_is_reverse and a.reference_name == a.next_reference_name
                     and abs(a.next_reference_start - a.reference_start) < 100000]#self.ms_window_size]
        return len(dlist)


    def refine_discordant_edge(self, e):
        # logging.debug("#TIME " + '%.3f\t'%clock() + " refine discordant edge " + str(e))
        v1min = max(0, (e.v1.pos - self.max_insert + self.read_length if e.v1.strand == 1 else e.v1.pos) - 1)
        v2min = max(0, (e.v2.pos - self.max_insert + self.read_length if e.v2.strand == 1 else e.v2.pos) - 1)
        v1max = min(e.v1.pos + self.max_insert - self.read_length if e.v1.strand == -1 else e.v1.pos, hg.chrLen[hg.chrNum(e.v1.chrom)]) - 1
        v2max = min(e.v2.pos + self.max_insert - self.read_length if e.v2.strand == -1 else e.v2.pos, hg.chrLen[hg.chrNum(e.v2.chrom)]) - 1
        d1list = [a for a in self.fetch(e.v1.chrom, v1min, v1max) if not a.is_unmapped]
        d2list = [a for a in self.fetch(e.v2.chrom, v2min, v2max) if not a.is_unmapped]
        d1Set = Set([(a.query_name, a.is_read1, a.is_reverse, a.is_secondary) for a in d1list])
        if e.v1.strand == e.v2.strand:
            d2Set = Set([(a.query_name, a.is_read1, not a.is_reverse, not a.is_secondary) for a in d2list])
        else:
            d2Set = Set([(a.query_name, a.is_read1, a.is_reverse, not a.is_secondary) for a in d2list])
        rSet = d1Set.intersection(d2Set)        
        if len(rSet) == 0:
            return (e, 0, [], None)
        multi_r = Set([])
        d1reads = {}
        d2reads = {}
        for a in d1list:
            if (a.query_name, a.is_read1, a.is_reverse, a.is_secondary) in d1reads:
                multi_r.add((a.query_name, a.is_read1, a.is_reverse, a.is_secondary))
            d1reads[(a.query_name, a.is_read1, a.is_reverse, a.is_secondary)] = a
        if e.v1.strand == e.v2.strand:
            for a in d2list:
                if (a.query_name, a.is_read1, not a.is_reverse, not a.is_secondary) in d2reads:
                    multi_r.add((a.query_name, a.is_read1, not a.is_reverse, not a.is_secondary))
                d2reads[(a.query_name, a.is_read1, not a.is_reverse, not a.is_secondary)] = a
        else:
            for a in d2list:
                if (a.query_name, a.is_read1, a.is_reverse, not a.is_secondary) in d2reads:
                    multi_r.add((a.query_name, a.is_read1, a.is_reverse, not a.is_secondary))
                d2reads[(a.query_name, a.is_read1, a.is_reverse, not a.is_secondary)] = a

        dpairs = defaultdict(lambda: [], {})
        for aa in rSet:
            if a.query_name in multi_r:
                continue
            a1 = d1reads[aa]
            a2 = d2reads[aa]
            a1clip_prefix = 0
            for a1c in a1.cigartuples:
                if a1c[0] == 5:
                    a1clip_prefix += a1c[1]
                else:
                    break
            a2clip_prefix = 0
            for a2c in a2.cigartuples:
                if a2c[0] == 5:
                    a2clip_prefix += a2c[1]
                else:
                    break
            a1clip_suffix = 0
            for a1c in a1.cigartuples[::-1]:
                if a1c[0] == 5:
                    a1clip_suffix += a1c[1]
                else:
                    break
            a2clip_suffix = 0
            for a2c in a2.cigartuples[::-1]:
                if a2c[0] == 5:
                    a2clip_suffix += a2c[1]
                else:
                    break

            if a1.is_reverse:
                r1 = (a1.infer_query_length() + a1clip_suffix - a1.query_alignment_end, a1.infer_query_length() + a1clip_suffix - a1.query_alignment_start - 1)
            else:
                r1 = (a1.query_alignment_start + a1clip_prefix, a1.query_alignment_end - 1 + a1clip_prefix)
            if a2.is_reverse:
                r2 = (a2.infer_query_length() + a2clip_suffix - a2.query_alignment_end, a2.infer_query_length() + a2clip_suffix - a2.query_alignment_start - 1)
            else:
                r2 = (a2.query_alignment_start + a2clip_prefix, a2.query_alignment_end - 1 + a2clip_prefix)

            if r1[0] <= r2[0] and r1[1] <= r2[1]:
                hom = r1[1] - r2[0] + 1
                prefix = True
            elif r1[0] >= r2[0] and r1[1] >= r2[1]:
                hom = r2[1] - r1[0] + 1
                prefix = False
            else:
                continue

            if ((e.v1.strand == 1) == (not a1.is_reverse)) != prefix:
                continue

            if hom > 0:
                # p1 = a1.reference_end - hom - 1 if e.v1.strand == 1 else a1.reference_start + hom
                # p2 = a2.reference_end - hom - 1 if e.v2.strand == 1 else a2.reference_start + hom
                p1 = a1.reference_end - 1 if e.v1.strand == 1 else a1.reference_start
                p2 = a2.reference_end - 1 if e.v2.strand == 1 else a2.reference_start
            else:
                p1 = a1.reference_end - 1 if e.v1.strand == 1 else a1.reference_start
                p2 = a2.reference_end - 1 if e.v2.strand == 1 else a2.reference_start
            if ((e.v1.chrom, e.v1.pos, e.v1.strand) != (e.v2.chrom, e.v2.pos, e.v2.strand)):
                dpairs[(hom, p1, p2)].append((a1, a2, r1, r2))
            elif (p1 >= p2):
                dpairs[(hom, p1, p2)].append((a1, a2, r1, r2))


        if len(dpairs) == 0:
            return (e, 0, [], None)
        max_s = max([len(s) for s in dpairs.values()])
        max_p = [p for p in dpairs.keys() if len(dpairs[p]) == max_s]

        if len(max_p) != 1:
            logging.debug("#TIME " + '%.3f\t'%clock() + " refine discordant edge max_p not 1 " + str(e) + " " + str(max_p))
            return (e, 0, [], None)
        hom = max_p[0][0]
        hom_seq = ''
        if dpairs[max_p[0]][0][0].is_secondary:
            vstrand = e.v2.strand
            a = dpairs[max_p[0]][0][1]
        else:
            vstrand = e.v1.strand
            a = dpairs[max_p[0]][0][0]
        if hom >= 0:
            if vstrand == 1:
                hom_seq = a.query_sequence[a.query_alignment_end - hom: a.query_alignment_end]
            else:
                hom_seq = a.query_sequence[a.query_alignment_start: a.query_alignment_start + hom]
        else:
            if vstrand == 1:
                hom_seq = a.query_sequence[a.query_alignment_end: a.query_alignment_end + abs(hom)]
            else:
                hom_seq = a.query_sequence[a.query_alignment_start - abs(hom): a.query_alignment_start]
        p1 = max_p[0][1]
        p2 = max_p[0][2]
        logging.debug("#TIME " + '%.3f\t'%clock() + " refine discordant edge found "
            + str(breakpoint_edge(breakpoint_vertex(e.v1.chrom, p1, e.v1.strand), breakpoint_vertex(e.v2.chrom, p2, e.v2.strand)))
            + " " + str(hom) + " " + str(len(dpairs[max_p[0]])) + " " + str(len(rSet)))
        return (breakpoint_edge(breakpoint_vertex(e.v1.chrom, p1, e.v1.strand), breakpoint_vertex(e.v2.chrom, p2, e.v2.strand), hom=hom, hom_seq=hom_seq), hom, dpairs[max_p[0]], hom_seq)

    def edge_has_high_mapq(self, read_list):
        bp1_mapq = max([rr[0].mapping_quality for rr in read_list])
        bp2_mapq = max([rr[1].mapping_quality for rr in read_list])
        logging.debug("#TIME %.3f\tbreakpoint_mapq: %d %d" % (clock(), bp1_mapq, bp2_mapq))
        if bp1_mapq < self.breakpoint_mapping_quality_cutoff:
            return False
        if bp2_mapq < self.breakpoint_mapping_quality_cutoff:
            return False
        return True

    def edge_has_high_entropy(self, read_list):
        try:
            bp1_entropy = max([stats.entropy(np.unique([x for x in rr[0].get_reference_sequence().upper() if x != 'N'], return_counts=True)[1]) for rr in read_list])
            bp2_entropy = max([stats.entropy(np.unique([x for x in rr[1].get_reference_sequence().upper() if x != 'N'], return_counts=True)[1]) for rr in read_list])
        except ValueError:
            # if the MD tag is missing from the BAM file (e.g. Isaac was used as the aligner, or some CRAM files), instead use the query sequence for entropy calc.
            bp1_entropy = max([stats.entropy(np.unique([x for x in rr[0].query_alignment_sequence.upper() if x != 'N'], return_counts=True)[1]) for rr in read_list])
            bp2_entropy = max([stats.entropy(np.unique([x for x in rr[1].query_alignment_sequence.upper() if x != 'N'], return_counts=True)[1]) for rr in read_list])

        logging.debug("#TIME %.3f\tbreakpoint_entropy: %.3f %.3f" % (clock(), bp1_entropy, bp2_entropy))
        if bp1_entropy < self.breakpoint_entropy_cutoff:
            return False
        if bp2_entropy < self.breakpoint_entropy_cutoff:
            return False
        return True

    def edge_passes_filters(self, read_list, e=None):
        logging.debug("#TIME %.3f\tedge_breakpoint_filter: %s" % (clock(), str(e)))
        if self.edge_has_high_mapq(read_list) and self.edge_has_high_entropy(read_list):
            return True
        return False

    def sa_tag_overlaps_primary(self, a):
        if not a.has_tag('SA'):
            return False
        t = a.get_tag('SA').split(',')
        if t[0] != a.reference_name:
            return False
        if (t[2] == '+') != a.is_reverse:
            return False
        if min(abs(int(t[1]) - a.reference_start), abs(int(t[1]) - a.reference_end)) > self.read_length:
            return False
        return True

    def sa_tag_mismatch_breakpoint(self, a, bp):
        if not a.has_tag('SA'):
            return False
        t = a.get_tag('SA').split(',')
        if t[0] != a.reference_name:
            return True
        if (t[2] == '+') != a.is_reverse:
            return True
        if bp.strand == -1 and (a.reference_start != bp.pos or int(t[1]) != bp.pos):
            return True
        if bp.strand == 1:
            if abs(a.reference_end - bp.pos) > 10:
                return True
            cigar_counts = [int(i) for i in re.findall(r'\d+', t[3])]
            cigar_op = [i for i in re.findall(r'\D', t[3])]
            sa_ref_len = sum([i[0] for i in zip(cigar_counts, cigar_op) if i[1] in 'MDNX'])
            if abs(int(t[1]) + sa_ref_len - bp.pos) > 10:
                return True
        return False

    def interval_discordant_edges(self, interval, filter_repeats=True, pair_support=-1, ms=None, amplicon_name=None):
        logging.debug("#TIME " + '%.3f\t'%clock() + " discordant edges " + str(interval))
        if pair_support == -1:
            pair_support = self.pair_support
        if type(interval) != hg.interval_list:
            ilist = hg.interval_list([interval])
        else:
            ilist = interval
        if (tuple([(i.chrom, i.start, i.end) for i in ilist]), filter_repeats, pair_support, not ms is None) in self.discordant_edge_calls:
            return self.discordant_edge_calls[(tuple([(i.chrom, i.start, i.end) for i in ilist]), filter_repeats, pair_support, not ms is None)]

        interval = ilist[0]
        dflist = []
        drlist = []
        for i in ilist:
            dflist += [a for a in self.fetch(i.chrom, max(1, i.start), i.end)
                       if not a.is_unmapped and not a.is_reverse and a.is_paired
                       and not a.is_proper_pair
                       and not a.mate_is_unmapped
                       and not a.is_secondary and a.reference_end is not None
                       and a.mapping_quality > self.mapping_quality_cutoff
                       and not (a.reference_name == a.next_reference_name
                                and a.mate_is_reverse
                                and abs(a.reference_start - a.next_reference_start) < self.max_insert)] # this section catches everted sequencing artifacts
            drlist += [a for a in self.fetch(i.chrom, max(1, i.start), i.end)
                  if not a.is_unmapped and a.is_reverse and a.is_paired
                       and not a.is_proper_pair
                       and not a.mate_is_unmapped
                       and not a.is_secondary and a.reference_end is not None
                       and a.mapping_quality > self.mapping_quality_cutoff
                       and not (a.reference_name == a.next_reference_name
                                and not a.mate_is_reverse
                                and abs(a.reference_start - a.next_reference_start) < self.max_insert)] # this section catches everted sequencing artifacts
        logging.debug("#TIME " + '%.3f\t'%clock() + " discordant edges: fetch discordant " + str(interval) + " " + str(len(dflist)) + " " + str(len(drlist)))
        # dflist = [a for a in dflist if not(a.reference_name == a.next_reference_name and a.mate_is_reverse and abs(a.template_length) < self.max_insert)]
        # drlist = [a for a in drlist if not(a.reference_name == a.next_reference_name and not a.mate_is_reverse and abs(a.template_length) < self.max_insert)]

        # dflist = [a for a in dflist if not(a.reference_name == a.next_reference_name and a.mate_is_reverse and abs(a.reference_start - a.next_reference_start) < self.max_insert)]
        # drlist = [a for a in drlist if not(a.reference_name == a.next_reference_name and not a.mate_is_reverse and abs(a.reference_start - a.next_reference_start) < self.max_insert)]


        logging.debug("#TIME " + '%.3f\t'%clock() + " discordant edges: discordant read pairs found: %s %s %s" % (str(interval) , len(dflist), len(drlist)))

        # perform biclustering for readpairs using union-find algorithm to give sets of connected read-pairs clist
        vlist = []
        vcount = 0
        vdict = {}
        for a in dflist + drlist:
            vlist.append((hg.absPos(a.reference_name, a.reference_start) * (-1 if a.is_reverse else 1), hg.absPos(a.next_reference_name, a.next_reference_start) * (-1 if a.mate_is_reverse else 1), a, vcount))
            vdict[vcount] = a
            vcount += 1
        # vlist = [(hg.absPos(a.reference_name, a.reference_start) * (-1 if a.is_reverse else 1), hg.absPos(a.next_reference_name, a.next_reference_start) * (-1 if a.mate_is_reverse else 1), a) for a in dflist + drlist]
        v0list = copy.copy(vlist)
        v0list.sort(key=lambda x: x[0])
        v1list = copy.copy(vlist)
        v1list.sort(key=lambda x: x[1])
        dlist = []
        v0listp = [v[0] for v in v0list]
        v1listp = [v[1] for v in v1list]
        plist = defaultdict(lambda: None, {})
        rlist = defaultdict(lambda: 0, {})
        nlist = defaultdict(lambda: 1, {})
        # identify edges with bisect and union-find algorithm
        # iii = 0
        for v in vlist:
            # iii += 1
            s0 = bisect.bisect_left(v0listp, v[0] - self.max_insert + self.read_length)
            e0 = bisect.bisect_right(v0listp, v[0] + self.max_insert - self.read_length)
            s1 = bisect.bisect_left(v1listp, v[1] - self.max_insert + self.read_length)
            e1 = bisect.bisect_right(v1listp, v[1] + self.max_insert - self.read_length)
            SS0 = [vv[3] for vv in v0list[s0:e0+1] if vv[3] > v[3]]
            SS1 = [vv[3] for vv in v1list[s1:e1+1] if vv[3] > v[3]]
            SS0.sort()
            SS1.sort()
            SS_intersect = []
            i0 = 0
            i1 = 0
            while True:
                if i0 == len(SS0) or i1 == len(SS1):
                    break
                if SS0[i0] == SS1[i1]:
                    SS_intersect.append(SS0[i0])
                    i0 += 1
                    i1 += 1
                elif SS0[i0] < SS1[i1]:
                    i0 += 1
                else:
                    i1 += 1
            if len(SS_intersect) >= pair_support:
                dlist.append(v[2])
            v1 = v[3]
            for v2 in SS_intersect:
                v1g = v1
                v2g = v2
                while plist[v1g] is not None:
                    v1g = plist[v1g]
                while plist[v2g] is not None:
                    v2g = plist[v2g]
                if v1g == v2g:
                    continue
                if rlist[v1g] > rlist[v2g]:
                    plist[v2g] = v1g
                    rlist[v1g] = max(rlist[v1g], rlist[v2g] + 1)
                    nlist[v1g] += nlist[v2g]
                else:
                    plist[v1g] = v2g
                    rlist[v2g] = max(rlist[v2g], rlist[v1g] + 1)
                    nlist[v2g] += nlist[v1g]
        clist = defaultdict(lambda: [], {})
        for v in plist:
            vg = v
            while plist[vg] is not None:
                vg = plist[vg]
            clist[vdict[vg]].append(vdict[v])

        mcdflist = []
        mcdrlist = []
        hgddict = {}
        for c in clist:
            if len(clist[c]) < pair_support:
                continue
            ml = clist[c]
            if filter_repeats:
                ml = [v for v in clist[c] if not hg.interval(v, bamfile=self.bamfile).filter_repeat() and v.mapping_quality > self.mapping_quality_cutoff]
                if len(ml) < pair_support:
                    continue
            hgl = hg.interval_list([])
            for v in ml:
                hgv = hg.interval(v, bamfile=self.bamfile)
                hgddict[hgv] = v
                hgl.append(hgv)
            hgl.sort()
            if c.is_reverse:
                mcdrlist.extend(hgl.merge_clusters(extend=self.max_insert - self.read_length))
            else:
                mcdflist.extend(hgl.merge_clusters(extend=self.max_insert - self.read_length))


        logging.debug("#TIME " + '%.3f\t'%clock() + " discordant edges: discordant clusters found: %s %d %d " % (str(interval), len(mcdflist), len(mcdrlist)))

        dnlist0 = []
        dnlist = []
        clist = hg.interval_list([c[0] for c in mcdflist + mcdrlist])
        clist.sort()
        ci = 0
        for c1 in mcdflist + mcdrlist:
            ci += 1
            neighbor_hglist = hg.interval_list([])
            for a1 in c1[1]:
                neighbor_hglist.append(hg.interval(hgddict[a1].next_reference_name, hgddict[a1].next_reference_start, hgddict[a1].next_reference_start))
            neighbor_hglist.sort()
            neighbor_hglist = hg.interval_list([a2[0] for a2 in neighbor_hglist.merge_clusters(extend=self.max_insert - self.read_length) if len(a2[1]) >= pair_support])
            for c2 in mcdflist + mcdrlist:
                if len(hg.interval_list([c2[0]]).intersection(neighbor_hglist, extend=self.max_insert)) == 0:
                    continue
                vl = []
                vlSet = Set([])
                vl1Set = Set([])
                vl2Set = Set([])
                for a1 in c1[1]:
                    for a2 in c2[1]:
                        aq1 = hgddict[a1]
                        aq2 = hgddict[a2]
                        if aq1.query_name == aq2.query_name and aq1.is_read1 != aq2.is_read1:
                            if aq1.reference_name == aq2.reference_name and abs(aq1.reference_start - aq2.reference_start) < self.read_length and abs(aq1.reference_end - aq2.reference_end) < self.read_length and aq1.is_reverse != aq2.is_reverse:
                                continue
                            if aq1.reference_name == aq2.reference_name and  aq1.is_reverse and not aq2.is_reverse and aq1.reference_start - aq2.reference_end + 1> 0 and aq1.reference_start - aq2.reference_end + 1< self.max_insert - 2 * self.read_length:
                                continue
                            if aq2.reference_name == aq1.reference_name and aq2.is_reverse and not aq1.is_reverse and aq2.reference_start - aq1.reference_end + 1 > 0 and aq2.reference_start - aq1.reference_end + 1 < self.max_insert - 2 * self.read_length:
                                continue
                            vl.append((aq1, aq2))
                            vlSet.add((aq1.reference_start, aq1.reference_end, aq2.reference_start, aq2.reference_end))
                            vl1Set.add((aq1.reference_start, aq1.reference_end))
                            vl2Set.add((aq2.reference_start, aq2.reference_end))
                if len(vl) == 0 or len([v for v in vl if v[1].reference_start*v[0].reference_start > 0]) == 0:
                    continue
                if not vl[0][0].is_reverse:
                    bp1 = breakpoint_vertex(c1[0].chrom, max([v[0].reference_end - 1 for v in vl if v[0].reference_start > 0]), 1)
                else:
                    bp1 = breakpoint_vertex(c1[0].chrom, min([v[0].reference_start for v in vl if v[0].reference_start > 0]), -1)
                if not vl[0][1].is_reverse:
                    bp2 = breakpoint_vertex(c2[0].chrom, max([v[1].reference_end - 1 for v in vl if v[1].reference_start > 0]), 1)
                else:
                    bp2 = breakpoint_vertex(c2[0].chrom, min([v[1].reference_start for v in vl if v[1].reference_start > 0]), -1)
                if ms is None:
                    ps = pair_support
                else:
                    ps = self.pair_support_count(bp1.chrom, bp1.pos, bp1.strand, ms)
                if len(vl) < ps or len(vl1Set) < pair_support or len(vl2Set) < pair_support:
                    continue

                if bp1.chrom == bp2.chrom and bp1.pos == bp2.pos and bp1.strand == bp2.strand and len(vl) < 2 * self.pair_support:
                    continue

                num_inverted = 0
                bp1c = None
                bp2c = None
                vl2 = []
                if bp1.chrom == bp2.chrom and bp1.strand == bp2.strand and abs(bp1.pos - bp2.pos) <= self.read_length:
                    non_inverted_reads = Set([])
                    multiple_non_inverted = False
                    if bp1.strand == 1:
                        for v in vl:
                            if v[0].reference_start == v[1].reference_start:
                                num_inverted += 1
                            elif self.sa_tag_overlaps_primary(v[0]):
                                num_inverted += 1
                            elif self.sa_tag_overlaps_primary(v[1]):
                                num_inverted += 1
                            else:
                                vl2.append(v)
                                if not multiple_non_inverted:
                                    non_inverted_reads.add(v[0].query_name)
                                    if len(non_inverted_reads) >= ps:
                                        multiple_non_inverted = True
                    else:
                        for v in vl:
                            if v[0].reference_end == v[1].reference_end:
                                num_inverted += 1
                            elif self.sa_tag_overlaps_primary(v[0]):
                                num_inverted += 1
                            elif self.sa_tag_overlaps_primary(v[1]):
                                num_inverted += 1
                            else:
                                vl2.append(v)
                                if not multiple_non_inverted:
                                    non_inverted_reads.add(v[0].query_name)
                                    if len(non_inverted_reads) >= ps:
                                        multiple_non_inverted = True
                    logging.debug("checking foldback2: " + str(bp1) + str(bp2) + " %s %s %d %d %d" % (bp1.strand, bp2.strand, len(vl), num_inverted, ps))

                    if len(vl2) < ps or (not multiple_non_inverted):
                        logging.debug("FOLDBACK: " + str(bp1) + str(bp2))
                        continue
                    vl = vl2
                    vl.sort(lambda x, y: x[0].reference_start - y[0].reference_start)
                    if bp1.strand == 1:
                        maxp = vl[0][0].reference_end - 1
                        maxn = 0
                        for v in vl[::-1]:
                            if len([v1 for v1 in vl if v1[0].reference_end <= v[0].reference_end and v1[0].reference_start > v[0].reference_end - 1 - self.max_insert + 2 * self.read_length]) > maxn:
                                maxp = v[0].reference_end
                                maxn = len([v1 for v1 in vl if v1[0].reference_end <= v[0].reference_end and v1[0].reference_end > v[0].reference_end - self.max_insert + 2 * self.read_length])
                        vl = [v for v in vl if v[0].reference_end - 1 <= maxp and v[0].reference_end - 1 > maxp -  self.max_insert + 2 * self.read_length]
                        if len(vl) < ps:
                            continue
                        bp1 = breakpoint_vertex(c1[0].chrom, max([v[0].reference_end - 1 for v in vl if v[0].reference_start > 0]), 1)
                        bp2 = breakpoint_vertex(c2[0].chrom, max([v[1].reference_end - 1 for v in vl if v[1].reference_start > 0]), 1)
                        if bp1.pos != bp2.pos:
                            bp1c = bp2
                            bp2c = bp1
                    else:
                        maxp = vl[-1][0].pos
                        maxn = 0
                        for v in vl:
                            if len([v1 for v1 in vl if v1[0].reference_start >= v[0].reference_start and v1[0].reference_start < v[0].reference_start + self.max_insert - 2 * self.read_length]) > maxn:
                                maxp = v[0].reference_start
                                maxn = len([v1 for v1 in vl if v1[0].reference_start >= v[0].reference_start and v1[0].reference_start < v[0].reference_start + self.max_insert - 2 * self.read_length])
                        vl = [v for v in vl if v[0].reference_start >= maxp and v[0].reference_start < maxp +  self.max_insert - 2 * self.read_length]
                        if len(vl) < ps:
                            continue
                        bp1 = breakpoint_vertex(c1[0].chrom, min([v[0].reference_start for v in vl if v[0].reference_start > 0]), -1)
                        bp2 = breakpoint_vertex(c2[0].chrom, min([v[1].reference_start for v in vl if v[1].reference_start > 0]), -1)
                        if bp1.pos != bp2.pos:
                            bp1c = bp2
                            bp2c = bp1
                bre_refine = self.refine_discordant_edge(breakpoint_edge(bp1, bp2))
                bre = bre_refine[0]

                if bp1.chrom == bp2.chrom and bp1.strand == bp2.strand and abs(bp1.pos - bp2.pos) <= self.read_length:
                    qname_exclude = set([])
                    for v in vl:
                        if (bp1.strand == 1 and max(v[0].reference_start, v[1].reference_start) > bre.v1.pos) or (bp1.strand == -1 and max(v[0].reference_end, v[1].reference_end) < bre.v1.pos):
                            qname_exclude.add(v[0].query_name)
                            continue
                        if (self.sa_tag_mismatch_breakpoint(v[0], bre.v1) or self.sa_tag_mismatch_breakpoint(v[0], bre.v1) or self.sa_tag_overlaps_primary(v[0]) or self.sa_tag_overlaps_primary(v[1])):
                            qname_exclude.add(v[0].query_name)
                            continue
                        if (bp1.strand == 1 and bre.v1.pos - v[0].reference_start + bre.v2.pos - v[1].reference_start > self.max_insert):
                            qname_exclude.add(v[0].query_name)
                            continue
                        if (bp2.strand == 1 and v[0].reference_end - bre.v1.pos + v[1].reference_end - bre.v2.pos > self.max_insert):
                            qname_exclude.add(v[0].query_name)
                            continue
                    vl = [v for v in vl if v[0].query_name not in qname_exclude]
                    if len(vl) < ps:
                        continue


                if bre.type() == 'everted' and abs(bre.v1.pos - bre.v2.pos) <= 30:
                    continue
                if bre.type() != 'concordant':
                    if self.edge_passes_filters(vl, bre):
                        dnlist0.append((bre, len(vl)))
                if bp1c is not None and bp2c is not None:
                    brec_refine = self.refine_discordant_edge(breakpoint_edge(bp1c, bp2c))
                    brec = brec_refine[0]
                    if brec.type() != 'concordant' and brec.v1.pos != brec.v2.pos:
                        if self.edge_passes_filters(vl, brec):
                            dnlist0.append((brec, len([(v[1], v[0]) for v in vl])))

        # remove local edges with no complementary edges and add warning if any found
        for bb1 in dnlist0:
            for bb2 in dnlist0:
                bre1 = bb1[0]
                bre2 = bb2[0]
                if bre1 == bre2 and (bre1.v1.chrom, bre1.v1.pos, bre1.v1.strand) != (bre1.v2.chrom, bre1.v2.pos, bre1.v2.strand):
                    continue
                if ((bre2.v2.chrom, bre2.v2.pos, bre2.v2.strand) == (bre1.v1.chrom, bre1.v1.pos, bre1.v1.strand) and
                    (bre2.v1.chrom, bre2.v1.pos, bre2.v1.strand) == (bre1.v2.chrom, bre1.v2.pos, bre1.v2.strand)) and bb1 not in dnlist:
                    dnlist.append(bb1)
                    continue
        if len(dnlist) != len(dnlist0):
            logging.warning("dnlists do not match " + str(len(dnlist0)) + " " + str(len(dnlist)))
            for bb1 in dnlist0:
                if bb1 not in dnlist:
                    logging.warning('dnlist0: ' + str(bb1[0]) + " " + str(bb1[1]))
            for bb1 in dnlist:
                if bb1 not in dnlist0:
                    logging.warning('dnlist: ' + str(bb1[0]) + " " + str(bb1[1]))

        logging.debug("#TIME " + '%.3f\t'%clock() + " discordant edges: local edges done " + str(interval) + " " + str(len(mcdflist)) + " " + str(len(mcdrlist)) + " " + str(len(dnlist)))
        self.get_mates_time = 0
        self.get_mates_num_calls = 0
        for c in mcdflist + mcdrlist:
            nlist = []
            if filter_repeats:
                if len(hg.interval_list([c[0]]).intersection(hg.conserved_regions)) > 0:
                    continue
            rep_content_time = 0
            intersection_time = 0
            nr_calls = 0
            for hga in c[1]:
                nmatelist = self.get_mates(hgddict[hga])
                if filter_repeats:
                    rpc = clock()
                    nmatelist = [a for a in nmatelist if not hg.interval(a, bamfile=self.bamfile).filter_repeat() and a.mapping_quality > self.mapping_quality_cutoff]
                    nr_calls += len(nmatelist)
                    rep_content_time += clock() - rpc
                ict = clock()
                nmatelist = [a for a in nmatelist if len(hg.interval_list([hg.interval(a, bamfile=self.bamfile)]).intersection(ilist)) == 0]
                intersection_time += clock() - ict
                nlist += nmatelist
            nflist = [n for n in nlist if not n.is_reverse]
            nrlist = [n for n in nlist if n.is_reverse]
            hgndict = {hg.interval(a, bamfile=self.bamfile):a for a in nflist + nrlist}
            hgnflist = hg.interval_list([hga for hga in hgndict if hga.strand==1])
            hgnrlist = hg.interval_list([hga for hga in hgndict if hga.strand==-1])
            hgnflist.sort()
            hgnrlist.sort()
            mcnflist = hgnflist.merge_clusters(self.max_insert - 2 * self.read_length)
            mcnrlist = hgnrlist.merge_clusters(self.max_insert - 2 * self.read_length)
            mcnflist = [m for m in mcnflist if len(m[1]) >= pair_support] 
            mcnrlist = [m for m in mcnrlist if len(m[1]) >= pair_support] 
            mcnlist = mcnflist + mcnrlist
            for cn in mcnlist:
                vl = []
                vlSet = Set([])
                vl1Set = Set([])
                vl2Set = Set([])
                if filter_repeats:
                    if len(hg.interval_list([cn[0]]).intersection(hg.conserved_regions)) > 0:
                        continue
                hgmi = 0
                for hgm in cn[1]:
                    hgmi += 1
                    if filter_repeats:
                        if hgm.filter_repeat() or hgndict[hgm].mapping_quality <= self.mapping_quality_cutoff:
                            continue
                    for a in self.get_mates(hgndict[hgm]):
                        if filter_repeats:
                            if hg.interval(a, bamfile=self.bamfile).filter_repeat() or a.mapping_quality <= self.mapping_quality_cutoff:
                                continue
                        if hg.interval(a, bamfile=self.bamfile).intersects(c[0]):
                            vl.append((a, hgndict[hgm]))
                            vlSet.add((a.reference_start, a.reference_end, hgndict[hgm].reference_start, hgndict[hgm].reference_end))
                            vl1Set.add((a.reference_start, a.reference_end))
                            vl2Set.add((hgndict[hgm].reference_start, hgndict[hgm].reference_end))
                            break
                if len(vl) == 0 or len([v for v in vl if v[1].reference_start*v[0].reference_start > 0]) == 0:
                    continue
                if not vl[0][0].is_reverse:
                    bp1 = breakpoint_vertex(vl[0][0].reference_name,
                                            max([v[0].reference_end - 1 for v in vl]), 1)
                else:
                    bp1 = breakpoint_vertex(vl[0][0].reference_name,
                                            min([v[0].reference_start for v in vl if v[0].reference_start > 0]), -1)
                if not vl[0][1].is_reverse:
                    bp2 = breakpoint_vertex(vl[0][1].reference_name,
                                            max([v[1].reference_end - 1 for v in vl]), 1)
                else:
                    bp2 = breakpoint_vertex(vl[0][1].reference_name,
                                            min([v[1].reference_start for v in vl if v[1].reference_start > 0]), -1)
                if ms is None:
                    ps = pair_support
                else:
                    ps = self.pair_support_count(bp1.chrom, bp1.pos, bp1.strand, ms)

                if len(vl) < ps or len(vl1Set) < pair_support or len(vl2Set) < pair_support:
                    continue
                num_inverted = 0
                non_inverted_reads = Set([])
                multiple_non_inverted = False
                if bp1.chrom == bp2.chrom and bp1.pos == bp2.pos and bp1.strand == bp2.strand:
                    if bp1.strand == 1:
                        for v in vl:
                            if v[0].reference_start == v[1].reference_start:
                                num_inverted += 1
                            elif not multiple_non_inverted:
                                non_inverted_reads.add(v[0].query_name)
                                if len(non_inverted_reads) >= ps:
                                    multiple_non_inverted = True
                    else:
                        for v in vl:
                            if v[0].reference_end == v[1].reference_end:
                                num_inverted += 1
                            elif not multiple_non_inverted:
                                non_inverted_reads.add(v[0].query_name)
                                if len(non_inverted_reads) >= ps:
                                    multiple_non_inverted = True
                    if len(vl) - num_inverted < ps or (not multiple_non_inverted):
                        continue
                bre_refine = self.refine_discordant_edge(breakpoint_edge(bp1, bp2))
                bre = bre_refine[0]
                if bre.type() != 'concordant':
                    if self.edge_passes_filters(vl, bre):
                        dnlist.append((bre, len(vl)))
        logging.debug("#TIME " + '%.3f\t'%clock() + " discordant edges: external edges done " + str(interval) + " " + str(self.get_mates_time) + " " +  str(self.get_mates_num_calls))
        dnlist.sort(key=lambda x: hg.absPos(x[0].v1.chrom, x[0].v1.pos) + 0.5 * x[0].v1.strand)
        for e in dnlist:
            logging.debug("#TIME %.3f\tdiscordant edges %s %s %s %s %d %f" % (clock(), e[0], e[1], e[0].type(), self.concordant_edge(e[0].v1)[0], len(self.concordant_edge(e[0].v1)[1]), hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos - e[0].v1.strand * self.max_insert).rep_content()))
        self.discordant_edge_calls[(tuple([(i.chrom, i.start, i.end) for i in ilist]),
                                    filter_repeats, pair_support, not ms is None)] = dnlist
        return dnlist


    def load_edges(self, edge_file):
        edge_lines = [line.strip().split() for line in open(edge_file)]
        edges = []
        for el in edge_lines:
            if el[2] == 'None':
                hom = None
                hom_seq = None
            else:
                hom = int(el[2])
                if hom != 0:
                    hom_seq = el[3]
                else:
                    hom_seq = ''
            e = breakpoint_edge(el[0], hom=hom, hom_seq=hom_seq)
            edges.append((e, int(el[1])))
        edges.sort(key=lambda x: hg.absPos(
                    x[0].v1.chrom, x[0].v1.pos) + 0.1 * x[0].v1.strand)
        return edges
        

    def get_sensitive_discordant_edges(self, ilist, msrlist, eilist=None, filter_repeats=True, pair_support=-1, ms_window_size0=10000, ms_window_size1=300, adaptive_counts=True, gcc=False, amplicon_name=None):
        if amplicon_name is not None and os.path.exists("%s_edges_cnseg.txt" % amplicon_name):
            return self.load_edges("%s_edges_cnseg.txt" % amplicon_name)
        if amplicon_name is not None and os.path.exists("%s_edges.txt" % amplicon_name):
            eilist = self.load_edges("%s_edges.txt" % amplicon_name)
        else:
            if eilist is None:
                if adaptive_counts:
                    eilist = self.interval_discordant_edges(
                        ilist, ms=msrlist, pair_support=pair_support)
                else:
                    eilist = self.interval_discordant_edges(ilist, pair_support=pair_support)
            eilist.sort(key=lambda x: hg.absPos(
                        x[0].v1.chrom, x[0].v1.pos) + 0.1 * x[0].v1.strand)
            if amplicon_name is not None:
                edge_file = open("%s_edges.txt" % amplicon_name, 'w')
                for e in eilist:
                    edge_file.write("%s\t%s\t%s\t%s\n" % (str(e[0]), e[1], e[0].hom, e[0].hom_seq))
                edge_file.close()
        eiSet = Set([(e[0].v1.chrom, e[0].v1.pos, e[0].v1.strand, e[0].v2.chrom, e[0].v2.pos, e[0].v2.strand) for e in eilist])
        for i, msr in zip(ilist, msrlist):
            elist = []
            for e in eilist:
                if e[0].v1.pos != -1 and hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos).intersects(i):
                    elist.append(e)
            ms_vlist = []
            msv_index = {}
            for msi in range(len(msr) - 1):
                if msr[msi + 1].info['cn'] < msr[msi].info['cn']:
                    msv = breakpoint_vertex(i.chrom, msr[msi].end, 1)
                else:
                    msv = breakpoint_vertex(i.chrom, msr[msi].end + 1, -1)
                ms_vlist.append(msv)
                msv_index[msv] = msi
            print "Meanshift", str(i), len(ms_vlist), ms_vlist
            sys.stdout.flush()
            for msv in ms_vlist:
                msi = msv_index[msv]
                if ('end_refined' in msr[msi].info) and msr[msi].info['end_refined']:
                    msve = [e for e in elist if e[0].v1.strand * (msr[msi].info['cn'] - msr[msi + 1].info['cn']) > 0 and abs(e[0].v1.pos - msr[msi].end) < self.max_insert + ms_window_size1]
                    if len(msve) == 0:
                        print "finesearch discordant edges", i.chrom, str(msr[msi]), str(msr[msi + 1])
                        efine = self.interval_discordant_edges(hg.interval(i.chrom, msv.pos - ms_window_size0-self.max_insert, msv.pos + ms_window_size1+self.max_insert), pair_support=2)
                        if len([e for e in efine if e[0].v1.strand * (msr[msi].info['cn'] - msr[msi + 1].info['cn']) > 0]) > 0:
                            if len([(e[1], e[0]) for e in efine if e[0].v1.strand * (msr[msi].info['cn'] - msr[msi + 1].info['cn']) > 0 and abs(e[0].v1.pos - msv.pos) < ms_window_size1]) > 0:
                                ebest = max([(e[1], e[0]) for e in efine if e[0].v1.strand * (
                                    msr[msi].info['cn'] - msr[msi + 1].info['cn']) > 0 and abs(e[0].v1.pos - msv.pos) < ms_window_size1])
                            else:
                                ebest = max([(e[1], e[0]) for e in efine if e[0].v1.strand * (
                                    msr[msi].info['cn'] - msr[msi + 1].info['cn']) > 0])
                            ebest = (ebest[1], ebest[0])
                            msve = [ebest]
                            print "finesearch discordant edge found", i.chrom, str(msr[msi]), str(msr[msi + 1]), str(ebest[0]), ebest[1]
                            if (ebest[0].v1.chrom, ebest[0].v1.pos, ebest[0].v1.strand, ebest[0].v2.chrom, ebest[0].v2.pos, ebest[0].v2.strand) not in eiSet:
                                elist.append(ebest)
                                eilist.append(ebest)
                                eiSet.add((ebest[0].v1.chrom, ebest[0].v1.pos, ebest[0].v1.strand, ebest[0].v2.chrom, ebest[0].v2.pos, ebest[0].v2.strand))
                                if len(hg.interval_list([hg.interval(ebest[0].v2.chrom, ebest[0].v2.pos, ebest[0].v2.pos)]).intersection(ilist)) > 0:
                                    if (ebest[0].v2.chrom, ebest[0].v2.pos, ebest[0].v2.strand, ebest[0].v1.chrom, ebest[0].v1.pos, ebest[0].v1.strand) not in eiSet:
                                        eilist.append((breakpoint_edge(ebest[0].v2, ebest[0].v1), ebest[1]))
                                        eiSet.add((ebest[0].v2.chrom, ebest[0].v2.pos, ebest[0].v2.strand, ebest[0].v1.chrom, ebest[0].v1.pos, ebest[0].v1.strand))
                                elist.sort(key=lambda x: hg.absPos(x[0].v1.chrom, x[0].v1.pos) + 0.1*x[0].v1.strand)
                                eilist.sort(key=lambda x: hg.absPos(x[0].v1.chrom, x[0].v1.pos) + 0.1*x[0].v1.strand)
                else:
                    print "msv end not refined", str(msr[msi]), str(msr[msi + 1])
                    msve = [e for e in elist if e[0].v1.strand * (msr[msi].info['cn'] - msr[msi + 1].info['cn']) > 0 and abs(
                        e[0].v1.pos - msr[msi].end) < self.max_insert + ms_window_size0]

        if amplicon_name is not None:
            edge_file = open("%s_edges_cnseg.txt" % amplicon_name, 'w')
            for e in eilist:
                edge_file.write("%s\t%s\t%s\t%s\n" %
                                (str(e[0]), e[1], e[0].hom, e[0].hom_seq))
            edge_file.close()
        return eilist


    def construct_segment(self, v):
        cpos = v.pos - v.strand * self.max_insert / 2
        cprevious = v.pos
        cflag = v.pos
        while abs(cflag - cprevious) < self.window_size:
            cprevious = cpos
            cpos = cpos - v.strand * self.max_insert / 2
            drange = [cpos, cpos + v.strand * self.max_insert]
            drange.sort()
            dlist = [a for a in self.fetch(v.chrom,
                                                   drange[0], drange[1])]
            if (len(dlist) * self.read_length <
                self.min_coverage * self.max_insert):
                continue
            cflag = cprevious
            if abs(cprevious - v.pos) > self.max_insert:
                v1 = breakpoint_vertex(v.chrom, cprevious, v.strand)
                discordant_neighbors = self.get_discordant_neighbors(v1)
                if len(discordant_neighbors) > 0:
                    return v1
            v2 = breakpoint_vertex(v.chrom, cpos, -1 * v.strand)
            discordant_neighbors = self.get_discordant_neighbors(v2)
            if len(discordant_neighbors) > 0:
                return v2
        return None


    # Methods to find all intervals in amplicon                 
    def interval_neighbors(self, i, ilist=[], rdlist=[], t=0, gcc=False):
        i2 = self.interval_extend(i)
        # i2 = i
        # i2 = self.interval_extend(i, ilist, rdlist)
        ms_window_size0 = 10000
        ms_window_size1 = 300
        logging.info("#TIME " + '%.3f\t'%clock() + " Calculating coverage meanshift segmentation")
        msrlist = [self.get_meanshift(i2, ms_window_size0, ms_window_size1, gcc)]
        logging.info("#TIME " + '%.3f\t'%clock() + " Detecting breakpoint edges")
        edges = self.interval_discordant_edges(i2, ms=msrlist)
        edges = [(e[1], e[0]) for e in edges]
        edges.sort(reverse=True)
        edges = [(e[1], e[0]) for e in edges]
        ei = 0
        logging.info("#TIME " + '%.3f\t'%clock() + " Selecting neighbors")
        neighbors = hg.interval_list([])
        while len(neighbors) < 10 and ei < len(edges):
            covered = False
            for i3 in ilist + neighbors:
                if i3.chrom == edges[ei][0].v2.chrom and edges[ei][0].v2.pos >= i3.start and edges[ei][0].v2.pos <= i3.end:
                    ei += 1
                    covered = True
                    break
            if covered:
                continue
            found_neighbor = False
            for i3 in rdlist:
                if i3.chrom == edges[ei][0].v2.chrom and edges[ei][0].v2.pos >= i3.start and edges[ei][0].v2.pos <= i3.end:
                    n = i3
                    n = hg.interval(i3.chrom, i3.start, i3.end)
                    found_neighbor = True
            if not found_neighbor:
                if edges[ei][0].v2.strand < 0:
                    n = self.interval_extend(hg.interval(edges[ei][0].v2.chrom, edges[ei][0].v2.pos, min(hg.chrLen[hg.chrNum(edges[ei][0].v2.chrom)] - 1, edges[ei][0].v2.pos + self.max_insert)))
                else:
                    n = self.interval_extend(hg.interval(edges[ei][0].v2.chrom, max(0, edges[ei][0].v2.pos - self.max_insert), edges[ei][0].v2.pos))
            if found_neighbor or n.size() > self.max_insert + 2:
                n.info = edges[ei][1]
                neighbors.append(n)
            ei += 1
        neighbors.sort()
        mc = neighbors.merge_clusters(extend=ms_window_size0)
        for c in mc:
            c[0].info = sum([c1.info for c1 in c[1]])
        nn = hg.interval_list([c[0] for c in mc])
        for e in nn:
            logging.debug("#TIME " + '%.3f\t'%clock() + " interval_neighbors: edges %s %s" % (str(i), str(e)))
        return nn


    def interval_hops(self, i=None, ilist=[], rdlist=[], gcc=False, explore=True):
        if type(i) == list or type(i) == hg.interval_list:
            i1list = i
            i = i[0]
        else:
            i1list = hg.interval_list([i])
        logging.debug("#TIME " + '%.3f\t'%clock() + " interval_hops: init " + str(i))
        ms_window_size0 = 10000
        i2list = hg.interval_list([])
        for i2 in i1list:
            ii = self.interval_extend(i2)
            logging.debug("#TIME " + '%.3f\t'%clock() + " interval_hops: interval extend " + str(i2) + " " + str(ii))
            i2list.append(ii)
        seen_list = hg.interval_list([])
        unseen_list = [(0,ii) for ii in i2list]
        heapq.heapify(unseen_list)
        clist = hg.interval_list(i2list)
        clist = hg.interval_list([ii[0] for ii in i2list.merge_clusters(extend=1)])
        while len(seen_list) < 10 and len(unseen_list) > 0:
            icc = heapq.heappop(unseen_list)
            ic = icc[1]
            if explore == False and len(hg.interval_list([ic]).intersection(i2list)) == 0:
                seen_list.append(ic)
                continue
            logging.debug("#TIME " + '%.3f\t'%clock() + " interval_hops: check rd " + str(i) + " " + str(ic) + " " + str(len(hg.interval_list([ic]).intersection(rdlist))))
            if len(hg.interval_list([ic]).intersection(i2list)) == 0 and len(hg.interval_list([ic]).intersection(rdlist)) > 0:
                seen_list.append(ic)
                continue
            logging.debug("#TIME " + '%.3f\t'%clock() + " interval_hops: search new " + str(i) + " " + str(ic))
            logging.info("#TIME " + '%.3f\t'%clock() + " Searching new neighbors for interval: " + str(ic))
            icn = self.interval_neighbors(ic, clist, rdlist=rdlist, gcc=gcc)
            logging.debug("#TIME " + '%.3f\t'%clock() + " interval_hops: neighbors " + str(i) + " " + str(ic) + " " + str(len(icn)))
            for ic2 in icn:
                logging.info("#TIME " + '%.3f\t'%clock() + " New neighbor: %s (weight=%d)" % (str(ic2), ic2.info))
                contained = False
                for i2 in clist:
                    if i2.contains(ic2):
                        contained = True
                if contained:
                    continue
                if ic2.size() < 2 * ms_window_size0 and len(self.interval_discordant_edges(ic2)) < 2:
                    continue
                if explore or len(hg.interval_list([ic]).intersection(i2list)) > 0:
                    heapq.heappush(unseen_list, (-ic2.info, ic2))
                clist.append(ic2)
            seen_list.append(ic)
        retlist = hg.interval_list(i2list + seen_list)
        retlist = [r[0] for r in retlist.merge_clusters(extend=1)]
        return retlist


    def interval_amplified(self, i, filter_conserved=True, filter_small=True):
        if len(hg.interval_list([i]).intersection(hg.conserved_regions) + hg.interval_list([i]).intersection(hg.centromere_list)) > 0:
            return False
        ms_window_size = 10000
        num_w = 0
        num_high = 0
        if filter_small and i.size() < 2 * ms_window_size and len(self.interval_discordant_edges(i)) < 2:
            return False 
        wc = self.window_coverage(i, ms_window_size, exact=False)
        mc = self.median_coverage()
        if self.span_coverage:
            arm_coverage = self.median_coverage(refi=i)
        else:
            arm_coverage = mc
        for w in wc:
            num_w += 1
            # if w[1] > mc[0] + 3 * mc[2]:
            if self.sensitivems == False:
                if mc[0] < arm_coverage[0] and w[1] > max(arm_coverage[0] + 3 * mc[2] * math.sqrt(arm_coverage[0] / mc[0]), arm_coverage[0] + 3.0 * mc[0] / 2.0):
                    num_high += 1
                elif mc[0] >= arm_coverage[0] and w[1] > max(mc[0] + 3 * mc[2], 5.0 * mc[0] / 2.0):
                    num_high += 1
            else:
                if mc[0] < arm_coverage[0] and w[1] > arm_coverage[0] + 3 * mc[2] * math.sqrt(arm_coverage[0] / mc[0]):
                    num_high += 1
                elif mc[0] >= arm_coverage[0] and w[1] > mc[0] + 3 * mc[2]:
                    num_high += 1
        # wc_high = len([w for w in wc if w[1] > mc[1] + 3 * mc[2]])
        if num_high > num_w / 5:
            return True
        elif filter_small == False and i.size() < 2 * ms_window_size and len(self.interval_discordant_edges(i)) >= 2:
            return True
        else:
            return False

    def interval_extend(self, i, strand=0, i0=None):
        ms_window_size = 10000
        extend_size = max(i.size() / ms_window_size, 1)
        max_window_size = 300000000
        if strand >= 0:
            extend_right = 1
            right_size = extend_size
        else:
            extend_right = -1
            right_size = 0
        if strand <= 0:
            extend_left = 1
            left_size = extend_size
        else:
            extend_left = -1
            left_size = 0
        ic = copy.copy(i)
        while ic.size() < max_window_size and (extend_left >= 0 or extend_right >= 0):
            if extend_right >= 0:
                if right_size < 1:
                    extend_right = -1
                elif ic.end + right_size * ms_window_size > hg.chrLen[hg.chrNum(ic.chrom)]:
                    if self.interval_amplified(hg.interval(ic.chrom, ic.end, hg.chrLen[hg.chrNum(ic.chrom)]), filter_small=False):
                        ic.end = hg.chrLen[hg.chrNum(ic.chrom)]
                        extend_right = -1
                    else:
                        extend_right = 0
                        right_size = right_size / 2
                elif self.interval_amplified(hg.interval(ic.chrom, ic.end, ic.end + right_size * ms_window_size), filter_small=False):
                    ic.end = ic.end + right_size * ms_window_size
                    if extend_right == 1:
                        right_size = 2 * right_size
                    else:
                        right_size = right_size / 2
                        if right_size < 1:
                            # ic.end = min(ic.end + ms_window_size, hg.chrLen[hg.chrNum(ic.chrom)])
                            extend_right = -1
                else:
                    extend_right = 0
                    right_size = right_size / 2
            if extend_left >= 0:
                if left_size < 1:
                    extend_left = -1
                elif ic.start - left_size * ms_window_size <= 1:
                    if self.interval_amplified(hg.interval(ic.chrom, 1, ic.start), filter_small=False):
                        ic.start = 1
                        extend_left = -1
                    else:
                        extend_left = 0
                        left_size = left_size / 2
                elif self.interval_amplified(hg.interval(ic.chrom, ic.start - left_size * ms_window_size, ic.start), filter_small=False):
                    ic.start = ic.start - left_size * ms_window_size
                    if extend_left == 1:
                        left_size = 2 * left_size
                    else:
                        left_size = left_size / 2
                        if left_size < 1:
                            # ic.start = max(ic.end - ms_window_size, 1)
                            extent_left = -1
                else:
                    extend_left = 0
                    left_size = left_size / 2
        if self.interval_amplified(hg.interval(ic.chrom, max(0, ic.end - 2 * ms_window_size), min(ic.end + 2 * ms_window_size, hg.chrLen[hg.chrNum(ic.chrom)])), filter_small=False):
            ic.end = min(ic.end + 10 * ms_window_size, hg.chrLen[hg.chrNum(ic.chrom)])
        if self.interval_amplified(hg.interval(ic.chrom, max(ic.start - 2 * ms_window_size, 0), min(ic.start + 2 * ms_window_size, hg.chrLen[hg.chrNum(ic.chrom)])), filter_small=False):
            ic.start = max(ic.start - 10 * ms_window_size, 0)
        if strand >= 0:
            ide = self.interval_discordant_edges(hg.interval(ic.chrom, ic.end + 1, min(hg.chrLen[hg.chrNum(ic.chrom)], ic.end + ms_window_size)))
            for e in ide:
                if e[0].v1.strand == 1:
                    ic.end = min(ic.end + 2 * ms_window_size, hg.chrLen[hg.chrNum(ic.chrom)])
                    break
        if strand <= 0:
            ide = self.interval_discordant_edges(hg.interval(ic.chrom, max(0, ic.start - ms_window_size), ic.start - 1))
            for e in ide:
                if e[0].v1.strand == -1:
                    ic.start = max(ic.start - 2 * ms_window_size, 0)
                    break
        # if ic.size() > ms_window_size:
        logging.debug('#TIME %.3f\t interval_extend: %s, %s, %s' % (clock(), str(i), strand, str(ic)))
        return ic


    # Method to create breakpoint graph, find network flow and cycle decomposition
    def interval_filter_vertices(self, ilist0, gcc=False, adaptive_counts=True, eilist=None, amplicon_name=None, runmode='FULL'):
        ms_window_size0 = 10000
        ms_window_size1 = 300
        ilist0.sort()
        ilist = hg.interval_list([a[0] for a in ilist0.merge_clusters()])

        # finesearch edges near refined meanshifts and add to eilist, create vertices corresponding to all  meanshifts and uncovered meanshifts
        all_msv = []
        msv_diff = {}
        all_msv_nocover = []
        logging.info("#TIME " + '%.3f\t'%(clock() - self.tstart) + " Calculating coverage meanshift segmentation")
        msrlist = [self.get_meanshift(i, ms_window_size0, ms_window_size1, gcc) for i in ilist]
        logging.info("#TIME " + '%.3f\t'%(clock() - self.tstart) + " Detecting breakpoint edges")
        sensitive_elist = self.get_sensitive_discordant_edges(ilist, msrlist, eilist, ms_window_size0=ms_window_size0, ms_window_size1=ms_window_size1, adaptive_counts=adaptive_counts, amplicon_name=amplicon_name)
        eilist = sensitive_elist
        logging.info("#TIME " + '%.3f\t'%(clock() - self.tstart) + " Building breakpoint graph")
        for i, msr in zip(ilist, msrlist):
            elist = []
            for e in eilist:
                if hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos).intersects(i):
                    elist.append(e)
            ms_vlist = []
            msv_index = {}
            for msi in range(len((msr)) - 1):
                if msr[msi + 1].info['cn'] < msr[msi].info['cn']:
                    msv = breakpoint_vertex(i.chrom, msr[msi].end, 1)
                else:
                    msv = breakpoint_vertex(i.chrom, msr[msi].end + 1, -1)
                msv_diff[msv] = msr[msi + 1].info['cn'] - msr[msi].info['cn']
                ms_vlist.append(msv)
                msv_index[msv] = msi
            all_msv.append(ms_vlist)
            print "Meanshift", str(i), len(ms_vlist), ms_vlist
            sys.stdout.flush()
            msve_match = {}
            for msv in ms_vlist:
                msi = msv_index[msv]
                if msr[msi].info['end_refined']:
                    msve = [e for e in elist if e[0].v1.strand * (msr[msi + 1].info['cn'] - msr[msi].info['cn']) < 0 and abs(e[0].v1.pos - msr[msi].end) < self.max_insert + ms_window_size1]
                else:
                    msve = [e for e in elist if e[0].v1.strand * (msr[msi + 1].info['cn'] - msr[msi].info['cn']) < 0 and abs(e[0].v1.pos - msr[msi].end) < self.max_insert + ms_window_size0]
                if len(msve) > 0:
                    msve_match[msv] = msve
            msv_nocover = [msv for msv in ms_vlist if msv not in msve_match]
            all_msv_nocover.append(msv_nocover)
            print "Meanshift no cover", str(i), msv_nocover

        # setup graph for flow optimization
        ngvlist_full = []
        elist_full = []
        ms_addlist = []
        kce = defaultdict(lambda: 0) # number of concordant reads
        koe = defaultdict(lambda: 0.0) # number of reads mapping outside the interval
        kbpe = defaultdict(lambda: 0.0) # number of discordant reads across breakpoint edge
        new_graph = breakpoint_graph()
        s = new_graph.new_vertex(ilist[0].chrom, -1, -1)
        for i, msr, ms_vlist, msv_nocover in zip(ilist, msrlist, all_msv, all_msv_nocover):
            ngvlist = []
            elist = []
            for e in eilist:
                if hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos).intersects(i):
                    elist.append(e)

            # add vertices to new_graph
            ei = 0
            nei = 0
            msi = 0
            if (len(elist) == 0 or elist[ei][0].v1.strand == 1 or elist[ei][0].v1.pos > i.start):
                if len(msv_nocover) == 0 or msv_nocover[msi].strand == 1 or msv_nocover[msi].pos > i.start:
                    nv = new_graph.new_vertex(i.chrom, i.start, -1)
                    ne = new_graph.new_edge(s, nv)
                    koe[ne] = len(self.interval_crossing_arcs(i.chrom, i.start, i.start + self.max_insert, -1, ilist))
                else: #len(ms_vlist) > 0 and ms_vlist[0].strand == -1 and ms_vlist[0].pos > i.start + self.max_insert
                    nv = new_graph.new_vertex(i.chrom, msv_nocover[msi].pos, msv_nocover[msi].strand)
                    ne = new_graph.new_edge(s, nv)
                    koe[ne] = len(self.interval_crossing_arcs(i.chrom, msv_nocover[msi].pos, msv_nocover[msi].strand + self.max_insert, -1, ilist))
                    ms_addlist.append(msv_nocover[msi])
                    msi += 1
            else:
                nv = new_graph.new_vertex(i.chrom, elist[0][0].v1.pos, -1)
                oecount = len(self.interval_crossing_arcs(nv.chrom, nv.pos, nv.pos + self.max_insert, -1, ilist))
                if oecount >= (self.pair_support if not adaptive_counts else self.pair_support_count(nv.chrom, nv.pos, -1, meanshift=msrlist, sensitivems=False)):
                    ne = new_graph.new_edge(s, nv)
                    koe[ne] = len(self.interval_crossing_arcs(nv.chrom, nv.pos, nv.pos + self.max_insert, -1, ilist))
                ei += 1
            ngvlist.append(nv)
            vc = breakpoint_vertex(ngvlist[0].chrom, ngvlist[0].pos, ngvlist[0].strand)
            while ei < len(elist) or msi < len(msv_nocover):
                vp = vc
                vc_type = 'edge'
                if msi >= len(msv_nocover):
                    vc = elist[ei][0].v1
                    ei += 1
                elif ei >= len(elist):
                    vc = msv_nocover[msi]
                    vc_type = 'meanshift'
                    ms_addlist.append(msv_nocover[msi])
                    msi += 1 
                elif elist[ei][0].v1.pos < msv_nocover[msi].pos:
                    vc = elist[ei][0].v1
                    ei += 1
                elif elist[ei][0].v1.pos == msv_nocover[msi].pos and elist[ei][0].v1.strand < msv_nocover[msi].strand:
                    vc = elist[ei][0].v1
                    ei += 1
                else:
                    vc = msv_nocover[msi]
                    vc_type = 'meanshift'
                    ms_addlist.append(msv_nocover[msi])
                    msi += 1
                if (vc.pos == vp.pos and vc.strand <= vp.strand) or (vc.pos == vp.pos + 1 and vc.strand < vp.strand):
                    continue
                logging.debug("#TIME " + '%.3f\t'%clock() + "interval_filter vertices new: " + str(vc) + " " + vc_type)
                if vc.strand == 1:
                    if ngvlist[nei].strand == 1:
                        nvc_prime = new_graph.new_vertex(ngvlist[nei].chrom, ngvlist[nei].pos+1, -1)
                        # oecount = len(self.interval_crossing_arcs(nvc_prime.chrom, nvc_prime.pos, nvc_prime.pos + self.max_insert, -1, ilist))
                        # if oecount >= (self.pair_support if not adaptive_counts else self.pair_support_count(nvc_prime.chrom, nvc_prime.pos, -1, meanshift=zip(ilist, msrlist, cnlist), sensitivems=False)):
                        #     ne = new_graph.new_edge(s, nvc_prime)
                        #     koe[ne] = oecount
                        ce = self.concordant_edge(vp)
                        nce = new_graph.new_edge(ngvlist[nei], nvc_prime)
                        kce[nce] = len(ce[1])
                        ngvlist.append(nvc_prime)
                        nei += 1
                    nv = new_graph.new_vertex(vc.chrom, vc.pos, 1)
                    if vc_type == 'meanshift':
                        oecount = len(self.interval_crossing_arcs(nv.chrom, max(0, nv.pos - 2 * self.max_insert), nv.pos + 2 * self.max_insert, 1, ilist))
                    else:
                        oecount = len(self.interval_crossing_arcs(nv.chrom, max(0, nv.pos - self.max_insert), nv.pos, 1, ilist))
                    # if vc_type == 'meanshift' or oecount >= (self.pair_support if not adaptive_counts else self.pair_support_count(nv.chrom, nv.pos, 1, meanshift=zip(ilist, msrlist, cnlist))):
                    if vc_type == 'meanshift':
                        ne = new_graph.new_edge(s, nv)
                        koe[ne] = oecount
                    ngvlist.append(nv)
                    nei += 1
                else:
                    logging.debug("#TIME " + '%.3f\t'%clock() + "interval_filter vertices: adding reverse edge = " + str(vc))
                    if ngvlist[nei].strand == 1 and not (ngvlist[nei].chrom == vc.chrom and ngvlist[nei].pos == vc.pos - 1):
                        nvc_prime = new_graph.new_vertex(ngvlist[nei].chrom, ngvlist[nei].pos+1, -1)
                        oecount = len(self.interval_crossing_arcs(nvc_prime.chrom, nvc_prime.pos, nvc_prime.pos + self.max_insert, -1, ilist))
                        # if oecount >= (self.pair_support if not adaptive_counts else self.pair_support_count(nvc_prime.chrom, nvc_prime.pos, -1, meanshift=zip(ilist, msrlist, cnlist), sensitivems=False)):
                        #     ne = new_graph.new_edge(s, nvc_prime)
                        #     koe[ne] = oecount
                        ce = self.concordant_edge(vp)
                        nce = new_graph.new_edge(ngvlist[nei], nvc_prime)
                        kce[nce] = len(ce[1])
                        ngvlist.append(nvc_prime)
                        nei += 1
                    if ngvlist[nei].strand == -1:
                        nvc_prime = new_graph.new_vertex(vc.chrom, vc.pos-1, 1)
                        oecount = len(self.interval_crossing_arcs(nvc_prime.chrom, max(0, nvc_prime.pos - self.max_insert), nvc_prime.pos, 1, ilist))
                        # if oecount >= (self.pair_support if not adaptive_counts else self.pair_support_count(nvc_prime.chrom, nvc_prime.pos, 1, meanshift=zip(ilist, msrlist, cnlist), sensitivems=False)):
                        #     ne = new_graph.new_edge(s, nvc_prime)
                        #     koe[ne] = oecount
                        ngvlist.append(nvc_prime)
                        nei += 1
                    nv = new_graph.new_vertex(vc.chrom, vc.pos, -1)
                    if vc_type == 'meanshift':
                        oecount = len(self.interval_crossing_arcs(nv.chrom, max(0, nv.pos - 2 * self.max_insert), nv.pos + 2 * self.max_insert, -1, ilist))
                    else:
                        oecount = len(self.interval_crossing_arcs(nv.chrom, nv.pos, nv.pos + self.max_insert, -1, ilist))
                    # if vc_type == 'meanshift' or oecount >= (self.pair_support if not adaptive_counts else self.pair_support_count(nv.chrom, nv.pos, -1, meanshift=zip(ilist, msrlist, cnlist), sensitivems=False)):
                    if vc_type == 'meanshift':
                        ne = new_graph.new_edge(s, nv)
                        koe[ne] = oecount
                    ce = self.concordant_edge(vc)
                    nce = new_graph.new_edge(ngvlist[nei], nv)
                    kce[nce] = len(ce[1])
                    ngvlist.append(nv)
                    nei += 1
                # ei += 1
            if ngvlist[nei].strand == -1:
                nv = new_graph.new_vertex(i.chrom, i.end, 1)
                ne = new_graph.new_edge(s, nv)
                koe[ne] = len(self.interval_crossing_arcs(nv.chrom, nv.pos - self.max_insert, nv.pos, 1, ilist))
                ngvlist.append(nv)
                nei += 1
            elif ngvlist[nei].strand == 1 and ngvlist[nei].pos < i.end:
                nvc_prime = new_graph.new_vertex(ngvlist[nei].chrom, ngvlist[nei].pos+1, -1)
                oecount = len(self.interval_crossing_arcs(nvc_prime.chrom, nvc_prime.pos, min(hg.chrLen[hg.chrNum(nvc_prime.chrom)], nvc_prime.pos + self.max_insert), -1, ilist))
                if oecount >= (self.pair_support if not adaptive_counts else self.pair_support_count(nvc_prime.chrom, nvc_prime.pos, -1, meanshift=msrlist, sensitivems=False)):
                    ne = new_graph.new_edge(s, nvc_prime)
                    koe[ne] = oecount
                ce = self.concordant_edge(vp)
                nce = new_graph.new_edge(ngvlist[nei], nvc_prime)
                kce[nce] = len(ce[1])
                ngvlist.append(nvc_prime)
                nei += 1
                nv = new_graph.new_vertex(i.chrom, i.end, 1)
                ne = new_graph.new_edge(s, nv)
                koe[ne] = len(self.interval_crossing_arcs(nv.chrom, nv.pos - self.max_insert, nv.pos, 1, ilist))
                ngvlist.append(nv)
                nei += 1
            ngvlist_full = ngvlist_full + ngvlist
            elist_full = elist_full + elist
            print "MSstats", len(ms_vlist), len(ms_addlist)
            for msa in ms_vlist:
                print "MSadd", str(msa), msv_diff[msa], self.foldup_count(msa.chrom, msa.pos, msa.strand), msa in ms_addlist#, self.pair_support_count(msa.chrom, msa.pos, msa.strand, ms, True)
        for e0 in elist_full:
            e = e0[0]
            if len(ilist.intersection([hg.interval(e.v2.chrom, e.v2.pos, e.v2.pos)])) > 0 and e.v1.pos >= e.v2.pos:
                ne = new_graph.add_edge(e)
                logging.debug("#TIME " + '%.3f\t'%clock() + "interval_filter vertices: added edge e = " + str(e))
                logging.debug("#TIME " + '%.3f\t'%clock() + "interval_filter vertices: added edge ne = " + str(ne) + " " + ne.edge_type)
                logging.debug("#TIME " + '%.3f\t'%clock() + "interval_filter vertices: added edge ne, v1.elist = " + str(ne.v1) + " " + ','.join(map(str, ne.v1.elist)))
                logging.debug("#TIME " + '%.3f\t'%clock() + "interval_filter vertices: added edge ne, v2.elist = " + str(ne.v2) + " " + ','.join(map(str, ne.v2.elist)))
                if ne is None:
                    raise ValueError("ne is None:" + str(e) + " " + str(len(e0[1])) + '\n'+','.join(map(str, new_graph.vs.values())))
                kbpe[ne] = e0[1]
            elif len(ilist.intersection([hg.interval(e.v2.chrom, e.v2.pos, e.v2.pos)])) == 0:
                ne = new_graph.add_edge(breakpoint_edge(breakpoint_vertex(s.chrom, s.pos, s.strand), e.v1))
                koe[ne] = e0[1]
        for nei in range(1,len(ngvlist_full)):
            if ngvlist_full[nei].strand == 1:
                new_graph.new_edge(ngvlist_full[nei-1], ngvlist_full[nei], edge_type='sequence')
            # else:
            #     new_graph.new_edge(ngvlist[nei-1], ngvlist[nei])
        for e in koe:
            koe[e] = max(0.0001, koe[e])
        # set up all constants
        logging.info("#TIME " + '%.3f\t'%(clock() - self.tstart) + " Optimizing graph copy number flow")
        C = self.median_coverage()[0] / 2
        print "C (haploid coverage) = ", C
        G = new_graph
        seqlist = [e for e in new_graph.es.values() if e.edge_type == 'sequence']
        n = len(seqlist)
        l = [abs(e.v2.pos - e.v1.pos)+1 for e in seqlist]
        k = [len([a for a in self.fetch(e.v1.chrom, e.v1.pos, e.v2.pos)]) for e in seqlist]
        # kgcc = [self.interval_coverage(hg.interval(i.chrom, e.v1.pos, e.v2.pos), gcc=True) * (e.v2.pos - e.v1.pos) / self.read_length for e in seqlist]
        # k = kgcc
        kcc = [self.interval_coverage(hg.interval(e.v1.chrom, e.v1.pos, e.v2.pos)) * (e.v2.pos - e.v1.pos) for e in seqlist]
        ke = {}
        ke.update(kbpe)
        ke.update(kce)
        ke.update(koe)
        K = [len([a for a in self.fetch(e.v1.chrom, e.v1.pos, e.v2.pos)]) * self.read_length/(abs(e.v2.pos - e.v1.pos) + 1.0) for e in seqlist]
        # edge read count kbpe defined above
        bplist = [e for e in new_graph.es.values() if (e.edge_type == 'discordant' or e.edge_type == 'breakpoint')]
        m = len(bplist)
        bpdict = {bplist[bpi]: bpi for bpi in range(len(bplist))}
        print "########## len bplist", len(bplist), ";   ################ kbpe, kce, koe = ", len(kbpe), len(kce), len(koe) 

        # set up problem size and variable types and constraint types
        numvar = n + m #r0 + r1-rn + m breakpoint edges + 2 dummy variables g0, g1 for separable objective + te (m values)
        numcon = 2 * n #flow for each sequence edge + flows for source and sink + ratio f0 to readdepth and one quadratic constraint + 2 constraints for dummy variables g0, g1 + te(1..m)
        bkx = [mosek.boundkey.lo] * (n + m)
        blx = [0.0] * (n + m)
        bux = [float('Inf')] * (n + m)
        bkc = [mosek.boundkey.fx] * (2 * n)
        blc = [0.0] * (2 * n)
        buc = [0.0] * (2 * n)
        asub = []
        aval = []
        qsubi = {}
        qsubj = {}
        qval = {}

        # setting up constraints
        for i in range(n):
            subarr = [i]
            valarr = [1.0]
            for e in seqlist[i].v1.elist:
                if e.edge_type == 'sequence':
                    continue
                if n + bpdict[e] in subarr:
                    j = subarr.index(n + bpdict[e])
                    valarr[j] += -1.0
                else:
                    subarr.append(n + bpdict[e])
                    valarr.append(-1.0)
            asub.append(np.array(subarr))
            aval.append(np.array(valarr))
            subarr = [i]
            valarr = [1.0]
            for e in seqlist[i].v2.elist:
                if e.edge_type == 'sequence':
                    continue
                if n + bpdict[e] in subarr:
                    j = subarr.index(n + bpdict[e])
                    valarr[j] += -1.0
                else:
                    subarr.append(n + bpdict[e])
                    valarr.append(-1.0)
            asub.append(np.array(subarr))
            aval.append(np.array(valarr))
        # setting up objective
        opro = [mosek.scopr.log] * (n + m)
        oprjo = range(n + m)
        oprfo = [-1 * ki for ki in k] + [-1 * ke[e] for e in bplist]
        oprgo = [C * li / self.read_length for li in l] + [(self.max_insert) * C / 2 / self.read_length for e in bplist]
        oprho = [0.0001] * (n + m)
        opcj = {cj:C * l[cj] / self.read_length for cj in range(len(l))}
        for e in bpdict:
            opcj[n + bpdict[e]] = self.max_insert * C / 2 / self.read_length
        
        def streamprinter(msg): 
            sys.stdout.write (msg) 
            sys.stdout.flush()
        env = mosek.Env()
        task = env.Task(0,0)
        task.set_Stream (mosek.streamtype.log, streamprinter)
        task.appendcons(numcon)
        task.appendvars(numvar)

        if mosek_major_version >= 9:
            for j in range(numvar):
                task.putvarbound(j, bkx[j], blx[j], bux[j])
            for i in range(numcon):
                task.putconbound(i, bkc[i], blc[i], buc[i])

        else:
            for j in range(numvar):
                task.putbound(mosek.accmode.var, j, bkx[j], blx[j], bux[j])
            for i in range(numcon):
                task.putbound(mosek.accmode.con, i, bkc[i], blc[i], buc[i])

        for i in range(numcon):
            task.putarow(i, asub[i], aval[i])
        # for i in qsubi:
        #    task.putqconk(2 * n + 2, qsubi[i], qsubj[i], qval[i])
        print "Problem set up done. Optimizing .."
        for cj in opcj:
            task.putcj(cj, opcj[cj])
        task.putobjsense(mosek.objsense.minimize)
        task.putSCeval(opro, oprjo, oprfo, oprgo, oprho)
        task.optimize()
        res = [ 0.0 ] * numvar
        print "Solution summary"
        task.solutionsummary(mosek.streamtype.log)
        task.getsolutionslice(mosek.soltype.itr, mosek.solitem.xx, 0, numvar, res)
        print ( "Solution is: %s" % res )
        
        wehc = {}
        for msv_ilist in zip(all_msv, ilist):
            slist = hg.interval_list([hg.interval('\t'.join(map(str, [sq[0].v1.chrom, sq[0].v1.pos, sq[0].v2.pos, sq[1]]))) for sq in zip(seqlist, res)])
            slist.sort()
            msl= [msv_ilist[1].start] + [v.pos for v in msv_ilist[0]] + [msv_ilist[1].end]
            mslist = hg.interval_list([hg.interval(msv_ilist[1].chrom, msl[i], msl[i + 1]) for i in range(len(msl) - 1)])
            for msi in mslist:
                if len(hg.interval_list([msi]).intersection(slist)) == 0:
                    print 'MSnointersection', str(msi), msl
                    for s in slist:
                        print str(s)
                    print '============================='
                    for s in seqlist:
                        print str(s)
                    exit()
                elif sum([ap[0].intersection(ap[1]).size() for ap in hg.interval_list([msi]).intersection(slist)]) == 0:
                    print 'MS0intersection', str(msi)
                    exit()

        edge_code = defaultdict(lambda:'discordant', {'concordant':'concordant', 'source':'source'})

        graph_logger.info("SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped")
        for si in range(n):
            graph_logger.info('sequence\t' + '\t'.join(map(str, [seqlist[si].v1,  seqlist[si].v2, res[si], K[si], seqlist[si].v2.pos - seqlist[si].v1.pos, k[si]])))
            wehc[seqlist[si]] = float(res[si])
        graph_logger.info("BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence")
        for bpi in range(m):
            # print  edge_code[bplist[bpi].type()], str(bplist[bpi]), res[n + bpi], ke[bplist[bpi]], bplist[bpi].kmer_homology()
            graph_logger.info('\t'.join(map(str, [edge_code[bplist[bpi].type()], bplist[bpi], res[n + bpi], ke[bplist[bpi]], bplist[bpi].hom, bplist[bpi].hom_seq])))
            wehc[bplist[bpi]] = float(res[n + bpi])
        lenlist = len(ilist)
        if len(ilist0) >= 10:
            lenlist = len(ilist0)
        all_msv_cat = reduce(lambda x, y: x+y, all_msv, [])
        oncolist = ','.join(Set([a[1].info['Name'] for a in ilist.intersection(hg.oncogene_list)]))+','
        istr = ','.join([i.chrom + ':' + str(i.start) + '-' + str(i.end) for i in ilist])
        summary_logger.info("TotalIntervalSize = " + str(sum([a.size() for a in ilist])))
        summary_logger.info("AmplifiedIntervalSize = " + str(sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5])))
        if len([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]) > 0:
            summary_logger.info("AverageAmplifiedCopyCount = " + str(sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos) for si in range(n) if res[si] >= 2.5]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos + 1 for si in range(n) if res[si] >= 2.5])))
        else:
            summary_logger.info("AverageAmplifiedCopyCount = 2")
        summary_logger.info("#Chromosomes = " + str(len(Set([i.chrom for i in ilist]))))
        summary_logger.info("#SeqenceEdges = " + str(n))
        summary_logger.info("#BreakpointEdges = " + str(len(kbpe)))
        summary_logger.info("#CoverageShifts = " + str(len(all_msv_cat)))
        summary_logger.info("#MeanshiftSegmentsCopyCount>5 = " + str(len([v for v in msv_diff.values() if v > 5])))
        summary_logger.info("#Foldbacks = " + str(len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1])))
        summary_logger.info("#CoverageShiftsWithBreakpointEdges = " + str(len([msa for msa in all_msv_cat if msa in ms_addlist])))


        # Summary, #intervals, t talsize, size>2.5, AvgCoverage>2.5, #chromosomes, #sequenceedges, #breakpointedges, #meanshiftbreaks, #meanshift>5, #msfoldbackedges, #msfoldbackedges, #mswithoutbreakpoint, oncogenes, representativestring, #bpedgeswithcommonkmers
        if len([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]) > 0:
            # print '\t'.join(map(str, ["Summary:", lenlist, sum([a.size() for a in ilist]), sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]), sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos) for si in range(n) if res[si] >= 2.5]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]), len(Set([i.chrom for i in ilist])), n, len(kbpe), len(all_msv_cat), len([v for v in msv_diff.values() if v > 5]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if msa in ms_addlist]), oncolist, istr, len([e for e in kbpe if e.kmer_homology()])]))
            print '\t'.join(map(str, ["Summary:", lenlist, sum([a.size() for a in ilist]), sum([seqlist[si].v2.pos - seqlist[si].v1.pos + 1 for si in range(n) if res[si] >= 2.5]), sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos + 1) for si in range(n) if res[si] >= 2.5]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos + 1 for si in range(n) if res[si] >= 2.5]), len(Set([i.chrom for i in ilist])), n, len(kbpe), len(all_msv_cat), len([v for v in msv_diff.values() if v > 5]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if msa in ms_addlist]), oncolist, istr]))
            for i in ilist:
                if sum([seqlist[si].v2.pos - seqlist[si].v1.pos + 1 for si in range(n) if res[si] >= 5 and hg.interval(seqlist[si].v1.chrom, seqlist[si].v1.pos, seqlist[si].v2.pos).intersects(i)]) == 0:
                    print "IntervalAmplifiedSize: ", i.chrom, i.start, i.end, 0, 2
                    continue
                print "IntervalAmplifiedSize: ", i.chrom, i.start, i.end, sum([seqlist[si].v2.pos - seqlist[si].v1.pos + 1 for si in range(n) if res[si] >= 5 and hg.interval(seqlist[si].v1.chrom, seqlist[si].v1.pos, seqlist[si].v2.pos).intersects(i)]), sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos + 1) for si in range(n) if res[si] >= 5 and hg.interval(seqlist[si].v1.chrom, seqlist[si].v1.pos, seqlist[si].v2.pos).intersects(i)]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos + 1 for si in range(n) if res[si] >= 5 and hg.interval(seqlist[si].v1.chrom, seqlist[si].v1.pos, seqlist[si].v2.pos).intersects(i)])
        else:
            # print '\t'.join(map(str, ["Summary:", lenlist, sum([a.size() for a in ilist]), sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]), sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos) for si in range(n)]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n)]), len(Set([i.chrom for i in ilist])), n, len(kbpe), len(all_msv_cat), len([v for v in msv_diff.values() if v > 5]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if msa in ms_addlist]), oncolist, istr, len([e for e in kbpe if e.kmer_homology()])]))
            print '\t'.join(map(str, ["Summary:", lenlist, sum([a.size() for a in ilist]), sum([seqlist[si].v2.pos - seqlist[si].v1.pos + 1 for si in range(n) if res[si] >= 2.5]), sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos + 1) for si in range(n)]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos + 1 for si in range(n)]), len(Set([i.chrom for i in ilist])), n, len(kbpe), len(all_msv_cat), len([v for v in msv_diff.values() if v > 5]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if msa in ms_addlist]), oncolist, istr]))

        if runmode == 'BPGRAPH':
            return
        logging.info("#TIME " + '%.3f\t'%(clock() - self.tstart) + " Plotting SV View")

        interval_index = 1
        for i in ilist:
            cycle_logger.info("Interval\t" + '\t'.join([str(interval_index), i.chrom, str(i.start), str(i.end)]))
            interval_index += 1

        new_graph.cycle_decomposition(wehc, s)


    # Plot coverage, meanshift copy count estimates and discordant edges in interval
    def plot_segmentation(self, ilist, amplicon_name, segments=[], scale_list=[], eilist=None, font='small'):
        fighsize = 12
        figvsize = 5
        if font == 'large':
            matplotlib.rcParams.update({'font.size': 18})
            figvsize = 5.85
        if font == 'all_amplicons':
            matplotlib.rcParams.update({'font.size': 48})
            figvsize = 5.21
            fighsize = 24
        fig = plt.figure(figsize=(fighsize,figvsize))
        plt.subplots_adjust(left=73/1000.0, right=1-73/1000.0, bottom=1/4.0, top=1-1/10.0)
        # dpi = 300
        if font == 'large':
            plt.subplots_adjust(left=73/1000.0, right=1-73/1000.0, bottom=2.1/5.85, top=90/100.0)
        if font == 'all_amplicons':
            plt.subplots_adjust(left=73/1000.0, right=1-73/1000.0, bottom=1/5.21, top=95/100.0)

        dpi = 1000.0/fighsize
        gs = gridspec.GridSpec(2, 1, height_ratios=[8,2])
        if font == 'all_amplicons':
            gs = gridspec.GridSpec(2, 1, height_ratios=[5,4])
        ax = fig.add_subplot(gs[0,0])
        if font == 'large':
            plt.title(os.path.basename(amplicon_name), fontsize=28)
        elif font != 'all_amplicons':
            plt.title(os.path.basename(amplicon_name))
        # if font == 'all_amplicons':
        #     plt.title(os.path.basename(amplicon_name), fontsize=56)
        ax2 = ax.twinx()
        ax2.set_ylabel("CN")
        ax3 = fig.add_subplot(gs[1,0], sharex=ax)
        ax.set_xlim(0, 1)
        ax.set_ylabel("Coverage")
        ax.yaxis.set_label_coords(-0.05, 0.25)
        ax2.yaxis.set_label_coords(1.05, 0.33)
        if font == 'all_amplicons':
            ax.set_ylabel("")
            ax2.set_ylabel("")
        for b in ilist.offset_breaks():
            ax.axvline(b[0], linestyle=b[1], color='k')
            ax3.axvline(b[0], linestyle=b[1], color='k')

        cx = []
        wc = []

        elist_dict = {}
        max_edge = 4
        scale_max_cov = 0
        scale_max_ms = 0
        msrlist = [self.get_meanshift(i) if i.size() > 50000 else self.meanshift_segmentation(i, window_size=300) for i in ilist]
        sensitive_elist = self.get_sensitive_discordant_edges(
            ilist, msrlist, eilist, ms_window_size0=10000, ms_window_size1=300, adaptive_counts=True, amplicon_name=amplicon_name)
        eilist = sensitive_elist


        for i, msr in zip(ilist, msrlist):
            de = [e for e in eilist if e[0].v1.pos != -1 and hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos).intersects(i)]  # self.interval_discordant_edges(i)
            elist_dict[i] = de
            elist_dict[i].sort(key=lambda x: hg.absPos(x[0].v1.chrom, x[0].v1.pos) + 0.1*x[0].v1.strand)
            for e in eilist:
                eposlist = []
                if e[0].v1.pos != -1:
                    eposlist.append(hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos))
                if e[0].v2.pos != -1:
                    eposlist.append(hg.interval(e[0].v2.chrom, e[0].v2.pos, e[0].v2.pos))
                if len(scale_list) == 0 or len(hg.interval_list(eposlist).intersection(scale_list)) > 0:
                    max_edge = max(max_edge, e[1])

        for i in ilist:
            if i.size() > 1000000:
                wc_i = [w for w in self.window_coverage(i, 10000, exact=False)]
            else:
                wc_i = [w for w in self.window_coverage(i, 100, exact=False)]
            cx += [((i.chrom, (c[0].start + c[0].end)/2), c[1]) for c in wc_i]
            wc += wc_i

        cx0 = [c for c in cx if ilist.xpos(c[0][0], c[0][1]) is not None]
        ax.bar([ilist.xpos(c[0][0], c[0][1]) for c in cx0], [c[1] for c in cx0], 0.0001, zorder=1, edgecolor='0.7', linewidth=0001)
        cmax = max([c[1] for c in wc])

        covl = []
        for i, msr in zip(ilist, msrlist):
            for seg in msr:
                avg_cov = np.average([c[1] for c in cx0 if c[0][0] == seg.chrom and c[0][1] >= seg.start and c[0][1] <= seg.end])
                if len(scale_list) == 0 or len(hg.interval_list([i]).intersection(scale_list)) > 0:
                    covl += [c[1] for c in cx0 if c[0][0] == seg.chrom and c[0][1] >= seg.start and c[0][1] <= seg.end]
                    scale_max_cov = max(scale_max_cov, avg_cov)
                    if seg.info['cn'] != float('inf'):
                        scale_max_ms = max(scale_max_ms, seg.info['cn'])

                    else:
                        scale_max_ms = max(scale_max_ms, 2000)

                ax2.plot((ilist.xpos(seg.chrom, max(i.start, seg.start)), ilist.xpos(seg.chrom, min(i.end, seg.end))), (seg.info['cn'], seg.info['cn']), linewidth=4, color='k')
        
        logging.debug("Max cov, max ms scales set to: " + str(scale_max_cov) + " " + str(scale_max_ms))
        covl.sort()
        if len(covl) > 0:
            m95cov = covl[-(len(covl)/20)]

        if scale_max_cov > 0 and m95cov > scale_max_cov:
            scale_max_ms = scale_max_ms * m95cov / scale_max_cov
            scale_max_cov = scale_max_cov * m95cov / scale_max_cov


        y_scale = 3.0
        # y_scale = 2.5
        if font == 'all_amplicons':
            y_scale = 2.5
        if scale_max_cov > 0:
            ax.set_ylim(0.1, y_scale * scale_max_cov)
        else:
            (ymin, ymax) = ax.get_ylim()
            ax.set_ylim(ymin, y_scale * ymax)
        if scale_max_ms > 0:
            (ymin, ymax) = (0.1, scale_max_ms)
            ax2.set_ylim(0.1, y_scale * ymax)
        else:
            (ymin, ymax) = ax2.get_ylim()
            ax2.set_ylim(0.1, y_scale * ymax)

        for i in ilist:
            for el in elist_dict[i]:
                e = el[0]
                if ilist.xpos(e.v2.chrom, e.v2.pos) is None and ilist.xpos(e.v1.chrom, e.v1.pos) is None:
                    continue
                elif ilist.xpos(e.v2.chrom, e.v2.pos) is None:
                    ax2.axvline(ilist.xpos(e.v1.chrom, e.v1.pos), color=ecolor[e.type()], linewidth=4.0 * min(1, float(el[1])/max_edge), alpha=0.5, zorder=10)
                    ax2.plot((ilist.xpos(e.v1.chrom, e.v1.pos), ilist.xpos(e.v1.chrom, e.v1.pos) - 0.01 * e.v1.strand), (0, 0), linewidth=8.0*min(1, float(el[1])/max_edge), color=ecolor[e.type()])
                elif ilist.xpos(e.v1.chrom, e.v1.pos) is None:
                    ax2.axvline(ilist.xpos(e.v2.chrom, e.v2.pos), color=ecolor[e.type()], linewidth=4.0 * min(1, float(el[1])/max_edge), alpha=0.5, zorder=10)
                    ax2.plot((ilist.xpos(e.v2.chrom, e.v2.pos), ilist.xpos(e.v2.chrom, e.v2.pos) - 0.01 * e.v2.strand), (0, 0), linewidth=8.0*min(1, float(el[1])/max_edge), color=ecolor[e.type()])
                else:
                    xmid = (ilist.xpos(e.v1.chrom, e.v1.pos) + ilist.xpos(e.v2.chrom, e.v2.pos)) / 2
                    xdia = abs(ilist.xpos(e.v2.chrom, e.v2.pos) - ilist.xpos(e.v1.chrom, e.v1.pos))
                    ydia = (1.0 + xdia) * 3 * ymax
                    pseudo_edge = breakpoint_edge(breakpoint_vertex(e.v1.chrom, hg.absPos(e.v1.chrom, e.v1.pos), e.v1.strand), breakpoint_vertex(e.v1.chrom, hg.absPos(e.v2.chrom, e.v2.pos), e.v2.strand))
                    ee = Arc((xmid, 0), xdia, ydia, fill=False, linewidth=4.0 * min(1, float(el[1])/max_edge), color=ecolor[pseudo_edge.type()], zorder=4, theta1=0.1, theta2=180)
                    ax2.add_patch(ee)
                    ax2.plot((ilist.xpos(e.v1.chrom, e.v1.pos), ilist.xpos(e.v1.chrom, e.v1.pos) - 0.01 * e.v1.strand), (0, 0), linewidth=8.0*min(1, float(el[1])/max_edge), color=ecolor[pseudo_edge.type()])
                    ax2.plot((ilist.xpos(e.v2.chrom, e.v2.pos), ilist.xpos(e.v2.chrom, e.v2.pos) - 0.01 * e.v2.strand), (0, 0), linewidth=8.0*min(1, float(el[1])/max_edge), color=ecolor[pseudo_edge.type()])
        ax2.axhline(2.0, alpha=0.8, linewidth=0.5, color='r')

        gparity = 0
        ry = 0.60
        ty = 0.65
        ogene_width = 4
        if font == 'large':
            ry = 0.85
            ry = 0.87
            ry = 0.77
            ogene_width = 12
        ogene_plotted = []
        for i in ilist:
            glist = hg.interval_list([i]).intersection(hg.oncogene_list)
            ogene_plotted += [g[1].info['Name'] for g in glist]
            for g in glist:
                if font == 'large':
                    ty = 0
                elif font == 'all_amplicons':
                    if gparity == 0:
                        ty = -0.1
                        ty = -0.07
                    else:
                        ty = 0.20
                        ty = 0.3
                else:
                    if gparity == 0:
                        ty = 0
                    else:
                        ty = 0.37
                if font == 'large':
                    ax3.plot([ilist.xpos(i.chrom, max(g[1].start, i.start)), ilist.xpos(i.chrom, min(g[1].end, i.end))], [ry, ry], 'r-',    linewidth=ogene_width)
                    ax3.text((ilist.xpos(i.chrom, max(g[1].start, i.start)) + ilist.xpos(i.chrom, min(g[1].end, i.end)))/2.0, ty, g[1].info['Name'], horizontalalignment='center', verticalalignment='bottom', fontsize=28, zorder=4)
                elif font == 'all_amplicons':
                    ogene_width = 36
                    ax3.plot([ilist.xpos(i.chrom, max(g[1].start, i.start)), ilist.xpos(i.chrom, min(g[1].end, i.end))], [0.85, 0.85], 'r-', linewidth=ogene_width)
                    ax3.text((ilist.xpos(i.chrom, max(g[1].start, i.start)) + ilist.xpos(i.chrom, min(g[1].end, i.end)))/2.0, -.05 + 0.37 * gparity, g[1].info['Name'], horizontalalignment='center', verticalalignment='bottom', fontsize=48, zorder=4)
                else:
                    ax3.plot([ilist.xpos(i.chrom, max(g[1].start, i.start)), ilist.xpos(i.chrom, min(g[1].end, i.end))], [ry, ry], 'r-',    linewidth=ogene_width)
                    ax3.text((ilist.xpos(i.chrom, max(g[1].start, i.start)) + ilist.xpos(i.chrom, min(g[1].end, i.end)))/2.0, ty, g[1].info['Name'], horizontalalignment='center', verticalalignment='bottom')
                gparity = (gparity + 1) % 2
            for s in segments:
                if not i.intersects(s):
                    continue
                ss = i.intersection(s)
                ax3.add_patch(Rectangle([ilist.xpos(i.chrom, max(ss.start, i.start)), 0.65], ilist.xpos(i.chrom, min(ss.end, i.end)) - ilist.xpos(i.chrom, max(ss.start, i.start)), 0.25, fc=chrcolor[s.info[1]], ec='k'))
                if font == 'large':
                    ax3.text((ilist.xpos(i.chrom, max(ss.start, i.start)) + ilist.xpos(i.chrom, min(ss.end, i.end)))/2.0, 0 , s.info[0], horizontalalignment='center', verticalalignment='bottom', fontsize=28)
                elif font == 'large' or font == 'all_amplicons':
                    ax3.text((ilist.xpos(i.chrom, max(ss.start, i.start)) + ilist.xpos(i.chrom, min(ss.end, i.end)))/2.0, 0 , s.info[0], horizontalalignment='center', verticalalignment='bottom', fontsize=48)
                else:
                    ax3.text((ilist.xpos(i.chrom, max(ss.start, i.start)) + ilist.xpos(i.chrom, min(ss.end, i.end)))/2.0, 0.2+int(s[0])%2*0.15, s.info[0], horizontalalignment='center', verticalalignment='top')
                # ax3.text((xpos(max(s[1].start, i.start)) + xpos(min(s[1].end, i.end)))/2.0, 0.2+0%2*0.15, s[0], horizontalalignment='center', verticalalignment='top')

        if font == 'large' or font == 'all_amplicons':
            axyticks = ax.get_yticks()
            ax.set_yticks([0, axyticks[-1]])
            ax2yticks = ax2.get_yticks()
            ax2.set_yticks([0, ax2yticks[-1]])


        ax.xaxis.set_visible(False)
        ax2.xaxis.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_linewidth(4)
        ax3.tick_params('x', length=0, which='major')
        ax3.tick_params('x', length=5, which='minor')
        if font == 'all_amplicons':
            previous_chrom = ''
            chrom_index = 1
            interval_poslist = []
            for i in ilist:
                if i.chrom != previous_chrom:
                    chrom_index = 1
                else:
                    chrom_index += 1
                previous_chrom = i.chrom
                imin = ilist.xpos(i.chrom, i.start)
                imax = ilist.xpos(i.chrom, i.end)
                segname = ''
                if imax - imin > 0.2:
                    segname = i.chrom + '.' + str(chrom_index)
                elif imax - imin > 0.05:
                    segname = i.chrom.strip('chr') + '.' + str(chrom_index)
                elif imax - imin > 0.02:
                    segname = i.chrom.strip('chr')
                interval_poslist.append((segname, (imax + imin) / 2))
            ax3.xaxis.set_major_locator(ticker.FixedLocator([c[1] for c in interval_poslist]))
            ax3.xaxis.set_major_formatter(ticker.FixedFormatter([c[0] for c in interval_poslist]))
        else:
            chrmin = {}
            chrmax = {}
            for i in ilist:
                if i.chrom not in chrmin:
                    chrmin[i.chrom] = ilist.xpos(i.chrom, i.start)
                    chrmax[i.chrom] = ilist.xpos(i.chrom, i.end)
                else:
                    chrmin[i.chrom] = min(ilist.xpos(i.chrom, i.start), chrmin[i.chrom])
                    chrmax[i.chrom] = max(ilist.xpos(i.chrom, i.end), chrmax[i.chrom])
            chrposlist = []
            for c in chrmin:
                chrposlist.append((c if chrmax[c] - chrmin[c] > 0.10 else c.strip('chr'), (chrmin[c] + chrmax[c]) / 2))
            ax3.xaxis.set_major_locator(ticker.FixedLocator([c[1] for c in chrposlist]))
            ax3.xaxis.set_major_formatter(ticker.FixedFormatter([c[0] for c in chrposlist]))
        xposlist = []
        if font != 'all_amplicons':
            for i in ilist:
                xposlist.append((str(i.start), ilist.xpos(i.chrom, i.start)))
                xposlist.append((str(i.end), ilist.xpos(i.chrom, i.end)))
            ax3.xaxis.set_minor_locator(ticker.FixedLocator([c[1] for c in xposlist]))
            ax3.xaxis.set_minor_formatter(ticker.FixedFormatter([c[0] for c in xposlist]))
            plt.setp(ax3.xaxis.get_minorticklabels(), rotation=90)
            ax3.tick_params(axis='x', which='minor', pad=15)
        # ax3.tick_params(axis='x', which='minor', pad=-5)
        ax3.yaxis.set_major_formatter(ticker.NullFormatter())
        ax3.set_ylim(0,1)

        # ax3.spines['bottom'].set_visible(False)
        # ax3.xaxis.set_visible(False)

        fig.subplots_adjust(hspace=0)
        try:
            fig.savefig(amplicon_name + '.png', dpi=dpi)
            fig.savefig(amplicon_name + '.pdf', dpi=dpi)
        except np.linalg.linalg.LinAlgError:
            logging.error("Numpy LinAlgError when forming amplicon plot! Cannot save " + amplicon_name + " image\n")

        plt.close()
