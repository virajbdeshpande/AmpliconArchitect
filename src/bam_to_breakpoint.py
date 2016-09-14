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


from breakpoint_graph import *
import hg19util as hg


summary_logger = logging.getLogger('summary')
graph_logger = logging.getLogger('graph')


class breakpoint_cluster:
    def __init__(self, edge, bamfile, max_insert):
        self.edge = edge


class bam_to_breakpoint():
    def __init__(self, bamfile, read_length=100, max_insert=400, insert_size=300,
                 window_size=10000, min_coverage=30, pair_support=5,
                 secondary_index=None, coverage_stats=None, coverage_windows=None):
        self.bamfile = bamfile
        self.window_size = window_size
        self.min_coverage = min_coverage
        self.max_insert = max_insert
        self.insert_size = insert_size
        self.read_length = read_length
        self.secondary_index=secondary_index
        self.gc_scale = defaultdict(lambda: 1.0) #self.gc_scaling()
        self.gc_set = False
        self.ms_window_size = 10000
        self.pair_support = pair_support
#        self.pair_support = max(3 * (self.median_coverage()[0] * (self.max_insert - self.read_length) / 2 / self.read_length), 3)
        for s in self.bamfile.header['SQ']:
            hg.update_chrLen([(c['SN'], c['LN']) for c in self.bamfile.header['SQ']])
        if coverage_stats is None:
            self.basic_stats_set = False
            print clock()
            self.median_coverage(window_list=coverage_windows)
            print clock()
        else:
            (wc_10000_median, wc_10000_avg, wc_10000_std, wc_300_median, wc_300_avg, wc_300_std, self.read_length, self.insert_size, self.insert_std, self.min_insert, self.max_insert, self.pair_support, self.percent_proper) = coverage_stats
            self.basic_stats = coverage_stats
            self.basic_stats_set = True

    def concordant_edge(self, v, ms=None):
        if v.strand == 1:
            dlist = [a for a in self.bamfile.fetch(v.chrom,
                                                   max(1, v.pos - self.max_insert),
                                                   v.pos)
                     if not a.is_unmapped and not a.is_reverse and
                     a.is_proper_pair and a.mpos >= v.pos]
            if len(dlist) > self.pair_support:
                v2= breakpoint_vertex(v.chrom,
                                         max(v.pos + 1, min([a.mpos for a in dlist])), -1)
                return (breakpoint_edge(v, v2), dlist)
        else:
            dlist = [a for a in self.bamfile.fetch(v.chrom,
                                                   max(1, v.pos - self.max_insert),
                                                   v.pos)
                     if not a.is_reverse and a.is_proper_pair and
                     a.mpos >= v.pos]
            if len(dlist) > self.pair_support:
                v2 = breakpoint_vertex(v.chrom,
                                         min(v.pos - 1, max([a.aend for a in dlist])), 1)
                return (breakpoint_edge(v, v2), dlist)

        return (None, dlist)

    def pair_support_count(self, chrom, position, strand, meanshift, foldup=False):
#        str(hg.interval(f[0][0][0].chrom, f[0][0][0].start, f[0][-1][0].end)), f[1], f[2]
        for fi in range(len(meanshift)):
            f = meanshift[fi]
            if f[0][0][0].chrom == chrom and f[0][0][0].start <= position and position <= f[0][-1][0].end:
                c = f[2]
                cd = c / 4
                if strand == -1:
                    if f[0][0][0].start >= position - self.ms_window_size / 2:
                        if fi > 0 and f[1] % 2 == 1:
                            cm = meanshift[fi - 1][2]
                            if c > cm:
                                cd = c - cm
                    if f[0][-1][0].end <= position + self.ms_window_size / 2:
                        if fi < len(meanshift) - 1 and f[1] % 4 >= 2:
                            cm = meanshift[fi - 1][2]
                            if cm > c and cm - c < cd:
                                cd = cm - c
                else:                                                             
                    if f[0][0][0].start >= position - self.ms_window_size / 2:
                        if fi > 0 and f[1] % 2 == 1:
                            cm = meanshift[fi - 1][2]
                            if cm > c and cm - c < cd:
                                cd = cm - c
                    if f[0][-1][0].end <= position + self.ms_window_size / 2:
                        if fi < len(meanshift) - 1 and f[1] % 4 >= 2:
                            cm = meanshift[fi - 1][2]
                            if cm < c and c - cm < cd:
                                cd = c - cm
                mc = self.median_coverage()
                cd = max(1, cd)
                pcount = max(((2/3.0) * (cd - 3 * math.sqrt(cd / mc[0]) * mc[2]) * (self.insert_size - self.read_length) / 2 / self.read_length), 2)
                pmincount = (mc[0] + 3 * mc[2]) * (self.insert_size - self.read_length) / self.read_length / 2
                if pcount < pmincount:
                    pcount = pmincount
                return pcount 
            
    def foldup_count(self, chrom, position, strand, cdiff=-1):
        interval = hg.interval(chrom, max(1, position - self.ms_window_size), min(hg.chrLen[hg.chrNum(chrom)], position + self.ms_window_size))
        if strand == 1:
            dlist = [a for a in self.bamfile.fetch(interval.chrom,
                                                   interval.start,
                                                   interval.end)
                     if not a.is_unmapped and not a.is_reverse and a.is_paired
                     and not a.is_proper_pair and not a.mate_is_unmapped
                     and not a.mate_is_reverse and a.mrnm==a.tid
                     and abs(a.mpos - a.pos) < 100000]#self.ms_window_size]
        else:
            dlist = [a for a in self.bamfile.fetch(interval.chrom,
                                                   interval.start,
                                                   interval.end)
                     if not a.is_unmapped and a.is_reverse and a.is_paired
                     and not a.is_proper_pair and not a.mate_is_unmapped
                     and a.mate_is_reverse and a.mrnm==a.tid
                     and abs(a.mpos - a.pos) < 100000]#self.ms_window_size]
        return len(dlist)

    def interval_discordant_edges(self, interval, filter_repeats=True, pair_support=-1, adaptive_counts=False, ms=None):
        logging.info("#TIME " + str(clock()) + " discordant edges" + str(interval))
        if pair_support == -1:
            pair_support = self.pair_support
#        if adaptive_counts == True:
#            raise "Warning: adaptive counts not supported"
#            adaptive_counts = False
        if type(interval) != hg.interval_list:
            ilist = hg.interval_list([interval])
        else:
            ilist = interval
        mslist = []
        if adaptive_counts:
            if ms is None:
                mslist = []
                for i in ilist:
                    mslist.append(self.meanshift_segmentation(i))
                ms = mslist[0]
            elif type(interval) != hg.interval_list:
                mslist = [ms]
            else:
                mslist = ms
                ms = mslist[0]
        interval = ilist[0]
        dflist = []
        drlist = []
        for i in ilist:
            dflist += [a for a in self.bamfile.fetch(i.chrom,
                                                max(1, i.start),
                                                i.end)
                  if not a.is_unmapped and not a.is_reverse and a.is_paired
                  and not a.is_proper_pair and not a.mate_is_unmapped]
            drlist += [a for a in self.bamfile.fetch(i.chrom,
                                                max(1, i.start),
                                                i.end)
                  if not a.is_unmapped and a.is_reverse and a.is_paired
                  and not a.is_proper_pair and not a.mate_is_unmapped]
        logging.info("#TIME " + str(clock()) + " discordant edges: fetch discordant " + str(interval))

        hgddict = {hg.interval(a, bamfile=self.bamfile): a for a in drlist + dflist}
        hgdflist = hg.interval_list([hga for hga in hgddict if hga.strand==1])
        hgdrlist = hg.interval_list([hga for hga in hgddict if hga.strand==-1])
        hgdflist.sort()
        hgdrlist.sort()
        mcdflist = hgdflist.merge_clusters(self.max_insert - 2 * self.read_length)
        mcdrlist = hgdrlist.merge_clusters(self.max_insert - 2 * self.read_length)
        mcdflist = hg.interval_list([l[0] for l in mcdflist if len(l[1]) >= pair_support])
        mcdrlist = hg.interval_list([l[0] for l in mcdrlist if len(l[1]) >= pair_support])
        dflist = [a for a in dflist if hg.interval_list([hg.interval(a, bamfile=self.bamfile)]).intersection(mcdflist) > 0]
        drlist = [a for a in drlist if hg.interval_list([hg.interval(a, bamfile=self.bamfile)]).intersection(mcdflist) > 0]
        if filter_repeats:
            dflist = [a for a in dflist if hg.interval(a, bamfile=self.bamfile).rep_content() <= 3]
            drlist = [a for a in drlist if hg.interval(a, bamfile=self.bamfile).rep_content() <= 3]
        hgddict = {hg.interval(a, bamfile=self.bamfile): a for a in drlist + dflist}
        hgdflist = hg.interval_list([hga for hga in hgddict if hga.strand==1])
        hgdrlist = hg.interval_list([hga for hga in hgddict if hga.strand==-1])
        hgdflist.sort()
        hgdrlist.sort()
        mcdflist = hgdflist.merge_clusters(self.max_insert - 2 * self.read_length)
        mcdrlist = hgdrlist.merge_clusters(self.max_insert - 2 * self.read_length)
        mcdflist = [l for l in mcdflist if len(l[1]) >= pair_support]
        mcdrlist = [l for l in mcdrlist if len(l[1]) >= pair_support]
        if adaptive_counts:
            mcdfadaptive = []
            mcdradaptive = []
            for mc in mcdflist:
                for (ms, i) in zip(mslist, ilist):
                    if hg.interval(mc[0].chrom, mc[0].end, mc[0].end).intersects(i) and len(mc[1]) > self.pair_support_count(mc[0].chrom, mc[0].end, 1, ms):
                        mcdfadaptive.append(mc)
                        break
            for mc in mcdrlist:
                for (ms, i) in zip(mslist, ilist):
                    if hg.interval(mc[0].chrom, mc[0].start, mc[0].start).intersects(i) and len(mc[1]) > self.pair_support_count(mc[0].chrom, mc[0].start, 1, ms):
                        mcdradaptive.append(mc)
                        break
            mcdflist = mcdfadaptive
            mcdrlist = mcdradaptive
        else:
            mcdflist = [mc for mc in mcdflist if len(mc[1]) >= pair_support]
            mcdrlist = [mc for mc in mcdrlist if len(mc[1]) >= pair_support]
        dnlist = []
        clist = hg.interval_list([c[0] for c in mcdflist + mcdrlist])
        clist.sort()
        logging.info("#TIME " + str(clock()) + " discordant edges: merge and sort discordant " + str(interval))
        for c1 in mcdflist + mcdrlist:
            for c2 in mcdflist + mcdrlist:
                vl = []
                for a1 in c1[1]:
                    for a2 in c2[1]:
                        aq1 = hgddict[a1]
                        aq2 = hgddict[a2]
                        if aq1.qname == aq2.qname and aq1.is_read1 != aq2.is_read1:
                            if abs(aq1.pos - aq2.pos) < self.read_length and abs(aq1.aend - aq2.aend) < self.read_length and aq1.is_reverse != aq2.is_reverse:
                                continue
                            if aq1.rname == aq2.rname and  aq1.is_reverse and not aq2.is_reverse and aq1.pos - aq2.aend > 0 and aq1.pos - aq2.aend < self.max_insert - 2 * self.read_length:
                                continue
                            if aq2.rname == aq1.rname and aq2.is_reverse and not aq1.is_reverse and aq2.pos - aq1.aend > 0 and aq2.pos - aq1.aend < self.max_insert - 2 * self.read_length:
                                continue
                            vl.append((aq1, aq2))
                if len(vl) == 0 or len([v for v in vl if v[1].pos*v[0].pos > 0]) == 0:
                    continue
                if not vl[0][0].is_reverse:
                    bp1 = breakpoint_vertex(c1[0].chrom, max([v[0].aend for v in vl if v[0].pos > 0]), 1)
                else:
                    bp1 = breakpoint_vertex(c1[0].chrom, min([v[0].pos for v in vl if v[0].pos > 0]), -1)
                if not vl[0][1].is_reverse:
                    bp2 = breakpoint_vertex(c2[0].chrom, max([v[1].aend for v in vl if v[1].pos > 0]), 1)
                else:
                    bp2 = breakpoint_vertex(c2[0].chrom, min([v[1].pos for v in vl if v[1].pos > 0]), -1)
                if not adaptive_counts:
                    ps = pair_support
                else:
                    ps = self.pair_support_count(bp1.chrom, bp1.pos, bp1.strand, ms)
                if len(vl) < ps:
                    continue
                num_inverted = 0
                bp1c = None
                bp2c = None
                if bp1.chrom == bp2.chrom and bp1.pos == bp2.pos and bp1.strand == bp2.strand:
                    if bp1.strand == 1:
                        for v in vl:
                            if v[0].pos == v[1].pos:
                                num_inverted += 1
                    else:
                        for v in vl:
                            if v[0].aend == v[1].aend:
                                num_inverted += 1
                    if len(vl) - num_inverted < ps:
                        continue
                    vl.sort(lambda x, y: x[0].pos - y[0].pos)
                    if bp1.strand == 1:
                        maxp = vl[0][0].aend
                        maxn = 0
                        for v in vl[::-1]:
                            if len([v1 for v1 in vl if v1[0].aend <= v[0].aend and v1[0].pos > v[0].aend - self.max_insert + 2 * self.read_length]) > maxn:
                                maxp = v[0].aend
                                maxn = len([v1 for v1 in vl if v1[0].aend <= v[0].aend and v1[0].aend > v[0].aend - self.max_insert + 2 * self.read_length])
                        vl = [v for v in vl if v[0].aend <= maxp and v[0].aend > maxp -  self.max_insert + 2 * self.read_length]
                        if len(vl) < ps:
                            continue
                        bp1 = breakpoint_vertex(c1[0].chrom, max([v[0].aend for v in vl if v[0].pos > 0]), 1)
                        bp2 = breakpoint_vertex(c2[0].chrom, max([v[1].aend for v in vl if v[1].pos > 0]), 1)
                        if bp1.pos != bp2.pos:
                            bp1c = bp2
                            bp2c = bp1
                    else:
                        maxp = vl[-1][0].pos
                        maxn = 0
                        for v in vl:
                            if len([v1 for v1 in vl if v1[0].pos >= v[0].pos and v1[0].pos < v[0].pos + self.max_insert - 2 * self.read_length]) > maxn:
                                maxp = v[0].pos
                                maxn = len([v1 for v1 in vl if v1[0].pos >= v[0].pos and v1[0].pos < v[0].pos + self.max_insert - 2 * self.read_length])
                        vl = [v for v in vl if v[0].pos >= maxp and v[0].pos < maxp +  self.max_insert - 2 * self.read_length]
                        if len(vl) < ps:
                            continue
                        bp1 = breakpoint_vertex(c1[0].chrom, min([v[0].pos for v in vl if v[0].pos > 0]), -1)
                        bp2 = breakpoint_vertex(c2[0].chrom, min([v[1].pos for v in vl if v[1].pos > 0]), -1)
                        if bp1.pos != bp2.pos:
                            bp1c = bp2
                            bp2c = bp1
#                    print '###', str(bp1), str(bp2), len(vl), num_inverted
                bre = breakpoint_edge(bp1, bp2)
                if bre.type() != 'concordant':
                    dnlist.append((bre, vl))
                if bp1c is not None and bp2c is not None:
                    brec = breakpoint_edge(bp1c, bp2c)
                    if brec.type() != 'concordant':
                        dnlist.append((brec, [(v[1], v[0]) for v in vl]))
        logging.info("#TIME " + str(clock()) + " discordant edges: local edges done " + str(interval) + " " + str(len(mcdflist)) + " " + str(len(mcdrlist)))
        self.get_mates_time = 0
        self.get_mates_num_calls = 0
        for c in mcdflist + mcdrlist:
            nlist = []
            if filter_repeats:
                if len(hg.interval_list([c[0]]).intersection(hg.conserved_regions)) > 0:
                    continue
            logging.info("#TIME " + str(clock()) + " discordant edges: global edges -1" + str(c[0]))
            rep_content_time = 0
            intersection_time = 0
            nr_calls = 0
            for hga in c[1]:
                nmatelist = self.get_mates(hgddict[hga])
                if filter_repeats:
                    rpc = clock()
                    nmatelist = [a for a in nmatelist if hg.interval(a, bamfile=self.bamfile).rep_content() <= 3]
                    nr_calls += len(nmatelist)
                    rep_content_time += clock() - rpc
                    ict = clock()
                    nmatelist = [a for a in nmatelist if len(hg.interval_list([hg.interval(a, bamfile=self.bamfile)]).intersection(ilist)) == 0]
                    intersection_time += clock() - ict
                nlist += nmatelist
            logging.info("#TIME " + str(clock()) + " discordant edges: global edges 0" + str(c[0]) + " " + str(rep_content_time) + " " + str(nr_calls))
            nflist = [n for n in nlist if not n.is_reverse]
            nrlist = [n for n in nlist if n.is_reverse]
            hgndict = {hg.interval(a, bamfile=self.bamfile):a for a in nflist + nrlist}
            hgnflist = hg.interval_list([hga for hga in hgndict if hga.strand==1])
            hgnrlist = hg.interval_list([hga for hga in hgndict if hga.strand==-1])
            hgnflist.sort()
            hgnrlist.sort()
            mcnflist = hgnflist.merge_clusters(self.max_insert - 2 * self.read_length)
            mcnrlist = hgnrlist.merge_clusters(self.max_insert - 2 * self.read_length)
            mcnflist = [m for m in mcnflist if len(m[1]) >= self.pair_support] 
            mcnrlist = [m for m in mcnrlist if len(m[1]) >= self.pair_support] 
            mcnlist = mcnflist + mcnrlist
            logging.info("#TIME " + str(clock()) + " discordant edges: global edges 1" + str(c[0]) + " " + str(len(mcnflist)) + " " + str(len(mcnrlist)) + " " + str(self.get_mates_time) + " " + str(self.get_mates_num_calls))
            for cn in mcnlist:
                vl = []
                if filter_repeats:
                    if len(hg.interval_list([cn[0]]).intersection(hg.conserved_regions)) > 0:
                        continue
                hgmi = 0
                for hgm in cn[1]:
                    hgmi += 1
                    if filter_repeats:
                        if hgm.rep_content() > 3:
                            continue
                    for a in self.get_mates(hgndict[hgm]):
                        if filter_repeats:
                            if hg.interval(a, bamfile=self.bamfile).rep_content() > 3:
                                continue
                        if hg.interval(a, bamfile=self.bamfile).intersects(c[0]):
                            vl.append((a, hgndict[hgm]))
                            break
                if len(vl) == 0 or len([v for v in vl if v[1].pos*v[0].pos > 0]) == 0:
                    continue
                if not vl[0][0].is_reverse:
                    bp1 = breakpoint_vertex(self.bamfile.getrname(vl[0][0].tid),
                                            max([v[0].aend for v in vl]), 1)
                else:
                    bp1 = breakpoint_vertex(self.bamfile.getrname(vl[0][0].tid),
                                            min([v[0].pos for v in vl if v[0].pos > 0]), -1)
                if not vl[0][1].is_reverse:
                    bp2 = breakpoint_vertex(self.bamfile.getrname(vl[0][1].tid),
                                            max([v[1].aend for v in vl]), 1)
                else:
                    bp2 = breakpoint_vertex(self.bamfile.getrname(vl[0][1].tid),
                                            min([v[1].pos for v in vl if v[1].pos > 0]), -1)
                if not adaptive_counts:
                    ps = pair_support
                else:
                    ps = self.pair_support_count(bp1.chrom, bp1.pos, bp1.strand, ms)
                if len(vl) < ps:
                    continue
                num_inverted = 0
                if bp1.chrom == bp2.chrom and bp1.pos == bp2.pos and bp1.strand == bp2.strand:
                    if bp1.strand == 1:
                        for v in vl:
                            if v[0].pos == v[1].pos:
                                num_inverted += 1
                    else:
                        for v in vl:
                            if v[0].aend == v[1].aend:
                                num_inverted += 1
#                    print '###', str(bp1), str(bp2), len(vl), num_inverted
                if len(vl) - num_inverted < ps:
                    continue
                bre = breakpoint_edge(bp1, bp2)
                if bre.type() != 'concordant':
                    dnlist.append((bre, vl))
            logging.info("#TIME " + str(clock()) + " discordant edges: global edges 2" + str(c[0]) + " " + str(len(mcnflist)) + " " + str(len(mcnrlist)) + " " + str(self.get_mates_time) + " " + str(self.get_mates_num_calls))
        logging.info("#TIME " + str(clock()) + " discordant edges: external edges done " + str(interval) + " " + str(self.get_mates_time) + " " +  str(self.get_mates_num_calls))
        dnlist.sort(key=lambda x: hg.absPos(x[0].v1.chrom, x[0].v1.pos) + 0.5 * x[0].v1.strand)
        print str(interval).strip(), len(dnlist)
        for e in dnlist:
            print '#', e[0], len(e[1]), e[0].type(), self.concordant_edge(e[0].v1, ms)[0], + len(self.concordant_edge(e[0].v1, ms)[1]), hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos - e[0].v1.strand * self.max_insert).rep_content()
#            if hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos - e[0].v1.strand * self.max_insert).rep_content() > 3 or (e[0].type() == 'everted' and (e[0].v1.pos - e[0].v2.pos) <30):
#                for a in e[1]:
#                    print '##', a[0].query_name, a[0].is_reverse, str(hg.interval(a[0], bamfile=self.bamfile)), '##', a[1].query_name, a[1].is_reverse, str(hg.interval(a[1], bamfile=self.bamfile)), hg.interval(a[0], bamfile=self.bamfile).rep_content()
        return dnlist

    def construct_segment(self, v):
        cpos = v.pos - v.strand * self.max_insert / 2
        cprevious = v.pos
        cflag = v.pos
        while abs(cflag - cprevious) < self.window_size:
            cprevious = cpos
            cpos = cpos - v.strand * self.max_insert / 2
            drange = [cpos, cpos + v.strand * self.max_insert]
            drange.sort()
            dlist = [a for a in self.bamfile.fetch(v.chrom,
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

    def get_mates(self, a):
        gmt = clock()
        self.get_mates_num_calls += 1
        try:
            miter = self.secondary_index.find(a.qname)
            retval = [m for m in miter if m.is_read1 != a.is_read1]
            self.get_mates_time += clock() - gmt
            return retval
        except:
#            print clock(), 'get_mates', str(a)
            retval = [a2 for a2 in self.bamfile.fetch(a.next_reference_name, a.next_reference_start, a.next_reference_start + 1) if a2.qname == a.qname]
#            retval = [self.bamfile.mate(a)]
#            print clock(), 'got_mates'
            self.get_mates_time += clock() - gmt
            return retval

    def interval_coverage(self, i, clip=False, gcc=False):
        if gcc:
            wc_raw = self.window_coverage(i)
            wc_corrected = 0
            j = 0
            for w in wc_raw:
                alist = [a for a in self.bamfile.fetch(w[0].chrom, w[0].start, w[0].end)]
                wc_corrected += w[0].size() * w[1] / self.gc_scaling()[int(w[0].gc_content() * 10)  / 10.0]
                if w[0].size() * w[1] / self.gc_scaling()[int(w[0].gc_content() * 10)  / 10.0]  > 10 * len(alist) * self.read_length:
                    print str(i).strip(), str(w[0]).strip(), wc_corrected, len(alist), w[1], self.gc_scaling()[int(w[0].gc_content() * 10)  / 10.0], w[0].gc_content(), w[0].sequence()
                    j += 1
                    if j > 100:
                        raise ValueError("j>100")
                        exit()
            return wc_corrected / i.size()
        s2 = i.start
        e2 = i.end
        if i.start < 0:
            s2 = 0
        if i.end > hg.chrLen[hg.chrNum(i.chrom)]:
            e2 = hg.chrLen[hg.chrNum(i.chrom)]
        if s2 >= e2:
            return 0
        alist = [a for a in self.bamfile.fetch(i.chrom, s2, e2)
                 if not a.is_unmapped]
#        if e2 - s2 >= window_size and clip == False and gcc == False:
#            return len(alist) * self.read_length / float(e2 - s2)
        sumb = 0
#        if clip == True or (clip == False and i.size() >= 100 * self.read_length):
#            return len(alist) * self.read_length / float(i.size())
        if clip == True:
            return sum([sum(a) for a in self.bamfile.count_coverage(i.chrom, s2, e2)]) / float(e2-s2)
        else:
            return len([a for a in alist if a.aend <= e2]) * self.read_length / float(e2-s2)
        for a in alist:
            ai = hg.interval(a, bamfile=self.bamfile).intersection(i)
            if ai is not None:
                sumb += ai.size()
        if sumb / float(i.size()) > 10 * len(alist) * self.read_length / float(i.size()):
            print str(i), sumb, len(alist)
            raise ValueError("isize exception")
            exit()
        return sumb / float(i.size())

    def window_coverage_stats(self, i, window_size=-1, gcc=False):
        if window_size == -1:
            window_size = self.max_insert - self.read_length
        j = range(i.start, i.end, window_size)
        jj = [hg.interval(i.chrom, k, k + window_size) for k in j]
        cc = [self.interval_coverage(k, gcc=gcc) for k in jj]
        dd = [abs(cc[j + 1] - cc[j]) for j in range(len(jj) - 1)]
        return (sum(cc)/len(cc), sum(dd)/len(dd))

    def window_coverage(self, i, window_size=-1, gcc=False):
#        print str(i)
        if window_size == -1:
            window_size = self.max_insert - self.read_length
        def win_breakup(i, window_size):
            for k in xrange(i.start, i.end, window_size):
                yield hg.interval(i.chrom, k, k + window_size - 1)
        for k in win_breakup(i, window_size):
            yield (k, self.interval_coverage(k, gcc=gcc))
#        return [(k, self.interval_coverage(k, gcc=gcc)) for k in jj]

    def median_coverage(self, window_size=-1, gcc=False, refi=-1, window_list=None):
        if (window_size == 10000 or window_size == -1) and self.basic_stats_set and refi == -1:
            return self.basic_stats
        if window_size == 300 and self.basic_stats_set and refi == -1:
            return self.basic_stats[3:6]

        num_iter = 1000
        iteri = 0
        sumchrLen = sum([l for l in hg.chrLen.values()])

        if not self.basic_stats_set:
            read_length = []
            insert_size = []
            window_list_index = 0
            while (window_list is not None and window_list_index < len(window_list)) or (window_list is None and iteri <= num_iter):
                if window_list is None:
                    newpos = random.random() * sumchrLen
                else:
                    cwindow = window_list[window_list_index]
                    window_list_index += 1
                    if cwindow.end - cwindow.start < 10000:
                        continue
                    newpos = hg.absPos(cwindow.chrom, ((cwindow.end + cwindow.start) / 2) - 5000)
                (c,p) = hg.chrPos(newpos)
                if c not in self.bamfile.references or p < 10000 or hg.chrLen[hg.chrNum(c)] < p + 10000 or len(hg.interval_list([hg.interval(c, p, p+10000)]).intersection(hg.conserved_regions, extend=10000)) > 0:
                    continue
                read_length += [a.infer_query_length(always=False) for a in self.bamfile.fetch(c, p, p+10000)]
                insert_size += [a.template_length for a in self.bamfile.fetch(c, p, p+10000) if a.is_proper_pair and not a.is_reverse and a.template_length < 10000 and a.template_length > 0]
                iteri += 1
            self.read_length = np.average(read_length)
            self.insert_size = np.average(insert_size)
            percent_proper = len(insert_size)*2.0/len(read_length)
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
                    newpos = random.random() * sumchrLen
                else:
                    cwindow = window_list[window_list_index]
                    window_list_index += 1
                    if cwindow.end - cwindow.start < 10000:
                        continue
                    newpos = hg.absPos(cwindow.chrom, ((cwindow.end + cwindow.start) / 2) - 5000)
                (c,p) = hg.chrPos(newpos)
                if c not in self.bamfile.references or p < ws or hg.chrLen[hg.chrNum(c)] < p + ws or len(hg.interval_list([hg.interval(c, p, p+ws)]).intersection(hg.conserved_regions, extend=ws)) > 0:
                    continue
                wc_ws.append(self.interval_coverage(hg.interval(c,p, p+ws), gcc=gcc))
                iteri += 1
            wc_ws.sort()
            wc_ws_median = np.median(wc_ws)
            wc_ws_filter = [c for c in wc_ws if c < 5 * wc_ws_median and c > 0]
            if len(wc_ws_filter) == 0:
                print len(wc_ws_filter), len(wc_ws), len([c for c in wc_ws if c > 0]), wc_ws_median
            wc_median.append(wc_ws_filter[len(wc_ws_filter) / 2])
            wc_avg.append(np.average(wc_ws_filter))
            wc_std.append(np.std(wc_ws_filter))

        if window_size not in [-1, 300, 10000]:
            return (wc_median[0], wc_avg[0], wc_std[0])
        
        (wc_10000_median, wc_10000_avg, wc_10000_std) = (wc_median[0], wc_avg[0], wc_std[0])
        (wc_300_median, wc_300_avg, wc_300_std) = (wc_median[1], wc_avg[1], wc_std[1])
        self.pair_support = max((max(wc_300_avg - 3*wc_300_std, wc_300_avg/10.0)  * (self.insert_size - self.read_length) / 2 / self.read_length)*percent_proper, 3)
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
        return rstats


    def gc_scaling(self):
        if self.gc_set:
            return self.gc_scale
        gc_total_rd = {i / 10.0:0 for i in range(11)}
        gc_num_windows = {i / 10.0:0 for i in range(11)}
        for ri in range(len(self.bamfile.references)):
#            print self.bamfile.references[ri]
#            print hg.chrLen
            if hg.chrNum(self.bamfile.references[ri]) not in hg.chrLen:
                continue
#            print self.bamfile.references[ri]
            wc = self.window_coverage(hg.interval(self.bamfile.references[ri], 0, self.bamfile.lengths[ri]))
            lwc = 0
            for w in wc:
                lwc += 1
                gc_total_rd[int(w[0].gc_content() * 10.0) / 10.0] += w[1]
                gc_num_windows[int(w[0].gc_content() * 10.0) / 10.0] += 1
#            print gc_num_windows, gc_total_rd, lwc
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

    def meanshift(self, i, window_size=-1, hb=2, cov=None, rd_global=-1, h0=-1, gcc=False, n=-1):
        if window_size == -1:
           window_size = self.max_insert - self.read_length
        if rd_global == -1:
            rd_global = self.median_coverage(window_size, gcc)[1]
        if h0 == -1:
            h0 = self.median_coverage(window_size, gcc)[2]
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
    #        j = range(i.start, i.end, window_size)
            jj = [hg.interval(i.chrom, k, k + window_size) for k in j]
            i2 = hg.interval(i.chrom, s2, e2)
            cov = [c for c in self.window_coverage(i2, window_size, gcc)]
            #cov = [self.interval_coverage(k) for k in jj]
#        print window_size, len(cov), str(cov[0][0]).strip(), cov[0][1], str(cov[1][0]).strip(), cov[1][1]
        def hr(wi):
            if cov[wi][1] < rd_global / 4.0:
                return h0 / 2.0
            else:
#                return math.sqrt(cov[wi][1] * self.read_length / window_size)
                return math.sqrt(cov[wi][1] / rd_global) * h0
        dfi = [(cov[wi][0], sum([wj * math.exp(-0.5 * wj ** 2 / hb ** 2) * math.exp(-0.5 * (cov[wi + wj][1] - cov[wi][1])** 2 / hr(wi)**2) for wj in range(-1 * n, n + 1)])) for wi in range(n, len(j) - n)]
#        print 'meanshift', str(i), len(cov), len(dfi), len([c for c in cov if c[0] is not None]), [(str(c[0][0]), c[0][1], c[1][1]) for c in zip([cc for cc in cov if cc[0] is not None], dfi) ]
        #[(interval,ms)]
        return dfi


    def meanshift_segmentation(self, i, window_size=-1, gcc=False, pvalue=0.01):
        if window_size == -1:
            window_size = 10000
        mc = self.median_coverage(window_size, gcc)
        rd_global = mc[0]
        h0 = mc[2]
        hb_profile = [2, 5, 10, 50, 100]
#        hb_profile = [2]
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
            e2 = hgl - hgl % window_size + e2 % window_size
            endskip = n - (e2 - i.end) / window_size
        i2 = hg.interval(i.chrom, s2, e2)
        cov = [c for c in self.window_coverage(i2, window_size, gcc)]
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
#                print 'THIS0', len(frozen), fi #, frozen[0][0][0].start
#                for ff in range(len(frozen)):
#                    print "THIS", ff, frozen[ff][0][0][0].start
                while msi < len(ms):
                    if fi >= 0 and fi < len(frozen) and ms[msi][0].start == frozen[fi][0][0][0].start:
#                        print "frozen segment", fi
                        if len(new_seg) > 0 and (frozen[fi][1] % 2 == 1 or (ms[msi][1] > 0 and ms[msi -1] <= 0)):
                            segs.append(new_seg)
                            new_seg = ms[msi: msi + len(frozen[fi][0])]
                        else:
                            new_seg += ms[msi: msi + len(frozen[fi][0])]
#                        segs.append(ms[msi: msi + len(frozen[fi][0])])
                        #len(segs[-1]), len(frozen[fi][0]), segs[-1][-1][0].end, frozen[fi][0][-1][0].end
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
#                if segs[si][0][0].start < 54857402 and segs[si][-1][0].end > 54857402:
#                    print (segs[si][0][0].start, segs[si][-1][0].end), (segs[si-1][0][0].start, segs[si-1][-1][0].end), (segs[si+1][0][0].start, segs[si+1][-1][0].end)
#                    print stats.ttest_ind([cc[1] for cc in cov[ci:ci + len(segs[si])]], [cs[1] for cs in cov[ci - len (segs[si - 1]):ci]], equal_var=False) 
#                    print abs(cp - c), 3 * math.sqrt(max(cp, c) / rd_global) * h0
#                    print abs(cn - c), 3 * math.sqrt(max(cn, c) / rd_global) * h0
#                    print  [cs[1] for cs in cov[ci - len (segs[si - 1]):ci]]
#                    print [cc[1] for cc in cov[ci:ci + len(segs[si])]]
#                    print  [cs[1] for cs in cov[ci + len (segs[si]):ci + len(segs[si]) + len(segs[si + 1])]]
                if si > 0:
                    if (len(segs[si]) < 15 or len(segs[si - 1]) < 15):
                        cp = cov2[ci - 1][1]
                        if abs(cp - c) > 3 * math.sqrt(max(cp, c) / rd_global) * h0:
#                        if abs(cp - c) > 2 * hr(c, window_size * len(segs[si])):
                            freeze |= 1
                    if stats.ttest_ind([cc[1] for cc in cov[ci:ci + len(segs[si])]], [cs[1] for cs in cov[ci - len (segs[si - 1]):ci]], equal_var=False)[1] < pvalue:
                        freeze |= 1
                if si < len(segs) - 1:
                    if (len(segs[si]) < 15 or len(segs[si + 1]) < 15):
                        cn = cov2[ci + len(segs[si])][1]
                        if abs(cn - c) > 3 * math.sqrt(max(cn, c) / rd_global) * h0:
#                        if abs(cn - c) > 2 * hr(c, window_size * len(segs[si])):
                            freeze |= 2
                    if stats.ttest_ind([cc[1] for cc in cov[ci:ci + len(segs[si])]], [cs[1] for cs in cov[ci + len (segs[si]):ci + len(segs[si]) + len(segs[si + 1])]], equal_var=False)[1] < pvalue:
                        freeze |= 2
#                if freeze > 0:
                frozen.append((segs[si], freeze, c, cov2[cpi:ci+len(segs[si])], cov[cpi:ci+len(segs[si])]))
                ci += len(segs[si])
                if freeze > 0:
                    cpi = ci
#        for f in frozen:
#            print str(hg.interval(f[0][0][0].chrom, f[0][0][0].start, f[0][-1][0].end)), f[1], f[2], str(i)
#        print '----...-------------...------------...-------------...-----------------...----------------------'
        #(list of windows[(windowinterval,ms)], left%2/right%4freeze, avg_coverage)

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
#        for a in shifts:
#           print a[0], a[1], a[2], len(a[3]), len(a[4]), str(i)
#        print '---------------------------------------------------------'

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
#                print [[stats.ttest_ind(s3i, s4i, equal_var=False)[1] for s3i in s3] for s4i in s4]
                if min([min([stats.ttest_ind(s3i, s4i, equal_var=False)[1] for s3i in s3]) for s4i in s4]) > pvalue:
#                    print shifts[shiftsi]
                    mergelist.append(shiftsi)
#            print mergelist
#            exit()
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
#            for a in shifts:
#                print a[0], a[1], a[2], len(a[3]), len(a[4])
#            print '---------------------------------------------------------'
        return shifts

    def meanshift_refined(self, i, window_size0=10000, window_size1=300, gcc=False):
        if hg.chrLen[hg.chrNum(i.chrom)] < 3 * window_size0:
            return self.meanshift_segmentation(i, window_size1, gcc)
        shifts0 = self.meanshift_segmentation(i, window_size0, gcc, pvalue=0.0027)
        shifts1 = reduce(lambda x,y: x+y, [self.meanshift_segmentation(hg.interval(i.chrom, s[0] - 3 * window_size0, s[0] + 3 * window_size0), window_size1, gcc, pvalue=0.05) for s in shifts0], [])

#        shifts1 = self.meanshift_segmentation(i, window_size1, gcc)

        matched_shifts = []
        for s0 in shifts0:
            S1 = [s1 for s1 in shifts1 if abs(s0[0] - s1[0]) < window_size0]
            S1_select = [s1 for s1 in S1 if (s0[2]-s0[1])*(s1[2] - s1[1]) > 0 and
                                            (s0[2] - s0[1])/(s1[2] - s1[1]) > 0.5 and
                                            (s0[2] - s0[1])/(s1[2] - s1[1]) < 2]
            bests1 = None
            bestscore = 0
            for s1 in S1_select:
                if bests1 is None:
                    bests1 = s1
                    bestscore = abs((s0[2] - s0[1]) - (s1[2] - s1[1]))
                elif abs((s0[2] - s0[1]) - (s1[2] - s1[1])) < bestscore:
                    bestscore = abs((s0[2] - s0[1]) - (s1[2] - s1[1]))
                    bests1 = s1
#            print 'meanshift refined:', (s0, S1_select, S1, bests1, bestscore)
            if bests1 is None:
                matched_shifts.append((s0[0], s0[1], s0[2], False))
                print  (s0[0], s0[1], s0[2], False)
            else:
                matched_shifts.append((bests1[0], s0[1], s0[2], True))
                print (bests1[0], s0[1], s0[2], True)
        return matched_shifts



    def interval_crossing_arcs(self, chrom, start, end, strand, ilist):
            if strand == -1:
                return [a for a in self.bamfile.fetch(chrom, max(0, start), min(end, hg.chrLen[hg.chrNum(chrom)]))
                        if not a.is_unmapped and a.is_reverse
                        and len(ilist.intersection([hg.interval(self.bamfile.getrname(a.mrnm), a.pnext, a.pnext)])) == 0]
            else:
                return [a for a in self.bamfile.fetch(chrom, max(0, start), min(end, hg.chrLen[hg.chrNum(chrom)]))
                        if not a.is_unmapped and not a.is_reverse 
                        and len(ilist.intersection([hg.interval(self.bamfile.getrname(a.mrnm), a.pnext, a.pnext)])) == 0]


    def interval_filter_vertices(self, ilist0):
        ms_window_size0 = 10000
        ms_window_size1 = 300
        ilist0.sort()
        ilist = hg.interval_list([a[0] for a in ilist0.merge_clusters()])

        #finesearch edges near refined meanshifts and add to eilist, create vertices corresponding to all  meanshifts and uncovered meanshifts
        all_msv = []
        msv_diff = {}
        msrlist = [self.meanshift_refined(i) for i in ilist]
        all_msv_nocover = []
        eilist = self.interval_discordant_edges(ilist, adaptive_counts=False, ms=msrlist)
        eilist.sort(key=lambda x: hg.absPos(x[0].v1.chrom, x[0].v1.pos) + 0.1*x[0].v1.strand)

        for i, msr in zip(ilist, msrlist):
            elist = []
            for e in eilist:
                if hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos).intersects(i):
                    elist.append(e)
            ms_vlist = []
            msv_index = {}
            for msi in range(len((msr))):
                if msr[msi][2] > msr[msi][1]:
                    msv = breakpoint_vertex(i.chrom, msr[msi][0], 1)
                else:
                    msv = breakpoint_vertex(i.chrom, msr[msi][0], -1)
                msv_diff[msv] = msr[msi][2] - msr[msi][1]
                ms_vlist.append(msv)
                msv_index[msv] = msi
            all_msv.append(ms_vlist)
            print "Meanshift", str(i), len(ms_vlist), ms_vlist
            sys.stdout.flush()
            msve_match = {}
            for msv in ms_vlist:
                ms = msr[msv_index[msv]]
                if ms[3]:
                    msve = [e for e in elist if e[0].v1.strand * (ms[1]-ms[2]) > 0 and abs(e[0].v1.pos - ms[0]) < self.max_insert + ms_window_size1]
                    if len(msve) == 0:
                        print "finesearch discordant edges", i.chrom, ms
                        efine = self.interval_discordant_edges(hg.interval(i.chrom, msv.pos - ms_window_size0-self.max_insert, msv.pos + ms_window_size1+self.max_insert), pair_support=2)
                        if len([e for e in efine if e[0].v1.strand * (ms[1]-ms[2]) > 0]) > 0:
                            if len([(e[1], e[0]) for e in efine if e[0].v1.strand * (ms[1]-ms[2]) > 0 and abs(e[0].v1.pos - msv.pos) < ms_window_size1]) > 0:
                                ebest = max([(e[1], e[0]) for e in efine if e[0].v1.strand * (ms[1]-ms[2]) > 0 and abs(e[0].v1.pos - msv.pos) < ms_window_size1])
                            else:
                                ebest = max([(e[1], e[0]) for e in efine if e[0].v1.strand * (ms[1]-ms[2]) > 0])
                            ebest = (ebest[1], ebest[0])
                            msve = [ebest]
                            print "finesearch discordant edge found", i.chrom, ms, str(ebest[0]), len(ebest[1])
                            if len(ebest[1]) < self.pair_support:
                                elist.append(ebest)
                                eilist.append(ebest)
                                if hg.interval(ebest[0].v2.chrom, e[0].v2.pos, e[0].v2.pos).intersects(i):
                                    elist.append((breakpoint_edge(ebest[0].v2, ebest[0].v1), [(hg.interval(self.get_mates(a[1])[0], bamfile=self.bamfile), self.get_mates(a[1])) for a in ebest[1]]))
                                if len(hg.interval_list([hg.interval(ebest[0].v2.chrom, e[0].v2.pos, e[0].v2.pos)]).intersection(ilist)) > 0:
                                    eilist.append((breakpoint_edge(ebest[0].v2, ebest[0].v1), [(hg.interval(self.get_mates(a[1])[0], bamfile=self.bamfile), self.get_mates(a[1])) for a in ebest[1]]))
                                elist.sort(key=lambda x: hg.absPos(x[0].v1.chrom, x[0].v1.pos) + 0.1*x[0].v1.strand)
                                eilist.sort(key=lambda x: hg.absPos(x[0].v1.chrom, x[0].v1.pos) + 0.1*x[0].v1.strand)
                else:
                    msve = [e for e in elist if e[0].v1.strand * (ms[1]-ms[2]) > 0 and abs(e[0].v1.pos - ms[0]) < self.max_insert + ms_window_size0]
                if len(msve) > 0:
                    msve_match[msv] = msve
            msv_nocover = [msv for msv in ms_vlist if msv not in msve_match]
            all_msv_nocover.append(msv_nocover)


        #setup graph for flow optimization
        ngvlist_full = []
        elist_full = []
        ms_addlist = []
        kce = defaultdict(lambda: 0) #number of concordant reads
        koe = defaultdict(lambda: 0.0) #number of reads mapping outside the interval
        kbpe = defaultdict(lambda: 0.0) #number of discordant reads across breakpoint edge
        new_graph = breakpoint_graph()
        s = new_graph.new_vertex(ilist[0].chrom, -1, -1)
        for i, msr, ms_vlist, msv_nocover in zip(ilist, msrlist, all_msv, all_msv_nocover):
            ngvlist = []
            elist = []
            for e in eilist:
                if hg.interval(e[0].v1.chrom, e[0].v1.pos, e[0].v1.pos).intersects(i):
                    elist.append(e)
            #add vertices to new_graph
            ei = 0
            nei = 0
            msi = 0
            if (len(elist) == 0 or elist[ei][0].v1.strand == 1 or elist[ei][0].v1.pos > i.start):
                if len(msv_nocover) == 0 or msv_nocover[msi].strand ==1 or msv_nocover[msi].pos > i.start:
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
                ne = new_graph.new_edge(s, nv)
                koe[ne] = len(self.interval_crossing_arcs(nv.chrom, nv.pos, nv.pos + self.max_insert, -1, ilist))
                ei += 1
            ngvlist.append(nv)
            vc = breakpoint_vertex(ngvlist[0].chrom, ngvlist[0].pos, ngvlist[0].strand)
            while ei < len(elist) or msi < len(msv_nocover):
                vp = vc
                if msi >= len(msv_nocover):
                    vc = elist[ei][0].v1
                    ei += 1
                elif ei >= len(elist):
                    vc = msv_nocover[msi]
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
                    ms_addlist.append(msv_nocover[msi])
                    msi += 1
                if (vc.pos == vp.pos and vc.strand <= vp.strand) or (vc.pos == vp.pos + 1 and vc.strand < vp.strand):
                    continue
                if vc.strand == 1:
                    if ngvlist[nei].strand == 1:
                        nvc_prime = new_graph.new_vertex(ngvlist[nei].chrom, ngvlist[nei].pos+1, -1)
                        oecount = len(self.interval_crossing_arcs(nvc_prime.chrom, nvc_prime.pos, nvc_prime.pos + self.max_insert, -1, ilist))
                        if oecount >= self.pair_support:
                            ne = new_graph.new_edge(s, nvc_prime)
                            koe[ne] = oecount
                        ce = self.concordant_edge(vp)
                        nce = new_graph.new_edge(ngvlist[nei], nvc_prime)
                        kce[nce] = len(ce[1])
                        ngvlist.append(nvc_prime)
                        nei += 1
                    nv = new_graph.new_vertex(vc.chrom, vc.pos, 1)
                    ne = new_graph.new_edge(s, nv)
                    koe[ne] = len(self.interval_crossing_arcs(nv.chrom, max(0, nv.pos - self.max_insert), nv.pos, 1, ilist))
                    ngvlist.append(nv)
                    nei += 1
                else:
                    if ngvlist[nei].strand == 1:# and ngvlist[nei].pos > self.max_insert - self.window_size:
                        nvc_prime = new_graph.new_vertex(ngvlist[nei].chrom, ngvlist[nei].pos+1, -1)
                        oecount = len(self.interval_crossing_arcs(nvc_prime.chrom, nvc_prime.pos, nvc_prime.pos + self.max_insert, -1, ilist))
                        if oecount >= self.pair_support:
                            ne = new_graph.new_edge(s, nvc_prime)
                            koe[ne] = oecount
                        ce = self.concordant_edge(vp)
                        nce = new_graph.new_edge(ngvlist[nei], nvc_prime)
                        kce[nce] = len(ce[1])
                        ngvlist.append(nvc_prime)
                        nei += 1
                    if ngvlist[nei].strand == -1:
                        nvc_prime = new_graph.new_vertex(vc.chrom, vc.pos-1, 1)
                        oecount = len(self.interval_crossing_arcs(nvc_prime.chrom, max(0, nvc_prime.pos - self.max_insert), nvc_prime.pos, 1, ilist))
                        if oecount >= self.pair_support:
                            ne = new_graph.new_edge(s, nvc_prime)
                            koe[ne] = oecount
                        ngvlist.append(nvc_prime)
                        nei += 1
                    nv = new_graph.new_vertex(vc.chrom, vc.pos, -1)
                    ne = new_graph.new_edge(s, nv)
                    koe[ne] = len(self.interval_crossing_arcs(nv.chrom, nv.pos, nv.pos + self.max_insert, -1, ilist))
                    ce = self.concordant_edge(vc)
                    nce = new_graph.new_edge(ngvlist[nei], nv)
                    kce[nce] = len(ce[1])
                    ngvlist.append(nv)
                    nei += 1
    #            ei += 1
            if ngvlist[nei].strand == -1:
                nv = new_graph.new_vertex(i.chrom, i.end, 1)
                ne = new_graph.new_edge(s, nv)
                koe[ne] = len(self.interval_crossing_arcs(nv.chrom, nv.pos - self.max_insert, nv.pos, 1, ilist))
                ngvlist.append(nv)
                nei += 1
            elif ngvlist[nei].strand == 1 and ngvlist[nei].pos < i.end:
                nvc_prime = new_graph.new_vertex(ngvlist[nei].chrom, ngvlist[nei].pos+1, -1)
                oecount = len(self.interval_crossing_arcs(nvc_prime.chrom, nvc_prime.pos, min(hg.chrLen[hg.chrNum(nvc_prime.chrom)], nvc_prime.pos + self.max_insert), -1, ilist))
                if oecount >= self.pair_support:
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
                if ne is None:
                    raise ValueError("ne is None:" + str(e) + str(len(e0[1])) + '\n'+','.join(map(str, new_graph.vs.values())))
                    exit()
                kbpe[ne] = len(e0[1])
        for nei in range(1,len(ngvlist_full)):
            if ngvlist_full[nei].strand == 1:
                new_graph.new_edge(ngvlist_full[nei-1], ngvlist_full[nei], edge_type='sequence')
#            else:
#                new_graph.new_edge(ngvlist[nei-1], ngvlist[nei])


        #set up all constants
        C = self.median_coverage()[0] / 2
        print "C (haploid coverage) = ", C
        G = new_graph
        seqlist = [e for e in new_graph.es.values() if e.edge_type == 'sequence']
        n = len(seqlist)
        l = [abs(e.v2.pos - e.v1.pos)+1 for e in seqlist]
        k = [len([a for a in self.bamfile.fetch(e.v1.chrom, e.v1.pos, e.v2.pos)]) for e in seqlist]
#        kgcc = [self.interval_coverage(hg.interval(i.chrom, e.v1.pos, e.v2.pos), gcc=True) * (e.v2.pos - e.v1.pos) / self.read_length for e in seqlist]
#        k = kgcc
        kcc = [self.interval_coverage(hg.interval(e.v1.chrom, e.v1.pos, e.v2.pos)) * (e.v2.pos - e.v1.pos) for e in seqlist]
        ke = {}
        ke.update(kbpe)
        ke.update(kce)
        ke.update(koe)
        K = [len([a for a in self.bamfile.fetch(i.chrom, e.v1.pos, e.v2.pos)]) * self.read_length/(abs(e.v2.pos - e.v1.pos) + 1.0) for e in seqlist]
        #edge read count kbpe defined above
        bplist = [e for e in new_graph.es.values() if e.edge_type == 'breakpoint']
        m = len(bplist)
        bpdict = {bplist[bpi]: bpi for bpi in range(len(bplist))}
        print "########## len bplist", len(bplist), ";   ################ kbpe, kce, koe = ", len(kbpe), len(kce), len(koe) 

        #set up problem size and variable types and constraint types
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

        #setting up constraints
        for i in range(n):
            subarr = [i]
            valarr = [1.0]
            for e in seqlist[i].v1.elist:
                if e.edge_type == 'sequence':
                    continue
                subarr.append(n + bpdict[e])
                valarr.append(-1.0)
            asub.append(np.array(subarr))
            aval.append(np.array(valarr))
            subarr = [i]
            valarr = [1.0]
            for e in seqlist[i].v2.elist:
                if e.edge_type == 'sequence':
                    continue
                subarr.append(n + bpdict[e])
                valarr.append(-1.0)
            asub.append(np.array(subarr))
            aval.append(np.array(valarr))

        #setting up objective
        opro = [mosek.scopr.log] * (n + m)
        oprjo = range(n + m)
        oprfo = [-1 * ki for ki in k] + [-1 * ke[e] for e in bplist]
        oprgo = [C * li / self.read_length for li in l] + [(self.max_insert) * C / 2 / self.read_length for e in bplist]
        oprho = [0.0] * (n + m)
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
        for j in range(numvar):
            task.putbound(mosek.accmode.var, j, bkx[j], blx[j], bux[j])
        for i in range(numcon):
            task.putbound(mosek.accmode.con, i, bkc[i], blc[i], buc[i])
        for i in range(numcon):
            task.putarow(i, asub[i], aval[i])
#        for i in qsubi:
#           task.putqconk(2 * n + 2, qsubi[i], qsubj[i], qval[i])
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
            mslc = [sum([ap[0].intersection(ap[1]).size() * float(ap[1].info[0]) for ap in hg.interval_list([msi]).intersection(slist)]) / sum([ap[0].intersection(ap[1]).size() for ap in hg.interval_list([msi]).intersection(slist)]) for msi in mslist]
            for ap in zip(mslist, mslc):
                print '\t'.join(map(str, ["MScount", str(ap[0]), ap[1]]))
        
        edge_code = defaultdict(lambda:'breakpoint', {'concordant':'concordant', 'source':'source'})

        graph_logger.info("Sequence edge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size")
        for si in range(n):
            graph_logger.info('sequence\t' + '\t'.join(map(str, [seqlist[si].v1,  seqlist[si].v2, res[si], K[si], seqlist[si].v2.pos - seqlist[si].v1.pos])))
            wehc[seqlist[si]] = float(res[si])
        graph_logger.info("Breakpoint edge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs")
        for bpi in range(m):
#            print  edge_code[bplist[bpi].type()], str(bplist[bpi]), res[n + bpi], ke[bplist[bpi]], bplist[bpi].kmer_homology()
            graph_logger.info('\t'.join(map(str, [edge_code[bplist[bpi].type()], bplist[bpi], res[n + bpi], ke[bplist[bpi]]])))
            wehc[bplist[bpi]] = float(res[n + bpi])
        lenlist = len(ilist)
        if len(ilist0) >= 10:
            lenlist = len(ilist0)
        all_msv_cat = reduce(lambda x, y: x+y, all_msv, [])
        oncolist = ','.join(Set([a[1].info['Name'] for a in ilist.intersection(hg.oncogene_list)]))+','
        istr = ','.join([i.chrom + ':' + str(i.start) + '-' + str(i.end) for i in ilist])
        summary_logger.info("#Intervals = " +  str(lenlist))
        summary_logger.info("TotalIntervalSize = " + str(sum([a.size() for a in ilist])))
        summary_logger.info("AmplifiedIntervalSize = " + str(sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5])))
        if len([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]) > 0:
            summary_logger.info("AverageAmplifiedCopyCount = " + str(sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos) for si in range(n) if res[si] >= 2.5]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5])))
        else:
            summary_logger.info("AverageAmplifiedCopyCount = 2")
        summary_logger.info("#Intervals = " + str(len(Set([i.chrom for i in ilist]))))
        summary_logger.info("#SeqenceEdges = " + str(n))
        summary_logger.info("#BreakpointEdges = " + str(len(kbpe)))
        summary_logger.info("#CoverageShifts = " + str(len(all_msv_cat)))
        summary_logger.info("#MeanshiftSegmentsCopyCount>5 = " + str(len([v for v in msv_diff.values() if v > 5])))
        summary_logger.info("#Foldbacks = " + str(len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1])))
        summary_logger.info("#CoverageShiftsWithBreakpointEdges = " + str(len([msa for msa in all_msv_cat if msa in ms_addlist])))
        summary_logger.info("OncogenesAmplified = " + str(oncolist))
        summary_logger.info("Intervals = " + str(istr))


        #Summary, #intervals, t talsize, size>2.5, AvgCoverage>2.5, #chromosomes, #sequenceedges, #breakpointedges, #meanshiftbreaks, #meanshift>5, #msfoldbackedges, #msfoldbackedges, #mswithoutbreakpoint, oncogenes, representativestring, #bpedgeswithcommonkmers
        if len([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]) > 0:
#            print '\t'.join(map(str, ["Summary:", lenlist, sum([a.size() for a in ilist]), sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]), sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos) for si in range(n) if res[si] >= 2.5]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]), len(Set([i.chrom for i in ilist])), n, len(kbpe), len(all_msv_cat), len([v for v in msv_diff.values() if v > 5]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if msa in ms_addlist]), oncolist, istr, len([e for e in kbpe if e.kmer_homology()])]))
            print '\t'.join(map(str, ["Summary:", lenlist, sum([a.size() for a in ilist]), sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]), sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos) for si in range(n) if res[si] >= 2.5]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]), len(Set([i.chrom for i in ilist])), n, len(kbpe), len(all_msv_cat), len([v for v in msv_diff.values() if v > 5]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if msa in ms_addlist]), oncolist, istr]))
        else:
#            print '\t'.join(map(str, ["Summary:", lenlist, sum([a.size() for a in ilist]), sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]), sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos) for si in range(n)]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n)]), len(Set([i.chrom for i in ilist])), n, len(kbpe), len(all_msv_cat), len([v for v in msv_diff.values() if v > 5]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if msa in ms_addlist]), oncolist, istr, len([e for e in kbpe if e.kmer_homology()])]))
            print '\t'.join(map(str, ["Summary:", lenlist, sum([a.size() for a in ilist]), sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n) if res[si] >= 2.5]), sum([res[si] * (seqlist[si].v2.pos - seqlist[si].v1.pos) for si in range(n)]) / sum([seqlist[si].v2.pos - seqlist[si].v1.pos for si in range(n)]), len(Set([i.chrom for i in ilist])), n, len(kbpe), len(all_msv_cat), len([v for v in msv_diff.values() if v > 5]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if self.foldup_count(msa.chrom, msa.pos, msa.strand) >= 1]), len([msa for msa in all_msv_cat if msa in ms_addlist]), oncolist, istr]))
        new_graph.cycle_decomposition(wehc, s)


    def interval_neighbors(self, i, ilist=[], rdlist=[], t=0):
        i2 = self.interval_extend(i)
##        i2 = self.interval_extend(i, ilist, rdlist)
        ms_window_size = 10000
        edges = self.interval_discordant_edges(i2)
        edges.sort(cmp=lambda x, y: len(x[1]) > len(y[1]))
        ei = 0
        neighbors = hg.interval_list([])
        print 'interval_neighbors found_discordant edges', str(i), str(i2)
        while len(neighbors) < 10 and ei < len(edges):
            covered = False
            for i3 in ilist + neighbors:
                if i3.chrom == edges[ei][0].v2.chrom and edges[ei][0].v2.pos >= i3.start and edges[ei][0].v2.pos <= i3.end:
                    ei += 1
                    covered = True
                    break
            if covered:
                continue
            if edges[ei][0].v2.strand < 0:
                n = self.interval_extend(hg.interval(edges[ei][0].v2.chrom, edges[ei][0].v2.pos, min(hg.chrLen[hg.chrNum(edges[ei][0].v2.chrom)] - 1, edges[ei][0].v2.pos + self.max_insert)))
            else:
                n = self.interval_extend(hg.interval(edges[ei][0].v2.chrom, max(1, edges[ei][0].v2.pos - self.max_insert), edges[ei][0].v2.pos))
            if n.size() > self.max_insert + 1:
                n.info = len(edges[ei][1])
                neighbors.append(n)
            ei += 1
        neighbors.sort()
        mc = neighbors.merge_clusters(extend=ms_window_size)
        for c in mc:
            c[0].info = sum([c1.info for c1 in c[1]])
        nn = hg.interval_list([c[0] for c in mc])
        for n in nn:
            print str(n)
        print "-----------------------------------------------------------------"
#        if t == 0:
#            for n in neighbors:
#                self.interval_neighbors(n, t=1)
        return nn

    def interval_hops(self, i, ilist=[], rdlist=[]):
        ms_window_size = 10000
        ii = self.interval_extend(i)
        ilist = hg.interval_list([ii])
        seen_list = hg.interval_list([])
        unseen_list = [(0,ii)]
        heapq.heapify(unseen_list)
        while len(ilist) < 10 and len(unseen_list) > 0:
            ic = heapq.heappop(unseen_list)[1]
            print "Interval hops", str(ic)
            icn = self.interval_neighbors(ic, ilist)
            for ic2 in icn:
                contained = False
                for i2 in ilist:
                    if i2.contains(ic2):
                        contained = True
                if contained:
                    continue
                if ic2.size() < 2 * ms_window_size and len(self.interval_discordant_edges(i)) < 2:
                    continue
                heapq.heappush(unseen_list, (-ic2.info, ic2))
                ilist.append(ic2)
            seen_list.append(ic)
            print "========================================================"
        ilist.sort()
        for ic in ilist:
            print "IntervalHops", str(ic)
        return ilist

    def interval_amplified(self, i, filter_conserved=True, filter_small=False):
        if len(hg.interval_list([i]).intersection(hg.conserved_regions)) > 0:
            return False
        ms_window_size = 10000
        num_w = 0
        num_high = 0
        if i.size() < 2 * ms_window_size and len(self.interval_discordant_edges(i)) < 2:
            return False 
        wc = self.window_coverage(i, ms_window_size)
        mc = self.median_coverage()
        for w in wc:
            num_w += 1
            if w[1] > mc[0] + 3 * mc[2]:
                num_high += 1
#        wc_high = len([w for w in wc if w[1] > mc[1] + 3 * mc[2]])
        if num_high > num_w / 5:
            return True
        else:
            return False

    def interval_extend(self, i, strand=0, i0=None):
        print "Interval extend start:", str(i), strand
        ms_window_size = 10000
        extend_size = max(i.size() / ms_window_size, 1)
        max_window_size = 3000000
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
#        print str(ic), ic.size(), extend_right, extend_left, right_size, left_size
        while ic.size() < max_window_size and (extend_left >= 0 or extend_right >= 0):
#            print str(ic), extend_right, extend_left, right_size, left_size
            if extend_right >= 0:
                if right_size < 1:
                    extend_right = -1
                elif ic.end + right_size * ms_window_size > hg.chrLen[hg.chrNum(ic.chrom)]:
                    if self.interval_amplified(hg.interval(ic.chrom, ic.end, hg.chrLen[hg.chrNum(ic.chrom)])):
                        ic.end = hg.chrLen[hg.chrNum(ic.chrom)]
                        extend_right = -1
                    else:
                        extend_right = 0
                        right_size = right_size / 2
                elif self.interval_amplified(hg.interval(ic.chrom, ic.end, ic.end + right_size * ms_window_size)):
                    ic.end = ic.end + right_size * ms_window_size
                    if extend_right == 1:
                        right_size = 2 * right_size
                    else:
                        right_size = right_size / 2
                        if right_size < 1:
#                            ic.end = min(ic.end + ms_window_size, hg.chrLen[hg.chrNum(ic.chrom)])
                            extend_right = -1
                else:
                    extend_right = 0
                    right_size = right_size / 2
            if extend_left >= 0:
                if left_size < 1:
                    extend_left = -1
                elif ic.start - left_size * ms_window_size <= 1:
                    if self.interval_amplified(hg.interval(ic.chrom, 1, ic.start)):
                        ic.start = 1
                        extend_left = -1
                    else:
                        extend_left = 0
                        left_size = left_size / 2
                elif self.interval_amplified(hg.interval(ic.chrom, ic.start - left_size * ms_window_size, ic.start)):
                    ic.start = ic.start - left_size * ms_window_size
                    if extend_left == 1:
                        left_size = 2 * left_size
                    else:
                        left_size = left_size / 2
                        if left_size < 1:
#                            ic.start = max(ic.end - ms_window_size, 1)
                            extent_left = -1
                else:
                    extend_left = 0
                    left_size = left_size / 2
        if self.interval_amplified(hg.interval(ic.chrom, ic.end - 2*ms_window_size, min(ic.end + 2 * ms_window_size, hg.chrLen[hg.chrNum(ic.chrom)]))):
            ic.end = min(ic.end + 10 * ms_window_size, hg.chrLen[hg.chrNum(ic.chrom)])
        if self.interval_amplified(hg.interval(ic.chrom, max(ic.start - 2 * ms_window_size, 0), ic.start + 2*ms_window_size)):
            ic.start = max(ic.start - 10 * ms_window_size, 0)
#        print str(ic)
        if strand >= 0:
            ide = self.interval_discordant_edges(hg.interval(ic.chrom, ic.end, min(hg.chrLen[hg.chrNum(ic.chrom)], ic.end + ms_window_size)))
            for e in ide:
                if e[0].v1.strand == -1:
                    ic.end = min(ic.end + 2 * ms_window_size, hg.chrLen[hg.chrNum(ic.chrom)])
                    break
        if strand <= 0:
            ide = self.interval_discordant_edges(hg.interval(ic.chrom, max(1, ic.start - ms_window_size), ic.start))
            for e in ide:
                if e[0].v1.strand == 1:
                    ic.start = max(ic.start - 2 * ms_window_size, 1)
                    break
#        if ic.size() > ms_window_size:
        print "Interval extend:", str(i), strand, str(ic)
        return ic

        
