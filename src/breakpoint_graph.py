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


from sets import Set
import heapq
import logging

from abstract_graph import *
import hg19util as hg

cycle_logger = logging.getLogger('cycle')

class breakpoint_vertex(abstract_vertex):
    def __init__(self, chrom='', pos=-2, strand=1, vid=0, graph=None):
        if pos == -2:
            vstring = chrom
            chrom = vstring[:vstring.find(':')]
            pos = int(vstring[vstring.find(':') + 1:-1])
            strand = 1 if vstring[-1] == '+' else -1
        if graph is not None and graph.has_vertex(chrom, pos, strand):
            raise Exception('Duplicate vertex added')
        abstract_vertex.__init__(self, vid, graph)
        self.chrom = chrom
        self.pos = pos
        self.strand = strand

    def __repr__(self):
        if self.strand == 1:
            return self.chrom + ':' + str(self.pos) + '+'
        else:
            return self.chrom + ':' + str(self.pos) + '-'

    def __hash__(self):
        return str(self).__hash__()

    def __gt__(self, y):
        return hg.absPos(self.chrom, self.pos) + 0.4 * self.strand > hg.absPos(y.chrom, y.pos) + 0.4 * y.strand



class breakpoint_edge(abstract_edge):
    def __init__(self, v1, v2=None, eid=0, graph=None, update_vertices=True, edge_type="discordant"):
        if type(v1) == str:
            estr = v1
            v1 = breakpoint_vertex(estr.split('>')[0][:-1])
            v2 = breakpoint_vertex(estr.split('>')[1])
        abstract_edge.__init__(self, v1, v2, eid, graph, update_vertices)
        if edge_type in ["concordant", "sequence"]:
            if v1.chrom != v2.chrom:
                raise Exception("Edge of type " + edge_type + " connects different chromosomes.")
        if edge_type in ["concordant", "sequence"]:
            if v1.strand == v2.strand:
                raise Exception("Edge of type " + edge_type + " connects same strand.")
        if edge_type == "concordant":
            if ((v1.strand == 1 and v1.pos + 1 != v2.pos) or
                (v2.strand == 1 and v2.pos + 1 != v1.pos)):
                raise Exception("Edge of type " + edge_type + " connects non-adjacent positions.")
        if edge_type == "sequence":
            if v1.strand == -1 and v1.pos > v2.pos:
                raise Exception("Start position for sequence edge greater than end position:" + str(v1) + '->' + str(v2))
            if v1.strand == 1 and v2.pos > v1.pos:
                raise Exception("Start position for sequence edge greater than end position")
        self.edge_type = edge_type

    def kmer_homology(self, k=10, span=100):
        seq1 = ''.join([a.capitalize() for a in hg.interval(self.v1.chrom, max(1,self.v1.pos - span), min(self.v1.pos + span, hg.chrLen[hg.chrNum(self.v1.chrom)]), self.v1.strand).sequence()])
        seq2 = ''.join([a.capitalize() for a in hg.interval(self.v2.chrom, max(1,self.v2.pos - span), min(self.v2.pos + span, hg.chrLen[hg.chrNum(self.v2.chrom)]), -1 * self.v2.strand).sequence()])
        kset1 = Set([seq1[i:i+10] for i in range(len(seq1) - k + 1)])
        kset2 = Set([seq2[i:i+10] for i in range(len(seq2) - k + 1)])
        return len(kset1.intersection(kset2))

    def type(self, min_insert=0, max_insert=500):
        if self.v1.pos == -1 or self.v2.pos == -1:
            return "source"
        elif self.v1.chrom != self.v2.chrom:
            return "interchromosomal"
        elif self.v1.pos <= self.v2.pos:
            vmin = self.v1
            vmax = self.v2
        else:
            vmin = self.v2
            vmax = self.v1
        if vmax.strand == 1 and vmin.strand == -1:
            return "everted"
        if vmax.strand == 1 and vmin.strand == 1:
            return "forward"
        if vmax.strand == -1 and vmin.strand == -1:
            return "reverse"
        if vmax.pos - vmin.pos > max_insert or vmax.pos - vmin.pos < min_insert:
            return "discordant"
        return "concordant"
               
    def __repr__(self):
        return str(self.v1) + '->' + str(self.v2)


class breakpoint_graph(abstract_graph):
    def __init__(self):
        abstract_graph.__init__(self)
        self.vhash = {}

    def has_vertex(self, chrom, pos, strand):
        vtemp = breakpoint_vertex(chrom, pos, strand)
        if vtemp.__hash__() in self.vhash:
            return self.vhash[vtemp.__hash__()]
        else:
            return None

    def new_vertex(self, chrom, pos, strand):
        v = self.has_vertex(chrom, pos, strand)
        if v is not None:
            return v
        v = breakpoint_vertex(chrom, pos, strand, graph=self)
        self.vhash[v.__hash__()] = v
        return v

    def new_edge(self, v1, v2, edge_type='breakpoint'):
        return breakpoint_edge(v1, v2, graph=self, edge_type=edge_type)

    def add_vertex(self, v):
        return self.new_vertex(v.chrom, v.pos, v.strand)

    def add_edge(self, e):
        if e.graph is not self:
            v1 = self.has_vertex(e.v1.chrom, e.v1.pos, e.v1.strand)
            v2 = self.has_vertex(e.v2.chrom, e.v2.pos, e.v2.strand)
            if v1 is None or v2 is None:
                return None
            return self.new_edge(v1, v2)
        return e

    def cycle_decomposition(self, w, s): #w is dict containing weights of edges and s is source vertex, this vertex has the exception of not having a sequence edge attached
        def thickest_cycle(hce, wehc):
#            print hce, wehc
            v1 = hce[1].v1
            a = [(-1 * hce[0], v1)]
            heapq.heapify(a)
            hdict = {v1: (hce[0], [hce[1]], None, Set([]))}
            seenSet = Set([])
            seenEdges = Set([])
            completed = False
            while len(a) > 0 and not completed:
#                print len(a), str(a[0]), str(hdict[a[0][1]])
                v1w, v1 = heapq.heappop(a)
#                if hce[1].v1.pos == 133027113:
#                    print "=============================================================="
#                    print 'here0', str(v1), v1w
                if v1 == hce[1].v1 and v1 in seenSet:
                    completed = True
                    break
#                if hce[1].v1.pos == 133027113:
#                    print 'here1'
                for e in v1.elist:
#                    if hce[1].v1.pos == 133027113:
#                        print 'here2', str(e)
                    if e.edge_type == 'sequence':
                        continue
                    else:
                        v2 =  e.neighbor(v1)
#                    if hce[1].v1.pos == 133027113:
#                        print 'here2', str(v2)
                    if v2 == s:
                        v3 = v2
                        if e in hdict[v1][3]:
                            nw = min(hdict[v1][0], wehc[e] / 2)
                        else:
                            nw = min(hdict[v1][0], wehc[e])
                        if not v3 in hdict or hdict[v3][2] is None or hdict[v3][0] < nw:
 #                           if hce[1].v1.pos == 133027113:
 #                               print 'here3', str(v3)
                            nhdict = hdict[v1][3].copy()
                            nhdict.add(e)
                            hdict[v3] = (nw, [e], v1, nhdict)
#                            print 'seen edges', e, v3, hdict[v3]
                            seenEdges.add(e)
                    else:
                        for e2 in v2.elist:
                            if e2.edge_type == 'sequence':
                                se = e2
                                v3 = e2.neighbor(v2)
                                break
                        if e in hdict[v1][3]:
#                            print 'e is seen', e, seenEdges
                            nw = min(hdict[v1][0], wehc[e] / 2)
                        elif se in hdict[v1][3]:
#                            print 'se is seen', se, seenEdges
                            nw = min(hdict[v1][0], wehc[e], wehc[se] / 2)
                        else:
                            nw = min(hdict[v1][0], wehc[e])
                        if not v3 in hdict or hdict[v3][2] is None or hdict[v3][0] < nw:
#                            if hce[1].v1.pos == 133027113:
#                                print 'here4', str(v3)
                            nhdict = hdict[v1][3].copy()
                            nhdict.add(e)
                            nhdict.add(se)
                            hdict[v3] = (nw, [e, se], v1, nhdict)
 #                           print 'seen edges', e, se, v3, hdict[v3]
                            seenEdges.add(e)
                            seenEdges.add(se)
                    if v3 in seenSet:
 #                       if hce[1].v1.pos == 133027113:
 #                           print 'here5', str(v3), str(hdict[v3][0]), str(hdict[v3][1]), str(hdict[v3][2])
                        continue
#                    if hce[1].v1.pos == 133027113:
#                        print 'here6', str(v3)
                    seenSet.add(v3)
                    heapq.heappush(a, (-1 * hdict[v3][0], v3))
            if len(a) == 0 and not completed:
                print "NOT COMPLETED", hce[1].v1
#            print hdict
            s2Set = Set([])
            tc = hdict[hce[1].v1][1]
            v2 = hdict[hce[1].v1][2]
            while v2 != hce[1].v1: #and not v2 in s2Set:
#                print hce[1].v1, v2, s2Set
                s2Set.add(v2)
                if v2 not in hdict:
                    print str(v2), str(hce[1].v1), str(tc)
                    for ee in hce[1].v1.elist:
                        print str(ee), wehc[ee]
                tc = hdict[v2][1] + tc
                v2 = hdict[v2][2]
                s2Set.add(v2)
 #               print v2, tc
            return tc, hdict[hce[1].v1][0]
        
        total_amplicon_content = sum([(e.v2.pos - e.v1.pos) * w[e] for e in w if e.edge_type == 'sequence'])
        print 'total_amplicon_content', total_amplicon_content
        amplicon_content_covered = 0
        w2 = w.copy()
        cycle_number = 1
        cycle_list = []
        while max(w2.values()) > 0:
            we = [(w2[e], e) for e in w2]
            we.sort()
            wer = we[::-1]
            we = wer
            wei = 0
            tcwmax  = -1
            tcmax = None
            tchwmax = -1
            tchmax = None
            tchw = -1
            while wei < len(we):# and (tcwmax == -1 or we[wei][0] >= tcwmax / 2.0):
#                if we[wei][1].edge_type == 'sequence':
#                    wei += 1
#                    continue
                tc, tcw = thickest_cycle(we[wei], w2)
                if len(tc) < 2:
                    print str(tc[0])
                    exit()
                if tcw > tcwmax:
                    tcmax = tc
                    tcwmax = tcw
#                sumlen = sum([abs(e.v1.pos - e.v2.pos) for e in tc if e.edge_type == 'sequence'])
#                if sumlen * tcw > tchwmax:
#                    tchwmax = sumlen * tcw
#                    tchmax = tc
#                    tchw = tcw
                wei += 1
            if tcwmax == -1:
                break
            tc = tcmax
            tcw = tcwmax
#            tc = tchmax
#            tcw = tchw
            for e in tc:
                if type(e) == tuple :
                    print str(e), str(tc)
            if -1 in [e.v1.pos for e in tc] + [e.v2.pos for e in tc]:
                csource = 0
                for ci in range(len(tc) - 1):
                    if -1 in [tc[ci].v1.pos, tc[ci].v2.pos] and -1 in [tc[ci +1].v1.pos, tc[ci + 1].v2.pos]:
                        csource = ci + 1
                        tc = tc[ci + 1:] + tc[0:ci + 1]
                        break
                if tc[0].v1 == tc[1].v1 or tc[0].v1 == tc[1].v2:
                    v2 = tc[0].v1
                    v1 = tc[0].v2
                else:
                    v2 = tc[0].v2
                    v1 = tc[0].v1
                for ci in range(len(tc)):
                    if tc[ci].v1.pos == v1.pos:
                        v2 = tc[ci].v2
                    else:
                        v2 = tc[ci].v1
                    if tc[ci].edge_type == 'sequence':
                        if v1.pos > v2.pos:
                            tc = tc[::-1]
                        break
                    v1 = v2
            else:
                if tc[0].v1 == tc[1].v1 or tc[0].v1 == tc[1].v2:
                    v2 = tc[0].v1
                    v1 = tc[0].v2
                else:
                    v2 = tc[0].v2
                    v1 = tc[0].v1
                for ci in range(len(tc)):
                    if tc[ci].v1.pos == v1.pos:
                        v2 = tc[ci].v2
                    else:
                        v2 = tc[ci].v1
                    if tc[ci].edge_type == 'sequence':
                        if v1.pos > v2.pos:
                            tc = tc[ci::-1] + tc[:ci:-1]
                        break
                    v1 = v2
                ci = 0
                while tc[ci].type() == 'concordant' or tc[ci-1].type() == 'concordant':
                    ci -= 1
                tc = tc[ci:] + tc[: ci]
                    
            if tcw == 0:
                print "tcw is 0"
                break
            print "Cycle ", cycle_number, ": Copy count = ",tcw, tc
            cycle_edge_list = []
            ci = 1
            v0 = None
            v0c = None
            if tc[0].v1 == tc[1].v1 or tc[0].v1 == tc[1].v2:
                v2 = tc[0].v1
                v1 = tc[0].v2
            else:
                v2 = tc[0].v2
                v1 = tc[0].v1
            if tc[0].edge_type == 'sequence':
                v0 = v1
                v0c = v2
            elif v1.pos == -1 or v2.pos == -1:
                print v1, "->", v2
                cycle_edge_list.append((v1,v2))
            v1 = v2
            while ci < len(tc):
                if tc[ci].v1.pos == v1.pos:
                    v2 = tc[ci].v2
                else:
                    v2 = tc[ci].v1
                if v1.pos == -1 or v2.pos == -1:
                    if v0 is not None:
                        print v0, "->", v0c
                        cycle_edge_list.append((v0,v0c))
                    print v1, "->", v2
                    cycle_edge_list.append((v1,v2))
                    v0 = None
                    v0c = None
                elif tc[ci].edge_type == 'sequence':
                    if v0 is None:
                        v0 = v1
                        v0c = v2
                    else:
                        v0c = v2
                elif tc[ci].type() != 'concordant':
                    if v0 is not None:
                        print v0, "->", v0c
                        cycle_edge_list.append((v0,v0c))
                        v0 = None
                        v0c = None
                v1 = v2
                ci += 1
            if v0 is not None:
                print v0, "->", v0c
                cycle_edge_list.append((v0,v0c))
            if amplicon_content_covered <= 0.9 * total_amplicon_content:
                cycle_list.append([cycle_number, tcw, tc, cycle_edge_list])
                acc = tcw * sum([abs(e[1].pos - e[0].pos) for e in cycle_edge_list if -1 not in [e[0].pos, e[1].pos]])
                amplicon_content_covered += acc
                print acc, amplicon_content_covered, total_amplicon_content, tcw, [abs(e[1].pos - e[0].pos) for e in cycle_edge_list]
            cycle_number += 1    
            #print tcw, tc
            for e in tc:
                w2[e] = w2[e] - tcw
#                if w2[e] == 0.0:
#                    w2.pop(e)

        segment_list = []
        for c in cycle_list:
            max_segment = c[3][0]
            max_orientation = '+'
            max_segi = 0
            segi = 0
            for e in c[3]:
                if (-1 in (max_segment[0].pos, max_segment[1].pos) and -1 not in (e[0].pos, e[1].pos)) or (abs(e[0].pos-e[1].pos) >= abs(max_segment[0].pos - max_segment[1].pos)):
                    max_segment = e
                    max_segi = segi
                    print str(e[0]), str(e[1]), e[0].pos, e[0].strand, e[1].pos, e[1].strand
                    if e[0].pos + 0.4*e[0].strand <= e[1].pos + 0.4*e[1].strand:
                        max_orientation = '+'
                    else:
                        max_orientation = '-'
                if e[0].pos + 0.4*e[0].strand <= e[1].pos + 0.4*e[1].strand:
                    if e not in segment_list:
                        segment_list.append(e)
                else:
                    if (e[1], e[0]) not in segment_list:
                        segment_list.append((e[1], e[0]))
                segi += 1
            if max_orientation == '+':
                c[3] = c[3][max_segi: ] + c[3][:max_segi]
            else:
                c[3] = [(e[1], e[0]) for e in c[3][:max_segi + 1][::-1]+c[3][max_segi + 1:][::-1]]

        segment_list.sort()
        segi = 1
        segment_index = {}
        for s in  [ss for ss in segment_list if ss[0].pos != -1 and ss[1].pos != -1]:
            segment_index[s] = segi
            segi += 1
        cycle_logger.info('List of cycle segments')
        for s in [ss for ss in segment_list if ss[0].pos == -1 or ss[1].pos == -1]:
            segment_index[s] = 0
        for s in [ss for ss in segment_list if ss[0].pos != -1 and ss[1].pos != -1]:
            cycle_logger.info('Segment\t' + '\t'.join([str(segment_index[s]), s[0].chrom, str(s[0].pos), str(s[1].pos)]))
        for c in cycle_list:
            seglist = []
            orientation_list = []
            for e in c[3]:
                if e in segment_index:
                    seglist.append(segment_index[e])
                    orientation_list.append('+')
                else:
                    seglist.append(segment_index[(e[1], e[0])])
                    orientation_list.append('-')
            cycle_logger.info("Cycle=" + str(c[0]) + ";Copy_count=" + str(c[1]) + ";Segments=" + ','.join([str(e[0])+str(e[1]) for e in zip(seglist, orientation_list)]))

        return None


    def __repr__(self):
        return '/n'.join(map(str, self.vs.values() + self.es.values())) + '\n'
