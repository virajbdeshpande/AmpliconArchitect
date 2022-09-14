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

import sys

from collections import defaultdict
import heapq
import logging

from abstract_graph import *
import ref_util as hg

cycle_logger = logging.getLogger('cycle')

class breakpoint_vertex(abstract_vertex):
    """Class representing breakpoint vertex derived from abstract_graph.abstract_vertex

    Attributes:
    chrom = chromosome name
    pos = 1-based chromosomal location
    strand = 1/-1 for forward/reverse strand
    vid = (optional)id of vertex
    graph = (optional) graph to which vertex belongs"""
    def __init__(self, chrom='', pos=-2, strand=1, vid=0, graph=None):
        """2 ways to initialize:
            1) chrom: breakpoint_vertex string in the format chrom:pos("+"/"-"")
            2) chrom, pos, strand: name(STR), pos (INT), strand("+"/"-"")"""
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
        """String format chrom:pos(+/-)"""
        if self.strand == 1:
            return self.chrom + ':' + str(self.pos) + '+'
        else:
            return self.chrom + ':' + str(self.pos) + '-'

    def __hash__(self):
        return str(self).__hash__()

    def __gt__(self, y):
        """Order vertices by absolute position (See hg19util.absPos) + strand"""
        return hg.absPos(self.chrom, self.pos) + 0.4 * self.strand > hg.absPos(y.chrom, y.pos) + 0.4 * y.strand


class breakpoint_edge(abstract_edge):
    """Class representing breakpoint edge derived from abstract_graph.abstract_edge

    Attributes:
    v1 = breakpoint_vertex 1 of the edge (recommend using v2 > v1)
    v2 = breakpoint_vertex 2 of the edge
    edge_type = "discordant"/"breakpoint" or "concordant" : genomic connectivity or source; "sequence": genomic interval
    eid = (optional) edge id
    graph = (optional) graph to which edge belongs"""
    def __init__(self, v1, v2=None, eid=0, graph=None, update_vertices=True, edge_type="discordant", hom=None, hom_seq=None):
        """2 ways to initialize:
        1) v1 = breakpoint_edge string in the format breakpoint_vertex1->breakpoint_vertex2
        2) v1,v2 = breakpoint_point_vertices
        Optional:
        eid: edge id
        graph: breakpoint_graph
        update_vertices: True if vertices edge should be added to vertex neighbor list
        edge_type: " = "discordant"/"breakpoint" or "concordant" : genomic connectivity or source; "sequence": genomic interval
        Required:
        If edge_type = "sequence": v1.chrom = v2.chrom, v1.pos > v2.pos else if equal v1.strand > v2.strand
        If edge_type = "concordant": v1.chrom = v2.chrom, |v1.pos - v2.pos| = 1 and the smaller has strand = 1 else -1"""
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
        self.hom = hom
        self.hom_seq = hom_seq

    def sequence(self, flank_size=-1):
        if self.edge_type == 'sequence':
            seq = hg.interval(self.v1.chrom, self.v1.pos, self.v2.pos).sequence()
            if flank_size > 0:
                seq = hg.interval(self.v1.chrom, self.v1.pos - flank_size + 1, self.v1.pos).sequence() + seq + hg.interval(self.v2.chrom, self.v2.pos, self.v2.pos + flank_size - 1).sequence()
        else:
            if self.hom == None:
                seq = 'N' * 20
            else:
                seq = self.hom_seq
            if flank_size == -1:
                flank_size = 1000
            if flank_size > 0:
                if self.hom is not None and self.hom > 0:
                    hom = self.hom
                else:
                    hom = 0
                if self.edge_type == 'source':
                    if self.v2.strand == -1:
                        right_seq = hg.interval(self.v2.chrom, self.v2.pos + hom, self.v2.pos + hom + flank_size - 1).sequence()
                        left_seq = ''
                    else:
                        left_seq = hg.interval(self.v2.chrom, self.v2.pos - hom - flank_size + 1, self.v2.pos - hom).sequence()
                        right_seq = ''
                elif self.v1.strand == 1:
                    left_seq = hg.interval(self.v1.chrom, self.v1.pos - hom - flank_size + 1, self.v1.pos - hom).sequence()
                    if self.v2.strand == -1:
                        right_seq = hg.interval(self.v2.chrom, self.v2.pos + hom, self.v2.pos + hom + flank_size - 1).sequence()
                    else:
                        right_seq = hg.interval(self.v2.chrom, self.v2.pos - hom - flank_size + 1, self.v2.pos - hom, strand=-1).sequence()
                else:
                    right_seq = hg.interval(self.v1.chrom, self.v1.pos + hom, self.v1.pos + hom + flank_size - 1).sequence()
                    if self.v2.strand == -1:
                        left_seq = hg.interval(self.v2.chrom, self.v2.pos + hom, self.v2.pos + hom + flank_size - 1, strand=-1).sequence()
                    else:
                        left_seq = hg.interval(self.v2.chrom, self.v2.pos - hom - flank_size + 1, self.v2.pos - hom).sequence()
            seq = left_seq + seq + right_seq
        return seq

    def kmer_homology(self, k=10, span=100):
        """Number of shared k-mers within "span" distance on either side of vertex positions"""
        seq1 = ''.join([a.capitalize() for a in hg.interval(self.v1.chrom, max(1,self.v1.pos - span), min(self.v1.pos + span, hg.chrLen[hg.chrNum(self.v1.chrom)]), self.v1.strand).sequence()])
        seq2 = ''.join([a.capitalize() for a in hg.interval(self.v2.chrom, max(1,self.v2.pos - span), min(self.v2.pos + span, hg.chrLen[hg.chrNum(self.v2.chrom)]), -1 * self.v2.strand).sequence()])
        kset1 = set([seq1[i:i+10] for i in range(len(seq1) - k + 1)])
        kset2 = set([seq2[i:i+10] for i in range(len(seq2) - k + 1)])
        return len(kset1.intersection(kset2))

    def type(self, min_insert=0, max_insert=500):
        """Determine type of "breakpoint"/"discordant edge
        Output values: 
        "source": Contains v.pos = -1, indicates end of linear contig.
        "interchromosomal": Different chromosomes.
        "everted": Forward strand of larger position connected to reverse strand of reverse, indicated by outward orientation of read-pairs, may suggest tandem duplication.
        "forward": Both vertex/paired-reads map to forward strand
        "reverse": Both vertex/paired-reads map to reverse strand
        "discordant": Alignment distance larger/smaller than max/min insert, may indicate deletion
        "concordant": Expected alignment length between min and max insert. NOTE: Different from edge_type
         """
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
        if vmax.pos == vmin.pos and vmax.strand != vmin.strand:
            return "everted"
        if vmax.strand == 1 and vmin.strand == 1:
            return "forward"
        if vmax.strand == -1 and vmin.strand == -1:
            return "reverse"
        if vmax.pos - vmin.pos > max_insert or vmax.pos - vmin.pos < min_insert:
            return "discordant"
        return "concordant"
               
    def __repr__(self):
        """breakpoint_vertex1->breakpoint_vertex2"""
        return str(self.v1) + '->' + str(self.v2)

    def __lt__(self, other):
        return min((self.v1.chrom, self.v1.pos), (self.v2.chrom, self.v2.pos)) < min((other.v1.chrom, self.v1.pos),
                                                                                     (other.v2.chrom, self.v2.pos))


class breakpoint_graph(abstract_graph):
    """Class representing breakpoint edge derived from abstract_graph.abstract_graph
    """
    def __init__(self, graphfile=None):
        """Creates an empty graph if no graphfile provided
        Loads graph from graph file in format defined in load_graphfile"""
        abstract_graph.__init__(self)
        self.vhash = {}
        if graphfile is not None:
            self.load_graphfile(graphfile)

    def has_vertex(self, chrom, pos, strand):
        vtemp = breakpoint_vertex(chrom, pos, strand)
        if vtemp.__hash__() in self.vhash:
            return self.vhash[vtemp.__hash__()]
        else:
            return None

    def new_vertex(self, chrom, pos, strand):
        """Create, add and return new breakpoint_vertex if similar vertex not already present"""
        v = self.has_vertex(chrom, pos, strand)
        if v is not None:
            return v
        v = breakpoint_vertex(chrom, pos, strand, graph=self)
        self.vhash[v.__hash__()] = v
        return v

    def new_edge(self, v1, v2, edge_type='discordant', hom=None, hom_seq=None):
        """Create, add and return breakpoint_edge to current graph. Recommend using "add_edge()". "new_edge()" may incorrectly add duplicate edges
        Arguments:
        v1,v2: breakpoint_vertex (These need to be vertices(objects) from current breakpoint graph)
        edge_type = "breakpoint"/"discordant"/"concordant"/"source"/"sequence" """
        return breakpoint_edge(v1, v2, graph=self, edge_type=edge_type, hom=hom, hom_seq=hom_seq)

    def add_vertex(self, v):
        """Create and add new vertex to graph if no similar vertex exists"""
        return self.new_vertex(v.chrom, v.pos, v.strand)

    def add_edge(self, e, edge_type='discordant'):
        """Add and return edge similar e to the graph. If e(object) already belongs to graph, return e.
        Checks if corresponding vertices already present else return None.
        If edge_type not defined, then inherits e.edge_type.
        """
        if e.edge_type is not None:
            edge_type = e.edge_type
        if e.graph is not self:
            v1 = self.has_vertex(e.v1.chrom, e.v1.pos, e.v1.strand)
            v2 = self.has_vertex(e.v2.chrom, e.v2.pos, e.v2.strand)
            if v1 is None or v2 is None:
                return None
            return self.new_edge(v1, v2, edge_type=edge_type, hom=e.hom, hom_seq=e.hom_seq)
        return e

    def load_graphfile(self, graphfile):
        """Load breakpoint_graph from file
        Format: edge_type  edge_string edge_copycount
        """
        graphfile_handle = open(graphfile)
        ll = [l.strip().split() for l in graphfile_handle]
        graphfile_handle.close()
        self.copy_count=defaultdict(lambda:0, {})
        for l in ll:
            if len(l) == 0:
                continue
            if l[0] == 'sequence':
                v1 = self.add_vertex(breakpoint_vertex(l[1]))
                v2 = self.add_vertex(breakpoint_vertex(l[2]))
                e = self.new_edge(v1, v2, edge_type='sequence')
                self.copy_count[e] = float(l[3])
            if l[0] == 'concordant':
                e = self.add_edge(breakpoint_edge(l[1], edge_type=l[0]))
                self.copy_count[e] = float(l[2])
            if l[0] == 'source' or l[0] == 'discordant' or l[0] == 'breakpoint':
                e = self.add_edge(breakpoint_edge(l[1], edge_type='discordant'))
                self.copy_count[e] = float(l[2])
        return

    def djikstra_distance(self, v1, v2, min_count=0):
        """Find shortest genomic path and distance between genomic locations (including strand) in the breakpoint graph with copy count = min_count.
        Return none if not found
        Return format:
        (distance, path, traversal_copy_count)
        distance: INT describing number of base-pairs in intermediate region
        path: list of alternating (sequence edge, strand(1/-1)) and (breakpoint edge, strand(1,-1)) such that sequence edges in first/last entries contain v1/v2
        """
        for e in self.es.values():
            if e.v1.chrom == v1.chrom and e.v1.pos <= v1.pos and v1.pos <= e.v2.pos:
                e1 = e
            if e.v1.chrom == v2.chrom and e.v1.pos <= v2.pos and v2.pos <= e.v2.pos:
                e2 = e
        if self.copy_count[e1] < min_count or self.copy_count[e2] < min_count:
            return None
        if v1.strand == v2.strand and e1 == e2 and (v2.pos - v1.pos) * v1.strand > 0:
            return (abs(v1.pos - v2.pos - 1), [(e1, v1.strand)], self.copy_count[e1])
        if v1.strand == 1:
            distance = e1.v2.pos - v1.pos
        else:
            distance = v1.pos - e1.v1.pos
        a = [(distance, [(e1, v1.strand)], self.copy_count[e1])]
        heapq.heapify(a)
        while len(a) > 0:
            d, path, cc = heapq.heappop(a)
            e, s = path[-1]
            if s == 1:
                e_new = e.v2.elist
                v = e.v2
            else:
                e_new = e.v1.elist
                v = e.v1
            e_new = [e_next for e_next in e_new if e_next.edge_type != 'sequence']
            e_search = []
            for en in e_new:
                min_c = min(cc, self.copy_count[en])
                if min_c < min_count:
                    continue
                if v == en.v1:
                    en_strand = 1
                    v_seq = en.v2
                else:
                    en_strand = -1
                    v_seq = en.v1
                if (en, en_strand) in path:
                    continue
                if (en, -1 * en_strand) in path:
                    min_c = min(min_c, self.copy_count[en] / 2.0)
                    if min_c < min_count:
                        continue
                en_seq, en_seqstrand = [(es, 1 if v_seq == es.v1 else -1) for es in v_seq.elist if es.edge_type == 'sequence'][0]
                min_c = min(min_c, self.copy_count[en_seq])
                if min_c < min_count:
                    continue
                if (en_seq, en_seqstrand) in path and not (en_seq == e1 and e1 == e2 and en_seqstrand == v1.strand):
                    continue
                if (en_seq, -1 * en_seqstrand) in path:
                    min_c = min(self.copy_count[en_seq] / 2.0, min_c)
                    if min_c < min_count:
                        continue
                if en_seq == e2 and v2.strand == en_seqstrand:
                    if v2.strand == 1:
                        dd = d + v2.pos - e2.v1.pos
                    else:
                        dd = d + e2.v2.pos - v2.pos
                    return (dd, path + [(en, en_strand), (en_seq, en_seqstrand)], min_c)
                heapq.heappush(a, (d + en_seq.v2.pos - en_seq.v1.pos + 1, path + [(en, en_strand), (en_seq, en_seqstrand)], min_c))
        return None

    def cycle_decomposition(self, w, s):

        """
        Decompose breakpoint_graph into 'simple' cycles.
        Simple cycles may contain a sequence edge atmost once along each strand.
        Reports maximum parsimonious cycles starting from thickest cycle until 80% of genomic content is covered.
        w is dict containing weights (counts) of edges
        s is source vertex, this vertex has the exception of not having a sequence edge attached"""
        def thickest_cycle(hce, wehc):
            # print hce, wehc
            v1 = hce[1].v1
            a = [(-1 * hce[0], v1)]
            heapq.heapify(a)
            hdict = {v1: (hce[0], [hce[1]], None, set())}
            seenSet = set()
            seenEdges = set()
            completed = False
            while len(a) > 0 and not completed:
                # print len(a), str(a[0]), str(hdict[a[0][1]])
                v1w, v1 = heapq.heappop(a)
                if v1 == hce[1].v1 and v1 in seenSet:
                    completed = True
                    break
                for e in v1.elist:
                    if e.edge_type == 'sequence':
                        continue
                    else:
                        v2 =  e.neighbor(v1)
                    if v2 == s:
                        v3 = v2
                        if e in hdict[v1][3]:
                            nw = min(hdict[v1][0], wehc[e] / 2)
                        else:
                            nw = min(hdict[v1][0], wehc[e])
                        if not v3 in hdict or hdict[v3][2] is None or hdict[v3][0] < nw:
                            nhdict = hdict[v1][3].copy()
                            nhdict.add(e)
                            hdict[v3] = (nw, [e], v1, nhdict)
                            seenEdges.add(e)
                    else:
                        for e2 in v2.elist:
                            if e2.edge_type == 'sequence':
                                se = e2
                                v3 = e2.neighbor(v2)
                                break
                        if e in hdict[v1][3]:
                            # print 'e is seen', e, seenEdges
                            nw = min(hdict[v1][0], wehc[e] / 2)
                        elif se in hdict[v1][3]:
                            # print 'se is seen', se, seenEdges
                            nw = min(hdict[v1][0], wehc[e], wehc[se] / 2)
                        else:
                            nw = min(hdict[v1][0], wehc[e])
                        if not v3 in hdict or hdict[v3][2] is None or hdict[v3][0] < nw:
                            nhdict = hdict[v1][3].copy()
                            nhdict.add(e)
                            nhdict.add(se)
                            hdict[v3] = (nw, [e, se], v1, nhdict)
                            # print 'seen edges', e, se, v3, hdict[v3]
                            seenEdges.add(e)
                            seenEdges.add(se)
                    if v3 in seenSet:
                        continue
                    seenSet.add(v3)
                    heapq.heappush(a, (-1 * hdict[v3][0], v3))
            if len(a) == 0 and not completed:
                print("NOT COMPLETED", hce[1].v1)
            s2Set = set()
            tc = hdict[hce[1].v1][1]
            v2 = hdict[hce[1].v1][2]
            while v2 != hce[1].v1: #and not v2 in s2Set:
                # print hce[1].v1, v2, s2Set
                s2Set.add(v2)
                if v2 not in hdict:
                    print(str(v2), str(hce[1].v1), str(tc))
                    for ee in hce[1].v1.elist:
                        print(str(ee), wehc[ee])
                tc = hdict[v2][1] + tc
                v2 = hdict[v2][2]
                s2Set.add(v2)
            #     print v2, tc
            return tc, hdict[hce[1].v1][0]

        total_amplicon_content = sum([(e.v2.pos - e.v1.pos) * w[e] for e in w if e.edge_type == 'sequence'])
        amplicon_content_covered = 0
        w2 = w.copy()
        cycle_number = 1
        cycle_list = []
        while max(w2.values()) > 0.1:
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

            # print "EEEEEEEEEEEEEE", len(w2)
            # for e in w2:
            #     print "EEEEEEEEEEEE", str(e), e.edge_type, w2[e]
            # print "EEEEEEEEE========================"
            while wei < len(we):# and (tcwmax == -1 or we[wei][0] >= tcwmax / 2.0):
                # if we[wei][1].edge_type == 'sequence':
                #     wei += 1
                #     continue
                if w2[we[wei][1]] < 0.1:
                    wei += 1
                    continue
                tc, tcw = thickest_cycle(we[wei], w2)
                if len(tc) < 2:
                    print(str(tc[0]))
                    exit()
                if tcw > tcwmax:
                    tcmax = tc
                    tcwmax = tcw
                # sumlen = sum([abs(e.v1.pos - e.v2.pos) for e in tc if e.edge_type == 'sequence'])
                # if sumlen * tcw > tchwmax:
                #     tchwmax = sumlen * tcw
                #     tchmax = tc
                #     tchw = tcw
                wei += 1
            if tcwmax == -1:
                break
            tc = tcmax
            tcw = tcwmax
            # tc = tchmax
            # tcw = tchw
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
                print("tcw is 0")
                break
            print("Cycle ", cycle_number, ": Copy count = ",tcw, tc)
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
                print(v1, "->", v2)
                cycle_edge_list.append((v1,v2))
            v1 = v2
            while ci < len(tc):
                if (tc[ci].v1.chrom, tc[ci].v1.pos, tc[ci].v1.strand) == (v1.chrom, v1.pos, v1.strand):
                    v2 = tc[ci].v2
                else:
                    v2 = tc[ci].v1
                if v1.pos == -1 or v2.pos == -1:
                    if v0 is not None:
                        print(v0, "->", v0c)
                        cycle_edge_list.append((v0,v0c))
                    print(v1, "->", v2)
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
                        print(v0, "->", v0c)
                        cycle_edge_list.append((v0,v0c))
                        v0 = None
                        v0c = None
                v1 = v2
                ci += 1
            if v0 is not None:
                print(v0, "->", v0c)
                cycle_edge_list.append((v0,v0c))
            if amplicon_content_covered <= 0.9 * total_amplicon_content or (tcw > 0.2 * cycle_list[0][1]):
                cycle_list.append([cycle_number, tcw, tc, cycle_edge_list])
                acc = tcw * sum([abs(e[1].pos - e[0].pos) for e in cycle_edge_list if -1 not in [e[0].pos, e[1].pos]])
                amplicon_content_covered += acc
            cycle_number += 1    
            # print tcw, tc
            for e in tc:
                w2[e] = w2[e] - tcw
            #     if w2[e] == 0.0:
            #         w2.pop(e)
            if amplicon_content_covered > total_amplicon_content:
                break

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


class graph_decomposition(object):
    """Class represents decomposition of a breakpoint_graph with balanced edge counts into cycles/walks
    Provides methods to merge and modify cycles into larger walks to represent architecture of complex rearrangements.
    """
    def __init__(self, segment_list=None, cycle_list=None, ilist=None, file=None, file_content=None):
        if file is not None or file_content is not None:
            self.segment_list = hg.interval_list([])
            self.segment_dict = {}
            self.cycle_dict = {}
            self.ilist = hg.interval_list([])

            if file_content:
                lines = file_content.split('\n')
            else:
                lines = str(open(file).read().decode()).split('\n')
            ll = [l.strip().split() for l in lines if len(l.strip()) > 0]
            for l in ll:
                if 'Segment' == l[0]:
                    s = hg.interval(l[2], int(l[3]), int(l[4]), info=[l[1]])
                    self.segment_dict[l[1]] = s
                    self.segment_list.append(s)
                elif 'Cycle=' in l[0]:
                    ls = l[0].split(';')
                    ci = ls[0].split('=')[1]
                    cn = float(ls[1].split('=')[1])
                    cl = []
                    for s in ls[2].split('=')[1].split(','):
                        if s[-1] == '+':
                            cl.append((s[:-1], 1))
                        else:
                            cl.append((s[:-1], -1))
                    self.cycle_dict[ci] = (ci, cn, cl)
                elif 'Interval' == l[0]:
                    self.ilist.append(hg.interval(l[2], int(l[3]), int(l[4]), info=[l[1]]))
        elif cycle_list is None:
            segment_set = hg.interval_list([hg.interval(ss[0], ss[1], ss[2]) for ss in {(s.chrom, s.start, s.end) for s in segment_list}])
            segment_set.sort()
            self.segment_list = segment_set
            self.segment_dict = {}
            seg_id = {}
            cl = []
            for s in enumerate(segment_set):
                self.segment_dict[str(s[0] + 1)] = s[1]
                seg_id[(s[1].chrom, s[1].start, s[1].end)] = str(s[0] + 1)
            for s in segment_list:
                cl.append((seg_id[(s.chrom, s.start, s.end)], s.strand))
            for ii in range(len(self.segment_list)):
                s = self.segment_list[ii]
                s.info = [seg_id[(s.chrom, s.start, s.end)]]
            self.cycle_dict = {'1':('1', 1, cl)}
            self.ilist = hg.interval_list([s[0] for s in segment_set.merge_clusters(extend=1)])
            for ii in range(len(self.ilist)):
                self.ilist[ii].info = [str(ii)]
        else:
            self.segment_list = segment_list
            self.segment_dict = {s.info[0]: s for s in segment_list}
            self.cycle_dict = {c[0]:c for c in cycle_list}
            if ilist is not None:
                self.ilist = ilist
            else:
                self.ilist = hg.interval_list([s[0] for s in segment_list.merge_clusters(extend=1)])
                for ii in range(len(self.ilist)):
                    self.ilist[ii].info = [str(ii)]

    def next_seg_id(self):
        mi = 0
        for i in self.segment_dict:
            if int(i) > mi:
                mi = int(i)
        return str(mi + 1)

    def next_cycle_id(self):
        mi = 1
        while str(mi) in self.cycle_dict:
            mi += 1
        return str(mi)

    def merge(self, c1, c2, si1, si2):
        cycle1 = self.cycle_dict[c1]
        cycle2 = self.cycle_dict[c2]
        # check if atmost 1 cycle has source vertex
        if '0' in [s[0] for s in cycle1[2]] and '0' in [s[0] for s in cycle2[2]]:
            raise Exception("Cannot merge 2 cycles with source vertices")
        # if cycle2 has source vertex, exchange c1,c2
        if '0' in [s[0] for s in cycle2[2]]:
            (c1, c2, si1, si2, cycle1, cycle2) = (c2, c1, si2, si1, cycle2, cycle1)
            if si1 == 0 or si1 == len(cycle1[2]) - 1:
                raise Exception("Cannot use source segment for merging")
        # check if segments overlap
        if not self.segment_dict[cycle1[2][si1][0]].intersects(self.segment_dict[cycle2[2][si2][0]]):
            raise Exception("Segments do not overlap" + str(self.segment_dict[cycle1[2][si1][0]]) + " " + str(self.segment_dict[cycle2[2][si2][0]]))
        # cnlist: (merged cn, cycle1cn, cycle2cn)
        if cycle1[1] == 0 or cycle2[1] == 0:
            raise Exception("Cycle copy numbers should be > 0 to merge")
        if cycle1[1] > cycle2[1]:
            cnlist = (cycle2[1], cycle1[1] - cycle2[1], 0.0)
        else:
            cnlist = (cycle1[1], 0.0, cycle2[1] - cycle1[1])
        seg1 = self.segment_dict[cycle1[2][si1][0]]
        seg2 = self.segment_dict[cycle2[2][si2][0]]
        seg1_found = False
        seg2_found = False
        for i in self.segment_list:
            if cycle1[2][si1][1] == 1 and (i.chrom, i.start, i.end) == (seg1.chrom, seg1.start, seg2.end):
                seg1_found = True
                ns1 = i.info[0]
                overlap1 = (ns1, cycle1[2][si1][1])
            elif cycle1[2][si1][1] == -1 and (i.chrom, i.start, i.end) == (seg1.chrom, seg2.start, seg1.end):
                seg1_found = True
                ns1 = i.info[0]
                overlap1 = (ns1, cycle1[2][si1][1])
            if cycle1[2][si1][1] == 1 and (i.chrom, i.start, i.end) == (seg1.chrom, seg2.start, seg1.end):
                seg2_found = True
                ns2 = i.info[0]
                overlap2 = (ns2, cycle1[2][si1][1])
            elif cycle1[2][si1][1] == -1 and (i.chrom, i.start, i.end) == (seg1.chrom, seg1.start, seg2.end):
                seg2_found = True
                ns2 = i.info[0]
                overlap2 = (ns2, cycle1[2][si1][1])
        if not seg1_found:
            ns1 = self.next_seg_id()
            overlap1 = (ns1, cycle1[2][si1][1])
            if cycle1[2][si1][1] == 1:
                self.segment_dict[ns1] = hg.interval(seg1.chrom, seg1.start, seg2.end, info=[ns1])
            else:
                self.segment_dict[ns1] = hg.interval(seg1.chrom, seg2.start, seg1.end, info=[ns1])
            self.segment_list.append(self.segment_dict[ns1])
        if not seg2_found:
            ns2 = self.next_seg_id()
            overlap2 = (ns2, cycle1[2][si1][1])
            if cycle1[2][si1][1] == 1:
                self.segment_dict[ns2] = hg.interval(seg1.chrom, seg2.start, seg1.end, info=[ns2])
            else:
                self.segment_dict[ns2] = hg.interval(seg1.chrom, seg1.start, seg2.end, info=[ns2])
            self.segment_list.append(self.segment_dict[ns2])
        cycle1_init = cycle1[2][:si1]
        if not cycle1[2][si1][1]:
            (overlap1, overlap2, ns1, ns2) = (overlap2, overlap1, ns2, ns1)
        if cycle1[2][si1][1] == cycle2[2][si2][1]:
            cycle2_span = cycle2[2][si2 + 1:] + cycle2[2][:si2]
        else:
            cycle2_span = [(s[0], -1 * s[1]) for s in cycle2[2][:si2][::-1] + cycle2[2][si2 + 1:][::-1]]
        cycle1_final = cycle1[2][si1 + 1:]
        mcycle = cycle1_init + [overlap1] + cycle2_span + [overlap2] + cycle1_final
        mcycle_id = self.next_cycle_id()
        self.cycle_dict[mcycle_id] = (mcycle_id, cnlist[0], mcycle)
        self.cycle_dict[c1] = (c1, cnlist[1], cycle1[2])
        self.cycle_dict[c2] = (c2, cnlist[2], cycle2[2])
        return

    def pivot(self, c1, si1, si2):
        cycle1 = self.cycle_dict[c1]
        # check if segments overlap
        if not self.segment_dict[cycle1[2][si1][0]].intersects(self.segment_dict[cycle1[2][si2][0]]):
            raise Exception("Segments do not overlap")
        # check if segments have opposite orientation
        if cycle1[2][si1][1] == cycle1[2][si2][1]:
            raise Exception("Segments should be in opposite orientation")
        seg1 = self.segment_dict[cycle1[2][si1][0]]
        seg2 = self.segment_dict[cycle1[2][si2][0]]
        seg1_found = False
        seg2_found = False
        for i in self.segment_list:
            if (i.chrom, i.start, i.end) == (seg1.chrom, seg1.start, seg2.end):
                seg1_found = True
                ns1 = i.info[0]
                overlap1 = (ns1, cycle1[2][si1][1])
            if (i.chrom, i.start, i.end) == (seg1.chrom, seg2.start, seg1.end):
                seg2_found = True
                ns2 = i.info[0]
                overlap2 = (ns2, cycle1[2][si2][1])
        if not seg1_found:
            ns1 = self.next_seg_id()
            overlap1 = (ns1, cycle1[2][si1][1])
            self.segment_dict[ns1] = hg.interval(seg1.chrom, seg1.start, seg2.end, info=[ns1])
            self.segment_list.append(self.segment_dict[ns1])
        if not seg2_found:
            ns2 = self.next_seg_id()
            overlap2 = (ns2, cycle1[2][si2][1])
            self.segment_dict[ns2] = hg.interval(seg1.chrom, seg2.start, seg1.end, info=[ns2])
            self.segment_list.append(self.segment_dict[ns2])
        cycle1_init = cycle1[2][:si1]
        if cycle1[2][si1][1] == -1:
            (overlap1, overlap2, ns1, ns2) = ((overlap2[0], -1 * overlap2[1]), (overlap1[0], -1 * overlap1[1]), ns2, ns1)
        cycle1_span = [(s[0], -1 * s[1]) for s in cycle1[2][si1 + 1:si2][::-1]]
        cycle1_final = cycle1[2][si2 + 1:]
        mcycle = cycle1_init + [overlap1] + cycle1_span + [overlap2] + cycle1_final
        mcycle_id = self.next_cycle_id()
        self.cycle_dict[mcycle_id] = (mcycle_id, cycle1[1], mcycle)
        self.cycle_dict[c1] = (c1, 0.0, cycle1[2])
        return

    def fasta_sequence(self, cycle_list=None, outfasta=None):
        if cycle_list is None:
            ccnlist = [(c[1], c[0]) for c in self.cycle_dict.values()]
            ccnlist.sort(reverse=True)
            print(ccnlist)
            cycle_list = [c[1] for c in ccnlist]
        fseq = ''
        if outfasta is not None:
            outfile = open(outfasta, 'w')
        for c in cycle_list:
            if outfasta is None:
                fseq += '>Cycle' + c + " Copy_count=" + str(self.cycle_dict[c][1]) + ";Segments=" + ','.join([seg[0] + ('+' if seg[1] == 1 else '-') for seg in self.cycle_dict[c][2]]) + '\n'
            else:
                outfile.write('>Cycle' + c + " Copy_count=" + str(self.cycle_dict[c][1]) + ";Segments=" + ','.join([seg[0] + ('+' if seg[1] == 1 else '-') for seg in self.cycle_dict[c][2]]) + '\n')
            for s in self.cycle_dict[c][2]:
                if s[0] == '0':
                    continue
                if s[1] == 1:
                    if outfasta is None:
                        fseq += self.segment_dict[s[0]].sequence(new_fa_file=self.fa_file)
                    else:
                        outfile.write(self.segment_dict[s[0]].sequence(new_fa_file=self.fa_file))
                else:
                    if outfasta is None:
                        fseq += hg.reverse_complement(self.segment_dict[s[0]].sequence(new_fa_file=self.fa_file))
                    else:
                        outfile.write(hg.reverse_complement(self.segment_dict[s[0]].sequence(new_fa_file=self.fa_file)))
            if outfasta is None:
                fseq += '\n'
            else:
                outfile.write('\n')
        if outfasta is not None:
            outfile.close()
        return fseq

    def __repr__(self):            
        s = ""
        for i in self.ilist:
            s += '\t'.join(["Interval", i.info[0], i.chrom, str(i.start), str(i.end)]) + '\n'
        for i in self.segment_list:
            s += '\t'.join(["Segment", i.info[0], i.chrom, str(i.start), str(i.end)]) + '\n'
        ccnlist = [(c[1], c[0]) for c in self.cycle_dict.values()]
        ccnlist.sort(reverse=True)
        for c in ccnlist:
            s += "Cycle=" + c[1] + ";Copy_count=" + str(c[0]) + ";Segments=" + ','.join([seg[0] + ('+' if seg[1] == 1 else '-') for seg in self.cycle_dict[c[1]][2]]) + '\n'
        return s
