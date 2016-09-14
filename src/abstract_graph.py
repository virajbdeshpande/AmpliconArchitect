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


# This file defines classes and methods for an abstract undirected graph, vertex and edge.


import logging

class abstract_vertex(object):
    def __init__(self, vid=0, graph=None):
        self.elist = []
        self.vid = vid  # vertexid
        self.graph = graph
        if self.vid == 0 and self.graph is not None:
                self.vid = self.graph.next_vid()
        if self.graph is not None:
                if vid in graph.vs:
                        raise Exception("Adding with duplicate vid")
                self.graph.include_vertex(self)

    def neighbors(self):
        return [e.v2 for e in self.oute]

    def __hash__(self):
        return self.vid

    def __repr__(self):
        return str(self.vid)


class abstract_edge(object):
    def __init__(self, v1, v2, eid=0, graph=None, update_vertices=True):
        self.v1, self.v2 = v1, v2
        self.eid = eid
        self.graph = graph
        if self.eid == 0 and self.graph is not None:
            self.eid = self.graph.next_eid()
        if self.graph is not None:
            if eid in self.graph.es:
                raise Exception("Adding edge with duplicate eid")
            self.graph.include_edge(self)
        if update_vertices:
            if v1.graph is not v2.graph:
                raise Exception("Adding edge between vertices of different graphs.")
            if graph is not None and v1.graph is not graph:
                raise Exception("Edge in different graph than vertex.")
            if graph is None and v1.graph is not None:
                graph = v1.graph
            v1.elist.append(self)
            v2.elist.append(self)

    def neighbor(self, v):
        if v == self.v1:
            return self.v2
        if v == self.v2:
            return self.v1
        raise Exception("Edge not connected to vertex")

    def __hash__(self):
        return self.eid

    def length(self):
        pass

    def __repr__(self):
        return str(self.v1) + '<->' + str(self.v2)


class abstract_graph(object):

    def __init__(self):
        self.es = {}  # key -->edges
        self.vs = {}  # key -->vertices
        #self.logger = logging.getLogger('Algae')
        self.max_eid = 1
        self.max_vid = 1

    def include_vertex(self, v):
        if v.vid in self.vs and self.vs[v.vid] is not v:
            raise "Adding vertex with duplicate vid"
        if v.graph is not None and v.graph is not self:
            raise "Adding vertex from another graph"
        if v.graph is None:
            v.graph = self
        self.vs[v.vid] = v

    def include_edge(self, e):
        if e.eid in self.es and self.es[e.eid] is not e:
            raise "Adding edge with duplicate eid"
        if e.graph is not None and e.graph is not self:
            raise "Adding edge from another graph"
        if e.graph is None:
            e.graph = self
        self.es[e.eid] = e

    def next_eid (self):
        while self.max_eid in self.es or -1 * self.max_eid in self.es:
            self.max_eid += 1
        return self.max_eid

    def next_vid (self):
        while self.max_vid in self.vs or -1 * self.max_vid in self.vs:
            self.max_vid += 1
        return self.max_vid
