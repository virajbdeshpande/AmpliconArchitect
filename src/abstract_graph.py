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


# This file defines classes and methods for an abstract undirected graph, vertex and edge.


import logging

class abstract_vertex(object):
	"""Class describing a graph vertex.
	Attributes:
	elist: List of abstract_edges
	vid: (optional) ID for the abstract_vertex
	graph: (optional) abstract_graph to which the vertex belongs"""
	def __init__(self, vid=0, graph=None):
		"""Initiate vertex with optional vid and graph"""
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
		"""Return list of vertices connected to abstract_vertex by a direct edge"""
		return [e.v2 for e in self.elist]

	def __hash__(self):
		"""Return hash based on vid to allow to efficiently check for presence of vid in graph, etc"""
		return self.vid

	def __repr__(self):
		"""Vertex is represented by vid"""
		return str(self.vid)


class abstract_edge(object):
	"""Class describing a graph edge.
	Attributes:
	v1, v2: Ordered pair of vertices connected by the edge
	eid: (optional) ID for the abstract_edge
	graph: (optional) abstract_graph to which the vertex belongs."""
	def __init__(self, v1, v2, eid=0, graph=None, update_vertices=True):
		"""Initiate edge
		Arguments: v1, v2, (optional)eid, (optional) graph.
		update_vertices: (optional True/False) to update vertices to include edge in v1.elist, v2.elist. (default=True)"""
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
		"""Given a vertex, return its neighbor along the edge"""
		if v == self.v1:
			return self.v2
		if v == self.v2:
			return self.v1
		raise Exception("Edge not connected to vertex")

	def __hash__(self):
		"""Return hash based on eid to allow to efficiently check for presence of eid in graph, etc"""
		return self.eid

	def length(self):
		"""Not implemented"""
		pass

	def __repr__(self):
		"""String representation of the form v1<->v2."""
		return str(self.v1) + '<->' + str(self.v2)


class abstract_graph(object):
	"""Class describing a graph.
	Attributes:
	vs: Dictionary from vid/key to vertex
	es: Dictionary from eid/key to edge
	max_vid: (internal) max_vid, used to assign vid for new vertex. Suggested to use function next_vid.
	max_eid: (internal) max_eid, used to assign eid for new edge. Suggested to use function next_eid."""
	def __init__(self):
		"""Initiate empty graph"""
		self.es = {}  # key -->edges
		self.vs = {}  # key -->vertices
		#self.logger = logging.getLogger('Algae')
		self.max_eid = 1
		self.max_vid = 1

	def include_vertex(self, v):
		"""Include orphan abstract_vertex in graph and update vertex.graph to point to self"""
		if v.vid in self.vs and self.vs[v.vid] is not v:
			raise "Adding vertex with duplicate vid"
		if v.graph is not None and v.graph is not self:
			raise "Adding vertex from another graph"
		if v.graph is None:
			v.graph = self
		self.vs[v.vid] = v

	def include_edge(self, e):
		"""Include orphan abstract_edge in graph and update edge.graph to point to self. Vertices should be updated separately"""
		if e.eid in self.es and self.es[e.eid] is not e:
			raise "Adding edge with duplicate eid"
		if e.graph is not None and e.graph is not self:
			raise "Adding edge from another graph"
		if e.graph is None:
			e.graph = self
		self.es[e.eid] = e

	def next_eid (self):
		"""Find the next eid available for assignment to new edge"""
		while self.max_eid in self.es or -1 * self.max_eid in self.es:
			self.max_eid += 1
		return self.max_eid

	def next_vid (self):
		"""Find the next vid available for assignment to new vertex"""
		while self.max_vid in self.vs or -1 * self.max_vid in self.vs:
			self.max_vid += 1
		return self.max_vid
