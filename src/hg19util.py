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


##This is a suite to load hg19 genome, genes, exons, repeat content and perform operations on this genome, compare variants

from bisect import bisect_left
from sets import Set
from collections import defaultdict
from time import clock
import pysam
import heapq
import copy
import os
import logging

try:
    DATA_REPO = os.environ["AA_DATA_REPO"]
except:
    logging.warning("#TIME " + '%.3f\t'%clock() + " Unable to set AA_DATA_REPO variable. Setting to working directory")
    DATA_REPO = '.'
if DATA_REPO == '.' or DATA_REPO == '':
    logging.warning("#TIME " + '%.3f\t'%clock() + " AA_DATA_REPO not set or empy. Setting to working directory")
    DATA_REPO = '.'

try:
    REF = [l.strip().split()[0] for l in open(DATA_REPO + "/reference.txt")][0]
except:
    logging.warning("#TIME " + '%.3f\t'%clock() + " Unable to find reference in $AA_DATA_REPO/reference.txt. Setting to working directory.")
    REF = ''

REF_files = defaultdict(lambda: '', {})
try:
    for l in open(DATA_REPO + '/' + REF + '/file_list.txt'):
        REF_files[l.strip().split()[0]] = l.strip().split()[1]
except:
    logging.warning("#TIME " + '%.3f\t'%clock() + " Unable to find reference in $AA_DATA_REPO/REF/file_list.txt. Setting to empty.")


class fake_fasta(object):
    def fetch(self, a=None, b=0, c=0):
        return ''.join(['N' for i in range(c - b + 1)])
try:
    fa_file = pysam.Fastafile(DATA_REPO + '/' + REF + '/' + REF_files['fa_file'])
except:
    logging.warning("#TIME " + '%.3f\t'%clock() + " Unable to open fasta file: \"" + DATA_REPO + '/' + REF + '/' + REF_files['fa_file'] + "\". Reference sequences will be set to N.")
    fa_file = fake_fasta()

chrLen_filename = DATA_REPO + '/' + REF + '/' + REF_files['chrLen_file']
duke35_filename = DATA_REPO + '/' + REF + '/' + REF_files['duke35_filename']
gene_filename = DATA_REPO + '/' + REF + '/' + REF_files['gene_filename']
exon_filename = DATA_REPO + '/' + REF + '/' + REF_files['exon_file']
oncogene_filename = DATA_REPO + '/' + REF + '/' + REF_files['oncogene_filename']
centromere_filename = DATA_REPO + '/' + REF + '/' + REF_files['centromere_filename']
conserved_regions_filename = DATA_REPO + '/' + REF + '/' + REF_files['conserved_regions_filename']
segdup_filename = DATA_REPO + '/' + REF + '/' + REF_files['segdup_filename']
complementary_nucleotide = defaultdict(lambda: 'N', {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n', 'N':'N'})
duke35 = []
duke35_exists = [True]

# Handling chromosome names, lengths, sorting, positions and addition of new chromosomes

chr_id = {}
chrName = {}
def chrNum(chrname, mode='append'):
    if chrname in chr_id:
        return chr_id[chrname]
    else:
        if mode == 'init':
            cnum = len(chr_id)
        else:
            cnum = 1000000 + len(chr_id)
        chr_id[chrname] = cnum
        chrName[cnum] = chrname
        return chr_id[chrname]

chrLen = defaultdict(lambda: 0, {})
try:
    for line in open(chrLen_filename):
        ll = line.strip().split()
        chrLen[chrNum(ll[0], mode='init')] = int(ll[1])
except:
    logging.warning("#TIME " + '%.3f\t'%clock() + " Unable to open chromosome lengths file: \"" + chrLen_filename + "\"")

chrOffset = {}
def absPos(chrname, pos=0):
    cnum = chrNum(chrname)
    if chrNum(chrname) not in chrOffset:
        chrkeys = chrName.keys()
        chrkeys.sort()
        sumlen = sum([chrLen[c] for c in chrLen if c in chrOffset])
        for i in range(len(chrkeys)):
            if chrkeys[i] not in chrOffset:
                chrOffset[chrkeys[i]] = sumlen
                sumlen += chrLen[chrkeys[i]]
            if cnum < chrkeys[i]:
                break
    return chrOffset[chrNum(chrname)] + pos 

for c in chrLen:
    ap = absPos(chrName[c]) 

def chrPos(abspos):
    for c in chrOffset:
        if chrOffset[c] < abspos and chrOffset[c] + chrLen[c] >= abspos:
            return (chrName[c], abspos - chrOffset[c])
    return None

def update_chrLen(len_list):
    for l in len_list:
        chrLen[chrNum(l[0])] = int(l[1])
    for l in len_list:
        cpos = absPos(l[0], 1)

def reverse_complement(seq):
    return ''.join([complementary_nucleotide[a] for a in seq][::-1])






class interval(object):
    def __init__(self, line, start=-1, end=-1, strand=1,
        file_format='', bamfile=None, info=''):
        self.info = ""
        self.file_format = file_format
        if type(line) == pysam.AlignedRead or type(line) == pysam.AlignedSegment:
            self.load_pysamread(line, bamfile)
        elif start == -1:
            self.load_line(line, file_format)
        elif end == -1:
            self.load_pos(line, start, start, strand)
        else:
            self.load_pos(line, start, end, strand)
        if len(info) > 0:
            self.info = info

    def load_line(self, line, file_format):
        if file_format == '':
            if len(line.strip().split()) == 1:
                self.chrom = line.split(':')[0]
                self.start = int(line.split(':')[1].split('-')[0])
                self.end = int(line.split(':')[1].split('-')[1])
                if self.start < self.end:
                    self.strand = 1
                else:
                    self.strand = -1
                return
            else:
                 file_format = 'bed'
        if file_format=='gff':
            ll = line.strip().split()
            self.chrom = ll[0]
            self.start, self.end = sorted([int(float(ll[3])), int(float(ll[4]))])
            if ll[6] =='+':
                self.strand = 1
            else:
                self.strand = -1
            self.info = {r[0: r.find('=')]: r[r.find('=') + 1: ]
                         for r in ll[8].strip().strip(';').split(';')}
            self.info['Variant'] = ll[5]
        elif file_format == 'bed':        
            ll = line.strip().split()
            self.chrom = ll[0]
            self.start, self.end = sorted([int(float(ll[1])), int(float(ll[2]))])
            if int(float(ll[2])) >= int(float(ll[1])):
                self.strand = 1
            else:
                self.strand = -1
            self.info = ll[3:]
        else:
            raise(Exception("Invalid interval format" + str(line)))

    def load_pos(self, chrom, start, end, strand):
        self.chrom = chrom 
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        if start > end:
            self.start = int(end)
            self.end = int(start)
            self.strand = -1 * strand

    def load_pysamread(self, line, bamfile):
        if bamfile is None:
            raise "Interval of pysam AlignedRead without bamfile"
        self.chrom = line.reference_name
        self.start = line.reference_start
        self.end = line.reference_end
        if line.is_reverse:
            self.strand = -1
        else:
            self.strand = 1

    def __gt__(self, y):
        if self.chrom != y.chrom:
            return chrNum(self.chrom) > chrNum(y.chrom)
        elif int(self.end) != int(y.end):
            return int(self.end) > int(y.end)
        else:
            return int(self.start) > int(y.start)

    def size(self):
        return self.end - self.start + 1

    def __str__(self):
        if len(str(self.info)) == 0:
            return '\t'.join(map(str, [self.chrom, self.start, self.end]))
        elif type(self.info) == list:
            return '\t'.join(map(str, [self.chrom, self.start, self.end] + self.info))
        elif type(self.info) == dict:
            return '\t'.join(map(str, [self.chrom, self.start, self.end] + [str(s) + '=' + str(self.info[s]) for s in self.info]))
        else:
            return '\t'.join(map(str, [self.chrom, self.start, self.end, self.info]))

    def gc_content(self):
        seq = fa_file.fetch(self.chrom, self.start, self.end)
        # if 'G' in seq:
        #     print seq, seq.count('G'), seq.count('C'), float(seq.count('G') + seq.count('C')) / len(seq)
        #     exit()
        if len(seq) == 0:
            return 0.5
        return float(seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')) / len(seq)

    def sequence(self, new_fa_file=None):
        if new_fa_file is not None:
            seq = new_fa_file.fetch(self.chrom, self.start, self.end)
        else:
            seq = fa_file.fetch(self.chrom, self.start, self.end)
        if self.strand == 1:
            return seq
        else:
            return ''.join([complementary_nucleotide[a] for a in seq][::-1])

    def intersects(self, n, extend=0, margin=0.0):
        if margin > 0.0:
            if self.intersects(interval(n.chrom, n.start, n.end - (1 - margin) * (n.end - n.start))) and self.intersects(interval(n.chrom, n.start + (1 - margin) * (n.end - n.start)), n.end):
                return True
            else:
                s = self
                if n.intersects(interval(s.chrom, s.start, s.end - (1 - margin) * (s.end - s.start))) and n.intersects(interval(s.chrom, s.start + (1 - margin) * (s.end - s.start)), s.end):
                    return True
            return False
        a = [self.chrom, max(0, self.start - extend), self.end + extend]
        b = [n.chrom, n.start, n.end]
        if (a[0] != b[0]):
            return False
        if (int(a[1])-int(b[1]))*(int(a[2])-int(b[1])) <= 0:
            return True
        if (int(a[1])-int(b[2]))*(int(a[2])-int(b[2])) <= 0:
            return True
        if (int(a[1])-int(b[1]))*(int(a[1])-int(b[2])) <= 0:
            return True
        if (int(a[2])-int(b[1]))*(int(a[2])-int(b[2])) <= 0:
            return True
        return False

    def intersection(self, y):
        if not self.intersects(y):
            return None
        return interval(self.chrom, max(self.start, y.start), min(self.end, y.end))

    def merge(self, y, extend=0):
        if not self.intersects(y, extend):
            return None
        return interval(self.chrom, min(self.start, y.start), max(self.end, y.end))

    def atomize(self, y):
        il = interval_list([self, y])
        il.sort()
        ilr = []
        if il[0].intersects(il[1]):
            ilint = il[0].intersection(il[1])
            if il[0].start < il[1].start:
                ilr.append((interval(il[0].chrom, il[0].start, ilint.start - 1), [il[0]]))
            elif il[1].start < il[0].start:
                ilr.append((interval(il[1].chrom, il[1].start, ilint.start - 1), [il[1]]))
            ilr.append((ilint, il))
            if il[0].end > il[1].end:
                ilr.append((interval(il[0].chrom, ilint.end + 1, il[0].end), [il[0]]))
            elif il[1].end > il[0].end:
                ilr.append((interval(il[1].chrom, ilint.end + 1, il[1].end), [il[1]]))
            return ilr
        else:
            return [(il[0], [il[0]]), (il[1], [il[1]])]

    def contains(self, x, y=-1, z=-1):
        if type(x) == interval:
            if self.intersects(x) and self.intersection(x).size() == x.size():
                return True
            else:
                return False
        if y != -1:
            if z == -1:
                z = y
            if self.intersects(interval(x, y, z)) and self.intersection(interval(x, y, z)).size() == interval(x, y, z).size():
                return True
        return False

    def rep_content(self):
        # logging.info("#TIME " + '%.3f\t'%clock() + " rep_content: init ")
        if self.chrom == 'chrM' or self.chrom == 'MT':
            return 5.0
        if self.chrom.strip('chr') not in map(str, range(1,23))+['X'+'Y']:
            return 1.0
        s34 = interval(self.chrom, self.start, max(self.start, self.end - 34))
        # logging.info("#TIME " + '%.3f\t'%clock() + " rep_content: to load duke ")
        if duke35_exists[0] and len(duke35) == 0:
            try:
                duke35file = open(duke35_filename)
                duke35.extend([l.strip() for l in duke35file])
                duke35file.close()
            except:
                logging.warning("#TIME " + '%.3f\t'%clock() + " rep_content: Unable to open mapability file \"" + duke35_filename + "\"." )
                duke35_exists[0] = False
                duke35.extend(["chr_Un  0   1   1"])
        # logging.info("#TIME " + '%.3f\t'%clock() + " rep_content: duke loaded")
        ictime = 0
        itime = 0
        hi = len(duke35) - 1
        lo = 0
        numiter = 0
        while hi - lo > 1:
            numiter += 1
            p = (hi + lo) / 2
            ctime = clock()
            m = interval(duke35[p])
            ictime += clock() - ctime
            ctime = clock()
            if s34.intersects(m) or m > s34:
                hi = p
            else:
                lo = p
            itime += clock() - ctime
        p = lo
        m = interval(duke35[p])
        sum_duke = 0
        len_duke = 0
        # logging.info("#TIME " + '%.3f\t'%clock() + " rep_content: found " + str(numiter) + " " + str(ictime) + " " + str(itime))
        while s34 > m or s34.intersects(m):
            if not s34.intersects(m):
                p += 1
                if p >= len(duke35) or p <= 0:
                    raise Exception('p index out of range: ' + str(p)+' '+str(lo)+' '+str(self)+' '+str(m) + ' '+str(interval(duke35[lo])))
                m = interval(duke35[p])
                continue
            repc = 5.0 if float(m.info[0]) == 0 else 1.0 / float(m.info[0])
            sum_duke += s34.intersection(m).size() * repc
            len_duke += s34.intersection(m).size()
            p += 1
            if p >= len(duke35):
                break
            m = interval(duke35[p])
        # logging.info("#TIME " + '%.3f\t'%clock() + " rep_content: done")
        # exit()
        if len_duke > 0:
            return sum_duke / len_duke
        else:
            return 1.0

    def num_unmasked(self):
        if self.chrom not in fa_file.references:
            return self.size()
        seq = fa_file.fetch(self.chrom, self.start, self.end)
        return len([c for c in seq if c in 'ACGT'])

    def segdup_uniqueness(self):
        sl = interval_list([self]).intersection(segdup_list)
        slsd = sum([self.intersection(i[1]).size() for i in sl])
        return float(self.size()) / (self.size() + slsd)
    
    def extend(self, extend_len=0):
        return interval(self.chrom, max(0, self.start - extend_len), min(self.end + extend_len, chrLen[chrNum(self.chrom)]), self.strand)


class interval_list(list, object):
    def __init__(self, ilist=None, file_format=None, sort=True):
        if ilist == None:
            ilist = []
        self.file_format = file_format
        if file_format in ['bed', 'gff']:
            self.bed_to_list(ilist)
        if file_format is None:
            list.__init__(self,ilist)
        if sort:
            self.sort()
        self.offset = None

    def bed_to_list(self, file_name):
        if file_name is not None:
            try:
                f = open(file_name)
                list.__init__(self, [interval(l, file_format=self.file_format)
                              for l in f if len(l.strip().split()) > 2
                              and l.strip()[0] != '#'])
                f.close()
            except:
                logging.warning("#TIME " + '%.3f\t'%clock() + " interval_list: Unable to open interval file \"" + file_name + "\"." )


    def merge_clusters(self, extend=0, margin=0.0):
        ml = []
        ci = None
        cl = []
        ai = 0
        for a in self[::-1]:
            ai += 1
            if ci is None or not a.intersects(ci, extend, margin):
                if ci is not None:
                    ml.append((ci, cl))
                ci = a
                cl = []
                if ai != sum([len(m[1]) for m in ml]) + 1:
                    print "divergent", ai, str(a)
                    exit()
            ci = ci.merge(a, extend)
            cl.append(a)
        if ci is not None:
            ml.append((ci,cl))
        return ml[::-1]

    def repeats(self, count=1):
        activeq = []
        if activeq is None:
            print "h1"
            exit()
        jinterval = None
        ilist = []
        for a in self[::-1]:
            while len(activeq) > 0 and not a.intersects(activeq[0][1]):
                heapq.heappop(activeq)
                if activeq is None:
                    print "h2"
                    exit()
            if len(activeq) < count and jinterval is not None:
                ilist.append((jinterval, copy.copy(aq)))
                if activeq is None:
                    print "h3"
                    exit()
                jinterval = None
            heapq.heappush(activeq, (-1 * a.start, a))
            if len(activeq) >= count:
                if jinterval is None:
                    jinterval = interval(a.chrom, activeq[0][1].start, a.end)
                    aq = copy.copy(activeq)
                else:
                    jinterval.start = min(jinterval.start, activeq[0][1].start)
                    heapq.heappush(aq, (-1 * a.start, a))
        if jinterval is not None:
            ilist.append((jinterval, copy.copy(aq)))
            jinterval = None
        return ilist[::-1]

    def intersection(self, l2, extend=0):
        si = len(self) - 1
        l2i = len(l2) - 1
        sj = len(self) - 1
        l2j = len(l2) - 1
        il = []
        while si >= 0:
            while l2i >= 0 and l2[l2i] > self[si] and not self[si].intersects(l2[l2i], extend=extend):
                l2i -= 1
            l2j = l2i
            while l2j >= 0 and self[si].intersects(l2[l2j], extend=extend):
                il.append((self[si], l2[l2j]))
                l2j -= 1
            si -= 1
        return il[::-1]

    def atomize(self, h2):
        i = 0
        j = 0
        atomlist = []
        if len(self) > 0:
            c1 = self[0]
        if len(h2) > 0:
            c2 = h2[0]
        c = None
        while i < len(self) or j < len(h2):
            # if c is not None:
            #     print "%%", i, j, str(c[0]), [str(aa) for aa in c[1]], [str(aa[0]) for aa in atomlist]
            # else:
            #     print "%%", i, j, [],  [str(aa[0]) for aa in atomlist]
            if c is not None:
                if i < len(self) and self[i] not in c[1] and (self[i].intersects(c[0], -1) or c[0] > self[i]):
                    atm = self[i].atomize(c[0])
                    atm = [(aa[0], [(lambda x: c[1][0] if x == c[0] else x)(aai) for aai in aa[1]]) for aa in atm]
                    # print "%i", [len(rr[1]) for rr in atm], [str(rr[0]) for rr in atm]
                    c = atm[-1]
                    i += 1
                    atomlist += atm[:-1]
                elif j < len(h2) and h2[j] not in c[1] and (h2[j].intersects(c[0], -1) or c[0] > h2[j]):
                    # print j, str(h2[j]), str(c[0]), c[0] > h2[j]
                    atm = c[0].atomize(h2[j])
                    atm = [(aa[0], [(lambda x: c[1][0] if x == c[0] else x)(aai) for aai in aa[1]]) for aa in atm]
                    # print "%j", [len(rr[1]) for rr in atm], [str(rr[0]) for rr in atm]
                    c = atm[-1]
                    j += 1
                    atomlist += atm[:-1]
                else:
                    atomlist.append(c)
                    # if i < len(self) and self[i] in c[1]:
                    # i += 1
                    # if j < len(h2) and h2[j] in c[1]:
                    # j += 1
                    c = None
            else:
                if i >= len(self):
                    atomlist.append((h2[j], [h2[j]]))
                    j += 1
                elif j >= len(h2):
                    atomlist.append((self[i], [self[i]]))
                    i += 1
                else:
                    atm = self[i].atomize(h2[j])
                    atomlist += atm[:-1]
                    c = atm[-1]
                    # if self[i] not in c[1]:
                    i += 1
                    # if h2[j] not in c[1]:
                    j += 1
        if c is not None:
            atomlist.append(c)
        return atomlist

    def get_repeat_content(self):
        try:
            duke35_file = open(duke35_filename)
            print "counting repeats", clock()
            self.sort()
            sum_duke = [0.0 for i in self]
            len_duke = [0.0 for i in self]
            lno = 0
            i = 0
            j = 0
            for line in duke35_file:
                lno += 1
                duke_int = interval(line)
                while not(duke_int.intersects(self[i])) and duke_int > self[i]:
                    i += 1
                if not duke_int.intersects(self[i]) and self[i] > duke_int:
                    continue
                j = i
                repc = 5.0 if float(duke_int.info[0]) == 0 else 1 / float(duke_int.info[0])
                while j < len(self) and self[j].intersects(duke_int):
                    sum_duke[j] += self[j].intersection(duke_int).size() * repc
                    len_duke[j] += self[j].intersection(duke_int).size()
                    j += 1
            duke35_file.close()
            return {self[i]:sum_duke[i] / len_duke[i] for i in range(len(interval_list))}
        except:
            logging.warning("#TIME " + '%.3f\t'%clock() + " get_repeat_content: Unable to open mapability file \"" + duke35_filename + "\"." )
            duke35_exists[0] = False
            duke35.extend(["chr_Un  0   1   1"])
            return {self[i]:1.0 for i in range(len(interval_list))}

    def offsets(self):
        if self.offset is not None:
            return self.offset
        gap = 0.1
        hratio = 0.8

        vlist = [i for i in self if chrNum(i.chrom) >= 100 and i.chrom[:3] != 'chr']
        hlist = [i for i in self if chrNum(i.chrom) < 100 or i.chrom[:3] == 'chr']
        v_count = len([i for i in self if chrNum(i.chrom) >= 100 and i.chrom[:3] != 'chr'])
        h_count = len(self) - v_count
        h_sum = sum([i.size() for i in hlist])
        v_sum = sum([i.size() for i in vlist])

        hK = len([i for i in hlist if i.size() < h_sum * gap / max(1, h_count)])
        hS = sum([i.size() for i in hlist if i.size > h_sum * gap / max(1, h_count)])
        min_hsize = hS / (max(1, h_count) / gap - hK)
        h_sum = hS + hK * min_hsize
        
        vK = len([i for i in vlist if i.size() < v_sum * gap / max(1, v_count)])
        vS = sum([i.size() for i in vlist if i.size > v_sum * gap / max(1, v_count)])
        min_vsize = vS / (max(1, v_count) / gap - vK)
        v_sum = vS + vK * min_vsize


        offset = {}

        h_start = 0
        hscale = 1 if v_count == 0 else hratio
        v_start = 0 if h_count == 0 else hratio
        vscale = 1 if h_count == 0 else (1 - hratio)

        hgap = gap / h_count if h_count > 0 else 0
        vgap = gap / v_count if v_count > 0 else 0
        hpos = h_start + (hgap / 2) * hscale
        vpos = v_start + (vgap / 2) * vscale
        for i in hlist:
            isize = max(i.size(), min_hsize)
            offset[i] = (hpos, hpos + ((1 - gap) * isize / h_sum) * hscale)
            hpos = hpos + ((1 - gap) * isize / h_sum + hgap) * hscale
        for i in vlist:
            isize = max(i.size(), min_vsize)
            offset[i] = (vpos, vpos + ((1 - gap) * isize / v_sum) * vscale)
            vpos = vpos + ((1 - gap) * isize / v_sum + vgap) * vscale
        self.offset = offset
        # for i in self:
        #     print str(i), offset[i], i.size(), hgap, h_sum, hscale, gap, hpos, vpos
        # exit()
        return offset

    def xpos(self, chrom, pos):
        offset = self.offsets()
        for i in self:
            if i.intersects(interval(chrom, max(0, pos - 1), pos)):
                o = offset[i]
                return (o[1] * (pos - i.start) + o[0] * (i.end - pos)) / (i.end - i.start)
        return None

    def offset_breaks(self):
        offset = self.offsets()
        gap = 0.1
        hratio = 0.8

        vlist = [i for i in self if chrNum(i.chrom) >= 100 and i.chrom[:3] != 'chr']
        hlist = [i for i in self if chrNum(i.chrom) < 100 or i.chrom[:3] == 'chr']
        v_count = len([i for i in self if chrNum(i.chrom) >= 100 and i.chrom[:3] != 'chr'])
        h_count = len(self) - v_count
        h_sum = sum([i.size() for i in hlist])
        v_sum = sum([i.size() for i in vlist])

        hscale = 1 if v_count == 0 else hratio
        vscale = 1 if h_count == 0 else (1 - hratio)

        hgap = gap / h_count if h_count > 0 else 0
        vgap = gap / v_count if v_count > 0 else 0

        breaks = []
        iprev = None
        for i in self:
            if iprev is None:
                iprev = i
                continue
            if i in hlist and iprev.chrom == i.chrom:
                breaks.append((offset[i][0] - hscale * hgap / 2, ':', i.chrom))
                print str(i), str(iprev), i in hlist, iprev.chrom == i.chrom
            elif i in hlist and iprev.chrom != i.chrom:
                breaks.append((offset[i][0] - hscale * hgap / 2, '--', i.chrom))
            elif i in vlist and iprev in hlist:
                breaks.append((offset[i][0] - vscale * vgap / 2, '-', i.chrom))
            elif i in vlist and i.chrom == iprev.chrom:
                breaks.append((offset[i][0] - vscale * vgap / 2, ':', i.chrom))
            else:
                breaks.append((offset[i][0] - vscale * vgap / 2, '--', i.chrom))

            iprev = i
        return breaks

    def __str__(self):
       return str(([str(i) for i in self]))



oncogene_list = interval_list(oncogene_filename, 'gff')
oncogene_list.sort()
gene_list = interval_list(gene_filename, 'gff')



exon_list = interval_list([])
def load_exons():
    if len(exon_list) > 0:
        return
    try:
        exon_file = open(exon_filename)
        exonFields = [interval(j, file_format='gff')
                      for j in exon_file.read().strip().split('\n')
                      if (len(j.strip()) > 0 and j.strip()[0] != '#' and
                          {r.split('=')[0]:r.split('=')[1]
                           for r in j.strip().split()[8].strip(';').split(';')
                          }['color'] == '000080')]
        exon_file.close()
        exon_list.extend((exonFields))
    except:
        logging.warning("#TIME " + '%.3f\t'%clock() + "unable to load exon file: \"" + exon_filename + "\"")

conserved_regions = interval_list(conserved_regions_filename, 'bed')
conserved_regions.sort()

centromere_list = interval_list(centromere_filename, 'bed')
centromere_list.sort()

segdup_list = interval_list(segdup_filename, 'bed')
segdup_list.sort()


