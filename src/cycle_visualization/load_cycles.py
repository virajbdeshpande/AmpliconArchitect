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


#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

sys.path.insert(1, os.path.join(sys.path[0], '/pedigree/projects/extrachromosome/src/'))
import hg19util as hg


parser = argparse.\
ArgumentParser(description="Cycles File")
parser.add_argument('--cycles', dest='cycles_file',
                    help="File listing cycles in amplicon", metavar='FILE',
                    action='store', type=str, nargs=1)
args = parser.parse_args()
cycles_file = args.cycles_file[0]
ll = [l.strip().split() for l in open(cycles_file) if len(l.strip()) > 0]

segments = hg.interval_list([hg.interval(l[2], int(l[3]), int(l[4]), info=[int(l[1])]) for l in ll if l[0]=='Segment'])
for s in segments:
    if s.chrom[:3] == 'chr':
        s.info.append('Human')
    else:
        s.info.append('Viral')
segments.sort()
segment_id_dict = {s.info[0]:s for s in segments}

cycles = []
for c in [l[0].split(';') for l in ll if 'Cycle=' in l[0]]:
    c_dict = {cc.split('=')[0]:cc.split('=')[1] for cc in c}
    new_dict = {}
    new_dict['Cycle'] = int(c_dict['Cycle'])
    new_dict['Copy_count'] = float(c_dict['Copy_count'])
    new_dict['Segments'] = [(int(s[:-1]), s[-1]) for s in c_dict['Segments'].split(',')]
    cycles.append(new_dict)
cycle_id_dict = {c['Cycle']:c for c in cycles}


segments_merge = segments.merge_clusters(extend=100000)
human_intervals = [i[0] for i in segments_merge if i[1][0].info[1] == 'Human']
human_segments = hg.interval_list([s for s in segments if s.info[1] == 'Human'])
viral_intervals = [i[0] for i in segments_merge if i[1][0].info[1] == 'Viral']
viral_segments = hg.interval_list([s for s in segments if s.info[1] == 'Viral'])
all_intervals = hg.interval_list([i[0] for i in segments_merge])

print str(all_intervals)
print [str(i) for i in human_intervals], [str(i) for i in viral_intervals]
print {i:str(segment_id_dict[i]) for i in segment_id_dict}
print cycle_id_dict




