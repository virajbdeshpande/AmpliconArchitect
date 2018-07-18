from collections import defaultdict
import pysam

import hg19util as hg

f = pysam.AlignmentFile("/pedigree2/projects/namphuon/data/SCC090/pacbio/merged.bam")

segs = defaultdict(lambda: [], {})
readlen = {}

refi = hg.interval_list([hg.interval(i) for i in f.references])

segi = 1

for l in f.fetch():
    ref = l.reference_name.split(':')[0]
    ref_start = int(l.reference_name.split(':')[1].split('-')[0]) + l.reference_start
    ref_end = int(l.reference_name.split(':')[1].split('-')[0]) + l.reference_end
    s = [segi, ref, ref_start, ref_end, l.query_name, l.query_alignment_start, l.query_alignment_end, -1 if l.is_reverse else 1]
    segs[l.query_name].append(s)
    segi += 1


for q in segs:
    segs[q].sort(key=lambda x: x[5])


c = open("pacbio_cycles.txt", 'w')
ri = 0
for r in refi:
    ri += 1
    c.write("Interval\t%s\t%s\t%s\t%s\n" % (ri, i.chrom, i.start, i.end))
c.write("List of cycle segments\n")
for q in segs:
    for s in segs[q]:
        c.write('Segment\t%s\n' % ('\t'.join([str(si) for si in s])))
ci = 0
for q in segs:
    ci += 1
    c.write('Cycle%s;Copy_count=1;Segments=%s\n' % (ci, ','.join([str(s[0]) + ('+') if s[7] == 1 else '-' for s in segs[q]])))

c.close()









