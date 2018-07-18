from collections import defaultdict
import pysam

import hg19util as hg

f = pysam.AlignmentFile("/pedigree2/projects/namphuon/data/SCC090/pacbio/merged.bam")

segs = defaultdict(lambda: [], {})
readlen = {}

refi = hg.interval_list([hg.interval(i) for i in f.references])

segi = 1

qi = 0
qindex = {}
qlist = []

for l in f.fetch():
    ref = l.reference_name.split(':')[0]
    ref_start = int(l.reference_name.split(':')[1].split('-')[0]) + l.reference_start
    ref_end = int(l.reference_name.split(':')[1].split('-')[0]) + l.reference_end
    qstart = l.query_alignment_start
    qend = l.query_alignment_end
    if l.query_name not in qindex:
        qindex[l.query_name] = qi
        qlist.append(l.query_name)
        qi += 1
    if l.is_reverse:
        qstart = l.infer_query_length() - l.query_alignment_end
        qend = l.infer_query_length() - l.query_alignment_start
    s = [segi, ref, ref_start, ref_end, l.query_name, qstart, qend, -1 if l.is_reverse else 1, l.infer_query_length()]
    segs[l.query_name].append(s)
    segi += 1


for q in segs:
    segs[q].sort(key=lambda x: x[5])


c = open("pacbio_cycles.txt", 'w')
ri = 0
for r in refi:
    ri += 1
    c.write("Interval\t%s\t%s\t%s\t%s\n" % (ri, r.chrom, r.start, r.end))
c.write("List of cycle segments\n")
for q in qlist:
    for s in segs[q]:
        c.write('Segment\t%s\n' % ('\t'.join([str(si) for si in s])))
ci = 0
for q in qlist:
    ci += 1
    c.write('Cycle%s;Copy_count=1;Segments=%s\n' % (ci, ','.join([str(s[0]) + ('+') if s[7] == 1 else '-' for s in segs[q]])))

c.close()









