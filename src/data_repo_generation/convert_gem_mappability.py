import sys
import os
from itertools import groupby
import re

# converts mappability from GEM to a limited range of mappabilities where 5 or more hits of a kmer result in zero mappability.
# write a bedgraph from the gem mappability format.
# used to generate mappability files for custom genomes in the AA data repo.
# alternate instructions for this task using more standard tools are here here: https://evodify.com/gem-mappability/

val_rep = {' ':'0', '!':'1', '"':'0.5', '#':'0.33', '$':'0.25', '%':'0'}

with open(sys.argv[1]) as infile, open(sys.argv[1] + ".bedgraph",'w') as outfile:
	for i in range(109):
		_ = infile.next()

	chars = ""
	for line in infile:
		if line.startswith("~chr"):
			if len(chars) > 0:
				groups = groupby(chars)
				result = [(label, sum(1 for _ in group)) for label, group in groups]
				ppos = 0
				for s,v in result:
					cend = ppos + v
					outfile.write("\t".join([currChrom,str(ppos),str(cend),val_rep[s]]) + "\n")
					ppos = cend

				chars = ""

			currChrom = line.rstrip()[1:]
			print(currChrom)

		else:
			rawchars = line.rstrip('\n\r')
			chars+=re.sub(r'[^!\"#$]', "%", rawchars)
			
sys.exit()
