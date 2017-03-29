#!/usr/bin/python
import sys
import os
import threading
import time
import subprocess
import math
import random
from Bio import SeqIO

start_time=time.time()
random.seed()
extensions = 50
# nthreads = int(sys.argv[1])
# fileDir = sys.argv[2]
nthreads = 23
fileDir = '.'
pStr = sys.argv[1]
if fileDir[-1] != "/":
	fileDir+="/"
ARTPath = "/pedigree2/projects/cancer_viral_integration/simulation/libs/art_bin_MountRainier/"

seqTech = "HS25"
readLen = 125
insert = 325

outpath = fileDir + "/../fastq/"
if not os.path.exists(outpath):
	os.makedirs(outpath)


class workerThread(threading.Thread):
	def __init__(self, threadID, name):
		threading.Thread.__init__(self)
		self.threadID = threadID
		self.name = name

	def run(self):
		print("Starting " + self.name)
		while fileStack:
			fileTuple = fileStack.pop()
			inputFName = fileTuple[0]
			inputD = fileTuple[1]
			inputFPre = os.path.basename(inputFName)

			#Tuple format
			#["coverage","copy_number","breakpoints","region_size"]
			print(self.name + " generating: " + inputFPre)
			print(str(len(fileStack)) + " items left to start")

			outF = outpath + os.path.splitext(inputFPre)[0].split("_extended")[0] + "_R"

			#Extend it 50 times, adjust coverage to be cov/50.0
			#coverage different
			currExt = 1.0
			if "ep" in inputFPre:
				currExt = float(extensions)

			threadCov = int(inputD["coverage"])*int(inputD["copy_number"])/(currExt)
			#print("InputFName: " + inputFName)
			#print("Coverage " + str(threadCov))
			#print("OutF: " + outF)
			#print("Reading: " + fileDir + inputFName)
			subprocess.call(ARTPath + "art_illumina -ss " + seqTech + " -f " + \
				str(threadCov) + " -d " + self.name + "_" + " -i " + fileDir + inputFName + \
				" -l " + str(readLen) + " -na -p -s 40 -o " + outF + " -m " + str(insert) + \
				" 2>&1 > /dev/null", shell=True)

			print("Finished " + self.name + " " + time.ctime(time.time()) + "\n Now zipping")
			subprocess.call("chmod 775 " + outF + "*",shell=True)
			subprocess.call("gzip -c " + outF + "1.fq -f > " + outF + "1.fastq.gz", shell=True)
			subprocess.call("gzip -c " + outF + "2.fq -f > " + outF + "2.fastq.gz", shell=True)
			print("Done zipping on " + self.name)
			
		return 1

uID_dict = {}
fList = os.listdir(fileDir)
#print fList[:5]
fileStack = []
paramList = ["coverage","copy_number","breakpoints","region_size"]
print("Prepping files\n")
for i in fList:
	if i.endswith(".fasta") and 'extended' not in i:
	# if i.endswith(".fasta") and i.startswith("cov"):
		print("prepping " + i)
		#three cases
		#episomalFasta
		if i.startswith("ep"):
			#read-the fasta
			with open(fileDir + i) as fasta_sequences:
				recordGen = SeqIO.parse(fasta_sequences, "fasta")
				seqObj = next(recordGen)
				seqID = seqObj.id
				seq = str(seqObj.seq)

			#re-write the FASTA
			concatSeq = reduce(lambda x,y: x+y,[seq]*extensions)
			outfname = i.rsplit(".fasta")[0] + "_extended.fasta"
			with open(fileDir + outfname,'w') as outfile:
				outfile.write(">" + seqID + "\n")
				outfile.write(concatSeq + "\n")
			paramsStr = i.rsplit("ep")[1].rsplit(".fasta")[0].rsplit("_")[1:]
			#GET RUN PARAMS AND PUT IT IN THE QUEUE



		elif i.startswith("bg"):			
			#background_integration_region
			#get run params and put it in the queue
			outfname = i
			paramsStr = i.rsplit("bg")[1].rsplit(".fasta")[0].rsplit("_")[1:]
			#put the copy number to 2
			paramsStr[1] = "1"
			#paramD = dict(zip(paramList,paramsStr))


for c in [1 4 16 32]:
	if i.endswith(".fasta") and i.startswith("cov"):
			#background_segs
			#run it at copy number 1 coverage 32
		outfname = i
		paramsStr = [str(c), "1", "0", str(10000 * 1000)]
		paramD = dict(zip(paramList,paramsStr))
		fileStack.append((outfname,paramD))

subprocess.call("chmod 775 " + fileDir + "*.fasta",shell=True)
print("Finished prepping, performing read gen\n")


threadL = []
for i in range(nthreads):
	threadL.append(workerThread(i,"t" + str(i)))
	threadL[i].start()

for t in threadL:
	t.join()

# subprocess.call("rm " + fileDir + "*_extended.fasta",shell=True)
end_time = time.time()
hours, rem = divmod(end_time-start_time, 3600)
minutes, seconds = divmod(rem, 60)
print("Finished running, elapsed time\n")
print("{:0>2} hrs, {:0>2} min, {:05.2f} sec".format(int(hours),int(minutes),seconds) + "\n")