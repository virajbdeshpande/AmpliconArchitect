import sys
import os
import bisect
from Bio import SeqIO
from Bio.Seq import Seq
from numpy import random as r 
from numpy import log2
import math



flankingLength = 2000

def conductSV(seq,segL,safeSeg,minEpisomeSize,SVprobs):
	revCompSafeSeg = str(Seq(safeSeg).reverse_complement())
	newSegL = segL
	safeCopies = 1
	origUniqueSegs = float(len(set(segL)))
	origLength = float(sum(len(s) for s in segL))

	#for i in range(int(math.ceil(log2(len(segL))+3))): #### MAKE THIS LINEAR
	i = 0
	while i < max(0,int(math.ceil(len(segL)/2.0))):
		safeCopies = segL.count(safeSeg) + segL.count(revCompSafeSeg)
		currUniqueSegs = len(set(segL));0
		meanSegsTogether = 1 + int(math.ceil(log2(len(newSegL))))
		totalLength = sum(len(s) for s in newSegL)
		#select runif strtpoint
		nonZeroLen = False
		while not nonZeroLen:
			strtP = r.random_integers(0,len(newSegL))
			#poisson seg nums
			numSegs = r.poisson(meanSegsTogether-1) + 1
			if numSegs >= len(newSegL):
				numSegs = len(newSegL)
			manipSegs = (newSegL + newSegL)[strtP:strtP+numSegs]
			loopedEndI = (strtP + numSegs) % len(newSegL)
			manipLen = sum(len(s) for s in manipSegs)
			compSubL = (newSegL + newSegL)[strtP+numSegs:len(newSegL)+strtP]
			startI = strtP
			endI = loopedEndI
			if startI != endI:
				nonZeroLen = True

		safeSegCount = manipSegs.count(safeSeg) + manipSegs.count(revCompSafeSeg)
		#deletion
		if r.random() < delProb and manipLen/origLength < 0.2 and currUniqueSegs/origUniqueSegs > 0.4:
			if safeCopies - safeSegCount >= 1 and float(len(manipSegs))/len(newSegL) <= 0.3:
				newSegL = compSubL
				i+=1
		
		#other SV
		else:
			#is duplicated
			dup = False
			if r.random() < dupProb and currUniqueSegs/origUniqueSegs < 1.5:
				dup = True
				i+=1

			if r.random() < invProb:
				revCompL = [""]*numSegs
				for ind,j in enumerate(manipSegs[::-1]):
					revCompL[ind] = str(Seq(j).reverse_complement())
				manipSegs = revCompL
				i+=1

			if dup:
				newSegL[endI:endI] = manipSegs
			else:
				newSegL = compSubL + manipSegs

			if r.random() < transLProb:
				#newSegL = newSegL[0:startI] + newSegL[endI:]
				newSegL = compSubL
				insP = r.random_integers(0,len(newSegL))
				newSegL[insP:insP] = manipSegs
				i+=1

	return newSegL			

#identify start chr
def transAbsToRelCoord(absPos,chrLengthD,chrStartInds):
	onStart = bisect.bisect_right(chrStartInds,absPos) - 1
	onChr = chrList[onStart]
	chrNormPos = absPos - chrStartInds[onStart]
	return chrNormPos,onChr


def restrictedExponential(mean, minV, maxV):
	if mean < minV or mean > maxV:
		sys.stderr.write("Improper args, distribution mean outside bounds\n")
		sys.exit()

	val = int(round(r.exponential(mean)))
	while val > maxV or val < minV:
		val = int(round(r.exponential(mean)))

	return val

#get list of valid chroms
def getRelChrs(CHROMS):
	chrSet = set()
	chromFile = open(CHROMS)
	for line in chromFile:
		chrSet.add(line.rstrip()[1:])
	chromFile.close()
	chrList = sorted(list(chrSet))
	return chrList

#takes a list of 2-tuples contaning numeric values and checks for overlap in the tuples
def checkOverlappingTuples(tupList,tupToAdd):
	for i in tupList:
		if tupToAdd[0] >= i[0] and tupToAdd[0] <= i[1]:
			return True
		elif tupToAdd[1] >= i[0] and tupToAdd[1] <= i[1]:
			return True

		return False



def readRefG(chrList,GENOME_BUILD_FASTA):
	#read human ref genome
	rgLength = 0
	chrLengthD = dict(zip(chrList,[0]*len(chrList)))
	seqD = {}
	sys.stdout.write("READING REF GENOME FASTA\n")
	fasta_sequences = SeqIO.parse(open(GENOME_BUILD_FASTA),'fasta')
	for fasta in fasta_sequences:
		name = fasta.id
		#sys.stdout.write("FINISHED READING " + name + "\n")
		if name in chrList:
			sequence = str(fasta.seq)
			seqD[name] = sequence
			rgLength+=len(sequence)
			chrLengthD[name] = len(sequence)

	#assign cumulative starting coordinates to all chromosomes for random window extraction
	chrStartInds = [0]*len(chrList)
	chrStartInds[0] = 1
	for ind,i in enumerate(chrList[1:]):
		chrStartInds[ind+1] = chrStartInds[ind] + chrLengthD[chrList[ind]]

	return rgLength, chrLengthD, chrStartInds, seqD

def readViG(VIRAL_SEQ_FILENAME):
	#read viral ref genome
	sys.stdout.write("READING VIRAL FASTA\n")
	fasta_sequences = SeqIO.parse(open(VIRAL_SEQ_FILENAME),'fasta')
	vSeq = fasta_sequences.next().seq.tostring()
	vSLen = len(vSeq)
	if vSLen < 5000:
		sys.stderr.write("VIRAL GENOME DANGEROUSLY SMALL - SIMULATION MAY CRASH\n") 

	#duplicate for circularization
	#vSeq = vSeq + vSeq
	return vSLen,vSeq