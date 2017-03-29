#!/usr/bin/env python
# Authors: Jens Luebeck, Viraj Deshpande (2016)
# For questions contact virajbdeshpande@gmail.com
import argparse
import sys
import os
import bisect
import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from numpy import random as r
import ViraGE_sim_funcs as vs

r.seed()
maxEpisomeSize = 3000000 #change to allow user specification
minEpisomeSize = 200000  #change to allow user specification
meanEpisomeSize = 800000
rep = sys.argv[1]
covV = [sys.argv[2]]
copyNV = [sys.argv[3]]
bpointV = [sys.argv[4]]
episSize = [sys.argv[5]]

outDir = "bm" + rep + "/sim_metadata/"



#set input files
VIRAL_SEQ_FILENAME="/home/jluebeck/viral_integration/data/ViralReferenceGenomes/hpv16.fasta"
VIRAL_NAME="hpv16ref_1"
GENOME_BUILD_FASTA="/home/jluebeck/ref_genomes/hg19.fa"
CHROMS="/home/jluebeck/ref_genomes/hg19_chr_names.txt"


TARGET_LOW_COMP_REGION = False
low_comp_on = TARGET_LOW_COMP_REGION





#print episome_size_param,num_breakpoint_param,uID

flankingLength = 2000
minSpace = 2000
numControl = 1000
controlLength = 10000

#default but allow spec otherwise to allow extrachromosomal fragments
allowExtraChromosomalFrags = False


SVprobs_dict = {}

SVprobs_dict['4'] = (0.0, 0.5, 0.6, 0.4)
SVprobs_dict['5'] = (0.25, 0.5, 0.6, 0.4)
SVprobs_dict['6'] = (0.75, 0.5, 0.6, 0.4)
SVprobs_dict['7'] = (0.5, 0.5, 0.6, 0.4)
SVprobs_dict['4l'] = (0.0, 0.5, 0.6, 0.4)
SVprobs_dict['5l'] = (0.25, 0.5, 0.6, 0.4)
SVprobs_dict['6l'] = (0.75, 0.5, 0.6, 0.4)
SVprobs_dict['7l'] = (0.5, 0.5, 0.6, 0.4)

SVprobs = SVprobs[rep]


# covV = ["1","4","16","32"]
# copyNV = ["4","16","32"]
# bpointV = ["0","4","8","16","32"]
# episSize = ["40000","160000","640000","2400000"]
# covV = ["32"]
# copyNV = ["32"]
# bpointV = ["0", "8", "32"]
# episSize = ["5000000","10000000"]






#select insertion site
def selectIntSite(low_comp_on):
	#dev future note mod to allow homology insertion basis
	built = False
	while not built:
		insertSite = r.random_integers(flankingLength,rgLength-flankingLength)

		if episome_size_param:
			leftFlank = vs.restrictedExponential(episome_size_param/2,episome_size_param/4,episome_size_param*3.0/4.0)
			rightFlank = episome_size_param - leftFlank
			safeWindow = max(leftFlank,rightFlank) 

		else:
			safeWindow = maxEpisomeSize/2 + 1 
			leftFlank = vs.restrictedExponential(meanEpisomeSize/2,minEpisomeSize/2,maxEpisomeSize/2)
			rightFlank = vs.restrictedExponential(meanEpisomeSize/2,minEpisomeSize/2,maxEpisomeSize/2)

		flankStart = insertSite - leftFlank
		flankEnd = insertSite + rightFlank

		#check if on different chrs
		onStart = bisect.bisect_right(chrStartInds,insertSite-safeWindow) - 1
		onEnd = bisect.bisect_right(chrStartInds,insertSite + safeWindow) - 1
		if onStart != onEnd:
			continue

		#identify start chr and respective chromosomal coords
		onChr = chrList[onStart]

		chrNormStart = flankStart - chrStartInds[onStart]
		chrNormEnd = flankEnd - chrStartInds[onStart]
		relIntPoint = insertSite - chrStartInds[onStart]

		#extract sequence
		fragSeq = seqD[onChr][chrNormStart:chrNormEnd]
		extSeq = seqD[onChr][chrNormStart-flankingLength:chrNormEnd+flankingLength]

		#disallow ambiguos bases
		#if float(fragSeq.count('N'))/len(fragSeq) > 0.00:
		#	continue
		if "N" in extSeq:
			continue

		totLow = float(sum(1 for c in fragSeq if c.islower()))
		lcPercent = totLow/len(fragSeq)

		if low_comp_on:
			#insert lies on lower
			if not lcPercent > 0.2:
				continue
			immediateWindowStart = relIntPoint - chrNormStart - 3000
			immediateWindowEnd = relIntPoint - chrNormStart + 3000
			if not fragSeq[immediateWindowStart:immediateWindowEnd].islower():
				continue

		logfile.write("Low complexity insert mode: " + str(low_comp_on) + "\nSegment low complexity content = " + str(lcPercent) + "\n")
		built = True

	return (flankStart,flankEnd),onChr,(chrNormStart,chrNormEnd), fragSeq, extSeq, relIntPoint

def checkBPointDist(bpointList,safeZone):
	for i in range(len(bpointList)-1):
		if bpointList[i+1] - bpointList[i] < minSpace:
			return False
		elif bpointList[i] >= safeZone[0] and bpointList[i] <= safeZone[1]:
			return False

	return True

def selectNBreakPoints(epiSeq,safeZone,vSLen,epiInsSite):
	#truncate the number of breakpoints if it episome size too small:
	finalNumBPoint = num_breakpoint_param
	if float(num_breakpoint_param)/episome_size_param > 1.0/5000:
		sys.stdout.write("While designing episome of length " + str(episome_size_param) + " with " + str(num_breakpoint_param)\
			+ " the spacing limit was violated. Limiting episome to have ")
		finalNumBPoint = min(int(episome_size_param/(minSpace*2.5)),num_breakpoint_param)
		if finalNumBPoint < 0:
			finalNumBPoint = 0
		sys.stdout.write(str(finalNumBPoint) + " breakpoints.\n")


	logfile.write("Episomal breakpoints\n")
	segL = []
	#breakL = []
	endFrag = epiSeq
	safeSeg = ""
	appropriateSpacing = False
	tries = 0
	while not appropriateSpacing:
		tries +=1
		breakSites = sorted(r.randint(0,len(epiSeq),finalNumBPoint))
		breakSites = [0] + breakSites + [len(epiSeq)]
		appropriateSpacing = checkBPointDist(breakSites,safeZone)

	sys.stdout.write("Took " + str(tries) + " attempts to get valid breakpoints\n")
	for i in range(len(breakSites)-1):
		segL.append(epiSeq[breakSites[i]:breakSites[i+1]])
		if breakSites[i] < safeZone[0] and breakSites[i+1] > safeZone[1]:
			safeSeg = epiSeq[breakSites[i]:breakSites[i+1]]


	if not safeSeg:
		sys.stderr.write("Safe zone error in breakpoint genertion. Quitting.\n")
		sys.exit()

	#breakL.append(len(epiSeq)-vSLen)
	#segL.append(endFrag)

	logfile.write(str(breakSites) + "\n")

	return breakSites,segL,safeSeg

def selectBreakPoints(epiSeq,safeZone,vSLen,epiInsSite):
	sys.stdout.write("Selecting breakpoints\n")
	logfile.write("Episomal breakpoints\n")
	segL = []
	breakL = [0]
	endFrag = epiSeq
	safeSeg = ""
	while len(endFrag) > 100000:
		breakSize = vs.restrictedExponential(5000,10000,110000)
		pointinSeq = breakL[-1] + breakSize
		if pointinSeq >= safeZone[0] and pointinSeq <= safeZone[1]:
			continue
		elif len(endFrag[breakSize:]) < 50000:
			continue
		elif pointinSeq <= safeZone[0] and pointinSeq <= safeZone[1]:
			safeSeg = endFrag[:breakSize]
		
		breakL.append(pointinSeq)
		segL.append(endFrag[:breakSize])
		endFrag = endFrag[breakSize:]

		#check if it covered the epiInsSite
	if not safeSeg:
		sys.stderr.write("Safe zone error in breakpoint genertion. Quitting.\n")
		sys.exit()

	breakL.append(len(epiSeq)-vSLen)
	segL.append(endFrag)

	logfile.write(str(breakL) + "\n")

	return breakL,segL,safeSeg

def selectUnifBreakPoints(epiSeq,safeZone):
	sys.stdout.write("Selecting breakpoints\n")
	logfile.write("Episomal breakpoints\n")
	segL = []
	endFrag = epiSeq
	safeSeg = ""

	if num_breakpoint_param == 0:
		breaks = [0,len(epiSeq)]

	elif len(epiSeq)/num_breakpoint_param < vSLen + 500:
		sys.stdout.write("Forcing interval spacing\n")
		validDistList = False
		while not validDistList:
			breaks = [0]
			for i in range(num_breakpoint_param-1):
				overShoot = 0
				if breaks[-1] > safeZone[0] - (len(epiSeq)/(num_breakpoint_param*4) + 500) and breaks[-1] < safeZone[0]:
					overShoot = vSLen
				breaks.append(breaks[-1] + overShoot + vs.restrictedExponential(len(epiSeq)/(num_breakpoint_param*4), len(epiSeq)/(num_breakpoint_param*5), max(vSLen+500,len(epiSeq)/(num_breakpoint_param)-500)))
			breaks  = breaks + [len(epiSeq)]
			diffs = [j-i for i, j in zip(breaks[:-1], breaks[1:])]
			caughtBad = False
			for i in breaks:
				if i > safeZone[0] - 50 and i < safeZone[1] + 50:
					caughtBad = True
					break
			if not caughtBad:
				validDistList = True

	else:
		validDistList = False
		while not validDistList:
			breaks = sorted(r.choice(xrange(1, len(epiSeq)), num_breakpoint_param-1, replace=False))
			breaks  = [0] + breaks + [len(epiSeq)]
			diffs = [j-i for i, j in zip(breaks[:-1], breaks[1:])]
			if not min(diffs) < min(1000,len(epiSeq)/(2.0*num_breakpoint_param)):
				caughtBad = False
				for i in breaks:
					if i > safeZone[0] - 50 and i < safeZone[1] + 50:
						caughtBad = True
						break
				if not caughtBad:
					validDistList = True

	segL = [epiSeq[x[0]:x[1]] for x in zip(breaks[:-1], breaks[1:])]

	logfile.write(str(breaks) + "\n")

	safeSeg = [epiSeq[x[0]:x[1]] for x in zip(breaks[:-1], breaks[1:]) if x[0] < safeZone[0] and x[1] > safeZone[1]][0]

	return breaks,segL,safeSeg


def breakPointListToSeqD(breakL,vSLen,epiInsSite,segL,onChr,relStart):
	bPointSeqD = {}
	for ind,i in enumerate(segL):
		containsVirus = 'False'
		breakP1 = breakL[ind]
		breakP2 = breakL[ind+1]
		if breakP1 <= epiInsSite and breakP2 >= epiInsSite:
			containsVirus = 'True'

		if breakP2 >= epiInsSite:
			breakP2-=vSLen
		if breakP1 >= epiInsSite:
			breakP1-=vSLen

		breakP1+=relStart
		breakP2+=relStart

		revCSeq = str(Seq(i).reverse_complement())
		bPointSeqD[i] = (onChr,str(breakP1),str(breakP2),'False',containsVirus)
		bPointSeqD[revCSeq] = (onChr,str(breakP1),str(breakP2),'True',containsVirus)

	return bPointSeqD

def annotateEpiL(epiSeqL,bPointSeqD,logfile,relIntSite,revCompVirus):
	outFields = "chrom","start","end","reversed","contains_virus"
	logfile.write("ITEMIZED EPISOME STRUCTURE\n")
	logfile.write("\t".join(outFields) + "\n")
	for i in epiSeqL:
		infoF = list(bPointSeqD[i])
		outLineD = dict(zip(outFields,infoF))
		if outLineD["contains_virus"] == "True":
			preStart = outLineD["start"]
			preEnd = str(relIntSite)
			logfile.write("\t".join([outLineD["chrom"],preStart,preEnd,outLineD["reversed"],"True"])+"\n")
			logfile.write("\t".join([VIRAL_NAME,"0",str(vSLen),str(bool(revCompVirus) != bool(outLineD["reversed"])),"True"])+"\n") 
			postStart = str(relIntSite+1)
			postEnd = outLineD["end"]
			logfile.write("\t".join([outLineD["chrom"],postStart,postEnd,outLineD["reversed"],"True"])+"\n")

		else:
			logfile.write("\t".join(infoF) + "\n")

def getRandHumanSegs(num_human_segs,priorWindows,windowSizes):
	windows = list(priorWindows)
	rgSeqs = []
	chrRelWindows = []
	extendSeq = []
	built = 0
	while len(windows) < num_human_segs + len(priorWindows):
		windowStart = r.random_integers(flankingLength,rgLength-flankingLength)
		windowEnd = windowStart + windowSizes[built]

		#check if overlapping with current
		if vs.checkOverlappingTuples(windows,(windowStart,windowEnd)):
			continue

		#identify start chr
		onStart = bisect.bisect_right(chrStartInds,windowStart) - 1
		onEnd = bisect.bisect_right(chrStartInds,windowEnd) - 1
		if onStart != onEnd:
			continue
		onChr = chrList[onStart]
		chrNormStart = windowStart - chrStartInds[onStart]
		chrNormEnd = windowEnd - chrStartInds[onStart]

		#extract sequence
		fragSeq = seqD[onChr][chrNormStart:chrNormEnd]
		extSeq = seqD[onChr][chrNormStart-2000:chrNormEnd+2000]

		#disallow ambiguos bases
		#if float(fragSeq.count('N'))/len(fragSeq) > 0.00:
		#	continue
		if "N" in extSeq:
			continue

		#passed all checks add to lists
		windows.append((windowStart,windowEnd))
		rgSeqs.append(fragSeq)
		chrRelWindows.append((onChr,chrNormStart,chrNormEnd))
		extendSeq.append(extSeq)
		built+=1

	return windows, rgSeqs, chrRelWindows, extendSeq

def doSim():
	global covSet
	sys.stdout.write("Parameter set " + str((covLevel,copyN,num_breakpoint_param,episome_size_param))+ "\n")
	#Select viral integration site
	absWindow,onChr,relWindow, fragSeq, extSeq, relIntSite = selectIntSite(low_comp_on)
	logfile.write("Host insertion point chosen:\n")
	logfile.write(str(onChr) + "\t" + str(relIntSite) + "\n")
	logfile.write("Host breakout window chosen:\n")
	logfile.write(str(onChr) + "\t" + str(relWindow) + "\n")
	logfile.write("Virus chosen " + VIRAL_SEQ_FILENAME.rsplit("/")[-1] + ":\n")
	logfile.write("Viral genome length: " + str(vSLen) + "\n")
	
	###ROTATE VIRUS 
	#rotate and insert virus at the assigned insertion point
	rotInd = r.randint(vSLen - 600) + 300
	vSeqRot = vSeq[rotInd:] + vSeq[:rotInd]
	logfile.write("Virus rotation pivot at [0 indexed]: " + str(rotInd) + "\n")
	revCompVirus = False
	if r.random() < 0.5:
		revCompVirus = True
		vSeqRot = str(Seq(vSeqRot).reverse_complement())

	logfile.write("Viral sequence reverse complemented: " + str(revCompVirus) + "\n")


	epiInsSite = relIntSite - relWindow[0]
	initEpi = fragSeq[:epiInsSite] + vSeqRot + fragSeq[epiInsSite:]

	#Assign breakpoints to initial episome
	safeZone = (epiInsSite,epiInsSite+len(vSeqRot))
	#breakL,segL,safeSeg = selectNBreakPoints(initEpi,safeZone,vSLen,epiInsSite)
	breakL,segL,safeSeg = selectUnifBreakPoints(initEpi,safeZone)
	bPointSeqD = breakPointListToSeqD(breakL,vSLen,epiInsSite,segL,onChr,relWindow[0])
	#parse functional regions of viral genome - dict of important parts (NOT IN THIS VERSION OF TOOL)

	#rearrange
	epiSeqL = vs.conductSV(initEpi,segL,safeSeg,minEpisomeSize,SVprobs)

	#annotate rearranged seq list
	annotateEpiL(epiSeqL,bPointSeqD,logfile,relIntSite,revCompVirus)

	episomal_sequence = "".join(epiSeqL)

	logfile.write("Total episomal sequence length = " + str(len(episomal_sequence)) + "\n")
	#write fasta of episomal sequence (and obfuscate start/end) & rg hosts
	episomalFa = open(outDir + "ep" + rep + "_" + uID + ".fasta",'w')
	episomalFa.write(">episomal_seq\n")
	episomalFa.write(episomal_sequence + "\n")
	episomalFa.close()

	#select and write other sample seqs and ids
	refFa = open(outDir + "bg" + rep + "_" + uID + ".fasta",'w')
	refCoords = open(outDir + "bg" + rep + "_" + uID + ".bed",'w')
	#refCoords.write("Chr\tStart\tEnd\n")
	c = onChr
	s = relWindow[0]-2000 + 1
	e = relWindow[1]+2000 +1 
	refFa.write(">inc_seg_" + c + "_" + str(s) + "_" + str(e) + "\n")
	refFa.write(extSeq + "\n")
	refCoords.write(c + "\t" + str(s) + "\t" + str(e) + "\n")

	# #generate sequences for checking variation in coverage
	# covSet = set()
	covSet.add(absWindow)
	# controlWindows, ctrlRGSeqs, ctrlRelWindows, ctrlExtendSeqs = getRandHumanSegs(numControl,covSet,[controlLength]*numControl)


	refCoords.close()
	refFa.close()





#get list of valid chroms
chrList = vs.getRelChrs(CHROMS)

#read host reference genome
rgLength, chrLengthD, chrStartInds, seqD  = vs.readRefG(chrList,GENOME_BUILD_FASTA)

#read viral reference genome
#basename = 
vSLen, vSeq = vs.readViG(VIRAL_SEQ_FILENAME)


#########PARAMETERS FOR BENCHMARKING
#select and write other sample seqs and ids
refFa = open(outDir + "cov" + rep + ".fasta",'w')
refCoords = open(outDir + "cov" + rep + ".bed",'w')
# refCoords.write("Chr\tStart\tEnd\n")
# c = onChr
# s = relWindow[0]-2000 + 1
# e = relWindow[1]+2000 +1 
# refFa.write(">inc_seg_" + c + "_" + str(s) + "_" + str(e) + "\n")
# refFa.write(extSeq + "\n")
# refCoords.write(c + "\t" + str(s) + "\t" + str(e) + "\n")

#generate sequences for checking variation in coverage
covSet = set()
#covSet.add(absWindow)

for covLevel in covV:
	for copyN in copyNV:
		for i in bpointV:
			num_breakpoint_param = int(i)
			for j in episSize:
				episome_size_param = int(j)
				uID_list = [covLevel,copyN,i,j]
				uID = "_".join(uID_list)
				logfile = open(outDir + "ViraGE_rep" + rep + "_log_" + uID + ".txt",'w')
				doSim()
				logfile.close()
				print ""



controlWindows, ctrlRGSeqs, ctrlRelWindows, ctrlExtendSeqs = getRandHumanSegs(numControl,covSet,[controlLength]*numControl)

#write seqs to fasta
for ind, i in enumerate(ctrlRelWindows):
	refFa.write(">ref_seq_" + i[0] + "_" + str(i[1]) + "_" + str(i[2]) + "\n" )
	refFa.write(ctrlExtendSeqs[ind] + "\n")
	refCoords.write(i[0] + "\t" + str(i[1]) + "\t" + str(i[2]) + "\n")

refCoords.close()
refFa.close()





