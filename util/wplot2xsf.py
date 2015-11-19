#!/usr/bin/env python
#
# wplot2xsf.py
# transforms output of wplot into XCrySDen XSF Format
# needs wploutout/psink/struct/wout/xyz file
# by Nikolaus Frohner <e0526087@student.tuwien.ac.at> 10/01/10

import sys
import os
import re
import getopt
from math import *

def usage():
	sys.stderr.write('[I] transforms output of wplot into XCrySDen XSF Format\n')
	sys.stderr.write('[I] Usage: ' + os.path.basename(sys.argv[0]) + ' [options] case wfidx\n')
	sys.stderr.write('[I] options:\n')
	sys.stderr.write('[I] -p or --nophase\tignore phase information\n')
	sys.stderr.write('[I] -s or --noshift\tdo not shift Wannier function centres\n')
	sys.exit()
	
try:
	opts, args = getopt.getopt(sys.argv[1:], "psh", ["nophase", "noshift", "help"])
except getopt.GetoptError, err:
	print "[E] "+str(err)
	usage()

doShift = True
doPhase = True
for o,a in opts:
	if o in ("-p", "--nophase"):
		doPhase = False
	elif o in ("-s", "--noshift"):
		doShift = False
	elif o in ("-h", "--help"):
		usage()
	
if len(args) != 2:
	usage()

###########
# globals #
###########
bohr2angstrom = 0.529177249
pi = 3.141592653589793
seedName = args[0]
wfIdx = args[1]
gridName = "MLWF"+str(wfIdx)
psinkFile = seedName + '_'+str(wfIdx)+'.psink'
psiargFile = seedName + '_'+str(wfIdx)+'.psiarg'
wplotoutFile = seedName + '.wplotout'
structFile = seedName + '.struct'
xyzFile = seedName + '_centres.xyz'
woutFile = seedName + '.wout'
atomMap = {}
centers = []

#############
# functions #
#############
def b2a(x): return x*bohr2angstrom
def a2b(x): return x/bohr2angstrom
	
def setAtomMap(file):
	i = 1
	for line in file:
		ind = line.find("Z: ")
		if ind != -1:
			atomMap[str(i)] = str(int(float(line[ind+3:].strip())))
			i += 1
	
def getNumberOfAtoms(file):
	foundIt = False
	for line in file:
		if line.find("SYMMETRY ADJUSTED POSITIONS OF THE BASIS ATOMS") != -1:
			foundIt = True
			break
	
	i=1
	if foundIt is True:
		for line in file:
			if i >= 2:
				break
			i += 1

	nAtoms = 0
	for line in file:
		if len(line) is 1:
			break
		nAtoms += 1
	return nAtoms

def printPsinkData(file,fPsiarg,waitLines=1):
	i=1
	for line in file:
		if i >= waitLines:
			break
		i = i+1
	
	# get grid points count
	i = 0
	gridCount = [0,0,0]
	
	for line in file:
		if i > 2:
			break
		gridCount[i] = int(re.split('\s+',line.strip())[0])
		i += 1
	
	xlim = gridCount[1]*gridCount[2]
	ylim = gridCount[2]
	zlim = gridCount[2]
	#sys.stderr.write('debug ' + str(gridCount[0]) + ' [options] <seed> <wfidx>\n')
	#sys.exit()
	grid = [[[0 for col in range(gridCount[0])] for row in range(gridCount[1])] for row2 in range(gridCount[2])]
	phase = [[[True for col in range(gridCount[0])] for row in range(gridCount[1])] for row2 in range(gridCount[2])]
	
	# we have to convert from row major to column major mode
	i = 0
	for line in file:
		#print line,
		tmpArr = re.split('\s+',line.strip())
		for elem in tmpArr:
			#print int(floor(i/xlim)) % gridCount[0], int(floor(i/ylim)) % gridCount[1], i % zlim
			grid[int(floor(i/xlim)) % gridCount[0]][int(floor(i/ylim)) % gridCount[1]][i % zlim] = elem
			i += 1

	if fPsiarg is not None:
		i = 0
		for line in fPsiarg:
			tmpArr = re.split('\s+',line.strip())
			for elem in tmpArr:
				if cos(float(elem)) <= 0:
					phase[int(floor(i/xlim)) % gridCount[0]][int(floor(i/ylim)) % gridCount[1]][i % zlim] = False
				i += 1
	
	c = 0
	for k in range(gridCount[0]):
		for j in range(gridCount[1]):
			for i in range(gridCount[2]):
				if phase[i][j][k] is False:
					sys.stdout.write(" -")	
				print grid[i][j][k],
				c += 1
				if c % 5 == 0:
					print
							
def getWfOffset(wfIdx):
	wfOffset = [0,0,0]
	
	if doShift is False:
		return wfOffset
	
	try:
		fCentres = open(xyzFile, 'r')
	except IOError:
		return wfOffset
	
	try:
		fWout = open(woutFile, 'r')
	except IOError:
		return wfOffset
	
	i = 1
	waitLines = 2 + int(wfIdx)
	for line in fCentres:
		if i >= waitLines:
			transCentre = re.split('\s+',line.strip())[1:4]
			break
		i += 1
	
	foundIt = False
	for line in fWout:
		if line.find("Final State") != -1:
			foundIt = True
			break
	
	i = 1
	waitLines = int(wfIdx)
	for line in fWout:
		if i >= waitLines:
			origCentre = re.split('[\s\,]+',line.strip())[6:9]
			break
		i += 1
	
	wfOffset = map(float.__sub__, map(float, transCentre), map(float, origCentre)) 
	return wfOffset

def printVectors(file,pattern,readLines=3,waitLines=0,printIndices=(2,3,4),preIndex=None,preMap=None,convB2A=True,wfOffset=None):
	lineCount = 0
	printVec = False

	for line in file:
		if lineCount == readLines:
			break
		
		if printVec is True:
			if waitLines is 0:
				lineCount += 1
				lineArr = re.split('\s+',line.strip())
				print " ",
				if preIndex is not None:
					if preMap is not None and lineArr[preIndex] in preMap.keys():
						print preMap[lineArr[preIndex]]+' ',
					else:
						print lineArr[preIndex],
				print ' ',
				for i in printIndices:
					if wfOffset is not None:
						print str(b2a(float(lineArr[i])+a2b(wfOffset[i-printIndices[0]])))+' ',
						continue
					if convB2A is True:
						print str(b2a(float(lineArr[i])))+' ',
					else: 
						print str(int(lineArr[i]))+' ',
				wfOffset=None
				print
			elif waitLines > 0:
				waitLines -= 1
		
		if line.find(pattern) != -1 and printVec is False:
			printVec = True

################
# preperations #
################
try:
	fStruct = open(structFile, 'r')
except IOError:
	print '[E] Could not find file ' + fStruct
	sys.exit()
	
try:
	fPsink = open(psinkFile, 'r')
except IOError:
	print '[E] Could not find file ' + psinkFile
	sys.exit()

try:
	fWplotout = open(wplotoutFile, 'r')
except IOError:
    print '[E] Could not find file ' + wplotoutFile
    sys.exit()

if doPhase is True:
	try:
		fPsiarg = open(psiargFile, 'r')
	except IOError:
		fPsiarg = None
else:
	fPsiarg = None

nAtoms = getNumberOfAtoms(fWplotout)
fWplotout.seek(0)
setAtomMap(fStruct)
wfOffset = getWfOffset(wfIdx)

###################
# XSF file output #
###################
print " CRYSTAL"
print
print " CONVVEC"
printVectors(fWplotout, "CONVENTIONAL UNIT CELL :", 3)
print
print " PRIMVEC"
printVectors(fWplotout, "PRIMITIVE UNIT CELL :", 3)
print
print " PRIMCOORD"
print " "+str(nAtoms)+" 1"
printVectors(fWplotout, "SYMMETRY ADJUSTED POSITIONS OF THE BASIS ATOMS", nAtoms, 2, (2,3,4), 0, atomMap)
print
print "BEGIN_BLOCK_DATAGRID_3D"
print "  my_datagrid"
print "  BEGIN_DATAGRID_3D_"+gridName
printVectors(fWplotout, "3D-NET OF GRID POINTS", 1, 1, (8,9,10), None, None, False)
printVectors(fWplotout, "PLOTTING AREA", 4, 4, (4,5,6), None, None, True, wfOffset)
printPsinkData(fPsink,fPsiarg)
print "  END_DATAGRID_3D"
print "END_BLOCK_DATAGRID_3D"
