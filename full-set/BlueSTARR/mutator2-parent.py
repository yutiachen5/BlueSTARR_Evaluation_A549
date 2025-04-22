#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import random
import ProgramName
import TempFilename
from FastaReader import FastaReader
from Pipe import Pipe
from Rex import Rex
rex=Rex()

# THESE SHOULD BE MOVED INTO A CONFIGURATION FILE:
TWO_BIT_TO_FA="twoBitToFa"
TWO_BIT="/hpc/group/igvf/A549/extra_GCs/hg38.2bit"
CHILD="/hpc/group/igvf/A549/full-set/BlueSTARR/mutator-child.py"

MAX_N=-1 # -1 = no maximum

# GLOBALS
ALPHA="ACGT"

# CLASS
class VarPosition:
    def __init__(self,text):
        if(not rex.find("(\S+):(\d+):ref=([A-Z]+):([A-Z]+),(.)",text)):
            raise Exception("Can't parse pos: "+text)
        self.chrom=rex[1]
        self.pos=int(rex[2])
        self.ref=rex[3]
        self.alleles=[rex[4],rex[5]]

# FUNCTIONS
def loadVARs(filename,MAX_N):
    n=0
    recs=[]
    with open(filename,"rt") as IN:
        for line in IN:
            # fields=line.rstrip().split()
            # if(len(fields)<2): continue
            # rec=CRE(fields[0])
            # for site in fields[1:]:
            rec=VarPosition(line)
            recs.append(rec)
            n+=1
            if(MAX_N>0 and n>=MAX_N): break
    return recs

def writeCoords(recs,seqLen,filename):
    halfLen=int(seqLen/2)
    COORD=open(filename,"wt")
    for i in range(len(recs)):
        rec=recs[i]
        chrom=rec.chrom; pos=rec.pos
        begin=pos-halfLen
        end=begin+seqLen
        print(chrom+":"+str(begin)+"-"+str(end),file=COORD) # len(interval) = seqLen

    COORD.close()

def makeFasta(coordFile,fastaFile):
    # twoBitToFa needs 0-based half-open coordinates
    cmd=TWO_BIT_TO_FA+" -noMask -seqList="+coordFile+" "+TWO_BIT+" "+fastaFile
    Pipe.run(cmd)

def runModel(recs,seqLen,fastaFile,inputFile,outputFile,JOBSIZE,model):
    halfLen=int(seqLen/2)
    OUT=open(outputFile,"wt")
    INPUTS=open(inputFile,"wt")
    reader=FastaReader(fastaFile)
    reader.doUppercase()
    varIndex=0
    n=0
    while(True):
        pair=reader.nextSequence()
        if(pair is None): break
        (defline,seq)=pair
        (actualInterval,attr)=FastaReader.parseDefline(defline)
        var=recs[varIndex]
        chrom=var.chrom;pos=var.pos; ref=var.ref; alleles=var.alleles

        for allele in alleles:
            if(seq[halfLen]!=ref[0]):
                raise Exception("Ref mismatch:"+seq[pos]+" vs "+ref[0])
            if len(ref) == len(allele): # SNV + MNV
                altSeq=seq[:halfLen]+allele+seq[(halfLen+len(allele)):]  
            elif len(ref) < len(allele): # INDEL 
                altSeq = seq[:halfLen]+allele+seq[(halfLen+len(ref)):]  # Insert `allele`, shift right
                altSeq = altSeq[:seqLen] # trimming at the end
            elif len(ref) > len(allele): # INDEL 
                missing_length=len(ref)-len(allele)
                altSeq=seq[:halfLen]+allele+seq[(halfLen+len(ref)):]+'N'*missing_length # padding at the end
            else:
                raise Exception("Unexpected mutation case")
            assert len(altSeq) == 2*halfLen
            print(actualInterval+"\tpos="+str(pos)+"\tref="+\
                    ref+"\t"+allele+"\t"+altSeq,file=INPUTS)
            n+=1
            if(n>=JOBSIZE):
                INPUTS.close()
                cmd=CHILD+" "+model+" "+inputFile
                print(cmd)
                pipe=Pipe(cmd)
                while(True):
                    line=pipe.readline()
                    if(line==None): break
                    if(rex.find("ref=",line)): print(line,file=OUT)
                pipe.close()
                n=0
                INPUTS=open(inputFile,"wt")
        varIndex+=1
    if(n>0):
        INPUTS.close()
        cmd=CHILD+" "+model+" "+inputFile
        pipe=Pipe(cmd)
        while(True):
            line=pipe.readline()
            if(line==None): break
            if(rex.find("ref=",line)): 
                print(line,file=OUT)
        pipe.close()
    reader.close()
    OUT.close(); INPUTS.close

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <model> <VARs.txt> <seq-len> <seqs-per-job> <out-file>\n")
(model,varFile,seqLen,JOBSIZE,outputFile)=sys.argv[1:]
seqLen=int(seqLen)
JOBSIZE=int(JOBSIZE)

# Make some temp files
coordFile=TempFilename.generate(".coords")
fastaFile=TempFilename.generate(".fasta")
inputFile=TempFilename.generate(".inputs")

print("Loading VARs",flush=True)
recs=loadVARs(varFile,MAX_N)
print("Writing coords file",flush=True)
localPos=writeCoords(recs,seqLen,coordFile)
print("Extracting genomic sequences",flush=True)
makeFasta(coordFile,fastaFile)
print("Running the model",flush=True)
runModel(recs,seqLen,fastaFile,inputFile,outputFile,JOBSIZE,model)

# Clean up
os.remove(coordFile); os.remove(fastaFile); os.remove(inputFile)



