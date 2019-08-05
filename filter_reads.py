#!/usr/bin/python

from __future__ import print_function
from Bio import SeqIO
from Bio.SeqIO import QualityIO
import sys
import re
import gzip
import commands

if sys.argv[1][-2:]=='gz': R1=SeqIO.parse(gzip.open(sys.argv[1]),'fastq')
else: R1=SeqIO.parse(open(sys.argv[1]),'fastq')
if  sys.argv[2][-2:]=='gz': R2=SeqIO.parse(gzip.open(sys.argv[2]),'fastq')
else: R2=SeqIO.parse(open(sys.argv[2]),'fastq')
out = gzip.open(sys.argv[3],"w+")

def filterSeq(read,maxn,minqual,minqual_maxbases):
	if (len(re.findall("N",str(read.seq)))/float(len(str(read.seq)))) < float(maxn):
		if 'Y' in read.description.split()[1]: return None
		bad_bases = 0
		for x in read.letter_annotations["phred_quality"]:
			if x < minqual:
				bad_bases += 1
		if bad_bases > minqual_maxbases:	
			return None
		return read
totreads = 0
passing_reads = 0
while True:
	try:
		read=R1.next()
		read2=R2.next()
		totreads += 1
	except StopIteration:
		break
        fil1=filterSeq(read,0.1,10,100)
        fil2=filterSeq(read2,0.1,10,100)
        if (fil1!=None) and (fil2!=None):
		sys.stdout.write(fil1.id + "\t" + str(fil1.seq) + "\t" + str(fil2.seq) + "\t" + QualityIO._get_sanger_quality_str(fil1) + "\t" + QualityIO._get_sanger_quality_str(fil2) + "\n")
		passing_reads += 1
	elif (fil1!=None):
		SeqIO.write(fil1,out,"fastq")
	elif (fil2!=None):
		SeqIO.write(fil2,out,"fastq")

sys.stderr.write("\t" + str(passing_reads) + " out of " + str(totreads) + " fragments passed the filtering" + "\n")
data=commands.getstatusoutput('date')
sys.stderr.write("2nd step: filtering out duplicated fragments at " + data[1] + "\n")

