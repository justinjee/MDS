#takes in an ALIGNED compressed fastq file and returns a nucchart (similar to fastq2nucchart)
#HOWEVER this is different in that it organizes sequences by barcode
#ALSO this keeps track of quality scores, not just number of reads
#each barcode keeps track of its own consensus sequence and how many reads made up the sequence

import numpy
import logging
import sys
from array import array

barlen = 14
bar = ''
bar2seq = {} #stores barcode:(combined quality score, consensus sequence)
#rpob
aligncheck = 'CTTCGGTTCCAGCCAGC'
alignshift = 4
endcheck = 'GACCCGTGAACG' #'TGCAGG' # to check for insertion/deletion from rev primer 'CCTGCA'
endalign = 103 #115 #the position endcheck should be in in subseq
#mrca
#aligncheck = 'GGCACAAGTGCCGGAA'
#alignshift = 0
#endcheck = 'GCGCCACCCA' # to check for insertion/deletion from rev primer 'CCTGCA'
#endalign = 107 #the position endcheck should be in in subseq 
ordshift = 33
qthresh=100

logging.basicConfig(filename='bc.log',level=logging.DEBUG)

infile = sys.argv[1] # 'compressed/pairedReads-2-3.txt'
outfile = sys.argv[2] #'compressed/barchart-2-3.txt'
auxfile = sys.argv[3] # output insertion/deletion counts

fastq = open(infile,'r')

linenum=0
lognum = 100000
for line in fastq:
	line = line.rstrip()
	tokens = line.split()
	sequence = tokens[1]
	quality = tokens[2]
	alignpos = sequence.find(aligncheck)
	if alignpos>0:
		bar = sequence[(alignpos-barlen-alignshift):(alignpos-alignshift)]
		subseq = sequence[(alignpos-alignshift):] #the part of the sequence after the barcode
	
		if bar in bar2seq:
			subseq = list(subseq)
			consensus = list(bar2seq[bar][1])
			for i in range(min([len(subseq),len(consensus),len(quality)])):
				if consensus[i]!=subseq[i]:
					consensus[i]='N'
				q = ord(quality[i+(alignpos-alignshift)])-ordshift
				if bar2seq[bar][0][i]<qthresh:
					bar2seq[bar][0][i]=bar2seq[bar][0][i]+q
			bar2seq[bar][1] = "".join(consensus)
		elif not 'N' in bar:
			bar2seq[bar] = [array('H',len(quality)*[0]),subseq]
			for i in range(min([len(quality),len(subseq)])):
				q = ord(quality[i+(alignpos-alignshift)])-ordshift
				bar2seq[bar][0][i]=q
	if linenum%lognum==0:
		logging.info('read line: '+str(linenum))
	linenum=linenum+1

logging.info('read file. now writing barchart')

bnmap = {'A':0,'C':1,'G':2,'T':3,'N':4}

barnum=0
readlen = 130
frameshiftCount = {} #maps frameshift number to count (ex: insertion of 1 occurs 3 times would be 1:3)
nucchart = numpy.zeros((len(bnmap),readlen))
for bar in bar2seq:
	sequence = bar2seq[bar][1]
	quality = bar2seq[bar][0]
	if max(quality)>=qthresh:
		endpos = sequence.rfind(endcheck)
		if endpos == endalign: 
			for i in range(len(sequence)):
				if quality[i]>qthresh:
					nucchart[bnmap[sequence[i]],i]=nucchart[bnmap[sequence[i]],i]+1
		elif not 'N' in sequence: #deal with insertions/deletions
			shift = endpos-endalign
			if shift in frameshiftCount:
				frameshiftCount[shift]=frameshiftCount[shift]+1
			else:
				frameshiftCount[shift]=1
	if barnum%lognum==0:
		logging.info('processing barcode '+str(barnum)+' of '+str(len(bar2seq)))
	barnum=barnum+1
numpy.savetxt(outfile,nucchart)
numpy.savetxt(auxfile,numpy.array(frameshiftCount.items()))

logging.info('done.')
