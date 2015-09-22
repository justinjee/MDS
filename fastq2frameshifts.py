#takes in an ALIGNED compressed fastq file and returns a nucchart (similar to fastq2nucchart)
#HOWEVER this is different in that it organizes sequences by barcode
#ALSO this keeps track of quality scores, not just number of reads
#each barcode keeps track of its own consensus sequence and how many reads made up the sequence

import numpy
import logging
import sys
from array import array

def locatediff(a,b):
#finds the position of the first difference between two strings
        u=zip(a,b)
        x=0
        for i in u:
            if not i[0]==i[1]:
               return x
            x=x+1
	return x

gene = sys.argv[3]

barlen = 14
bar = ''
bar2seq = {} #stores barcode:(combined quality score, consensus sequence)
aligncheck = ''
alignshift = 0
endcheck = ''
endalign = 0
target = ''
if gene=='rpob':
	aligncheck = 'CTTCGGTTCCAGCCAGC'
	alignshift = 4
	endcheck = 'GACCCGTGAACG' #'TGCAGG' # to check for insertion/deletion from rev primer 'CCTGCA'
	endalign = 103 #115 #the position endcheck should be in in subseq
	target = 'AGTTCTTCGGTTCCAGCCAGCTGTCTCAGTTTATGGACCAGAACAACCCGCTGTCTGAGATTACGCACAAACGTCGTATCTCCGCACTCGGCCCAGGCGGTCT'
else:
	aligncheck = 'GGCACAAGTGCCGGAA'
	alignshift = 0
	endcheck = 'GCGCCACCCA' # to check for insertion/deletion from rev primer 'CCTGCA'
	endalign = 107 #the position endcheck should be in in subseq 
	target = 'GGCACAAGTGCCGGAAGTGAACTCGGCGCTGGTGTCGATCAATCCGCAAAACGGTGCCGTTATGGCGCTGGTCGGTGGCTTTGATTTCAATCAGAGCAAGTTTAACC'
ordshift = 33
qthresh=3

logging.basicConfig(filename='bc.log',level=logging.DEBUG)

infiles = sys.argv[1] # 'compressed/pairedReads-2-3.txt'
outfile = sys.argv[2] #'compressed/barchart-2-3.txt'

for infile in infiles.split('..'):
	fastq = open(infile,'r')
	linenum=0
	lognum = 100000
	for line in fastq:
		line = line.rstrip()
		tokens = line.split()
		sequence = tokens[1]
		alignpos = sequence.find(aligncheck)
		if alignpos>0:
			bar = sequence[(alignpos-barlen-alignshift):(alignpos-alignshift)]
			subseq = sequence[(alignpos-alignshift):] #the part of the sequence after the barcode
	
			#only save in situations where there is a frameshift
			endpos = subseq.rfind(endcheck)
			if not endpos==endalign:
				if bar in bar2seq:
					subseq = list(subseq)
					consensus = list(bar2seq[bar][1])
					for i in range(min([len(subseq),len(consensus)])):
						if consensus[i]!=subseq[i]:
							consensus[i]='N'
					bar2seq[bar][0]=bar2seq[bar][0]+1
					bar2seq[bar][1] = "".join(consensus)
				elif not 'N' in bar:
					bar2seq[bar] = [1,subseq]
		if linenum%lognum==0:
			logging.info('read line: '+str(linenum))
		linenum=linenum+1

	logging.info('read file. now writing barchart')

bnmap = {'A':0,'C':1,'G':2,'T':3,'N':4}

barnum=0
readlen = 130

of = open(outfile,'w')
seqrad = 5

indels = {}

for bar in bar2seq:
	sequence = bar2seq[bar][1]
	quality = bar2seq[bar][0]
	if quality>=qthresh:
		endpos = sequence.rfind(endcheck)
		shift = endpos-endalign
		#of.write(str(shift)+'\t'+bar+'\t'+sequence+'\n')
		shiftpos = locatediff(target,sequence)
		preseq = sequence[max(0,shiftpos-seqrad):shiftpos]
		indelseq = ''
		if shift>0:
			indelseq = sequence[shiftpos:(shiftpos+shift)]
		else:
			indelseq = target[shiftpos:(shiftpos-shift)]
		endshiftpos = max(0,shift)+shiftpos
		postseq = sequence[endshiftpos:(endshiftpos+seqrad)]
		k = str(shift)+'\t'+str(shiftpos)+'\t'+preseq+' '+indelseq+' '+postseq
		if k in indels:
			indels[k]=indels[k]+1
		else:
			indels[k]=1
for k in indels:
	of.write(k+'\t'+str(indels[k])+'\n')
 

logging.info('done.')
