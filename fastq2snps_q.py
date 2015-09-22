#takes in an ALIGNED compressed fastq file and returns a nucchart (similar to fastq2nucchart)
#HOWEVER this is different in that it organizes sequences by barcode
#ALSO this keeps track of quality scores, not just number of reads
#each barcode keeps track of its own consensus sequence and how many reads made up the sequence

import numpy
import logging
import sys
from array import array

def seqcompare(a,b):
	u=zip(a,b)
	x=[]
	for i in u: 
	    if i[0]==i[1]:
	        x.append('*') 
	    else: 
	        x.append(i[1])
	return "".join(x) 


barlen = 14
bar = ''
bar2seq = {} #stores barcode:(combined quality score, consensus sequence)
aligncheck = ''
alignshift = 0
endcheck = ''
endalign = 0
target = ''

gene = sys.argv[3]

if gene=='rpob':
	aligncheck = 'CTTCGGTTCCAGCCAGC'
	alignshift = 4
	endcheck = 'GACCCGTGAACG' #'TGCAGG' # to check for insertion/deletion from rev primer 'CCTGCA'
	endalign = 103 #115 #the position endcheck should be in in subseq
	target = 'AGTTCTTCGGTTCCAGCCAGCTGTCTCAGTTTATGGACCAGAACAACCCGCTGTCTGAGATTACGCACAAACGTCGTATCTCCGCACTCGGCCCAGGCGGTC'
elif gene=='mrca':
	aligncheck = 'GGCACAAGTGCCGGAA'
	alignshift = 0
	endcheck = 'GCGCCACCCA' # to check for insertion/deletion from rev primer 'CCTGCA'
	endalign = 107 #the position endcheck should be in in subseq 
	target = 'GGCACAAGTGCCGGAAGTGAACTCGGCGCTGGTGTCGATCAATCCGCAAAACGGTGCCGTTATGGCGCTGGTCGGTGGCTTTGATTTCAATCAGAGCAAGTTTAACC'
ordshift = 33
qthresh=90

logging.basicConfig(filename='bc.log',level=logging.DEBUG)

infiles = (sys.argv[1]).split('..') # 'compressed/pairedReads-2-3.txt'
outfile = sys.argv[2] #'compressed/barchart-2-3.txt'

for infile in infiles:
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
	
			#only save in situations where there is a SNP and no frameshift
			endpos = subseq.rfind(endcheck)
			if endpos==endalign and not target in subseq:
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

of = open(outfile,'w')

for bar in bar2seq:
	sequence = bar2seq[bar][1]
	sequence = sequence[0:endalign]
	quality = bar2seq[bar][0]
	if max(quality)>=qthresh and not 'N' in sequence:
		snpseq = seqcompare(target,sequence)
		c = endalign-snpseq.count('*')
		of.write(str(c)+'\t'+bar+'\t'+sequence+'\t'+snpseq+'\n')

logging.info('done.')
