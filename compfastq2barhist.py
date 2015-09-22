#reads a fastq file from Illumina
#1) Checks to see if sequencing primer annealed to the right spot by checking pad
#2) Gets barcode sequence
#3) Increments number of reads of a given barcode by 1
#note: barcodes and sequences stored as num (see bpaux.py for conversions) to save space
#example: python compfastq2barchart.py test.fastq, output.txt, rpob

import time
import logging
import sys

fnames = (sys.argv[1]).split('..')
hname = sys.argv[2] #'hist.txt'
gene = sys.argv[3]

logging.basicConfig(filename='mbd.log',level=logging.DEBUG)
logfreq = 100000
starttime = time.time()
barlen = 14
bar = ''
bardict = {} 
aligncheck = ''
alignshift = 0
if gene=='rpob':
	aligncheck = 'TTCGGTTCCAGCCAGC'
	alignshift = 5 
elif gene=='mrca':
	aligncheck = 'GGCACAAGTGCCGGAA'
	alignshift = 0
linenum = 0

for fname in fnames:
	sfile = open(fname,'r')
	for line in sfile:
		if linenum%logfreq==0:
			logging.info('processing line '+str(linenum)+'. time elapsed: '+str(time.time()-starttime))
		line = line.rstrip()
		tokens = line.split()
		seq = tokens[1]
		alignpos = seq.find(aligncheck)
		bar = seq[(alignpos-barlen-alignshift):(alignpos-alignshift)]
		if alignpos>0 and not ('N' in bar):
			if bar in bardict:
				bardict[bar] = bardict[bar]+1
			else:
				bardict[bar] = 1
		linenum = linenum+1

	logging.info('done reading file')

#outfile = open(bname,'w')
countdict = {}
for keys in bardict:
#	outfile.write(str(keys)+"\t"+str(bardict[keys])+"\n")
	if bardict[keys] in countdict:
		countdict[bardict[keys]]=countdict[bardict[keys]]+1
	else:
		countdict[bardict[keys]]=1

logging.info('done storing bardict. Total time: '+str(time.time()-starttime)+'s')

outfile = open(hname,'w')
for keys in countdict:
	outfile.write(str(keys)+"\t"+str(countdict[keys])+"\n")
