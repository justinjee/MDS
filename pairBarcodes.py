#code based on compfastq2barhist.py
#takes 2 compressed fastq files and compares their barcodes
#returns a new file that has only consensus barcodes

import time
import logging
import sys

f1name = 'compressed/fwd/compressed-0-1.fastq'
f2name = 'compressed/revrc/compressed-6-7.fastq'
bname = 'compressed/pairedBarList-0-1.txt'

logging.basicConfig(filename='pair.log',level=logging.DEBUG)
logfreq = 100000
starttime = time.time()
barlen = 13
bar = ''
id2bar = {} 
aligncheck = 'CTTCGGTTCCAGCCAGC'
alignshift = 3 
linenum = 0
idCol = 0
seqCol = 1
qualCol = 2
unionConvention='-' #sign we'll put in front of barcode sequences to indicate they were in both files

file1 = open(f1name,'r')
for line in file1:
	if linenum%logfreq==0:
		logging.info('processing line '+str(linenum)+'. time elapsed: '+str(time.time()-starttime))
	line = line.rstrip()
	tokens = line.split()
	alignpos = line.find(aligncheck)

	sequence = tokens[seqCol]
	bar = sequence[(alignpos-barlen-alignshift):alignpos]
	if not ('N' in bar):
		id2bar[tokens[idCol]] = bar
	linenum = linenum+1

logging.info('done reading file 1')

file2 = open(f2name,'r')
linenum = 0
for line in file2:
        if linenum%logfreq==0:
                logging.info('processing line '+str(linenum)+'. time elapsed: '+str(time.time()-starttime))
        line = line.rstrip()
        tokens = line.split()
        alignpos = line.find(aligncheck)

        sequence = tokens[seqCol]
        bar = sequence[(alignpos-barlen-alignshift):alignpos]
        if tokens[idCol] in id2bar and bar == id2bar[tokens[idCol]]:
                id2bar[tokens[idCol]] = unionConvention+bar
        linenum = linenum+1


outfile = open(bname,'w')
for keys in id2bar:
	if len(id2bar[keys])>0 and unionConvention in id2bar[keys]:
		bar = id2bar[keys].split(unionConvention)
		bar = bar[1]
		outfile.write(keys+"\t"+bar+"\n")

logging.info('Total time: '+str(time.time()-starttime)+'s')

