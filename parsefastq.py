#takes fastq files from N lanes and stores only the ones with barcodes we will use in a single file
import time
import logging
import bpaux
import sys

starttime = time.time()

flag=''
rc=False
baseline=0
oftag = sys.argv[1]
readno=1

if oftag=='mrcaf':
	flag = 'GGCACAAGTGCCGGAA' #mrca fwd
	baseline=15
elif oftag=='mrcar':
	flag = 'AGTGCCTGGGTGGCGC' #mrca rev
	baseline=0
	rc=True
	readno=2
elif oftag=='rpobf':
	flag = 'CTTCGGTTCCAGCCAGC'
	baseline=18
elif oftag=='rpobr':
	flag = 'CCTGCACGTTCACGGGT'
	baseline=1
	rc=True
	readno=2

infiles = ['hiseq/2015-04-10-A/lane1_NoIndex_L001_R'+str(readno)+'_001.fastq','hiseq/2015-04-10-A/lane2_NoIndex_L002_R'+str(readno)+'_001.fastq']

seqnum = 1
qualnum = 3
blocknum = 4
lognum = 500000

outfiles = []
#we will create a separate file for each phase
nphases = 8
for i in range(0,nphases):
	fout = open('compressed/'+oftag+'/compressed-'+str(i)+'.fastq','w')
	outfiles.append(fout)

logging.basicConfig(filename='compress.log',level=logging.DEBUG)
readno = 0
for f in infiles:
	linenum = 0
	phase = -1
	fin = open(f,'r')
	savequal = False
	for line in fin:
		if linenum%lognum==0:
                        logging.info('processed line '+str(linenum)+'. time: '+str(time.time()-starttime))
		if linenum%blocknum==seqnum:
			readno = readno+1
			line = line.rstrip()
			loc = line.find(flag)
			phase = loc-baseline
			if phase>=0 and phase<nphases:
				if rc:
					line = bpaux.revcomp(line)
				outfiles[phase].write(str(readno)+"\t"+line+"\t")
				savequal = True
		if linenum%blocknum==qualnum and savequal:
			line = line.rstrip()
			if rc:
				line = line[::-1]
			outfiles[phase].write(line+"\n")
			savequal=False
		linenum=linenum+1
	logging.info("completed reading "+f)
logging.info("compressed.fastq generated. time="+str(time.time()-starttime))
