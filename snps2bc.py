import sys
import numpy as np
infile = sys.argv[1]
countfile = sys.argv[2]
outfile = sys.argv[3]
gene = sys.argv[4]
rthresh = 3
if len(sys.argv)>4:
 rthresh=int(sys.argv[5])

f = open(infile,'r')

target = ''
thresh=0 #to prevent double-frameshift mutations from affecting the outcome
if gene=='rpob':
	target = 'AGTTCTTCGGTTCCAGCCAGCTGTCTCAGTTTATGGACCAGAACAACCCGCTGTCTGAGATTACGCACAAACGTCGTATCTCCGCACTCGGCCCAGGCGGTC' #rpob
	thresh=2
else:
	target = 'GGCACAAGTGCCGGAAGTGAACTCGGCGCTGGTGTCGATCAATCCGCAAAACGGTGCCGTTATGGCGCTGGTCGGTGGCTTTGATTTCAATCAGAGCAAGTTTAACC' #mrca
	thresh=1 #2 for rpob

bnmap = {'A':0,'C':1,'G':2,'T':3,'N':4}
bc = np.zeros((len(bnmap),len(target)))

counts = np.loadtxt(countfile)
total = sum(counts[counts[:,0]>=rthresh,1])

mutcount = 0
for line in f:
 mutcount = mutcount+1
 line = line.rstrip()
 tokens = line.split()
 seq = tokens[2]
 nreads = int(tokens[-1])
 nsnps = int(tokens[0])
 if nsnps<=thresh and nreads>=rthresh:
  for i in range(min(len(target),len(seq))):
   bc[bnmap[seq[i]],i]=bc[bnmap[seq[i]],i]+1
 
for i in range(len(target)):
 bc[bnmap[target[i]],i]=bc[bnmap[target[i]],i]+total-mutcount

np.savetxt(outfile,bc)
