from collections import defaultdict
import sys

f = open(sys.argv[1],'r')

thresh = 3
if len(sys.argv)>2:
    thresh = int(sys.argv[2])

fs2count = defaultdict(int)
for line in f:
    if not 'N' in line:
        tokens = (line.rstrip()).split()
        fs = int(tokens[0])
        q = int(tokens[-2])
        n = int(tokens[-1])
        if q>=thresh:
            fs2count[fs]=fs2count[fs]+n

for key in sorted(fs2count):
    print str(key)+'\t'+str(fs2count[key])
