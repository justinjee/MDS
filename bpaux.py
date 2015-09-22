#includes functions for converting strings of base pairs to integers and vice versa

def base2num(bases):
#converts bases into a base4 number, with the base on the left having smallest value
	bnmap = {'A':0,'C':1,'G':2,'T':3} #,'N':4}
	n = 0
	for i in range(0,len(bases)):
		n = n+bnmap[bases[i]]*(4**i)
	return n

def getBaseAtPos(n,pos):
#grabs the letter (base) at position pos from numeric code n
	nbmap = {0:'A',1:'C',2:'G',3:'T'} #,4:'N'}
	n = int(n/(4**pos))
	return nbmap[n%4]

def num2base(num,seqlen):
	bases = ''
	for i in range(0,seqlen):
		bases = bases+getBaseAtPos(num,i)
	return bases

def revcomp(seq):
	seq=seq.replace('A','X')
	seq=seq.replace('T','A')
	seq=seq.replace('X','T')
	seq=seq.replace('C','X')
	seq=seq.replace('G','C')
	seq=seq.replace('X','G')
	return seq[::-1]
