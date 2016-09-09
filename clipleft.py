

import sys
from Bio import SeqIO


f=open(sys.argv[1],'r')

g=open(sys.argv[1]+".clip",'w')

fmt = sys.argv[3]

l = int(sys.argv[2])

for i in SeqIO.parse(f,fmt):
	
	SeqIO.write(i[l:],g,fmt)