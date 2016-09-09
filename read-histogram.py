from Bio import SeqIO
import sys
import re
import glob
import statistics


f = open(sys.argv[1],'r') #open(sys.arv[1],'r')

fmt = sys.argv[2]

bins = int(sys.argv[3])

g=open(sys.argv[1]+".hist",'w')

c=0
p=[]
for x in SeqIO.parse(f,fmt):
	c=c+1
	
	p.append(int(len(x.seq)))
		
print p[:100]
mean1 =statistics.mean(p)
stdev1=statistics.stdev(p)
h1=max(p)
l1=min(p)

step=bins#int(h1/bins)

for i in range(0,h1,step):

	t = sum(1 for v in p if v >= i and v < (i+step))

	print i,i+step,t

print "file\tNum\tmean\tstdev\thigh\tlow"
print "%s\t%d\t%d\t%d\t%d\t%d\t" %(i,c,mean1,stdev1,h1,l1)
	


		






