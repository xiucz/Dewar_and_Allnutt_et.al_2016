

import sys
import re
import glob

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))


#n.b. taxa in final column as in utax output.. take the sum of replicate otus in a table, take the otus in sort order

f=open(sys.argv[1],'r')

g=open(sys.argv[2],'w')

taxa={}
tot=0
data={}
otus=[]
for line in f:
	if line[0]=="#":
		g.write(line)
	else:
		otus.append(line.split("\t")[0])
f.seek(0)
otus.sort(key=tokenize)


for line in f:
	if line[0]<>"#":
		data[line.split("\t")[0]]=line.split("\t")[1:]
		data[line.split("\t")[0]][-1]=data[line.split("\t")[0]][-1].rstrip("\n")
for i in otus:

	taxonomy=data[i][-1]
	
	if taxonomy in taxa.keys():
	
		c=0
		for p in taxa[taxonomy][1:]:
			c=c+1
			taxa[taxonomy][c]=float(taxa[taxonomy][c])+float(data[i][c-1])
	else:
		data[i].insert(0,i)
		taxa[taxonomy]=data[i][:-1]
		
		#print taxa[taxonomy]
		tot=tot+1
			
print "total unique taxa=", tot
	
for i in taxa.keys():

	g.write("\t".join(str(a) for a in taxa[i])+"\t"+i+"\n")
			
			