#!/usr/bin/python
from Bio import SeqIO
import sys
import os
import re
import glob
#from StringIO import StringIO # Python 2


	#usage convert.py infile informat outformat
digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))



folder = sys.argv[1] # directory containing only files to be parsed
informat = str(sys.argv[2])
outformat = str(sys.argv[3])
filelist=glob.glob(folder+"/*")

filelist.sort(key=tokenize)
print filelist


try:
	num = str(sys.argv[4])
except:
	num=0
for inputfile in filelist:
		
	f = open(inputfile,'rb')
	outputfile = inputfile[:len(inputfile)-len(informat)] + outformat
	print inputfile,outputfile
	
	g = open(outputfile,'w')
	if num<>0:
		count = SeqIO.convert(f[0:num], informat,g, outformat)
	else:
		count = SeqIO.convert(f, informat,g, outformat)

	print("Converted %i records" % count)
