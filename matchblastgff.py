#!/usr/bin/env python
# -*- coding: utf-8 -*-
from collections import defaultdict
from pprint import pprint
import argparse
import sys
import re

parser = argparse.ArgumentParser(description='A script for matching blast hit coordinates with gff file gene coordinates')
parser.add_argument('-b','--blast', help='blast (tabular) file name',required=True)
parser.add_argument('-g','--gff',help='gff file name', required=True)
args = parser.parse_args()
bdi = defaultdict(list)
with open(args.blast) as f:
	for line in f:
		if line.startswith( 'Query_' ):
			parts = line.split("\t")
			parts[-1] = parts[-1].rstrip()
			bdi[parts[1]].append(parts[2:])

gdi = defaultdict(list)
edi = defaultdict(list)
with open(args.gff) as f:
	for line in f:
		if line.startswith( 'Pe' ):
			parts = line.split("\t")
			if parts[2] == "gene":
				mf = re.match(r'ID=(.*);N', parts[8])	
				parts[8] = mf.group(1)
				gdi[parts[0]].append(parts[2:9])	
			elif parts[2] == "exon":
				mf = re.match(r'ID=(.*):exon', parts[8])	
				parts[8] = mf.group(1)
				edi[parts[0]].append(parts[2:9])	
for key in bdi.keys():
	if key in gdi:
		for i in range(0,len(bdi[key])):
			bleft = int(bdi[key][i][6])
			bright = int(bdi[key][i][7])
			for j in range(0,len(gdi[key])):
				gleft = int(gdi[key][j][1])
				gright = int(gdi[key][j][2])
				if ((bright - bleft < 0) and (gright - gleft > 0)):
					bright, bleft = bleft, bright
				if ((gright > bleft > gleft) or (gleft > bright > gright)):
					print key, bdi[key][i], gdi[key][j]
		for i in range(0,len(bdi[key])):
			bleft = int(bdi[key][i][6])
			bright = int(bdi[key][i][7])
			for j in range(0,len(edi[key])):
				eleft = int(edi[key][j][1])
				eright = int(edi[key][j][2])
				if ((bright - bleft < 0) and (eright - eleft > 0)):
					bright, bleft = bleft, bright
				if ((eright > bleft > eleft) or (eleft > bright > eright)):
					print key, bdi[key][i], edi[key][j]
				
