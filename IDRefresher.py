#!/usr/bin/env python
# -*- coding: utf-8 -*-http://localhost:8888/tree
from collections import defaultdict
from pprint import pprint
import argparse
import sys
import re

parser = argparse.ArgumentParser(description='A script for getting new IDs from PreviousVersionID')
parser.add_argument('-l','--list', help='previous ID list file',required=True)
parser.add_argument('-g','--gff',help='gff file name', required=True)
args = parser.parse_args()

with open(args.list) as f:
	oldids = file.read(f).splitlines()
oldids.sort()

gdi = defaultdict(list)
with open(args.gff) as f:
	for line in f:
		parts = line.split("\t")
		try:
			if parts[2] == "gene":
				particles = parts[8].split(";")
				particles[0]=particles[0][3:]
				for num in particles:
					if num.startswith('PreviousVersionID'):
						mf = re.match(r'PreviousVersionID=(.*)', num)
						if mf.group(1) in oldids:
							gdi[mf.group(1)].append(particles[0])
		except IndexError:
			pass					
for k, v in gdi.iteritems():
	print '%s: %s' % (k, v)
if set(oldids).difference(gdi.keys()):
        print "Warning: some of the IDs were not found in the list of PreviousIDs."
        print set(oldids).difference(gdi.keys())
