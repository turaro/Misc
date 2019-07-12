#!/usr/bin/env python
# -*- coding: utf-8 -*-http://localhost:8888/tree
from collections import defaultdict
from pprint import pprint
import argparse
import sys
import re

parser = argparse.ArgumentParser(description='A script for parsing SnpEff output')
parser.add_argument('-g','--genes', help='gene table file name',required=True)
args = parser.parse_args()
cat = defaultdict(list)

with open(args.genes) as f:
	for line in f:
		if line.startswith('#GeneName'):
			title = line.split("\t")
			title[-1] = title[-1].rstrip()
		elif line.startswith('Peaxi'):
			parts = line.split("\t")
			parts[-1] = parts[-1].rstrip()
			for i in range(0,len(parts)):
				try:
					if int(parts[i])>0:
						cat[i].append(parts[0])
				except ValueError:
					pass
for key in range(0,len(title)):
	try:
		cat[title[key]] = cat.pop(key)
	except KeyError:
					pass
del cat["variants_impact_LOW"]		#deleting all generic/high amount of genes categories
del cat["variants_impact_MODIFIER"]
del cat["variants_impact_MODERATE"]
del cat["variants_effect_downstream_gene_variant"]
del cat["variants_effect_intron_variant"]
del cat["variants_effect_upstream_gene_variant"]
del cat["variants_effect_3_prime_UTR_variant"]
del cat["variants_effect_synonymous_variant"]
del cat["variants_effect_missense_variant"]
for k, v in cat.items():
	print('%s: %s' % (k, v))
