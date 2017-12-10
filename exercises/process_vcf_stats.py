#!/usr/bin/env python

import sys
import re

infile = sys.argv[1]
outfile = sys.argv[2]
hapstart = sys.argv[3]
hapend = sys.argv[4]
startcol = 9

with open(infile, 'r') as inhandle, open(outfile, 'w') as outhandle:
    line = inhandle.readline()
    while not line.startswith('#CHROM'):
        line = inhandle.readline()
    sampleids = line.strip().split()[startcol:]
    heterozygosites = [0] * len(sampleids)
    haplotypes = [''] * (2*len(sampleids))
    nsnps = 0
    for line in inhandle:
        fields = line.strip().split()
        mobj = re.search('AC=(\d+);.*;DP=(\d+)', fields[7])
        if not mobj:
            print("No match for:\n", fields[7], file=sys.stderr)
        elif int(mobj.group(1)) > 10 and int(mobj.group(2)) > 1000:     
            nsnps += 1
            pos = fields[1]
            rechap = True if pos >= hapstart and pos <= hapend else False 
            bases = (fields[3], fields[4])
            for ind in range(0, len(sampleids)):
                call = fields[ind + startcol].split(':')[0]
                gt = tuple(map(int, call.split('/') ) )
                if gt[0] != gt[1]:
                    heterozygosites[ind] += 1
                if rechap:
                    haplotypes[2*ind] += bases[ gt[0] ]
                    haplotypes[2*ind+1] += bases[ gt[1] ]
    heterozygosites = [ x / nsnps for x in heterozygosites ]
    for ind,het in zip(sampleids, heterozygosites):
        print(ind, round(het, 3), sep="\t", file=outhandle)
    hap_counts = {}
    for hap in haplotypes:
        hap_counts[hap] = hap_counts.get(hap, 0) + 1
    for hap,count in hap_counts.items():
        print(hap, count, sep="\t", file=outhandle)


