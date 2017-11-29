#!/usr/bin/env python

import sys

infile = sys.argv[1]
outfile = sys.argv[2]
startcol = 9

with open(infile, 'r') as inhandle, open(outfile, 'w') as outhandle:
    line = inhandle.readline()
    while line.split()[0] != "#CHROM":
        line = inhandle.readline()
    sampleids = line.strip().split()[startcol:]
    heterozygosities = [0] * len(sampleids)
    nsnps = 0
    for line in inhandle:
        nsnps += 1
        fields = line.strip().split()
        genotypes = fields[startcol:]
        for ind in range(len(genotypes)):
            gt = genotypes[ind]
            call = gt.split('|')
            if call[0] != call[1]:
                heterozygosities[ind] += 1
    mean_hets = [ hets / nsnps for hets in heterozygosities ]
    for ind in range(len(mean_hets)):
        print(sampleids[ind],": ", round(mean_hets[ind], 3), sep="")




