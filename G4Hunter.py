#!/usr/bin/env python3
# Reformats G4Hunter's predictions (just for chrM for now).

import os
import pdb

cwd = os.path.dirname(__file__)

g4HunterFile = ('../git/G4-hunter/output/'
    + 'Results_Mitochondria_NC_012920_1/'
    + 'Mitochondria_NC_012920_1-Merged.txt')
outputFile = cwd + '/G4Hunter.bed'
outputLocationOnly = cwd + '/G4Hunter_locationOnly.bed'

chrom = 'chrM'
with open(g4HunterFile) as f, \
      open(outputFile, 'w') as of, \
      open(outputLocationOnly, 'w') as locationOnlyFile:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            continue
        if line.startswith('Start'):
            continue
        fields = line.split('\t')
        score1 = float(fields[4].strip())
        score = abs(score1)
        if score < 0.5:
            continue
        strand = '+' if score1 > 0 else '-'
        r = [chrom, fields[0].strip(), fields[1].strip(),
            '.', str(score), strand]
        of.write('\t'.join(r) + '\n')
        r = [chrom, fields[0].strip(), fields[1].strip(),
            '.', '1', strand]
        locationOnlyFile.write('\t'.join(r) + '\n')

