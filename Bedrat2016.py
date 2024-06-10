#!/usr/bin/env python3
# Slightly reformats Bedrat 2016 data.

import os
import pdb

import pandas

cwd = os.path.dirname(__file__)

g4_all = pandas.read_csv(cwd + '/Bedrat2016_TableS2.csv')
g4_all.seqnames = 'chrM'
# just using location of G4s (not their scores)
g4_all['score'] = 1
# renaming to avoid spaces
g4_all['name'] = g4_all.Type

g4 = g4_all[ g4_all.Type=='G4' ]
g4 = g4[[ 'seqnames', 'start', 'end', 'name', 'score', 'strand' ]]
g4.to_csv(cwd + '/Bedrat2016_G4.bed', sep='\t', header=False, index=False)

g4_G4_UG4 = g4_all[ [t in ['G4', 'UG4'] for t in g4_all.Type] ]

g4_G4_UG4 = g4_G4_UG4[[
  'seqnames', 'start', 'end', 'name', 'score', 'strand' ]]
g4_G4_UG4.to_csv(
  cwd + '/Bedrat2016_G4_UG4.bed', sep='\t', header=False, index=False)

