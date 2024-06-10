#!/usr/bin/env python3
# Downloads Mao 2018 data.

import subprocess

outputDir = '../pub/hg19/ChIP/G4/'
bedFile = 'GSE107690_K562_High_confidence_peaks.bed'

# download the file
subprocess.run(['curl',
  'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE107690&format=file&file=' + bedFile + '.gz'],
  stdout = open(outputDir + bedFile + '.gz', 'w'))
# uncompress the file
# ??? using a pipe might be cleaner
subprocess.run(['gunzip', '--force', outputDir + bedFile + '.gz'])

