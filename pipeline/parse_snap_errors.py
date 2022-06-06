#!/usr/bin/env python

#################
# Python script to parse log for and remove SNAP models flagged as errors
# Author: Kasey Pham
# Date updated: 12/10/2021
# Based on Daren Card's pipeline for genome annotation using MAKER
# Part of pipeline for semi-automated annotation of NitFix project long read genomes
# Usage: python parse_snap_errors.py [TAXON] [ERROR_FILE]
#################

# Import modules
import sys
import os

# Read arguments
ERROR_FILE = "error.log"
TAXON = sys.argv[1]
if len(sys.argv) > 2:
    ERROR_FILE = sys.argv[2]

error_in = open(ERROR_FILE, "r")
error_text = error_in.readlines()
error_in.close()

# if there are flagged models, there will be two or more lines in error file
if len(error_text) > 1:
    # remove summary line at the end of the error file
    flag_lines = error_text[0:(len(error_text) - 1)]
    # go through each flag and extract model number
    for line in flag_lines:
        error_model = line.split()[1].strip()
        os.system("grep -vwE {model} {taxon}_maker01.zff.length50_aed0.25.ann > {taxon}_maker01.zff.length50_aed0.25.noerr.ann".format(model = error_model, taxon = TAXON))
