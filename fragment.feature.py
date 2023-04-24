# coding:utf8
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import re
import json
import pysam
import random
import numpy as np
import gzip
import math
import pprint
#import beeprint
from collections import Counter


def isSoftClipped(cigar):
    # Op BAM Description
    # M 0 alignment match (can be a sequence match or mismatch)
    # I 1 insertion to the reference
    # D 2 deletion from the reference
    # N 3 skipped region from the reference
    # S 4 soft clipping (clipped sequences present in SEQ)
    # H 5 hard clipping (clipped sequences NOT present in SEQ)
    # P 6 padding (silent deletion from padded reference)
    # = 7 sequence match
    # X 8 sequence mismatch
    for (op, count) in cigar:
        if op in [4, 5, 6]:
            return True
    return False


def readIterator(filename):
    if os.path.exists(filename) and (os.path.exists(
            filename.replace(".bam", ".bai"))
                                     or os.path.exists(filename + ".bai")):
        input_file = pysam.Samfile(filename, "rb")
        for read in input_file.fetch():
            yield read
        input_file.close()


def calculatePerBase(filenames, sample, sample_type):
    maxLength = 500
    minLength = 50
    total = 0
    result = {}
    for read in readIterator(filenames):
        if read.mapq < 5:
            continue
        if read.is_duplicate or read.is_qcfail or read.is_unmapped:
            continue
        if isSoftClipped(read.cigar):
            continue
        if read.is_paired:
            if read.mate_is_unmapped:
                continue
            if read.rnext != read.tid:
                continue
            if read.is_read1:
                if read.isize == 0:
                    continue
                rstart = min(read.pos, read.pnext) + \
                    1
                rend = rstart + abs(read.isize) - 1
                rmid = (int(rend) + int(rstart)) / 2
                rlength = rend - rstart + 1
                if minLength <= rlength <= maxLength:
                    total += 1
                    if rlength in result:
                        result[rlength] += 1
                    else:
                        result[rlength] = 1
    for k in result.keys():
        vaf = format(float(result[k]) / float(total), '.10f')
        print '\t'.join([sample, sample_type, str(k), str(vaf)])
    return


def main():
    filename = sys.argv[1]  # BAM file path
    sample = sys.argv[2]  # Sample id
    sample_type = sys.argv[3]  # Sample type
    vplots = calculatePerBase(filename, sample, sample_type)


if __name__ == "__main__":
    main()
