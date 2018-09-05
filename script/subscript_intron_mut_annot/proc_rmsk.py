#! /usr/bin/env python

import sys, gzip

rmsk_file = sys.argv[1]

with gzip.open(rmsk_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        F[5] = F[5].replace("chr", "")
        print '\t'.join([F[5], F[6], F[7], F[9], F[10], F[11], F[12]])


