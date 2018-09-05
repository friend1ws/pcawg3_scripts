#! /usr/bin/env python

import sys, gzip

gene_bed_file = sys.argv[1]
cancer_gene_file = sys.argv[2]

gene2cginfo ={}
with open(cancer_gene_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        gene2cginfo[F[0]] = '\t'.join(F[1:])


cg2chr = {}
cg2start = {}
cg2end = {}
with gzip.open(gene_bed_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[3] in gene2cginfo:
            cg2chr[F[3]] = F[0]
            cg2start[F[3]] = F[1] if F[3] not in cg2start else str(min(int(F[1]), int(cg2start[F[3]])))
            cg2end[F[3]] = F[2] if F[3] not in cg2start else str(max(int(F[2]), int(cg2start[F[3]])))


for cg in cg2chr:
    print '\t'.join([cg2chr[cg], cg2start[cg], cg2end[cg], cg, gene2cginfo[cg]])





