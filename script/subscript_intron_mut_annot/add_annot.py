#! /usr/bin/env python

import sys

input_file = sys.argv[1]
gene_file = sys.argv[2]
intron_file = sys.argv[3]
cg_file = sys.argv[4]
rmsk_file = sys.argv[5]


mut2gene = {}
with open(gene_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mut = '\t'.join(F[:5])
        if mut not in mut2gene:
            mut2gene[mut] = [F[8] + ',' + F[10]]
        else:
            mut2gene[mut].append(F[8] + ',' + F[10])

mut2int_dist = {}
with open(intron_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mut = '\t'.join(F[:5])
        int_dist = min(int(F[1]) - int(F[6]) + 1, int(F[7]) - int(F[2]))

        if mut not in mut2int_dist:
            mut2int_dist[mut] = int_dist
        elif int_dist < mut2int_dist[mut]:
            mut2int_dist[mut] = int_dist




mut2cg = {}
with open(cg_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mut = '\t'.join(F[:5])
        mut2cg[mut] = ','.join(F[8:]) if mut not in mut2cg else mut2cg[mut] + '|' + ','.join(F[8:])


mut2rmsk = {}
with open(rmsk_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mut = '\t'.join(F[:5])
        mut2rmsk[mut] = ','.join(F[8:]) if mut not in mut2rmsk else mut2rmsk[mut] + '|' + ','.join(F[8:])


with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mut = '\t'.join(F[:5])
        gene = '|'.join(list(set(mut2gene[mut]))) if mut in mut2gene else "---"
        cg = mut2cg[mut] if mut in mut2cg else "---"
        rmsk = mut2rmsk[mut] if mut in mut2rmsk else "---"
        int_dist = str(mut2int_dist[mut]) if mut in mut2int_dist else "0"

        print '\t'.join(F) + '\t' + gene + '\t' + int_dist + '\t' +  cg + '\t' + rmsk

