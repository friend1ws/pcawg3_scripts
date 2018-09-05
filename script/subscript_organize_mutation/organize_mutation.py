#! /usr/bin/env python

import sys, gzip

snv_file = sys.argv[1]
indel_file = sys.argv[2]
summary_file = sys.argv[3]
output_dir = sys.argv[4]


aliquote_id2donor_id = {}

with open(summary_file, 'r') as hin:

    header = hin.readline().rstrip('\n').split('\t')
    header2ind = {}
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.strip('\n').split('\t')
        for aliquot_id in F[header2ind["tumor_wgs_aliquot_id"]].split(','):
            aliquote_id2donor_id[aliquot_id] = F[header2ind["icgc_donor_id"]]


donor_id2file_handle = {} 
with gzip.open(snv_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[1] not in aliquote_id2donor_id: continue
        donor_id = aliquote_id2donor_id[F[1]]
        if donor_id not in donor_id2file_handle:
            donor_id2file_handle[donor_id] = open(output_dir + "/" + donor_id + ".mutation.txt", 'w')

        print >> donor_id2file_handle[donor_id], '\t'.join(F[5:10])


with open(indel_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[1] not in aliquote_id2donor_id: continue
        donor_id = aliquote_id2donor_id[F[1]]
        if donor_id not in donor_id2file_handle:
            donor_id2file_handle[donor_id] = open(output_dir + "/" + donor_id + ".mutation.txt", 'w')

        print >> donor_id2file_handle[donor_id], '\t'.join(F[5:10])


for donor_id in donor_id2file_handle:
    donor_id2file_handle[donor_id].close()





