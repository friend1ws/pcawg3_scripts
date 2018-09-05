#! /usr/bin/env python

import sys

savnet_result_file = sys.argv[1]

key2exists = {}
with open(savnet_result_file, 'r') as hin:

    header = hin.readline().rstrip('\n').split('\t')
    header2ind = {}
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind["Mutation_Type"]] not in ["splicing acceptor creation", "splicing donor creation"]: continue

        key = F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Mutation_Key"]]
        if key in key2exists: continue
        key2exists[key] = 1

        tchr, tpos, tref, talt = F[header2ind["Mutation_Key"]].split(',')

        if len(tref) == 1 and len(talt) == 1:
            tstart, tend = tpos, tpos
        elif len(tref) > 1:
            tstart, tend = str(int(tpos) + 1), str(int(tpos) + len(tref))
            tref, talt = tref[1:], "-"
        elif len(talt) > 1:
            tstart, tend = tpos, tpos 
            tref, talt = "-", talt[1:]
 
        print '\t'.join([tchr, tstart, tend, tref, talt, F[header2ind["Cancer_Type"]], \
                         F[header2ind["Sample_Name"]], F[header2ind["Splicing_Key"]], \
                         F[header2ind["Mutation_Type"]].replace("splicing ", "").replace(" creation", "")])

