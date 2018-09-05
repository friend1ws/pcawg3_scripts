#! /usr/bin/env python

import sys, os, glob, subprocess

input_file = sys.argv[1]
output_file = sys.argv[2]

# SJ_control = "../../170108/control/SJ/output/control_2_3.bed.gz"
SJ_control = "../data/control/SJ/output/control_2_3.bed.gz"

hout = open(output_file, 'w')

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        if int(F[10]) < 100: continue

        half_exonizaiton_junction = F[7]

        SJ_out_tab_file = "../data/SJ/" + F[5] + "/" + F[6] + ".SJ.out.tab"
        # SJ_out_tab_file = "/home/eva/rawdata/icgc_pan_rna/output/" + F[5] + "/star/" + F[6] + "/" + F[6] + ".SJ.out.tab"
        # if not os.path.exists(SJ_out_tab_file):
        #     SJ_out_tab_file = "/home/eva/rawdata/icgc_pan_rna/single/output/" + F[5] + "/star/" + F[6] + "/" + F[6] + ".SJ.out.tab"

        output_tmp_file = output_file + ".exonization.tmp.txt"
        open(output_tmp_file,  'w').close()

        subprocess.call(["junc_utils", "exonization_pair", half_exonizaiton_junction, SJ_out_tab_file, output_tmp_file, \
                         "--grc", "--control_file", SJ_control, "--read_num_thres", "1"])

        opposite_junction = "---\t---"
        with open(output_tmp_file, 'r') as hin2:
            for line2 in hin2:
                F2 = line2.rstrip('\n').split('\t')
                opposite_junction = F2[0] + ":" + F2[1] + "-" + F2[2] + "\t" + F2[6]

        print >> hout, '\t'.join(F) + '\t' + opposite_junction
        # print '\t'.join(F) + '\t' + opposite_junction

        subprocess.call(["rm", "-rf", output_tmp_file])

