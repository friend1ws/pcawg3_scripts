#! /usr/bin/env python

import sys, re
import pysam
from Bio import pairwise2

input_file = sys.argv[1]

rmsk_tb = pysam.TabixFile("../db/rmsk.bed.gz")

reference = "/home/w3varann/database/GRCh37/GRCh37.fa"

AluJ0 = """ggccgggcgcggtggctcacgcctgtaatcccagcactttgggaggccgaggcgggaggattgcttgagcc
caggagttcgagaccagcctgggcaacatagcgagaccccgtctctacaaaaaatacaaaaattagccggg
cgtggtggcgcgcgcctgtagtcccagctactcgggaggctgaggcaggaggatcgcttgagcccaggagt
tcgaggctgcagtgagctatgatcgcgccactgcactccagcctgggcgacagagcgagaccctgtctca"""

AluJ0 = AluJ0.replace('\n', '').upper()


def get_seq(reference, chr, start, end):

    seq = ""
    for item in pysam.faidx(reference, chr + ":" + str(start) + "-" + str(end)):
        # if item[0] == ">": continue
        seq = seq + item.rstrip('\n')
    seq = seq.replace('>', '')
    seq = seq.replace(chr + ":" + str(start) + "-" + str(end), '')

    return seq

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))



with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[12] == "---": continue

        genes = F[9].split("|")
        rmsks = F[12].split("|")

        # check SINE sense or anti-sense
        gene_dirs = []
        for elm in genes:
            elms = elm.split(',')
            gene_dirs.append(elms[1])
        gene_dir = ','.join(list(set(gene_dirs)))

        if gene_dir not in ["+", "-", "+,-", "-,+"]:
            print >> sys.stderr, "Something is wrong!"
            sys.exit(1)

        sine_flag, sine_s_flag, sine_a_flag = False, False, False
        if rmsks[0] != "---":
            for elm in rmsks:
                elms = elm.split(',')
                if elms[2] == "SINE":
                    sine_flag = True
                    if gene_dir in ["+", "-"]:
                        if gene_dir == elms[0]:
                            sine_s_flag = True
                        else:
                            sine_a_flag = True

        if gene_dir not in ["+", "-"]: continue

        if sine_a_flag == False: continue


        # get new donor and acceptor sites
        new_donor_pos, new_acceptor_pos = "---", "---"
        junc_match = re.match(r'([^ \t\n\r\f\v,]+)\:(\d+)\-(\d+)', F[7])
        if junc_match is None:
            print "Something is wrong"
            sys.exit(1)
        junc_chr, junc_start, junc_end = junc_match.group(1), int(junc_match.group(2)), int(junc_match.group(3))

        if F[8] == "donor":
            new_donor_pos = junc_start - 1 if gene_dir == "+" else junc_end + 1
        else:
            new_acceptor_pos = junc_end + 1 if gene_dir == "+" else junc_start - 1

        if F[13] != "---":
            junc_match = re.match(r'([^ \t\n\r\f\v,]+)\:(\d+)\-(\d+)', F[13])
            if junc_match is None:
                print "Something is wrong"
                sys.exit(1)
            junc_chr, junc_start, junc_end = junc_match.group(1), int(junc_match.group(2)), int(junc_match.group(3))

            if F[8] == "donor":
                new_acceptor_pos = junc_end + 1 if gene_dir == "+" else junc_start - 1
            else:
                new_donor_pos = junc_start - 1 if gene_dir == "+" else junc_end + 1



        # get the records for control junction data for the current position
        tabixErrorFlag = 0
        try:
            records = rmsk_tb.fetch(F[0], int(F[1]) - 10, int(F[1]) + 10)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1

        rep_start, rep_end, rep_strand, rep_name, rep_class, rep_family = "---", "---", "---", "---", "---", "---"
        new_donor_pos_rel, new_acceptor_pos_rel = "---", "---"
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                if record[5] != "SINE": continue
                rep_chr, rep_start, rep_end, rep_strand, rep_name, rep_class, rep_family = record[0], int(record[1]) + 1, int(record[2]), record[3], record[4], record[5], record[6]
        
                if rep_strand == '+':

                    if new_donor_pos != "---":
                        if new_donor_pos < rep_start:
                            new_donor_pos_rel = - float("inf")
                        elif new_donor_pos > rep_end:   
                            new_donor_pos_rel = float("inf")
                        else: 
                            new_donor_pos_rel = new_donor_pos - rep_start + 1

                    if new_acceptor_pos != "---":
                        if new_acceptor_pos < rep_start:
                            new_acceptor_pos_rel = - float("inf")
                        elif new_acceptor_pos > rep_end:
                            new_acceptor_pos_rel = float("inf")
                        else:
                            new_acceptor_pos_rel = new_acceptor_pos - rep_start + 1

                else:

                    if new_donor_pos != "---":
                        if new_donor_pos < rep_start:
                            new_donor_pos_rel = float("inf")
                        elif new_donor_pos > rep_end:   
                            new_donor_pos_rel = - float("inf")
                        else: 
                            new_donor_pos_rel = rep_end - new_donor_pos + 1
 
                    if new_acceptor_pos != "---":
                        if new_acceptor_pos < rep_start:
                            new_acceptor_pos_rel = float("inf")
                        elif new_acceptor_pos > rep_end:
                            new_acceptor_pos_rel = - float("inf")
                        else:
                            new_acceptor_pos_rel = rep_end - new_acceptor_pos + 1


                rep_seq = get_seq(reference, rep_chr, rep_start, rep_end)
                if rep_strand == '-': rep_seq = reverse_complement(rep_seq)

                alignments = pairwise2.align.localxs(rep_seq, AluJ0, -1, -1)
 
 
                mod_new_donor_pos_rel, mod_new_acceptor_pos_rel = new_donor_pos_rel, new_acceptor_pos_rel
                query_seq_pos = 0
                reference_seq_pos = 0
                for i in range(len(alignments[0][0])):
                    if alignments[0][0][i] != "-": query_seq_pos = query_seq_pos + 1
                    if alignments[0][1][i] != "-": reference_seq_pos = reference_seq_pos + 1

                    if query_seq_pos == new_donor_pos_rel: mod_new_donor_pos_rel = reference_seq_pos
                    if query_seq_pos == new_acceptor_pos_rel: mod_new_acceptor_pos_rel = reference_seq_pos

                if F[1] == "103246947":
                    print alignments[0][0]
                    print alignments[0][1]

                if new_donor_pos != "---":
                    exon_type = "donor-acceptor" if new_acceptor_pos != "---" else "donor-only"
                elif new_acceptor_pos != "---":
                    exon_type = "acceptor-only"


                print ','.join(F[:5]) + '\t' + F[5] + '\t' + F[6] + '\t' + '\t'.join(F[8:13]) + '\t' + \
                        rep_chr + ':' + str(rep_start) + '-' + str(rep_end) + '\t' + str(new_donor_pos) + '\t' + str(new_acceptor_pos) + '\t' + str(new_donor_pos_rel) + '\t' + str(new_acceptor_pos_rel) + '\t' + \
                        str(mod_new_donor_pos_rel) + '\t' + str(mod_new_acceptor_pos_rel) + '\t' + str(exon_type)

                # print '\t'.join(F) + '\t' + str(mod_new_donor_pos_rel) + '\t' + str(mod_new_acceptor_pos_rel)



# 7   56094771    56094771    G   A   BRCA-US DO4635  7:56088925-56094772 donor   PSPH,-  4850    --- +,AluSp,SINE,Alu    7:56094892-56099621 2
# 9   2821298 2821298 G   C   BRCA-US DO1287  9:2820099-2821298   donor   PUM3,-  1201    --- +,AluY,SINE,Alu 9:2821437-2823780   3
# 16  29705436    29705436    C   A   BRCA-US DO2593  16:29705434-29705984    donor   QPRT,+  548 --- -,AluSx1,SINE,Alu   16:29690532-29705330    4



