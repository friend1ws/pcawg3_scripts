#! /usr/bin/env python

deep_thres = 100

def count_annot_file(input_file, annot_start_column):

    mut_count = 0
    cg1_count = 0
    cg2_count = 0
    cg3_count = 0
    sine_count = 0
    sine_s_count = 0
    sine_a_count = 0

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            int_dist = int(F[annot_start_column + 1])
            if int_dist < deep_thres: continue
            mut_count = mut_count + 1

            genes = F[annot_start_column + 0].split("|")
            cgs = F[annot_start_column + 2].split("|")
            rmsks = F[annot_start_column + 3].split("|")

            gene_dirs = [] 
            for elm in genes:
                elms = elm.split(',')
                gene_dirs.append(elms[1])
            gene_dir = ','.join(list(set(gene_dirs)))
            
            if gene_dir not in ["+", "-", "+,-", "-,+"]:
                print >> sys.stderr, "Something is wrong!"
                sys.exit(1)


            cg1_flag, cg2_flag, cg3_flag = False, False, False
            if cgs[0] != "---":
                for elm in cgs:
                    elms = elm.split(',')
                    if elms[1] != "---": cg1_flag = True
                    if elms[2] != "---": cg2_flag = True
                    if elms[4] != "---": cg3_flag = True

            if cg1_flag == True: cg1_count = cg1_count + 1
            if cg2_flag == True: cg2_count = cg2_count + 1
            if cg3_flag == True: cg3_count = cg3_count + 1


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

            if sine_flag == True: sine_count = sine_count + 1
            if sine_s_flag == True: sine_s_count = sine_s_count + 1
            if sine_a_flag == True: sine_a_count = sine_a_count + 1


    return '\t'.join([str(x) for x in [mut_count, cg1_count, cg2_count, cg3_count, sine_count, sine_s_count, sine_a_count]])

