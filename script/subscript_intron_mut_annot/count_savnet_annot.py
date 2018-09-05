#! /usr/bin/env python

import sys, glob, os
from count_annot_file import count_annot_file

savnet_mut_annot_file = sys.argv[1]

deep_thres = 100

print '\t'.join(["Cancer_Type", "Sample_Name", "Deep_Mut_Count", "CG1_Count", "CG2_Count", "CG3_Count", "SINE_Count", "SINE_Sense_Count", "SINE_Antisense_Count"])
print "savnet\tsavnet\t" + count_annot_file(savnet_mut_annot_file, 9)
 
