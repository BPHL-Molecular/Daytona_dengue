#!/usr/bin/env python

import glob
import re
filelist=glob.glob("./output/assembly_qc_pass/*.consensus.fa")

for afile in filelist:
   paths = afile.strip().split("/")
   file2=open(afile)
   content2 = file2.readlines()
   line1=re.sub('>SER\d+_','>',content2[0])
   with open("./output/assembly_qc_pass/SER_del/"+paths[-1], "w") as f:
      f.write(line1)
      f.write(content2[1])