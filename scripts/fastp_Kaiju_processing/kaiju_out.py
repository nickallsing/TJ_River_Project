#Creates a list for R1 reads
import os
import re
kaiju_out_files=[]
files = [f for f in os.listdir('../../data/fastp_Kaiju_processing/output/fastp_out')]
for f in files:
    if "R1" in f:
        kaiju_out_files.append('../kaiju_out/'+f[:-3])
    else:
        pass
kaiju_out_files.sort()
kaiju_out = ",".join(kaiju_out_files)
print(kaiju_out)
