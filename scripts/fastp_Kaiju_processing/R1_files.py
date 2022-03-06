#Creates a list for R1 reads
import os
import re
qc_R1_files=[]
files = [f for f in os.listdir('../../data/fastp_Kaiju_processing/output/fastp_out')]
for f in files:
    if "R1" in f:
        qc_R1_files.append(f)
    else:
        pass
qc_R1_files.sort()
R1_files = ",".join(qc_R1_files)
print(R1_files)
