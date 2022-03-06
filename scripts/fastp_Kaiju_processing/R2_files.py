#Creates a list for R2 reads
import os
import re
qc_R2_files=[]
files = [f for f in os.listdir('../../data/fastp_Kaiju_processing/output/fastp_out')]
for f in files:
    if "R2" in f:
        qc_R2_files.append(f)
    else:
        pass
qc_R2_files.sort()
R2_files = ",".join(qc_R2_files)
print(R2_files)
