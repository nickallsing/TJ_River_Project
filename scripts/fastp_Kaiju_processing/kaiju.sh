#!/bin/bash
max_cores=30
R1=$(python R1_files.py)
R2=$(python R2_files.py)
OUT=$(python kaiju_out.py)
cd ../../data/fastp_Kaiju_processing/output/fastp_out
kaiju-multi -z $max_cores -t ../../input/nodes.dmp -f ../../input/kaiju_db_nr_euk.fmi -i $R1 -j $R2 -o $OUT