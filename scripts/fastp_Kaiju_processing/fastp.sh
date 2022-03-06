#!/bin/bash
max_cores=32
cd ../../data/fastp_Kaiju_processing/input
ls *fastq.gz | parallel -j$max_cores --max-args=2 fastp -h {1.}.html -j {1.}.json -i {1} -I {2} --out1 {1.}.out.gz --out2 {2.}.out.gz
mv *.html ../output/html_out/
mv *.json ../output/html_out
mv *.out.gz ../output/fastp_out
