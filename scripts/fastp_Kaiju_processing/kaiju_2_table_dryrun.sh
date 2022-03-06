#!/bin/bash
max_cores=32
cd ../../data/fastp_Kaiju_processing/output/kaiju_out
ls *.out | parallel --dryrun -j$max_cores kaiju2table -t ../../input/nodes.dmp -n ../../input/names.dmp -r species -l domain,superkingdom,phylum,class,order,family,genus,species -o {.}.tsv {}
