#!/bin/bash

for i in "AA" "AC" "AG" "AT" "CA" "CC" "CG" "CT" "GA" "GC" "GG" "GT" "TA" "TC" "TG" "TT"
do
   echo ${i}
   grep -c ${i} ../seq/SRR2155482.fasta_0.fasta
done
