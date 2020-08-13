#!/usr/bin/env bash
#declare -i K
#declare -i T
K=13
rm execute_results.csv
touch execute_results.csv
echo K,THREAD,READ_TIME,CONSTRUCTION_TIME,PROCESSING_TIME,REMAIN_PROCESSING_TIME,WRITING_TIME,EXTERNAL_TIME,MEMORY_USED > execute_results.csv
while true; do
  T=4
  ./cfrk-cpu input_files/seq1.fasta ../results/teste.txt $K 4
  while true; do
      /usr/bin/time --format=",%e,%M" -o execute_results.csv -a ./cfrk-cpu input_files/seq1.fasta ../results/teste.txt $K $T >> execute_results.csv
      T=$(( T - 1 ))
      if (( T < 1 )); then
        break
      fi
    done
  K=$(( K - 1 ))
  if (( K < 1 )); then
    break
  fi
done
