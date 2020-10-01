#!/usr/bin/env bash
S=$1
F=$2
K=$3
T=$4
rm output_files/$F.txt
rm output_files/$F.csv
touch $F
echo K,THREAD,READ_TIME,CONSTRUCTION_TIME,PROCESSING_TIME,REMAIN_PROCESSING_TIME,WRITING_TIME,EXTERNAL_TIME,MEMORY_USED > $F.csv
while true; do
  ../bin/$S input_files/seq1.fasta output_files/$F.txt $K $T
  while true; do
      /usr/bin/time --format=",%e,%M" -o $F -a ../bin/$S input_files/seq1.fasta output_files/$F.txt $K $T >> $F
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
