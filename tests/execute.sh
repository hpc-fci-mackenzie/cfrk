#!/usr/bin/env bash

for arg in "$@"; do
  index=$(echo "$arg" | cut -f1 -d=)
  val=$(echo "$arg" | cut -f2 -d=)
  case $index in
    script) S=$val;;
    output-file) F=$val;;
    K) K=$val;;
    threads) T=$val;;
  esac
done

rm output_files/"$F"_*.txt
rm output_files/"$F".csv
touch output_files/"$F".csv
echo K,THREAD,READ_TIME,CONSTRUCTION_TIME,PROCESSING_TIME,REMAIN_PROCESSING_TIME,WRITING_TIME,EXTERNAL_TIME,MEMORY_USED > output_files/"$F".csv
while true; do
  IT=$T
  ../bin/"$S" input_files/seq1.fasta output_files/"$F".txt "$K" "$IT" > /dev/null 2>&1
  while true; do
      /usr/bin/time --format=",%e,%M" -o output_files/"$F".csv -a ../bin/"$S" input_files/seq1.fasta output_files/"$F"_"$K".txt $K $IT >> output_files/"$F".csv
      IT=$(( IT - 1 ))
      if (( IT < 1 )); then
        break
      fi
    done
  K=$(( K - 1 ))
  if (( K < 1 )); then
    break
  fi
done
