# Automated Test for CFRK

# NOTA: Linux nega a permissão ao rodar esse script!
#       Executar este comando antes de rodar o bash:
#
#          sudo chmod u+x ./test.sh
#
#!/usr/bin/env bash

echo "Automated test for CFRK!"
cd ../bin

# ./cfrk seq1.fasta out.cfrk 2 12 8192

diff out.cfrk out-seq1.cfrk
ret=$?

if [[ $ret -eq 0 ]]; then
	echo "Successful test with seq1.fasta!"
else
	exit 1
fi

# ./cfrk seq2.fasta out.cfrk 2 12 8192

# diff out.cfrk out-seq2.cfrk
# ret=$?

# if [[ $ret -eq 0 ]]; then
# 	echo "Successful test with seq2.fasta!"
# else
# 	exit 1
# fi
