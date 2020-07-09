# Automated Test for CFRK

# NOTA: Linux nega a permissão ao rodar esse script!
#       Executar este comando antes de rodar o bash:
#
#          sudo chmod u+x ./test.sh
#

#!/usr/bin/env bash -e

echo "Automated test for CFRK!"

../bin/./cfrk ../sample/seq1.fasta out.cfrk 2 12 8192

diff out.cfrk out-seq1.cfrk

../bin/./cfrk ../sample/seq2.fasta out.cfrk 2 12 8192

diff out.cfrk out-seq2.cfrk
