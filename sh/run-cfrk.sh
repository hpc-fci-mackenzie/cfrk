# NOTA: Ubuntu nega a permissão ao rodar esse script!
#       Executar este comando antes de rodar o bash:
#          sudo chmod u+x ./run-cfrk.sh

#!/usr/bin/env bash
echo "Rodando CFRK CPU version!"
cd ../bin
# arg1 = dataset name
# arg2 = outfile name
# arg3 = k
# arg4 = nt
./cfrk-cpu.out teste.fasta result.cfrk 2 12 1
# ./crfk-cpu.out