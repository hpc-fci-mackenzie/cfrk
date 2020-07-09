import subprocess

command = '/home/mferreira/repos/cfrk/bin/cfrk-cpu.out'
input_file = 'seq1.fasta'
output_file = 'cpu_seq1_output.txt'

with open('execute_results.csv', 'w') as f:
    f.write('K,THREAD,READ_TIME,CONSTRUCTION_TIME,PROCESSING_TIME,REMAIN_PROCESSING_TIME\n')
    out_str = ''
    for k in range(1, 300):
        for t in range(1, 5):
            # output_file = f'cpu_SRR11510555_output_k{k}_t{t}.txt'
            process = subprocess.Popen([command, input_file, output_file, str(k), str(t)],
                             stdout= subprocess.PIPE)
            out, err = process.communicate()
            if err: print(err)
            # f.write(out.decode("utf-8"))
            out_str += out.decode("utf-8")
    f.write(out_str)
    f.close()
print("Finished")
