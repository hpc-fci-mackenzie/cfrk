import subprocess


def warm_up(command, input_file, output_file, k, max_threads):
    process = subprocess.Popen([command, input_file, output_file, str(k), str(max_threads)],
                               stdout=subprocess.PIPE)
    subprocess.Popen(['rm', output_file])


def run_tests(command, input_file, output_file, k, max_threads):
    out_str = ''
    for t in range(1, (max_threads + 1)):
        process = subprocess.Popen([command, input_file, output_file + f'_{t}.txt', str(k), str(t)],
                                   stdout=subprocess.PIPE)
        out, err = process.communicate()
        if err:
            print(err)
        out_str += out.decode("utf-8")
    return out_str


def main(args):
    with open('execute_results.csv', 'w') as f:
        f.write('K,THREAD,READ_TIME,CONSTRUCTION_TIME,PROCESSING_TIME,REMAIN_PROCESSING_TIME\n')
        command, input_file, max_k, max_threads = args
        for k in range(1, int(max_k) + 1):
            warm_up(command, input_file, f'warm_up_result.txt', int(k), int(max_threads))
            f.write(run_tests(command, input_file, f'../results/{k}k', int(k), int(max_threads)))
        f.close()


if __name__ == '__main__':
    import sys
    print(main(sys.argv[1:]))