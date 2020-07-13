import numpy as np


def freq_matrix(filename, k):
    f = open(filename)
    f = f.readlines()
    m_count = np.zeros((len(f), 4**k))
    for j in range(len(f)):
        l = f[j].split(" ")
        for i in range(4**k):
            index, value = l[i].split(':')
            m_count[j][int(index)] = int(value)
    return m_count


def main(args):
    file_1, file_2, k = args
    result_1, result_2 = freq_matrix(file_1, int(k)), freq_matrix(file_2, int(k))
    diff = result_1 - result_2
    return len(diff.nonzero()[0]) == 0


if __name__ == '__main__':
    import sys
    print(main(sys.argv[1:]))
