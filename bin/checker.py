import numpy as np

def FreqMatrix(filename, k):
    f = open(filename)
    f = f.readlines()
    m_count = np.zeros((len(f), 4**k))
    for j in range(len(f)):
        l = f[j].split(" ")
        for i in range(4**k):
            index, value = l[i].split(':')
            m_count[j][int(index)] = int(value)
    return m_count

for k in range(1, 7):

    matheus = FreqMatrix(f"./cpu_SRR11510555_output_k{k}_t1.txt", k)

    fabricio = FreqMatrix(f"./gpu_SRR115510555_1_teste{k}.txt", k)

# print(fabricio)
# print(matheus)
    diff = fabricio-matheus
    print(len(diff.nonzero()[0]) == 0)
    # print(len(diff.nonzero()[0]))
