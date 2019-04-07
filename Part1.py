import random
import numpy as np
from scipy.optimize import fsolve
import os
import pathlib


ICPC = [1, 1.5, 2]
ML = [6, 7]
SL = 500
SC = [5, 20]


def generate_seq(seq_count, seq_len):
    sequences = []
    vfunc = np.vectorize(seq_func)

    for i in range(seq_count):
        seq = np.random.uniform(0, 4, seq_len)
        seq = vfunc(seq)
        line = ''.join(seq)
        sequences.append(line)
    return sequences


def seq_func(x):
    E = ['A', 'T', 'C', 'G']
    x = int(x)
    return E[x]


def generate_motif(motif_len, icpc):
    E = ['A', 'T', 'C', 'G']
    if icpc == 2:
        result = []
        for i in range(motif_len):
            motif = [0, 0, 0, 0]
            one = np.random.choice(E, 1)
            motif[E.index(one)] = 1
            result.append(motif)
        return result

    motif = []
    for i in range(motif_len):
        w0 = random.randint(1, 5)
        w1 = random.randint(1, 5)
        w0 = w0 / 100
        w1 = w1 / 100
        rest = 1 - w0 - w1
        equ = w0 * np.log2(4 * w0) + w1 * np.log2(4 * w1) - icpc
        func = lambda x: x * np.log2(4 * x) + (rest - x) * np.log2(4 * (rest - x)) + equ
        w2 = fsolve(func, 1e-9)
        w3 = rest - w2[0]
        line = [w0, w1, w2[0], w3]
        np.random.shuffle(line)
        motif.append(line)
    return motif


def generate_site(motif, seq_count):
    sites = []
    E = ['A', 'T', 'C', 'G']
    for i in range(seq_count):
        site = []
        for row in motif:
            letter = np.random.choice(E, 1, p=row)
            site.append(letter[0])
        line = ''.join(site)
        sites.append(line)
    return sites


def plant(sites, sequences):
    result = []
    seq_len = len(sequences[0])
    mot_len = len(sites[0])
    for index, seq in enumerate(sequences):
        s = sites[index]
        location = random.randint(0, seq_len - mot_len)
        line = ''.join((seq[:location], s, seq[location + mot_len:]))
        result.append(line)
    return result


def write_file(sequences, sites, motif, motif_len, count):
    os.mkdir('benchmark/dataset'+ str(count))

    # write sequences to 'sequences.fa'
    file = open('benchmark/dataset' + str(count) + '/sequences.fa', 'w+')
    for row in sequences:
        file.write(row + '\n')
    file.close()

    # write sites to 'sites.txt'
    file = open('benchmark/dataset' + str(count) + '/sites.txt', 'w+')
    for row in sites:
        file.write(row + '\n')
    file.close()

    # write motif to 'motif.txt'
    file = open('benchmark/dataset' + str(count) + '/motif.txt', 'w+')
    header = 'MOTIF{0}    {1}\n'.format(str(count), str(motif_len))
    file.write(header)
    for row in motif:
        line = ' '.join(str(x) for x in row)
        file.write(line + '\n')
    file.close()

    # write motif_len to 'motiflength.txt'
    file = open('benchmark/dataset' + str(count) + '/motiflength.txt', 'w+')
    file.write(str(motif_len))
    file.close()

# # directory for all the benchmarks
# os.mkdir('benchmark', 0o777)

count = 0
for i in range(10):
    # generate combination (a)
    ml = 8
    sc = 10
    for icpc in ICPC:
        sequences = generate_seq(sc, SL)
        motif = generate_motif(ml, icpc)
        sites = generate_site(motif, sc)
        result = plant(sites, sequences)
        write_file(sequences, sites, motif, ml, count)
        print('finish dataset ', count)
        count += 1

    # generate combination (b)
    icpc = 2
    sc = 10
    for ml in ML:
        sequences = generate_seq(sc, SL)
        motif = generate_motif(ml, icpc)
        sites = generate_site(motif, sc)
        result = plant(sites, sequences)
        write_file(sequences, sites, motif, ml, count)
        print('finish dataset ', count)
        count += 1

    # generate combination (c)
    ml = 8
    icpc = 2
    for sc in SC:
        sequences = generate_seq(sc, SL)
        motif = generate_motif(ml, icpc)
        sites = generate_site(motif, sc)
        result = plant(sites, sequences)
        write_file(sequences, sites, motif, ml, count)
        print('finish dataset ', count)
        count += 1