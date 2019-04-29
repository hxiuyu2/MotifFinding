import numpy as np


SEQ_LEN = 500
ALPHA = ['A', 'T', 'C', 'G']


# EM approach: MEME algorithm
def meme(sequences, motif_len):
    # loop to find a relatively good start point
    max_score = 0
    best_p = []
    walk_list = range(0, SEQ_LEN - motif_len)
    for l in walk_list:
        init_p = np.zeros(shape=(4, motif_len + 1)) + 0.17
        init_p[0][0] = 0.25
        init_p[1][0] = 0.25
        init_p[2][0] = 0.25
        init_p[3][0] = 0.25

        for j in range(motif_len):
            character = ALPHA.index(sequences[0][j + l])
            init_p[character][j + 1] = 0.5
        z_t = update_z(sequences, init_p, motif_len)
        p_matrix = m_step(sequences, z_t, motif_len)

        # calculate score
        cur_score = 0
        for i in range(motif_len):
            for j in range(4):
                cur_score += p_matrix[j][i + 1] * np.log2(p_matrix[j][i + 1] / p_matrix[j][0])

        if cur_score > max_score:
            max_score = cur_score
            best_p = p_matrix

    # assign the best start point
    p_matrix = best_p
    prev_ind = [0] * (motif_len + 1)

    # loop until converge
    while True:
        z_t = e_step(sequences, motif_len, p_matrix)
        p_next = m_step(sequences, z_t, motif_len)

        # if change in p < e or sites do not change
        diff = np.sum(abs(np.array(p_matrix[:][1:]) - np.array(p_next[:][1:])))
        ind = np.argmax(p_next, axis=0)
        diff_ind = np.sum(abs(np.array(ind[:][1:]) - np.array(prev_ind[:][1:])))
        if diff < 0.1 or diff_ind == 0:
            break
        else:
            p_matrix = p_next
            prev_ind = ind
    z_t = e_step(sequences, motif_len, p_matrix)
    return z_t, p_next


def e_step(sequences, motif_len, p_matrix):
    z_t = np.zeros((len(sequences), SEQ_LEN-motif_len+1))
    for i in range(len(sequences)):
        for j in range(SEQ_LEN-motif_len+1):
            pr = get_prob(sequences[i], p_matrix, j, motif_len)
            z_t[i][j] = pr
    row_sums = np.sum(z_t, axis=1)
    z_t = np.divide(z_t, row_sums[:, np.newaxis])
    return z_t


def m_step(sequences, z_t, motif_len):
    p_t = np.zeros(shape=(4, motif_len+1))

    # calculate n(c, k), k>0
    n_ck = np.zeros(shape=(4, motif_len+1))
    for i in range(len(sequences)):
        for j in range(SEQ_LEN - motif_len + 1):
            n_ck[ALPHA.index(sequences[i][j])][0] += 1
            for k in range(motif_len):
                n_ck[ALPHA.index(sequences[i][j + k])][k + 1] += z_t[i][j]

    # n(c, k) when k = 0
    for i in range(4):
        n_ck[i][0] -= np.sum(n_ck[i][1:])

    # calculate p(c, k)
    for i in range(motif_len + 1):
        sum_n = np.sum([num[i] + 1 for num in n_ck])
        for j in range(4):
            p_t[j][i] = (n_ck[j][i] + 1) / sum_n
    return p_t


def get_prob(sequence, p_matrix, start, motif_len):
    result = 1
    for i in range(SEQ_LEN):
        if i < start or i >= start + motif_len:
            result *= p_matrix[ALPHA.index(sequence[i])][0]
        else:
            result *= p_matrix[ALPHA.index(sequence[i])][i - start + 1]
    return result


def update_z(sequences, p_matrix, motif_len):
    z_t = np.zeros((len(sequences), SEQ_LEN - motif_len + 1))
    for i in range(len(sequences)):
        for j in range(SEQ_LEN - motif_len + 1):
            result = 1
            for k in range(motif_len):
                result *= p_matrix[ALPHA.index(sequences[i][k+j])][k]
            z_t[i][j] = result
    row_sums = np.sum(z_t, axis=1)
    z_t = np.divide(z_t, row_sums[:, np.newaxis])
    return z_t


# test for benchmarks
import os
import time

# os.mkdir('MEME', 0o777)

for i in range(70):
    seq_file = open('benchmark/dataset{}/sequence.fa'.format(str(i)))
    seq_list = []
    for line in seq_file:
        if line.startswith('>'):
            pass
        else:
            seq_list.append(line.strip())
    ml_file = open('benchmark/dataset{}/motiflength.txt'.format(str(i)))
    ml = int(ml_file.readline().strip())

    # calculate runtime
    start = time.time()
    z, result = meme(seq_list, ml)
    end = time.time()

    os.mkdir('output/dataset{}'.format(str(i)), 0o777)

    # write sites to file
    pred_site = open('output/dataset{}/predictedsites.txt'.format(str(i)), 'w+')
    for z_i in z:
        pred_site.write(str(np.argmax(z_i)))
        pred_site.write('\n')
    pred_site.close()

    # write motif to file
    pred_motif = open('output/dataset{}/predictedmotif.txt'.format(str(i)), 'w+')
    pred_motif.write('>motif{}    {}\n'.format(str(i), ml))
    for j in range(1, ml+1):
        pred_motif.write('{} {} {} {}\n'.format(str(result[0][j]), str(result[1][j]), str(result[2][j]), str(result[3][j])))
    pred_motif.write('<')
    pred_motif.close()

    print('for dataset {}, runtime: {}sec'.format(str(i), str(end - start)))


