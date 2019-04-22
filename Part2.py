import numpy as np
from copy import copy, deepcopy

LOOP_TIME = 300
SEQ_LEN = 500
ALPHA = ['A', 'T', 'C', 'G']


# def gibbs_sampler(sequences, motif_len):
#     seq_cnt = len(sequences)
#     sites = np.random.randint(0, SEQ_LEN - motif_len, seq_cnt).tolist()
#     for i in range(LOOP_TIME):
#         for j in range(seq_cnt):
#             temp_seq = sequences[:j]+sequences[j+1:]
#             temp_st = sites[:j]+sites[j+1:]
#             temp_motif = get_motif(temp_seq, temp_st, motif_len)
#             cur_prob = get_p(temp_motif, sequences[j], motif_len)
#             sites[j] = np.random.choice(SEQ_LEN - motif_len, 1, p=cur_prob)[0]
#     return sites, get_motif(sequences, sites, motif_len)
#
#
# def get_motif(sequences, sites, motif_len):
#     seq_cnt = len(sequences)
#     motif= np.zeros([4, motif_len]) + 0.1/seq_cnt
#     motif = motif.tolist()
#     for i in range(seq_cnt):
#         plant = sequences[i][sites[i]:sites[i]+motif_len]
#         for j in range(motif_len):
#             pos = ALPHA.index(plant[j])
#             motif[pos][j] += 1/seq_cnt
#     return motif
#
#
# def get_p(motif, sequence, motif_len):
#     p_a = np.prod([motif[ALPHA.index(a)][i] for i, a in enumerate(sequence[0:motif_len])])
#     prob = [p_a]
#     for i in range(SEQ_LEN - motif_len - 1):
#         next_prob = prob[-1]* motif[ALPHA.index(sequence[i+motif_len])][-1] / motif[ALPHA.index(sequence[i])][0]
#         prob.append(next_prob)
#     prob = prob / np.sum(prob)
#     return prob
#
#
# seq = []
# sequence_txt = open('benchmark/dataset5/sequences.fa')
# for line in sequence_txt:
#     seq.append(line.strip())
# sequence_txt.close()
#
# sts = []
# sites_txt = open('benchmark/dataset5/sites.txt')
# for line in sites_txt:
#     sts.append(int(line.strip()))
#
# ml = 8
# print(sts)
# ALPHA = ['A','T','C','G']
# for i in range(5):
#     print(seq[i][sts[i]:sts[i]+ml])
#

# EM approach
# MEME algorithm
def meme(sequences, motif_len):
    # loop to find start points
    max_score = 0
    best_p = []
    for k in range(len(sequences)):
        rand_list = np.random.choice(range(SEQ_LEN-motif_len), 20)
        for l in rand_list:
            init_p = np.zeros(shape=(4, motif_len+1)) + 0.17
            init_p[0][0] = 0.25
            init_p[1][0] = 0.25
            init_p[2][0] = 0.25
            init_p[3][0] = 0.25
            for j in range(motif_len):
                character = ALPHA.index(sequences[k][j+l])
                init_p[character][j+1] = 0.5
            z_t = e_step(sequences, motif_len, init_p)
            p_matrix = m_step(sequences, z_t, motif_len)

            # calculate score
            cur_score = 0
            for i in range(motif_len):
                for j in range(4):
                    cur_score += p_matrix[j][i+1]*np.log2(p_matrix[j][i+1]/p_matrix[j][0])

            if cur_score > max_score:
                max_score = cur_score
                best_p = deepcopy(p_matrix)

    # assign p w/ max likelihood as start point
    p_matrix = deepcopy(best_p)

    # loop until converge
    while True:
        z_t = e_step(sequences, motif_len, p_matrix)
        p_next = m_step(sequences, z_t, motif_len)

        # if change in p < e
        diff = np.sum(abs(np.array(p_matrix) - np.array(p_next)))
        print(diff)
        if diff < 0.5:
            break
        else:
            p_matrix = deepcopy(p_next)
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
        for j in range(SEQ_LEN-motif_len):
            p_t[ALPHA.index(sequences[i][j])][0] += 1/SEQ_LEN/len(sequences)
            for k in range(motif_len):
                n_ck[ALPHA.index(sequences[i][j+k-1])][k+1] += z_t[i][j]

    # calculate p(c, k)
    for i in range(1, motif_len+1):
        for j in range(4):
            sum_n = np.sum([num[i] for num in n_ck])
            p_t[j][i] = n_ck[j][i] / sum_n
    return p_t


def get_prob(sequence, p_matrix, start, motif_len):
    result = 1
    for i in range(SEQ_LEN):
        if i < start or i > start + motif_len:
            result *= p_matrix[ALPHA.index(sequence[i])][0]
        else:
            result *= p_matrix[ALPHA.index(sequence[i])][i-start]
    return result


seq_file = open('benchmark/dataset5/sequences.fa')
seq_list = []
for line in seq_file:
    seq_list.append(line.strip())
ml = 8
# GAAATTGC
# 22
# 6
# 158
# 437
# 115
import time
start = time.time()
z, result = meme(seq_list, ml)
end = time.time()
for z_i in z:
    print(np.argmax(z_i))
indexs = np.argmax(result, axis=0)
print([ALPHA[i] for i in indexs[1:]])
print('runtime: ',end - start)