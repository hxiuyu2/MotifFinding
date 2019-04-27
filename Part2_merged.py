import numpy as np
import random
import copy
import os
import time

SEQ_LEN = 500
LOOP_NUM = 100000
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


# Gibbs Sampler without F
def generate_prob(pro_matrix,k_mer):
    dic = {'A':0,'T':1,'C':2,'G':3}
    p = 1
    for i in range(len(k_mer)):
        p = p * pro_matrix[dic[k_mer[i]]][i]
    return p


def background_prob(fa_Seq,sites,motif_len):
    dic = {'A':0,'T':1,'C':2,'G':3}
    background_c = np.zeros(4)
    for i in range(len(fa_Seq)):
        for j in range(len(fa_Seq[i])):
            if j not in list(range(sites[i],sites[i]+motif_len)):
                background_c[dic[fa_Seq[i][j]]] +=1
    background_p = background_c/sum(background_c)
    return background_p


def gibbs_without_F_until_converge(fa_Seq, motif_len):
    fa_Num = len(fa_Seq)
    sites = [random.randint(0, (len(fa_Seq[0]) - motif_len)) for i in range(fa_Num)]
    k_mers = []
    for i in range(fa_Num):
        k_mers.append(list(fa_Seq[i][sites[i]:sites[i] + motif_len]))
    f = [float('inf') for i in range(motif_len)]
    A = []
    B = []
    # start iteration
    hide_index = 0
    equal_count = 0
    step_count = 0
    while equal_count < 3 and step_count<LOOP_NUM:
        step_count += 1
        hide_same = True
        while hide_same:
            new_hide_index = random.randint(0, fa_Num - 1)  # randonly hide a row
            if new_hide_index != hide_index:
                hide_same = False
                hide_index = new_hide_index
        temp_motif_matrix = k_mers[:hide_index] + k_mers[hide_index + 1:]  # generate temperate motif matrix
        psudo_count_matrix = [[1 for i in range(motif_len)] for j in range(4)]  # generate psudo_count_matrix
        for i in range(motif_len):
            for j in range(fa_Num - 1):
                if temp_motif_matrix[j][i] == 'A':
                    psudo_count_matrix[0][i] += 1
                elif temp_motif_matrix[j][i] == 'T':
                    psudo_count_matrix[1][i] += 1
                elif temp_motif_matrix[j][i] == 'C':
                    psudo_count_matrix[2][i] += 1
                elif temp_motif_matrix[j][i] == 'G':
                    psudo_count_matrix[3][i] += 1
        psudo_count_matrix = np.array(psudo_count_matrix)  # pro_matrixpsudo_count_matrix
        pro_matrix = psudo_count_matrix / (fa_Num - 1 + 4)  # generate pro_matrix
        prob = []
        hiden_string = fa_Seq[hide_index]
        for i in range(len(fa_Seq[0]) - motif_len + 1):
            temp_k_mer = hiden_string[i:i + motif_len]
            prob.append(generate_prob(pro_matrix, temp_k_mer))
        # site temperature update
        updated_site = np.random.choice(len(fa_Seq[0]) - motif_len + 1, 1, p=prob / sum(prob))[0]
        temp_sites = copy.deepcopy(sites)
        temp_sites[hide_index] = updated_site  # update site(temp_site)
        # k_mers temperature update
        temp_k_mers = copy.deepcopy(k_mers)
        temp_k_mers[hide_index] = list(
            fa_Seq[hide_index][sites[hide_index]:sites[hide_index] + motif_len])  # 新的k_mers(temp_k_mers)
        if temp_sites == sites:
            equal_count += 1
        else:
            equal_count = 0
        sites = copy.deepcopy(temp_sites)
        k_mers = copy.deepcopy(temp_k_mers)
        if equal_count ==3 or step_count == (LOOP_NUM-1):
            A.append(sites)
            B.append(k_mers)
    B = B[-1]
    M = [[0 for i in range(motif_len)] for j in range(4)]  # 生成psudo_count_matrix
    for i in range(motif_len):
        for j in range(fa_Num):
            if B[j][i] == 'A':
                M[0][i] += 1
            elif B[j][i] == 'T':
                M[1][i] += 1
            elif B[j][i] == 'C':
                M[2][i] += 1
            elif B[j][i] == 'G':
                M[3][i] += 1
    M = np.array(M)  # pro_matrixpsudo_count_matrix
    return (M / fa_Num).T


# test for benchmarks
for i in range(70):
    seq_file = open('benchmark/dataset{}/sequences.fa'.format(str(i)))
    seq_list = []
    for line in seq_file:
        seq_list.append(line.strip())
    ml_file = open('benchmark/dataset{}/motiflength.txt'.format(str(i)))
    ml = int(ml_file.readline().strip())

    # calculate runtime
    start = time.time()
    z, _ = meme(seq_list, ml)
    result = gibbs_without_F_until_converge(seq_list, ml)
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
    pred_motif.write('motif{}    {}\n'.format(str(i), ml))
    for j in range(ml):
        pred_motif.write('{} {} {} {}\n'.format(str(result[j][0]), str(result[j][1]), str(result[j][2]), str(result[j][3])))
    pred_motif.close()

    print('for dataset {}, runtime: {}sec'.format(str(i), str(end - start)))