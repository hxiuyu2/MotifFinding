import random
import numpy as np
import math
import copy
import os
import timeit

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

def gibbs_without_F_unitil_converge(file_number):
    fa_path = "D:/IntroToDataMining-CS412/final project/benchmark/dataset"+str(file_number)+"/sequences.fa"
    fa_in = open(fa_path,"r")
    fa_Seq = []
    fa_Num = 0
    for line in fa_in.readlines():
        line = line.rstrip()
        fa_Num = fa_Num + 1
        fa_Seq.append(line)
    # read motiflen
    motif_len_path = "D:/IntroToDataMining-CS412/final project/benchmark/dataset"+str(file_number)+"/motiflength.txt"
    motiflen_in = open(motif_len_path,"r")
    motif_len = int(motiflen_in.readline())
    # set initial value
    ## step1 random pick initial site
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
    while equal_count < 3 and step_count<200000:
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
        if equal_count ==3 or step_count == (200000-1):
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
    # write into file:
    os.mkdir('Gibbs_No_F/dataset' + str(file_number))
    ## write sites to 'sites.txt'
    file = open('Gibbs_No_F/dataset' + str(file_number) + '/sites.txt', 'w+')
    for row in A[-1]:
        file.write(str(row) + '\n')
    file.close()
    ## write motif to 'motif.txt'
    file = open('Gibbs_No_F/dataset' + str(file_number) + '/motif.txt', 'w+')
    header = 'MOTIF{0}    {1}\n'.format(str(file_number), str(motif_len))
    file.write(header)
    for row in (M/fa_Num).T:
        line = ' '.join(str(x) for x in row)
        file.write(line + '\n')
    file.close()

def gibbs_without_F_not_to_converge(file_number):
    fa_path = "D:/IntroToDataMining-CS412/final project/benchmark/dataset"+str(file_number)+"/sequences.fa"
    fa_in = open(fa_path,"r")
    fa_Seq = []
    fa_Num = 0
    for line in fa_in.readlines():
        line = line.rstrip()
        fa_Num = fa_Num + 1
        fa_Seq.append(line)
    # read motiflen
    motif_len_path = "D:/IntroToDataMining-CS412/final project/benchmark/dataset"+str(file_number)+"/motiflength.txt"
    motiflen_in = open(motif_len_path,"r")
    motif_len = int(motiflen_in.readline())
    # set initial value
    ## step1 random pick initial site
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
    iteration_n = 150000
    for itera in range(iteration_n):
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
        if itera == (150000 - 1):
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
    # write into file:
    os.mkdir('Gibbs_No_F/dataset' + str(file_number))
    ## write sites to 'sites.txt'
    file = open('Gibbs_No_F/dataset' + str(file_number) + '/sites.txt', 'w+')
    for row in A[-1]:
        file.write(str(row) + '\n')
    file.close()
    ## write motif to 'motif.txt'
    file = open('Gibbs_No_F/dataset' + str(file_number) + '/motif.txt', 'w+')
    header = 'MOTIF{0}    {1}\n'.format(str(file_number), str(motif_len))
    file.write(header)
    for row in (M/fa_Num).T:
        line = ' '.join(str(x) for x in row)
        file.write(line + '\n')
    file.close()





def loop_gibbs_without_F():
    for i in range(70):
        start = timeit.default_timer()

        if i % 7 in [0,1] :
            gibbs_without_F_not_to_converge(i)
        else:
            gibbs_without_F_unitil_converge(i)
        print(str(i)+"th dataset")

        stop = timeit.default_timer()
        print('Time: ', stop - start)
loop_gibbs_without_F()


