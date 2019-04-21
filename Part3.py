from scipy import *
import numpy as np
import matplotlib.pyplot as plt
import time
import Part2 as p2


def compute_dkl():
    res = [[]]*7
    for i in range(70):
        with open("benchmark/dataset{}/motif.txt".format(str(i)),"r") as file1, open ("prediction/dataset{}/predictedmotif.txt".format(str(i)),"r") as file2:
            next(file1)
            next(file2)
            line1 = file1.readline()
            line2 = file2.readline()
            sum = 0
            while line1 and line2:
                P = line1.split()
                Q = line2.split()
                for j in range(4):

                    if float(Q[j]) == 0.0 or float(P[j]) == 0.0:
                        temp = 0.0
                    else:
                        temp = float(P[j]) * log(float(P[j]) / float(Q[j]))
                    sum += temp

                line1 = file1.readline()
                line2 = file2.readline()
            res[i%7].append(sum)
    return res


def get_position():
    result = [[]]*7
    for i in range(70):
        # read motif length from 'motiflength.txt'
        motif_len_txt = open('benchmark/dataset{}/motiflength.txt'.format(str(i)))
        motif_len = int(motif_len_txt.readline().strip())
        motif_len_txt.close()

        # get number of overlap position
        predicted = open('prediction/dataset{}/predictedsites.txt'.format(str(i)))
        true_site = open('benchmark/dataset{}/sites.txt'.format(str(i)))
        for k in range(motif_len):
            pred_val = int(predicted.readline().strip())
            true_val = int(true_site.readline().strip())
            positions = motif_len - abs(pred_val - true_val)
            result[i%7].append(positions)
        predicted.close()
        true_site.close()
    return result


def get_overlap():
    result = [[]] * 7
    for i in range(70):
        # get number of overlap position
        predicted = open('prediction/dataset{}/predictedmotif.txt.txt'.format(str(i)))
        true_motif = open('benchmark/dataset{}/motif.txt'.format(str(i)))

        # get motif length from the first line of motif.txt
        motif_len = int(true_motif.readline().strip().split()[1])
        match_cnt = 0
        for k in range(motif_len):
            pred_val = predicted.readline().strip()
            true_val = true_motif.readline().strip()

            # convert to float list
            pred_val = np.array([float(i) for i in pred_val])
            true_val = np.array([float(i) for i in true_val])

            # get site[k]
            pred_alpha = np.argmax(pred_val)
            true_alpha = np.argmax(true_val)
            if pred_alpha == true_alpha:
                match_cnt += 1
        result[i%7].append(match_cnt)
        predicted.close()
        true_motif.close()
    return result


def get_runtime():
    start = time.time()
    # TODO:
    # code for part 2
    end = time.time()
    return end - start


def avg_std(metrics):
    means = np.mean(metrics, axis=1)
    stds = np.std(metrics, axis=1)
    return means, stds


def draw(y, filename):
    plt.clf()
    plt.plot(range(len(y)), y)
    plt.savefig(filename)
