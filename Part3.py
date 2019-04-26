from scipy import *
import numpy as np
import matplotlib.pyplot as plt
import time
#import Part2 as p2


def compute_dkl():
    res = []
    for i in range(70):
        with open("benchmark/dataset{}/motif.txt".format(str(i)),"r") as file1, open ("MEME/dataset{}/predictedmotif.txt".format(str(i)),"r") as file2:
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
            res.append(sum)
    return res


def compute_position():
    result = []
    for i in range(70):
        with open("benchmark/dataset{}/motif.txt".format(str(i)), "r") as file1, open(
                "MEME/dataset{}/predictedsites.txt".format(str(i)), "r") as file3, open(
            "benchmark/dataset{}/sites.txt".format(str(i)), "r") as file2:
            num = int(file1.readline().split()[1])

            line1 = file2.readline()
            line2 = file3.readline()
            sum = 0

            while line1 and line2:
                P = line1.split()
                Q = line2.split()
                sum += num - max(abs(int(P[0]) - int(Q[0])), num)
                #
                # if int(P[0]) < int(Q[0]):
                #     temp = int(P[0]) + num - int(Q[0])
                # else:
                #     temp = int(Q[0]) + num - int(P[0])
                # if temp <= 0:
                #     temp = 0
                #
                # temp += temp
                line1 = file2.readline()
                line2 = file3.readline()
            # result.append(temp)
            result.append(sum)
    return result



def comput_overlap():
    result = []
    for i in range(70):
        with open("benchmark/dataset{}/sequences.fa".format(str(i)), "r") as file1, open(
                "MEME/dataset{}/predictedsites.txt".format(str(i)), "r") as file3, open(
            "benchmark/dataset{}/sites.txt".format(str(i)), "r") as file2, open(
            "benchmark/dataset{}/motiflength.txt".format(str(i)),"r") as file4:
            length = int(file4.readline().split()[0])
            seq = file1.readline()
            line1 = file2.readline()
            line2 = file3.readline()
            res_count=0
            while line1 and line2:
                P = line1.split()
                Q = line2.split()
                original = seq[int(P[0]):int(P[0])+length-1]
                predict = seq[int(Q[0]):int(Q[0])+length-1]
                count = 0

                # index1=0
                # index2=0
                # temp_str=""
                #print(original,predict)
                # while index1 < len(original) and index2 < len(predict):
                #     if original[index1] == predict[index2]:
                #         count+=1
                #     index1+=1
                #     index2+=1

                for pred_char, true_char in zip(predict, original):
                    if pred_char == true_char:
                        count += 1

                if count >= length / 2:
                    res_count+=1

                seq = file1.readline()
                line1 = file2.readline()
                line2 = file3.readline()
            result.append(res_count)
    return result


def seven_list(List):
    res = []
    for i in range(0,7):
        temp = []
        for j in range(i,len(List),7):
            temp.append(List[j])
        res.append(temp)

    return res


# def get_runtime():
#     start = time.time()
#     # TODO:
#     # code for part 2
#     end = time.time()
#     return end - start


def avg_std(metrics):
    means = np.mean(metrics, axis=1)
    stds = np.std(metrics, axis=1)
    return means, stds



def draw(y, filename):
    plt.clf()
    plt.plot(range(len(y)), y)
    plt.savefig(filename)



if __name__ == '__main__':

    dkl = compute_dkl()
    print(len(dkl))
    print(dkl)
    dkl_final = seven_list(dkl)

    position = compute_position()
    print(len(position))
    print(position)
    position_final = seven_list(position)

    overlap = comput_overlap()
    print(len(overlap))
    print(overlap)
    overlap_final = seven_list(overlap)


    print(dkl_final)
    print(position_final)
    print(overlap_final)

    print('average + std for dkl: ', avg_std(dkl_final))
    print('average + std for position: ', avg_std(position_final))
    print('average + std for overlap: ', avg_std(overlap_final))


