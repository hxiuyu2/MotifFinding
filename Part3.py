from scipy import *
import numpy as np
import matplotlib.pyplot as plt


# algorithm = 'Gibbs_No_F'

def compute_dkl(algorithm):
    res = []
    for i in range(70):
        with open("benchmark/dataset{}/motif.txt".format(str(i)),"r") as file1, open (algorithm+"/dataset{}/predictedmotif.txt".format(str(i)),"r") as file2:
            line1 = file1.readline()
            file2.readline()
            ML = int(line1.split()[1])
            line1 = file1.readline()
            line2 = file2.readline()
            sum = 0
            while line1 and line2:
                if line1.startswith('<') or line2.startswith('<'):
                    break
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
            res.append(sum/ML)
    return res


def compute_position(algorithm):
    result = []
    for i in range(70):
        with open("benchmark/dataset{}/motif.txt".format(str(i)), "r") as file1, open(
                algorithm+"/dataset{}/predictedsites.txt".format(str(i)), "r") as file3, open(
            "benchmark/dataset{}/sites.txt".format(str(i)), "r") as file2:
            num = int(file1.readline().split()[1])

            line1 = file2.readline()
            line2 = file3.readline()
            sum = 0

            while line1 and line2:
                P = line1.split()
                Q = line2.split()
                sum += num - min(abs(int(P[0]) - int(Q[0])), num)
                line1 = file2.readline()
                line2 = file3.readline()
            result.append(sum)
    return result



def comput_overlap(algorithm):
    result = []
    for i in range(70):
        with open("benchmark/dataset{}/sequences.fa".format(str(i)), "r") as file1, open(
                algorithm+"/dataset{}/predictedsites.txt".format(str(i)), "r") as file3, open(
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


def avg_std(metrics):
    means = np.mean(metrics, axis=1)
    stds = np.std(metrics, axis=1)
    return means, stds



def draw(y, filename):
    plt.clf()
    # plt.plot(range(len(y[0])), y[0], label='Gibbs_F_2')
    plt.plot(range(len(y)), y, label='MEME')
    # plt.plot(range(len(y[2])), y[2], label='Gibbs_No_F')
    # plt.plot(range(len(y)), y, label='merged')
    # plt.plot(range(len(y[1])), y[1], label='Gibbs_No_F')
    plt.legend(loc='upper left')
    plt.xticks(np.arange(7), ['ML=8,SC=10,ICPC=1','ML=8,SC=10,ICPC=1.5','ML=8,SC=10,ICPC=2','ML=6,SC=10,ICPC=2','ML=7,SC=10,ICPC=2','ML=8,SC=5,ICPC=2','ML=8,SC=20,ICPC=2'], rotation=90)
    plt.tight_layout()
    plt.savefig(filename)


def get_runtime():
    all_runtime = []

    # runtime_file = open('Gibbs_F_2/runtime.txt')
    # gibbs_f = np.zeros((10,7))
    # cnt = 0
    # for line in runtime_file:
    #     if line.startswith('Time'):
    #         t = line.split()[1]
    #         t = float(t)
    #         gibbs_f[int(cnt / 7)][cnt % 7] = t
    #         cnt += 1
    # all_runtime.append(gibbs_f.tolist())
    # runtime_file.close()

    runtime_file = open('MEME/runtime.txt')
    meme_time = np.zeros((7,10))
    cnt = 0
    for line in runtime_file:
        t = line.split(': ')[1][:-4]
        t = float(t)
        meme_time[cnt % 7][int(cnt / 7)] = t
        cnt += 1
    runtime_file.close()
    all_runtime.append(meme_time.tolist())

    # runtime_file = open('Gibbs_No_F/runtime.txt')
    # gibbs_time = np.zeros((10,7))
    # cnt = 0
    # for line in runtime_file:
    #     if line.startswith('Time'):
    #         t = line.split()[1]
    #         t = float(t)
    #         gibbs_time[int(cnt / 7)][cnt % 7] = t
    #         cnt += 1
    # all_runtime.append(gibbs_time.tolist())
    # runtime_file.close()
    # return all_runtime
    return meme_time


if __name__ == '__main__':
    kld = []
    pos = []
    site = []

    for algo in ['MEME']:
        dkl = compute_dkl(algo)
        dkl_final = seven_list(dkl)
        # dkl_final = np.array(dkl).reshape((7,10))

        position = compute_position(algo)
        position_final = seven_list(position)
        # position_final = np.array(position).reshape((7,10))

        overlap = comput_overlap(algo)
        overlap_final = seven_list(overlap)
        # overlap_final = np.array(overlap).reshape((7,10))

        avg, _ = avg_std(dkl_final)
        kld.append(avg)
        avg, _ = avg_std(position_final)
        pos.append(avg)
        avg, _ = avg_std(overlap_final)
        site.append(avg)

    draw(kld[0], 'KLDf.png')
    draw(pos[0], 'POSf.png')
    draw(site[0], 'SITEf.png')
    # time, _ = avg_std(get_runtime())
    # draw(time, 'TIME.png')
    # print(time)
    #
    # print(kld)
    # print(pos)
    # print(site)
    # print(time)
    #
    # avg, std = avg_std(kld)
    # print('average for KLD is {}, standard error for KLD is {}'.format(avg, std))
    # avg, std = avg_std(pos)
    # print('average for overlap position is {}, standard error for overlap position is {}'.format(avg, std))
    # avg, std = avg_std(site)
    # print('average for overlap site is {}, standard error for overlap site is {}'.format(avg, std))
    # avg, std = avg_std(time)
    # print('average for runtime is {}, standard error for runtime is {}'.format(np.mean(time), np.std(time)))
