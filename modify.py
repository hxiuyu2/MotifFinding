path = 'project_submission/output/gibbs_no_f/loop80000'
for i in range(70):
    file = open(path+'/dataset{}/predictedmotif.txt'.format(str(i)), 'r')
    content = file.read()
    file.close()

    file = open(path + '/dataset{}/predictedmotif.txt'.format(str(i)), 'w')
    file.seek(0,0)
    content = '>'+content+'<'
    file.write(content)
    file.close()

    # file = open(path+'/dataset{}/predictedmotif.txt'.format(str(i)), 'a+')
    # file.write('<')
    # file.close()
