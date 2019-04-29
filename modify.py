for i in range(70):
    file = open('benchmark/dataset{}/motif.txt'.format(str(i)), 'r+')
    content = file.read()
    file.close()
