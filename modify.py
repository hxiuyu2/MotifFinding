# for i in range(70):
#     file = open('benchmark/dataset{}/motif.txt'.format(str(i)), 'wr+')
#     ind = 0
#     contents = []
#     for line in file:
#         line = '>sequence'+str(ind)+'\n'+line
#         contents.append(line)
#         ind += 1
#     contents = ''.join(contents)
#     file.write(contents)
#     file.close()
for i in range(70):
    file = open('benchmark/dataset{}/sequences.fa'.format(str(i)), 'r')
    ind = 0
    contents = []
    for line in file:
        line = '>sequence'+str(ind)+'\n'+line
        contents.append(line)
        ind += 1
    file.close()
    file = open('benchmark/dataset{}/sequences.fa'.format(str(i)), 'w')
    contents = ''.join(contents)
    file.write(contents)
    file.close()
