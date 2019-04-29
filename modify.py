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
    file = open('dataset/dataset{}/sites.txt'.format(str(i)), 'r')
    contents = file.read()
    file.close()

    contents = contents.split(',')
    contents = '\n'.join(contents)

    file = open('dataset/dataset{}/sites.txt'.format(str(i)), 'w')
    file.write(contents)
    file.close()
