handle1 = open("pcc_nodes.txt",'r') #list of nodes from cyto
genes1 = set()
for line in handle1:
    genes1.add(line)

genes2 = set()
handle2 = open("reg_nodes.txt",'r') #list of nodes from cyto
for line in handle2:
    genes2.add(line)

print(len(genes1 & genes2))
