handle1 = open("mi_edgelist_nodups.txt",'r')
seen1 = set()
for line in handle1:
    element = line.split()
    genes1 = '-'.join(sorted(element[0:2]))
    if genes1 not in seen1:
        seen1.add(genes1)

handle2 = open("pcc_edgelist_nodups_quote.txt",'r')
seen2 = set()
for line in handle2:
    element = line.split()
    genes2 = '-'.join(sorted(element[0:2]))
    if genes2 not in seen2:
        seen2.add(genes2)

w = open("edges_both.txt",'w')
for item in list(seen1 & seen2):
    w.write("%s\n" % item)
w.close()
