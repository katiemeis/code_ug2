"""Converts probe IDs to names after reformatting aracne results. Input file is reformatted aracne results (adjconvert.py) and output file is ready for Cytoscape"""


import os
import sys
import getopt

if __name__ == "__main__":

    filename = ''
    outfilename = ''

    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "i:o:")
    except getopt.GetoptError:
        print "Problem getting options\n"
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-i":
            filename = arg
        elif opt == "-o":
            outfilename = arg
        else:
            print "Don't know %s \n" % opt
            sys.exit(2)

    probes = {}

    for line in open("GeneList.txt", 'r'):

        item = line.split()
        key, value = item [0], item [2]
        probes[key] = value

    with open(filename, 'r') as r, open(outfilename, 'a') as w:

        for line in r:
            field = line.split()
            probe1 = field[0]
            probe2 = field[1]
            MI = field[2]

            try:
                name1 = probes[field[0]]
                new1 = field[0].replace(probe1, name1)
            except KeyError:
                new1 = field[0]

            try:
                name2 = probes[field[1]]
                new2 = field[1].replace(probe2, name2)
            except KeyError:
                new2 = field[1]

            w.write(new1 + ' ' + new2 + ' ' + MI + '\n')
