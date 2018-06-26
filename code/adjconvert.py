import csv
import os
import sys
import getopt

class adjFile:

    def __init__(self, filename):
    "converts output adjacency matrix from aracne into table for Cytoscape"
		
        with open(filename, "r") as f, open(outfilename, "w") as of:
            reader = csv.reader(f, delimiter='\t')
            writer = csv.writer(of, delimiter=' ')
            for lines in reader:
                key = lines.pop(0)
                for i in range(0, len(lines), 2):
                    writer.writerow([key, lines[i], lines[i+1]])


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
        if opt == "-o":
            outfilename = arg

    adj = adjFile(filename)
