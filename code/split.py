import os
import sys
import getopt

class splitfile:

    def __init__(self, filename):

        filehandle = open(filename, "r")
        table = filehandle.read()
        lines = table.splitlines()
        outputfile = open(outfile, "w")
        tmpfile = os.tmpnam()

        for line in lines[1:]:
            part = line.split("\t")
            probe = part[0]
            os.system("~/aracne/source/ARACNE/aracne2 -H ~/aracne/source/ARACNE/ -i {0} -h {1} -o {2}".format(filename, probe, tmpfile))
            os.system("cat {0} >> {1}".format(tmpfile, outfile))
            os.remove(tmpfile)

        os.system("grep -v '>' {} > {}".format(outfile, finalout))


if __name__ == "__main__":

    filename = ''
    outfile = ''
    finalout = ''

    try:    
        (opts, args) = getopt.getopt(sys.argv[1:], "i:o:f:")
    except getopt.GetoptError:
        print "Problem getting options\n"
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-i":
            filename = arg
        if opt == "-o":
            outfile = arg
	if opt == "-f":
            finalout = arg

    split = splitfile(filename)
