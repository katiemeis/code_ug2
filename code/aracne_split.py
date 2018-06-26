import os
import sys
import getopt

class splitfile:

    def __init__(self, filename):
    "splits input file by row and runs each separately"

        try:
            filehandle = open(filename, "r")
            try:
                table = filehandle.read()
        except IOError:
            print "File %s doesn't exist\n" % filename
            sys.exit(2)

        lines = table.splitlines()
        outputfile = open("outfile", "w")
        tmpfile = os.tmpnam()
        os.system("cut -f1 -d '\t' {} > {}".format(filename, tmpfile))

        with open (tmpfile) as f:
            n = [x.strip() for x in f.readlines()]

        for i in n[1:]:
            os.system("~/aracne/source/ARACNE/aracne2 -H ~/aracne/source/ARACNE/ -i {0} -h {1} -o {2}".format(filename, i, tmpfile))
            os.system("cat {0} >> {1}".format(tmpfile, "outfile"))
            os.remove(tmpfile)

        os.system("grep -v '>' {} > {}".format("outfile", finalout))

if __name__ == "__main__":

    filename = ''
    outfile = ''

    try:    
        (opts, args) = getopt.getopt(sys.argv[1:], "i:o:")
    except getopt.GetoptError:
        print "Problem getting options\n"
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-i":
            filename = arg
        elif opt == "-o":
            finalout = arg
        else:
            print "Not an option"

    split = splitfile(filename)
