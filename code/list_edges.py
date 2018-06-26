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



    with open(filename, 'r') as handle, open(outfilename, 'a') as w:

        seen = set()
        for line in handle:
            element = line.split()
            genes = '-'.join(sorted(element[0:2]))
            if genes not in seen:
                seen.add(genes)
        w.write('\n'.join(seen))