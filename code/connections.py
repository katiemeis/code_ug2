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


    hemgenes = ['Gata1', 'Csf1r', 'Csf2ra', 'Csf2rb', 'Il3ra',  'Gata2', 'Spi1', 'Cebpa', 'Cebpd', 'Gfi1', 'Gfi1b', 'Epor', 'Klf12', 'Klf1', 'Klf11', 'Klf16', 'Klf17', 'Klf15', 'Klf13', 'Klf10', 'Klf14']

    with open(filename, 'r') as handle, open(outfilename, 'a') as w:

        for line in handle:
            element = line.split()
            keygene = element[0]
            
            if keygene in hemgenes:
                w.write(line)
