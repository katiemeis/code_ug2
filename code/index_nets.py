

import os
import sys
import getopt
import csv

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

    node_index = {}
    latitude_dict = {}
    longitude_dict = {}

    with open("Swz_data2.txt", 'r') as r:
        for line in r:
            item = line.split()
            node_index[item[1]] = item[0]
            latitude_dict[item[1]] = item[5]
            longitude_dict[item[1]] = item[6]

    with open(filename, 'r') as r, open(outfilename, 'a') as w:
        w.write("%s\t%s\t%s\n" % ("DataIndex", "Latitude", "Longitude"))
        for line in r:
            field = line.split()
            barcode = field[1]
            index = node_index[barcode]
            lat = latitude_dict[barcode]
            lon = longitude_dict[barcode]
            w.write("%s\t%s\t%s\n" % (index, lat, lon))
