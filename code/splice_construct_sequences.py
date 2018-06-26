""" Read in the fasta file containing the CECRs, including the
promoter, and a fasta file of the intervening vector sequence,
and splice sequences like so: CECR+Intervening+Promoter. """

import os
import sys
import getopt
import re
from Bio import SeqIO

if __name__ == "__main__":

	format = "fasta"
	pattern = ""

	try:	
		(opts, args) = getopt.getopt(sys.argv[1:], "f:i:I:o:p:")
	except getopt.GetoptError:
		print "Problem getting options\n"
		sys.exit(2)

	for opt, arg in opts:
		if opt == "-f":
			format = arg
		elif opt == "-i":
			inputfile = os.path.abspath(arg)
		elif opt == "-I":
			interveningseqfile = os.path.abspath(arg)
		elif opt == "-o":
			outputfile = os.path.abspath(arg)
		elif opt == "-p":
			promoter = arg
		else:
			print "Don't know %s \n" % opt
			sys.exit(2)

	if len(args) > 0:
		print """Need no arguments.\n"""
		sys.exit(2)

	interveningsequence_rec = list(SeqIO.parse(\
								interveningseqfile,\
								format))[0]
								
	input_dict = SeqIO.index(inputfile, format)
	print "Total records: %d" % len(input_dict)
	promoter_rec = input_dict[promoter]

	output_recs = []

	for input_record in SeqIO.parse(inputfile, format):

		if input_record.id == promoter:
			output_record = interveningsequence_rec + \
							promoter_rec

			print "Input: %d Output: %d" % \
									(len(promoter_rec),\
									len(output_record))

			output_record.id = promoter_rec.id
			output_record.name = promoter_rec.name
			output_record.description = promoter_rec.description

			output_recs.append(output_record)

		else:

			output_record = input_record + \
							interveningsequence_rec + \
							promoter_rec

			print "Input: %d Output: %d" % \
									(len(input_record),\
									len(output_record))

			output_record.id = input_record.id
			output_record.name = input_record.name
			output_record.description = input_record.description

			output_recs.append(output_record)

	
	SeqIO.write(output_recs, outputfile, format)
