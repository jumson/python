#!/usr/bin/python
import binascii
import csv 
from optparse import OptionParser

## Setup the option parser
parser = OptionParser()
parser.add_option("-i", "--infile", dest="infile",
                  help="CSV filename output from Saleae Analyzer", metavar="INFILE")
parser.add_option("-o", "--outfile", dest="outfile",
                  help="Binary output file for the bytes (it will overwrite)", metavar="OUTFILE")       
parser.add_option("-b", "--block_size", dest="block_size", type="int", default=512,
                  help="This script will write bytes to outfile in this size chunks. This might help with gigantic files and memory issues", metavar="BLOCK_SIZE")   

## Putting the options into option
(options, args) = parser.parse_args()

# check for required stuff
if not options.infile or not options.outfile:
    print("\t You must specify an infile and an outfile (outfile will be overwritten)")
    exit(-1)

ASCII_byte = ""
BINARY_byte = b""

# https://www.geeksforgeeks.org/working-csv-files-python/
csvfile = open(options.infile, 'r')
binfile = open(options.outfile, 'w')

try:
    # creating a csv reader object 
    csvreader = csv.reader(csvfile) 
    # throw away first row (headers)
    csvreader.next() 
    # throw away next three rows (three null bytes?)
    csvreader.next() 
    csvreader.next() 
    csvreader.next() 
    # extracting each data row one by one 
    for row in csvreader: 
        # as strings: '(0xAD)'
        a = row[2].split()[5]
        # appending the slice of just the byte value
        ASCII_byte = a[3:5]
        # making it into an actual byte
        BINARY_byte += binascii.unhexlify(ASCII_byte)
        # write when it is big enough 
        if len(BINARY_byte) > options.block_size:
            binfile.write(BINARY_byte)
            BINARY_byte = b""

    binfile.write(BINARY_byte)
    csvfile.close()
    binfile.close()

except:
    csvfile.close()
    binfile.close()

