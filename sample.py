"""Sample script to test the test harness
Doubles the input text """

from optparse import OptionParser
import os

def perform_double(inpath, outpath):
    infile = open(inpath, 'rt')
    outfile = open(outpath, 'wt')
    outfile.write(infile.read()*2)


if __name__ == '__main__':
    parser = OptionParser(usage='%prog [options] inputfile outputfile')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        parser.error('Must specify input and output files')
    if len(args) > 2:
        parser.error('Too many arguments')
    inpath, outpath = args
    if not os.path.exists(inpath):
        parser.error('Input file does not exist')
    perform_double(inpath, outpath)