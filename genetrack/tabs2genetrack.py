"""
File transformer. Takes a tabuolar file (bed or gff) and 
transforms them to a file used as genetrack input.
Optional parameters may be used to shift the 5' ends with 
a specified amount. This is useful if the file corresponds to data 
with fixed fragment widths you can move each fragment to the 
center location. The program may be invoked in multiple ways. 
As a standalone script::

    python tabs2genetrack.py

As a python module::

    python -m genetrack.scripts.tabs2genetrack

Or in other python scripts::

>>>
>>> from genetrack.scripts import tabs2genetrack
>>> tabs2genetrack.transform( inpname, outname, format='bed', shift=0)
>>>

Run the script with no parameters to see the options that it takes.

**Note1**: The script will invoke the system `sort` command to sort the file 
that is substantially faster under Unix than Windows.

"""
import os, sys, csv, shutil, re, time
from itertools import *

def path_join(*args):
    "Builds absolute path"
    return os.path.abspath(os.path.join(*args))

curr_dir = os.path.dirname( __file__ )
TEMP_DATA_DIR = path_join( curr_dir, 'tempdir' )

# create the temporary data directory if not present
if not os.path.isdir( TEMP_DATA_DIR ):
    os.mkdir( TEMP_DATA_DIR )

def tempdata(*args):
    "Generates paths to temporary data"
    return path_join(TEMP_DATA_DIR, *args)

def commify(n):
    """
Formats numbers with commas

>>> commify(10000)
'10,000'
"""
    n = str(n)
    while True:
        (n, count) = re.subn(r'^([-+]?\d+)(\d{3})', r'\1,\2', n)
        if count == 0:
            break
    return n

class Timer(object):
    """
A timer object for display elapsed times.

>>> timer = Timer()
>>> timer.format(30)
'30.00 seconds'
>>> timer.format(320)
'5.3 minutes'
>>> timer.format(3200)
'53.3 minutes'
>>> timer.format(30500)
'8.5 hours'
"""
    def __init__(self):
        self.start()

    def start(self):
        self.start_time = time.time()

    def format(self, value):
        min1 = 60.0
        hour = 60 * min1
        if value < 60:
            return '%4.2f seconds' % value
        elif value < hour:
            return '%3.1f minutes' % (value/min1)
        else:
            return '%3.1f hours' % (value/hour)

    def report(self):
        elapsed = self.stop()
        return self.format( elapsed )

    def stop(self):
        elapsed = time.time() - self.start_time
        self.start()
        return elapsed

def consolidate( inpname, outname, format):
    """
    Consolidates an input file. 
    Merges multiple indentical indices into one line
    """
    fp = open(outname, 'wt')

    # recover the original basename
    basename = os.path.basename(inpname).replace('.sorted', '')

    # create a few information headers
    fp.write("#\n# created with tabs2genetrack\n")
    fp.write("# source: %s, format %s\n#\n" % (basename, format) )
    fp.write("chrom\tindex\tforward\treverse\tvalue\n")
    
    # the main reader
    reader = csv.reader(open(inpname, 'rb'), delimiter='\t')

    # will duplicate some code to avoid an extra conditional internally
    chrom, index, fwd, rev, val = reader.next()
    lastindex, fwd, rev, val = int(index), int(fwd), int(rev), int(val)
    collect = [ chrom, lastindex, fwd, rev, val ]

    for row in reader:
        chrom, index, fwd, rev, val = row
        index, fwd, rev, val = int(index), int(fwd), int(rev), int(val)
        if index == lastindex:
            collect[2] += fwd
            collect[3] += rev
            collect[4] += val
        else:
            # write out the collected data
            fp.write( '%s\t%s\t%s\t%s\t%s\n' % tuple(collect) )
            collect = [ chrom, index, fwd, rev, val ]
            lastindex = index

    fp.close()

def transform(inpname, outname, format, shift=0, index=False, options=None):
    """
    Transforms reads stored in bedfile to a genetrack input file.
    Requires at least 6 bed columns to access the strand.
    """

    # detect file formats
    if format == "BED":
        CHROM, START, END, STRAND = 0, 1, 2, 5
    elif format == "GFF":
        CHROM, START, END, STRAND = 0, 3, 4, 6
    else:
        raise Exception('Invalid file format' % format)

    # two sanity checks, one day someone will thank me
    if format == 'BED' and inpname.endswith('gff'):
        raise Exception('BED format on a gff file?')
    if format == 'GFF' and inpname.endswith('bed'):
        raise Exception('GFF format on a bed file?')

    # find the basename of the outputname
    basename = os.path.basename(outname)
   
    # two files store intermediate results
    flat   = tempdata( '%s.flat' % basename )
    sorted = tempdata( '%s.sorted' % basename )

    # check for track information on first line, 
    # much faster this way than conditional checking on each line
    fp = file(inpname, 'rU')
    first = fp.readline()
    fp.close()

    # create the reader
    reader = csv.reader(file(inpname, 'rU'), delimiter='\t')

    # skip if trackline exists
    if first.startswith == 'track':
        reader.next()

    # unwind the comments
    list(takewhile(lambda x: x[0].startswith('#'), reader))

    # copious timing info for those who enjoy these
    timer, full = Timer(), Timer()

    print("parsing '%s'" % inpname)
    print("output to '%s'" % outname)

    # create the unsorted output file and apply corrections
    print("unsorted flat file '%s'" % flat)

    fp = file(flat, 'wt')
    for linec, row in enumerate(reader):
        try:
            chrom, start, end, strand = row[CHROM], row[START], row[END], row[STRAND]
        except Exception, exc:
            first = row[0][0]
            # may be hitting the end of the file with other comments
            if  first == '>':
                break # hit the sequence content of the gff file
            elif first == '#':
                continue # hit upon some comments
            else:
                print(row)
                raise Exception(exc) 

        if strand == '+':
            # on forward strand, 5' is at start
            idx = int(start) + shift
            fwd, rev, val = 1, 0, 1
        elif strand == '-':
            # on reverse strand, 5' is at end
            idx = int(end) - shift
            fwd, rev, val = 0, 1, 1
        else:
            # no strand specified, generate interval centers
            idx = (int(start)+int(end))/2
            fwd, rev, val = 0, 0, 1

        # it is essential be able to sort the index as a string! 
        fp.write('%s\t%012d\t%s\t%s\t%s\n' % (chrom, idx, fwd, rev, val))

    fp.close()
    linet = commify(linec)
    print("parsing %s lines finished in %s" % (linet, timer.report()))

    # if it is producing coverage then it will expand reads into full intervaals

    # now let the sorting commence
    cmd = "sort %s > %s" % (flat, sorted)
    print("sorting into '%s'" % sorted)
    os.system(cmd)
    print("sorting finished in %s" % timer.report() )

    print("consolidating into '%s'" % outname)
    consolidate( sorted, outname, format=format)
    print("consolidate finished in %s" % timer.report() )
    print("output saved to '%s'" % outname)
    print("full conversion finished in %s" % full.report() )

    # attempting to cleanup the remaining files
    for name in (flat, sorted):
        print("removing temporary file '%s'" % name )
        os.remove(name)

def option_parser():
    "The option parser may be constructed in other tools invoking this script"
    import optparse

    usage = "usage: %prog -i inputfile -o outputfile -f format"

    parser = optparse.OptionParser(usage=usage)

    # setting the input file name
    parser.add_option(
        '-i', '--input', action="store", 
        dest="inpname", type='str', default=None,
        help="the input file name (required)"
    )

    # setting the output file name
    parser.add_option(
        '-o', '--output', action="store", 
        dest="outname", type='str', default=None,
        help="output file name (required)"
    )

    # file formats
    parser.add_option(
        '-f', '--format', action="store", 
        dest="format", type="str", default='', 
        help="input file format, bed or gff (required)",
    )
    
    # correction shift added in 5' direction for start/end coordinates
    parser.add_option(
        '-s', '--shift', action="store", 
        dest="shift", type="int", default=0, 
        help="shift for the 5' end on each strand (default=0)",
    )

    # verbosity can be 0,1 and 2 (increasing verbosity)
    parser.add_option(
        '-v', '--verbosity', action="store", 
        dest="verbosity", type="int", default=1, 
        help="sets the verbosity (0, 1) (default=1)",
    )

    parser.add_option("-x", "--index",
        action="store_true", dest="index", default=False,
        help="also creates an hdf index for the file")

    parser.add_option(
        '-w', '--workdir', action="store", 
        dest="workdir", type='str', default=None,
        help="work directory (optional)"
    )

    return parser

if __name__ == '__main__':
    
    parser = option_parser()

    options, args = parser.parse_args()

    # uppercase the format
    options.format = options.format.upper()

    if options.format not in ('BED', 'GFF'):
        sys.stdout = sys.stderr
        parser.print_help()
        sys.exit(-1)

    # missing file names
    if not (options.inpname and options.outname and options.format):
        parser.print_help()
        sys.exit(-1)
    else:
        transform(inpname=options.inpname, outname=options.outname,\
            format=options.format, shift=options.shift, index=options.index, options=options)
