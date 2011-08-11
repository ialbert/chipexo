# idxtogff.py
#
# Converts reads in .idx format to .gff
#
# By Pindi Albert, 2011
#
# Input: .idx format reads
# Format: tab-separated chromosome (chr##), index, + reads, - reads
#
# Output: .gff format reads
# Format: standard gff, one line per read
#
# Run with no arguments or -h for usage and command line options

import csv, os, logging
from optparse import OptionParser

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def get_output_path(input_path, options):
    directory, fname = os.path.split(input_path)
    
    fname = ''.join(fname.split('.')[:-1]) # Strip extension (will be re-added as appropriate)
    
    output_dir = os.path.join(directory, 'idxtogff')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    return os.path.join(output_dir, '%s.gff' % fname)


def process_file(path, options):
    logging.info('Processing file "%s"' % path)
    
    output_path = get_output_path(path, options)
    
    reader = csv.reader(open(path,'rU'), delimiter='\t')
    writer = csv.writer(open(output_path, 'wt'), delimiter='\t')
    
    for line in reader:
        if len(line) < 4:
            continue
        try:
            chr, index, forward, reverse = line[0], int(line[1]), int(line[2]), int(line[3])
        except ValueError:
            continue
        for i in range(forward):
            writer.writerow((chr, 'idxtogff', '.', index, index+1, 1, '+', '.', '.'))
        for i in range(reverse):
            writer.writerow((chr, 'idxtogff', '.', index, index+1, 1, '-', '.', '.'))
        
        
def run():   
    parser = OptionParser(usage='%prog [options] input_paths')
    (options, args) = parser.parse_args()
        
        
    if not args:
        parser.print_help()
        sys.exit(1)
        
        
    for path in args:
        if not os.path.exists(path):
            parser.error('Path %s does not exist.' % path)
        if os.path.isdir(path):
            for fname in os.listdir(path):
                fpath = os.path.join(path, fname)
                if os.path.isfile(fpath) and not fname.startswith('.') and not fname.endswith('.gff'):
                    process_file(fpath, options)
        else:
            process_file(path, options)
            
            
run()
    
