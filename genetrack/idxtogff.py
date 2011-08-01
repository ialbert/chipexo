import csv, os, logging
from optparse import OptionParser

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def process_file(path):
    logging.info('Processing file "%s"' % path)
    dir, fname = os.path.split(path)
    comps = fname.split('.')
    out_path = os.path.join(dir, '%s.gff' % '.'.join(comps[:-1]))
    r = csv.reader(open(path, 'rU'), delimiter='\t')
    out = csv.writer(open(out_path, 'wt'), delimiter='\t')
    for line in r:
        if len(line) < 4:
            continue
        try:
            chr, index, forward, reverse = line[0], int(line[1]), int(line[2]), int(line[3])
        except ValueError:
            continue
        for i in range(forward):
            out.writerow((chr, 'idxtogff', '.', index, index+1, 1, '+', '.', '.'))
        for i in range(reverse):
            out.writerow((chr, 'idxtogff', '.', index, index+1, 1, '-', '.', '.'))
        
        
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
                    process_file(fpath)
        else:
            process_file(path)
            
            
run()
    
