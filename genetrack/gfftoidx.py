import csv, os, logging, sys
from optparse import OptionParser, IndentedHelpFormatter

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

from genetrack import ChromosomeManager

def get_output_path(input_path, options):
    directory, fname = os.path.split(input_path)
    
    fname = ''.join(fname.split('.')[:-1]) # Strip extension (will be re-added as appropriate)
    
    output_dir = os.path.join(directory, 'gfftoidx')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    return os.path.join(output_dir, '%s.idx' % fname)


def process_file(path, options):
    logging.info('Processing file "%s"' % path)
    
    output_path = get_output_path(path, options)
    
    reader = csv.reader(open(path,'rU'), delimiter='\t')
    writer = csv.writer(open(output_path, 'wt'), delimiter='\t')
    writer.writerow(['chrom', 'index', 'forward', 'reverse'])
    
    manager = ChromosomeManager(reader)

    while not manager.done:
        cname = manager.chromosome_name()
        logging.info('Processing chromosome %s' % cname)
        data = manager.load_chromosome()
        for read in data:
            writer.writerow([cname] + read)
        
usage = '''
input_paths may be:
- a file or list of files to run on
- a directory or list of directories to run on all files in them
- "." to run in the current directory
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description
        
def run():   
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-v', action='store_true', dest='verbose', help='Verbose mode: displays debug messages')
    parser.add_option('-q', action='store_true', dest='quiet', help='Quiet mode: suppresses all non-error messages')
    (options, args) = parser.parse_args()
    
    
    if options.verbose:
        logging.getLogger().setLevel(logging.DEBUG) # Show all info/debug messages
    if options.quiet:
        logging.getLogger().setLevel(logging.ERROR) # Silence all non-error messages
        
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
            
     

if __name__ == '__main__':       
    run()
    
