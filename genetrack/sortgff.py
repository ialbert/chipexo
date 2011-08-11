# gfftoidx.py
#
# Sort a .gff file into order suitable for input to genetrack.py
#
# By Pindi Albert, 2011
#
# Input: .gff format reads
# Format: standard gff, score interpreted as read count
#
# Output: sorted reads
# Format: same
# Order: grouped by chromosome and sorted by start index within chromosome group
#
# Run with no arguments or -h for usage and command line options

import csv, os, logging, sys, subprocess
from optparse import OptionParser, IndentedHelpFormatter

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def get_output_path(input_path, options):
    directory, fname = os.path.split(input_path)
    
    fname = ''.join(fname.split('.')[:-1]) # Strip extension (will be re-added as appropriate)
    
    output_dir = os.path.join(directory, 'sortgff')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    return os.path.join(output_dir, '%s.gff' % fname)




def process_file(path, options):
    
    if not path.endswith('.gff'):
        logging.error('File "%s" is not .gff' % path)
        return
    
    logging.info('Processing file "%s"' % path)
    
    output_path = get_output_path(path, options)
    
    command = 'sort -k 1,1 -k 4,4n "%s" > "%s"' % (path, output_path)
    logging.info('Executing "%s"' % command)
    p = subprocess.Popen(command, shell=True)
    p.wait()
        
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
    
