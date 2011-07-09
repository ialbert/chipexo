from optparse import OptionParser, IndentedHelpFormatter
import csv, os, sys, subprocess, shutil, threading, logging

from genetrack import get_output_path

def chromosome_list(fpath):
    file = open(fpath, 'rU')
    chromosomes = set()
    r = csv.reader(file, delimiter='\t')
    for line in r:
        chr = line[0]
        if chr not in chromosomes:
            chromosomes.add(chr)
    return chromosomes

class ProcessFileThread(threading.Thread):
    semaphore = threading.Semaphore(1)
    def __init__(self, command):
        threading.Thread.__init__(self)
        self.command = command
    def run(self):
        ProcessFileThread.semaphore.acquire()
        p = subprocess.Popen(self.command, shell=True)
        p.wait()
        ProcessFileThread.semaphore.release()


def process_file(fpath, exclusion, sigma):
    chroms = chromosome_list(fpath)
    outputs = []
    threads = []
    # Process each chromosome
    for chrom in chroms:
        c = 'python genetrack.py %s -s %d -e %d -c %s' % (fpath, sigma, exclusion, chrom)
        t = ProcessFileThread(c)
        t.start()
        threads.append(t)
        output_path = get_output_path(fpath, sigma, exclusion, chrom)
        outputs.append(output_path)
    # Wait for completion
    for thread in threads:
        thread.join()
    # Merge together output files
    real_output = open(get_output_path(fpath, sigma, exclusion, ''), 'wt')
    for output in outputs:
        shutil.copyfileobj(open(output, 'rt'), real_output)
        os.unlink(output)
    real_output.close
    
           
usage = '''
input_paths may be:
- a file or list of files to run on
- a directory or list of directories to run on all files in them
- "." to run in the current directory

example usages:
python genetrack.py -s 10 /path/to/a/file.txt path/to/another/file.txt
python genetrack.py -s 5 -e 50 /path/to/a/data/directory/
python genetrack.py .
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

def run():   
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-s', action='store', type='int', dest='sigma', default=5,
                      help='Sigma to use when smoothing reads to call peaks. Default 5.')
    parser.add_option('-e', action='store', type='int', dest='exclusion', default=20,
                      help='Exclusion zone around each peak that prevents others from being called. Default 20.')
    parser.add_option('-p', action='store', type='int', dest='processes', default=1,
                      help='Number of processes to run concurrently')
    (options, args) = parser.parse_args()
    
    ProcessFileThread.semaphore = threading.Semaphore(options.processes)

    if not args:
        parser.print_help()
        sys.exit(1)
        
    for path in args:
        if not os.path.exists(path):
            parser.error('Path %s does not exist.' % path)
        if os.path.isdir(path):
            for fname in os.listdir(path):
                fpath = os.path.join(path, fname)
                if os.path.isfile(fpath) and not fname.startswith('.'): 
                    process_file(fpath, options.exclusion, options.sigma)
        else:
            process_file(path, options.exclusion, options.sigma)
            
if __name__ == '__main__':
    logging.error('Multiprocess script is currently not working')