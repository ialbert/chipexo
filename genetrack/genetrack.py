from optparse import OptionParser, IndentedHelpFormatter
import csv, logging, numpy, math, bisect, sys, os

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

WIDTH = 100

class Peak(object):
    def __init__(self, index):
        self.index = index
        self.start = index - 10
        self.end = index + 10
        self.value = 0
        self.deleted = False
        self.safe = False
    def __repr__(self):
        return '[%d] %d' % (self.index, self.value)

def parse_reads(reader):
    chromosomes = {}
    reader.next()
    for line in reader:
        cname, index, forward, reverse = line[:4]
        if cname not in chromosomes:
            chromosomes[cname] = []
        chromosomes[cname].append([int(index), int(forward), int(reverse)])
    return chromosomes

def make_keys(data):
    return [read[0] for read in data]
    
def make_peak_keys(peaks):
    return [peak.index for peak in peaks]

def get_window(data, start, end, keys):
    ''' Returns all reads from the data set with index between the two indexes'''
    start_index = bisect.bisect_left(keys, start)
    end_index = bisect.bisect_right(keys, end)
    return data[start_index:end_index]
    
def get_index(value, keys):
    ''' Returns the index of the value in the keys using bisect '''
    return bisect.bisect_left(keys, value)

def allocate_array(data, width=0):
    ''' Allocated a new array with the dimensions required to fit all reads in the
    argument. The new array is totally empty.'''
    hi = max([item[0] for item in data])
    return numpy.zeros(hi+width*2, numpy.float)
    
def normal_array(width, sigma, normalize=True):
    ''' Returns an array of the normal distribution of the specified width '''
    log2, sigma2 = math.log(2), float(sigma)**2
    
    def normal_func(x):
        return math.exp( -x * x / ( 2 * sigma2 ))
        values = map( func, range(-width, width) )
        
    # width is the half of the distribution
    values = map( normal_func, range(-width, width) )
    values = numpy.array( values, numpy.float )

    # normalization
    if normalize:
        values = 1.0/math.sqrt(2 * numpy.pi * sigma2) * values 

    return values

def call_peaks(array, data, keys, direction, width=WIDTH, exclusion=0):
    peaks = []
    def find_peaks():
        # Go through the array and call each peak
        results = (array > numpy.roll(array, 1)) & (array > numpy.roll(array, -1))
        indexes = numpy.where(results)
        for index in indexes[0]:
            peaks.append(Peak(index))
        #history = [array[0], array[1]] # Last 2 values
        #for i, value in enumerate(array[2:]):
        #    if history[1] > history[0] and history[1] > value:
        #        # history[1] is a peak
        #        peaks.append(Peak(i+1-WIDTH, history[1]))
        #    history[0] = history[1] # Shift history values
        #    history[1] = value
    find_peaks()
        
    def calculate_reads():
        # Calculate the number of reads in each peak
        for peak in peaks:
            reads = get_window(data, peak.start, peak.end, keys)
            peak.value = sum([read[direction] for read in reads])
    calculate_reads()
        
    before = len(peaks)
        
    def perform_exclusion():
        # Process the exclusion zone
        peak_keys = make_peak_keys(peaks)
        peaks_by_value = peaks[:]
        peaks_by_value.sort(key=lambda peak: -peak.value)
        for peak in peaks_by_value:
            peak.safe = True
            window = get_window(peaks, peak.index-exclusion, peak.index+exclusion, peak_keys)
            for excluded in window:
                if excluded.safe:
                    continue
                i = get_index(excluded.index, peak_keys)
                del peak_keys[i]
                del peaks[i]
    perform_exclusion()
            
    after = len(peaks)
    logging.debug('%d of %d peaks (%d%%) survived exclusion' % (after, before, after*100/before))
            
    return peaks
    
def process_file(path, sigma, exclusion):
    
    logging.info('Processing file "%s"' % path)
    
    directory, fname = os.path.split(path)
    
    if fname.startswith('INPUT'):
        fname = fname[5:].strip('_') # Strip "INPUT_" from the file if present
    fname = ''.join(fname.split('.')[:-1]) # Strip extension (will be re-added as appropriate)

    output_dir = os.path.join(directory, 'genetrack_s%de%d' % (sigma, exclusion))
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    output_path = os.path.join(output_dir, fname+'O_%s.txt' % fname)
    
    chromosomes = parse_reads(csv.reader(open(path,'rU'), delimiter='\t'))
    writer = csv.writer(open(output_path, 'wt'), delimiter='\t')
    writer.writerow(('chrom', 'strand', 'start', 'end', 'value'))
    
    for cname, data in chromosomes.items():
        logging.debug('Processing chromosome %s' % cname)
        data = chromosomes[cname]
        keys = make_keys(data)
        # Create the arrays that hold the sum of the normals
        forward_array = allocate_array(data, WIDTH)
        reverse_array = allocate_array(data, WIDTH)
        normal = normal_array(WIDTH, sigma)
        
        # Add each read's normal to the array
        for read in data:
            index, forward, reverse = read
            # Add the normals to the appropriate regions
            forward_array[index:index+WIDTH*2] += normal * forward
            reverse_array[index:index+WIDTH*2] += normal * reverse
            
        logging.debug('Calling forward strand')
        forward_peaks = call_peaks(forward_array, data, keys, 1, exclusion=exclusion)
        logging.debug('Calling reverse strand')
        reverse_peaks = call_peaks(reverse_array, data, keys, 2, exclusion=exclusion)
    
        
        for peak in forward_peaks:
            writer.writerow((cname, '+', peak.start, peak.end, peak.value))
        for peak in reverse_peaks:
            writer.writerow((cname, '-', peak.start, peak.end, peak.value))   
    

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
                if os.path.isfile(fpath) and not fname.startswith('.'): 
                    process_file(fpath, options.sigma, options.exclusion)
        else:
            process_file(path, options.sigma, options.exclusion)
            
if __name__ == '__main__':
    #run()
    import cProfile
    cProfile.run('run()', 'profile.bin')
            
