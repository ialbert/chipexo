from optparse import OptionParser, IndentedHelpFormatter
import csv, logging, numpy, math, bisect, sys, os, copy

from chrtrans import convert_data

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

WIDTH = 100

class Peak(object):
    def __init__(self, index, pos_width, neg_width):
        self.index = index
        self.start = index - neg_width
        self.end = index + pos_width
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

class CSVReaderWrapper(object):
    ''' A wrapper around a CSV reader that can rewind one line '''
    def __init__(self, reader):
        last = None
        rewound = False
    def rewind(self):
        rewound = True
    def next(self):
        if rewound:
            return last
        else:
            return aaa
       
def is_int(i):
    try:
        int(i)
        return True
    except ValueError:
        return False
       
def is_valid(line):
    if len(line) != 5:
        return False
    try:
        [int(i) for i in line[1:]]
        return True
    except ValueError:
        return False
       
def next_valid(reader):
    LINE = reader.next()
    s = 0
    while not is_valid(LINE):
        LINE = reader.next()
        s += 1
    if s > 0:
        logging.info('Skipped %d line(s) of file' % s)
    return LINE
       
def parse_line(line):
    cname, index, forward, reverse = line[:4]
    return [int(index), int(forward), int(reverse)]

LINE = None
current_chromosome = None
STOP = False


def chromosome_iterator(reader):
    global LINE, STOP, current_chromosome
    chromosomes = []
    def reads_iterator():
        global LINE, STOP, current_chromosome
        while True:
            try:
                LINE = next_valid(reader)
            except StopIteration:
                STOP = True
                raise
            if LINE[0] != current_chromosome:
                current_chromosome = LINE[0]
                return
            yield parse_line(LINE)
    LINE = next_valid(reader)
    current_chromosome = LINE[0]
    while not STOP:
        yield current_chromosome, reads_iterator()
                


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

def call_peaks(array, data, keys, direction, options):
    peaks = []
    def find_peaks():
        # Go through the array and call each peak
        results = (array > numpy.roll(array, 1)) & (array > numpy.roll(array, -1))
        indexes = numpy.where(results)
        for index in indexes[0]:
            pos = options.down_width or options.exclusion // 2
            neg = options.up_width or options.exclusion // 2
            if direction == 2: # Reverse strand
                pos, neg = neg, pos # Swap positive and negative widths
            peaks.append(Peak(int(index)-WIDTH, pos, neg))
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
            window = get_window(peaks, peak.index-options.exclusion, peak.index+options.exclusion, peak_keys)
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
    
def process_chromosome(cname, data, writer, options):
    
    logging.debug('Processing chromosome %s' % cname)
    keys = make_keys(data)
    # Create the arrays that hold the sum of the normals
    forward_array = allocate_array(data, WIDTH)
    reverse_array = allocate_array(data, WIDTH)
    normal = normal_array(WIDTH, options.sigma)
    
    
    def populate_array():
        # Add each read's normal to the array
        for read in data:
            index, forward, reverse = read
            # Add the normals to the appropriate regions
            if forward:
                forward_array[index:index+WIDTH*2] += normal * forward
            if reverse:
                reverse_array[index:index+WIDTH*2] += normal * reverse
    populate_array()
        
    logging.debug('Calling forward strand')
    forward_peaks = call_peaks(forward_array, data, keys, 1, options)
    logging.debug('Calling reverse strand')
    reverse_peaks = call_peaks(reverse_array, data, keys, 2, options)

    # Convert chromosome name in preparation for writing our
    cname = convert_data(cname, 'zeropad', 'numeric')
    
    
    for peak in forward_peaks:
        writer.writerow((cname, '+', peak.start, peak.end, peak.value))
    for peak in reverse_peaks:
        writer.writerow((cname, '-', peak.start, peak.end, peak.value))   
    
    
    
def get_output_path(input_path, options):
    directory, fname = os.path.split(input_path)
    
    if fname.startswith('INPUT'):
        fname = fname[5:].strip('_') # Strip "INPUT_" from the file if present
    fname = ''.join(fname.split('.')[:-1]) # Strip extension (will be re-added as appropriate)

    output_dir = os.path.join(directory, 'genetrack_s%de%d' % (options.sigma, options.exclusion))
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if options.chromosome:
        fname = options.chromosome + '_' + fname
    return os.path.join(output_dir, '%s_s%de%d.txt' % (fname, options.sigma, options.exclusion))
    
    
def process_file(path, options):
    
    global WIDTH
    WIDTH = options.sigma * 5
    
    logging.info('Processing file "%s" with s=%d, e=%d' % (path, options.sigma, options.exclusion))
    
    output_path = get_output_path(path, options)
    
    reader = csv.reader(open(path,'rU'), delimiter='\t')
    #chromosomes = parse_reads(reader)
    writer = csv.writer(open(output_path, 'wt'), delimiter='\t')
    writer.writerow(('chrom', 'strand', 'start', 'end', 'value'))
    
    for cname, data in chromosome_iterator(reader):
        if not options.chromosome or options.chromosome == cname:
            process_chromosome(cname, list(data), writer, options)
        else:
            list(data) # Consume iterator even if not used
    

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
    parser.add_option('-u', action='store', type='int', dest='up_width', default=0,
                      help='Upstream width of called peaks. Default uses half exclusion zone.')
    parser.add_option('-d', action='store', type='int', dest='down_width', default=0,
                      help='Downstream width of called peaks. Default uses half exclusion zone.')
    parser.add_option('-c', action='store', type='string', dest='chromosome', default='',
                      help='Chromosome (ex chr11) to limit to. Default process all.')
    parser.add_option('-f', action='store', type='string', dest='config_file', default='',
                      help='Optional file to load sigma and exclusion parameters per input file.')
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
        
    CONFIG = {}
    if options.config_file:
        logging.info('Loading configuration file "%s"' % options.config_file)
        reader = csv.reader(open(options.config_file, 'rt'), delimiter='\t')
        for line in reader:
            if len(line) == 3 and is_int(line[1]) and is_int(line[2]):
                CONFIG[line[0]] = {'sigma':int(line[1]), 'exclusion':int(line[2])}
        logging.info('Loaded configuration settings for %d files.' % len(CONFIG))
        
    for path in args:
        if not os.path.exists(path):
            parser.error('Path %s does not exist.' % path)
        if os.path.isdir(path):
            files = []
            for fname in os.listdir(path):
                fpath = os.path.join(path, fname)
                if os.path.isfile(fpath) and not fname.startswith('.'): 
                    files.append(fpath)
        else:
            files = [path]
        for fpath in files:
            dir, fname = os.path.split(fpath)
            if fname in CONFIG:
                current_options = copy.deepcopy(options)
                current_options.sigma = CONFIG[fname]['sigma']
                current_options.exclusion = CONFIG[fname]['exclusion']
            else:
                if options.config_file:
                    logging.warning('File "%s" not found in config file' % fname)
                current_options = options
            process_file(fpath, current_options)
            
if __name__ == '__main__':
    #reader = csv.reader(open('data/INPUT_genetrack_Reb1_rep2.idx', 'rU'), delimiter='\t')
    #it = chromosome_iterator(reader)
    #for cname, data in it:
    #    print cname, len(list(data))
    run()
    #import cProfile
    #cProfile.run('run()', 'profilev6.bin')
            
