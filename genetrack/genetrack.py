"""
genetrack.py

Peak calling script

By Pindi Albert, 2011

Run with no arguments or -h for usage and command line options

Input: either idx, bed or gff format. Format will be autodetected.

Output: Called peaks in either gff

.txt format: tab-separated chromosome, strand, start, end, read count
.gff format: standard gff 9 fields
.bed format: 6 column bed with strand information
"""

from optparse import OptionParser, IndentedHelpFormatter
import csv, logging, numpy, math, bisect, sys, os, copy

logging.basicConfig(format='%(levelname)s:%(message)s')

WIDTH = 100

def gff_row(cname, start, end, score, source, type='.', strand='.', phase='.', attrs={}):
    global readsize
    if strand != '+':
        start += readsize
        end += readsize

    attrs['ID']= (start+end)/2
    return (cname, source, type, start, end, score, strand, phase, gff_attrs(attrs))

def gff_attrs(d):
    if not d:
        return '.'
    return ';'.join('%s=%s' % item for item in d.items())

class InvalidFileError(Exception):
    pass

class Peak(object):
    def __init__(self, index, pos_width, neg_width, height):
        self.index = index
        self.start = index - neg_width
        self.end = index + pos_width
        self.value = 0
        self.readcount = 0
        self.height = height
        self.deleted = False
        self.safe = False
    def __repr__(self):
        return '[%d] %d' % (self.index, self.value)

def is_int(i):
    try:
        int(i)
        return True
    except ValueError:
        return False
        
        
class ChromosomeManager(object):
    ''' Manages a CSV reader of an index file to only load one chrom at a time '''
    def __init__(self, reader):
        self.done = False
        self.reader = reader
        self.processed_chromosomes = []
        self.current_index = 0
        self.next_valid()
        
    def next(self):
        self.line = self.reader.next()
    
    def detect_format(self, line):
        pass
    
    def is_valid(self, line):
        global readsize
        if len(line) not in [4, 5, 6, 9]:
            return False
        
        # try genetrack index format
        try:
            [int(i) for i in line[1:]]
            self.format = 'idx'
            readsize = 0
            self.parse_line = self.parse_idx_line
            logging.debug('IDX format detected')
            return True
        except ValueError:
            # try next format
            pass
    
        
        # trying BED format
        try:
            start = int(line[1])
            end   = int(line[2])
            readsize = end - start
            self.format = 'bed'
            self.parse_line = self.parse_bed_line
            logging.debug('BED format detected')        
            return True
        except ValueError:
            # try next format
            pass
    
        # trying the GFF format
        try:
            start = int(line[3])
            end   = int(line[4])
            readsize = end - start + 1
            self.format = 'gff'
            self.parse_line = self.parse_gff_line
            logging.debug('GFF format detected')
            return True
        except ValueError:
            # try next format
            pass
        
        return False
        
    def next_valid(self):
        ''' Advance to the next valid line in the reader '''
        self.line = self.reader.next()
        s = 0
        while not self.is_valid(self.line):
            self.line = self.reader.next()
            s += 1
        if s > 0:
            logging.info('Skipped initial %d line(s) of file' % s)
            
    def parse_gff_line(self, line):
        return [int(line[3]), line[6], line[5]]

    def parse_idx_line(self, line):
        return [int(line[1]), int(line[2]), int(line[3])]
    
    def parse_bed_line(self, line):
        # turn it into one based interval
        return [int(line[1]) + 1, line[5], line[4]]
        
    def chromosome_name(self):
        ''' Return the name of the chromosome about to be loaded '''
        return self.line[0]
        
    def load_chromosome(self, collect_data=True):
        ''' Load the current chromosome into an array and return it '''
        cname = self.chromosome_name()
        if cname in self.processed_chromosomes:
            logging.error('File is not grouped by chromosome')
            raise InvalidFileError
        self.data = []

        while self.line[0] == cname:
            
            if collect_data:
                read = self.parse_line(self.line)
                if read[0] < self.current_index:
                    logging.error('Reads in chromosome %s are not sorted by index. (At index %d)' % (cname, self.current_index))
                    raise InvalidFileError
                self.current_index = read[0]
                self.add_read(read)
            try:
                self.next()
            except StopIteration:
                self.done = True
                break
        self.processed_chromosomes.append(cname)
        self.current_index = 0
        data = self.data
        del self.data # Don't retain reference anymore to save memory
        return data
    
    def add_read(self, read):
        if self.format == 'idx':
            self.data.append(read)
        else:
            index, strand, value = read
            '''
            if value == '' or value == '.':
                value = 1
            else:
                #value = int(value)
                value = 1
            '''
            value = 1
            if not self.data:
                self.data.append([index, 0, 0])
                current_read = self.data[-1]
            if self.data[-1][0] == index:
                current_read = self.data[-1]
            elif self.data[-1][0] < index:
                self.data.append([index, 0, 0])
                current_read = self.data[-1]
            else:
                logging.error('Reads in chromosome %s are not sorted by index. (At index %d)' % (self.chromosome_name(), index))
                raise InvalidFileError
            if strand == '+':
                current_read[1] += value
            elif strand == '-':
                current_read[2] += value
            else:
                logging.error('Strand "%s" at chromosome "%s" index %d is not valid.' % (strand, self.chromosome_name(), index))
                raise InvalidFileError
        
    
    def skip_chromosome(self):
        ''' Skip the current chromosome, discarding data '''
        self.load_chromosome(collect_data=False)
    
            
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

def get_range(data):
    lo = min([item[0] for item in data])
    hi = max([item[0] for item in data])
    return lo, hi

def get_chunks(lo, hi, size, overlap=500):
    ''' Divides a range into chunks of maximum size size. Returns a list of 2-tuples
    (slice_range, process_range), each a 2-tuple (start, end). process_range has zero overlap
    and should be given to process_chromosome as-is, and slice_range is overlapped and should be used to
    slice the data (using get_window) to be given to process_chromosome. '''
    chunks = []
    for start_index in range(lo, hi, size):
        process_start = start_index
        process_end = min(start_index + size, hi) # Don't go over upper bound
        slice_start = max(process_start - overlap, lo) # Don't go under lower bound
        slice_end = min(process_end + overlap, hi) # Don't go over upper bound
        chunks.append(((slice_start, slice_end), (process_start, process_end)))
    return chunks
    
def allocate_array(data, width):
    ''' Allocates a new array with the dimensions required to fit all reads in the
    argument. The new array is totally empty. Returns the array and the shift (number to add to
    a read index to get the position in the array it should be at).'''
    lo, hi = get_range(data)
    rng = hi - lo
    shift = width - lo
    return numpy.zeros(rng+width*2, numpy.float), shift
    
def normal_array(width, sigma, normalize=True):
    ''' Returns an array of the normal distribution of the specified width '''
    sigma2 = float(sigma)**2
    
    def normal_func(x):
        return math.exp( -x * x / ( 2 * sigma2 ))
        
    # width is the half of the distribution
    values = map( normal_func, range(-width, width) )
    values = numpy.array( values, numpy.float )

    # normalization
    #if normalize:
    #    values = 1.0/math.sqrt(2 * numpy.pi * sigma2) * values 

    return values

def call_peaks(array, shift, data, keys, direction, options):
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
            peaks.append(Peak(int(index)-shift, pos, neg, array[index]))
    find_peaks()
        
    def calculate_reads():
        # Calculate the number of reads in each peak
        for peak in peaks:
            reads = get_window(data, peak.start, peak.end, keys)
            peak.value = peak.height
            peak.readcount = sum([read[direction] for read in reads])
            indexes = [r for read in reads for r in [read[0]] * read[direction]] # Flat list of indexes with frequency
            peak.stddev = numpy.std(indexes)
    calculate_reads()
        
    before = len(peaks)
        
    def perform_exclusion():
        # Process the exclusion zone
        peak_keys = make_peak_keys(peaks)
        peaks_by_value = peaks[:]
        peaks_by_value.sort(key=lambda peak: -peak.value)
        for peak in peaks_by_value:
            peak.safe = True
            window = get_window(peaks, peak.index-options.exclusion//2, peak.index+options.exclusion//2, peak_keys)
            for excluded in window:
                if excluded.safe:
                    continue
                i = get_index(excluded.index, peak_keys)
                del peak_keys[i]
                del peaks[i]
    perform_exclusion()
            
    after = len(peaks)
    if before != 0:
        logging.debug('%d of %d peaks (%d%%) survived exclusion' % (after, before, after*100/before))
            
    return peaks
    
def process_chromosome(cname, data, process_bounds, options):
    '''
    Process a chromosome. Takes the chromosome name, list of reads,
    the bounds (2-tuple) to write results in, and options.
    '''
    if data:
        logging.info('Processing chromosome %s indexes %d-%d' % (cname, process_bounds[0], process_bounds[1]))
    else:
        logging.info('Skipping chromosome %s indexes %d-%d because no reads within these bounds' % (cname, process_bounds[0], process_bounds[1]))
        return
    keys = make_keys(data)
    # Create the arrays that hold the sum of the normals
    forward_array, forward_shift = allocate_array(data, WIDTH)
    reverse_array, reverse_shift = allocate_array(data, WIDTH)
    normal = normal_array(WIDTH, options.sigma)
    
    
    def populate_array():
        # Add each read's normal to the array
        for read in data:
            index, forward, reverse = read
            # Add the normals to the appropriate regions
            if forward:
                forward_array[index+forward_shift-WIDTH:index+forward_shift+WIDTH] += normal * forward
            if reverse:
                reverse_array[index+reverse_shift-WIDTH:index+reverse_shift+WIDTH] += normal * reverse
    populate_array()

    if options.bedgraph:
        logging.debug('Generating bedgraph files')
        f = open('forward.bedgraph','at')
        indexes = numpy.where(forward_array >= 0.5)[0]
        for index in indexes:
            f.write('%s\t%d\t%d\t%s\n' % (cname, index - forward_shift, index - forward_shift + 1, forward_array[index]))
        f.close()

        f = open('reverse.bedgraph','at')
        indexes = numpy.where(reverse_array > 0.5)[0]
        for index in indexes:
            pos = index - reverse_shift + readsize
            f.write('%s\t%d\t%d\t%s\n' % (cname, pos, pos + 1, reverse_array[index]))
        f.close()

    logging.debug('Calling forward strand')
    forward_peaks = call_peaks(forward_array, forward_shift, data, keys, 1, options)
    logging.debug('Calling reverse strand')
    reverse_peaks = call_peaks(reverse_array, reverse_shift, data, keys, 2, options)

    def write(cname, strand, peak):
        start = max(peak.start, 1)
        end = peak.end
        value = peak.value
        stddev = peak.stddev
        if value > options.filter:
            if options.format == 'gff':
                writerow(gff_row(cname=cname, source='genetrack', start=start, end=end,
                                        score=value, strand=strand, attrs={'stddev':stddev, 'height':peak.height, 'readcount':peak.readcount}))
            else:
                writerow((cname, strand, start, end, value))
    
    for peak in forward_peaks:
        if process_bounds[0] < peak.index < process_bounds[1]:
            write(cname, '+', peak)
    for peak in reverse_peaks:
        if process_bounds[0] < peak.index < process_bounds[1]:
            write(cname, '-', peak)

def writerow( elems ):
    print "\t".join( map(str, elems) )
    
def process_file(path, options):
    
    global WIDTH
    WIDTH = options.sigma * 4
    
    logging.info('Processing file "%s" with s=%d, e=%d' % (path, options.sigma, options.exclusion))

    if options.bedgraph:
        f = open('forward.bedgraph','wt')
        f.write("track type=bedGraph color=200,0,0\n")
        f.close()
            
        f = open('reverse.bedgraph','wt')
        f.write("track type=bedGraph color=0,0,200\n")
        f.close()

    if path == '-':
        reader = csv.reader(sys.stdin, delimiter='\t')
    elif os.path.exists(path):
        reader = csv.reader(open(path,'rt'), delimiter='\t')
    else:
        logging.error('Path "%s" does not exist.' % path)
        return

    if options.format == 'idx':
        writerow(('chrom', 'strand', 'start', 'end', 'value'))
    else:
        writerow(['##gff-version 3'])
    
    manager = ChromosomeManager(reader)
     
    while not manager.done:
        cname = manager.chromosome_name()
        if not options.chromosome or options.chromosome == cname: # Should we process this chromosome?
            logging.info('Loading chromosome %s' % cname)
            data = manager.load_chromosome()
            if not data:
                continue
            keys = make_keys(data)
            lo, hi = get_range(data)
            for chunk in get_chunks(lo, hi, size=options.chunk_size * 10 ** 6, overlap=WIDTH):
                (slice_start, slice_end), process_bounds = chunk
                window = get_window(data, slice_start, slice_end, keys)
                process_chromosome(cname, window, process_bounds, options)
            if options.chromosome: # A specific chromosome was specified. We're done, so terminate
                break
            #process_chromosome(cname, list(data), writer, options)
        else:
            logging.info('Skipping chromosome %s' % cname)
            manager.skip_chromosome()
    
usage = '''
input_paths may be:

    - a file to run on
    - "-" to run on standard input

example usage:

    python genetrack.py -s 10 /path/to/a/file.txt
    python genetrack.py -s 5 -e 50 -
'''.lstrip()
 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

def run():   
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-s', action='store', type='int', dest='sigma', default=5,
                      help='Sigma to use when smoothing reads to call peaks. Default %default')
    parser.add_option('-e', action='store', type='int', dest='exclusion', default=20,
                      help='Exclusion zone around each peak that prevents others from being called. Default %default.')
    parser.add_option('-u', action='store', type='int', dest='up_width', default=0,
                      help='Upstream width of called peaks. Default uses half exclusion zone.')
    parser.add_option('-d', action='store', type='int', dest='down_width', default=0,
                      help='Downstream width of called peaks. Default uses half exclusion zone.')
    parser.add_option('-F', action='store', type='int', dest='filter', default='3',
                      help='Absolute read filter; outputs only peaks with larger peak height. Default %default. ')
    parser.add_option('-c', action='store', type='string', dest='chromosome', default='',
                      help='Chromosome (ex chr11) to limit to. Default process all.')
    parser.add_option('-k', action='store', type='int', dest='chunk_size', default=10,
                      help='Size, in millions of base pairs, to chunk each chromosome into when processing. Each 1 million size uses approximately 20MB of memory. Default %default.')
    parser.add_option('-o', action='store', type='string', dest='format', default='gff',
                      help='Output format for called peaks. Valid formats are gff and txt. Default %default.')
    parser.add_option('-b', action='store_true', dest='bedgraph', help='Output bed graph tracks.')

    parser.add_option('-v', action='store_true', dest='verbose', help='Verbose mode: displays debug messages.')
 
    (options, args) = parser.parse_args()
    
    if options.verbose:
        logging.getLogger().setLevel(logging.DEBUG) # Show all info/debug messages
    else:
        logging.getLogger().setLevel(logging.ERROR)

    if options.format not in ['gff','txt']:
        parser.error('%s is not a valid format. Use -h option for a list of valid methods.' % options.format)
                
    if not args:
        parser.print_help()
        sys.exit(1)

    if len(args) > 1:
        parser.error('Only one file can be processed at a time.')
        
    path = args[0]

    process_file(path, options)

            
if __name__ == '__main__':
    run()
    #import cProfile
    #cProfile.run('run()', 'profilev8.bin')
            
