import csv, logging, numpy, math, bisect

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

WIDTH = 100

class Peak(object):
    def __init__(self, index, value):
        self.index = index
        self.start = index - 10
        self.end = index + 10
        self.value = value
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

def call_peaks(array, keys, direction, width=WIDTH, exclusion=0):
    # Go through the array and call each peak
    peaks = []
    history = [array[0], array[1]] # Last 2 values
    for i, value in enumerate(array[2:]):
        if history[1] > history[0] and history[1] > value:
            # history[1] is a peak
            peaks.append(Peak(i+1-WIDTH, history[1]))
        history[0] = history[1] # Shift history values
        history[1] = value
        
    # Calculate the number of reads in each peak
    for peak in peaks:
        reads = get_window(data, peak.start, peak.end, keys)
        peak.value = sum([read[direction] for read in reads])
        
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
    return peaks
    
        
chromosomes = parse_reads(csv.reader(open('data/INPUT_genetrack_Reb1_rep2.idx','rU'), delimiter='\t'))
writer = csv.writer(open('data/OUTPUT.txt', 'wt'), delimiter='\t')
writer.writerow(('chrom', 'strand', 'start', 'end', 'value'))

for cname, data in chromosomes.items():
    cname = 'chr01'
    data = chromosomes[cname]
    keys = make_keys(data)
    # Create the arrays that hold the sum of the normals
    forward_array = allocate_array(data, WIDTH)
    reverse_array = allocate_array(data, WIDTH)
    normal = normal_array(WIDTH, 10)
    
    # Add each read's normal to the array
    for read in data:
        index, forward, reverse = read
        # Add the normals to the appropriate regions
        forward_array[index:index+WIDTH*2] += normal * forward
        reverse_array[index:index+WIDTH*2] += normal * reverse
        
    forward_peaks = call_peaks(forward_array, keys, 1)
    reverse_peaks = call_peaks(reverse_array, keys, 2)

    
    for peak in forward_peaks:
        writer.writerow((cname, '+', peak.start, peak.end, peak.value))
    for peak in reverse_peaks:
        writer.writerow((cname, '-', peak.start, peak.end, peak.value))   

    print forward_peaks[:100]
    break
