import csv, logging, numpy, math

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

WIDTH = 100

def parse_reads(reader):
    chromosomes = {}
    reader.next()
    for line in reader:
        cname, index, forward, reverse = line[:4]
        if cname not in chromosomes:
            chromosomes[cname] = []
        chromosomes[cname].append([int(index), int(forward), int(reverse)])
    return chromosomes

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
        
chromosomes = parse_reads(csv.reader(open('data/INPUT_genetrack_Reb1_rep2.idx','rU'), delimiter='\t'))

for chromosome, data in chromosomes.items():
    data = chromosomes['chr01']
    # Create the arrays that hold the sum of the normals
    forward_array = allocate_array(data, WIDTH)
    #reverse_array = allocate_array(data, WIDTH)
    base_normal = normal_array(WIDTH, 10)
    
    # Add each read's normal to the array
    for read in data:
        index, forward, reverse = read
        normal = base_normal * forward
        forward_array[index:index+WIDTH*2] += normal # Add the normal to the appropriate region    
    
    peaks = []
    history = [forward_array[0], forward_array[1]] # Last 2 values
    # Go through the array and call each peak
    for i, value in enumerate(forward_array[2:]):
        if history[1] > history[0] and history[1] > value:
            # History[1] is a peak
            peaks.append([i+1-WIDTH, history[1]])
        history[0] = history[1]
        history[1] = value
    
    print peaks[:100]
    break
