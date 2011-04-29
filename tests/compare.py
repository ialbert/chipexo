from optparse import OptionParser
import csv

def print_graph(graph):
	print ' '.join(['%3d' % x for x in graph]) + '   [%d]' % len(graph)

def compare(pathA, pathB, column=0, bins=20):
    dataA = []
    dataB = []
    for line in csv.reader(open(pathA,'rt'), delimiter='\t'):
        dataA.append(float(line[column]))
    for line in csv.reader(open(pathB,'rt'), delimiter='\t'):
        dataB.append(float(line[column]))
    lower = min(dataA+dataB)
    upper = max(dataA+dataB)
    width = (upper-lower)/bins
    graphA = [0] * bins
    graphB = [0] * bins
    for i in range(bins):
        graphA[i] = len([x for x in dataA if lower+width*i < x < lower+width*(i+1)])
        graphB[i] = len([x for x in dataB if lower+width*i < x < lower+width*(i+1)])
    print [int(lower+i*width) for i in range(bins)]
    print_graph(graphA)
    print_graph(graphB)
    print graphA == graphB
        
if __name__ == '__main__':
    parser = OptionParser(usage='%prog [options] fileA fileB')
    parser.add_option('-c', action='store', type='int', dest='columns', default=0)
    parser.add_option('-b', action='store', type='int', dest='bins', default=20)
    (options, args) = parser.parse_args()
    if len(args) == 0:
        parser.error('Must specify files')
    if len(args) > 2:
        parser.error('Too many arguments')
    if len(args) == 1:
        compare(args[0], args[0], options.columns, options.bins)
    else:
        pathA, pathB = args
        compare(pathA, pathB, options.columns, options.bins)
    