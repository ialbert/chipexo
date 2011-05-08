from optparse import OptionParser, IndentedHelpFormatter
import re, sys, os

ROMAN = ['0', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X',
         'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'XVII', 'XVIII', 'XIX', 'XX',
         'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI', 'XXVII', 'XXVIII', 'XXIX', 'XXX']


def noop(data):
    return data

def zeropad_to_numeric(data):
    return re.sub(r'chr0(\d)', r'chr\1', data)
    
def numeric_to_zeropad(data):
    return re.sub(r'chr(\d([^\d]|$))', r'chr0\1', data)
    
def roman_to_numeric(data):
    for num, roman in enumerate(ROMAN):
        r = re.compile('chr%s([^IVX]|$)' % roman, flags=re.IGNORECASE)
        data = r.sub(r'chr%d\1' % num, data)
    return data

def numeric_to_roman(data):
    for num, roman in enumerate(ROMAN):
        data = re.sub('chr%d([^\d]|$)' % num, r'chr%s\1' % roman,  data)
    return data

FORMATS = ['zeropad', 'numeric', 'roman']
IN_CONVERT = {'zeropad':zeropad_to_numeric, 'roman':roman_to_numeric, 'numeric':noop}
OUT_CONVERT = {'zeropad':numeric_to_zeropad, 'roman':numeric_to_roman, 'numeric':noop}

def conversion_functions(in_fmt, out_fmt):
    ''' Returns the proper list of functions to apply to perform a conversion '''
    return [IN_CONVERT[in_fmt], OUT_CONVERT[out_fmt]]
    
def autodetect_format(data):
    if re.search('chr0\d', data):
        return 'zeropad'
    if re.search('chr[IVXivx]', data):
        return 'roman'
    return 'numeric'
    
def process_file(path, in_fmt, out_fmt):
    f = open(path, 'rt')
    data = f.read()
    if in_fmt == 'autodetect':
        in_fmt = autodetect_format(data)
    f.close()
    f = open(path, 'wt')
    for fn in conversion_functions(in_fmt, out_fmt):
        data = fn(data)
    f.write(data)
    f.close()
    
    
usage = '''
input_paths may be:
- a file or list of files to run on
- a directory or list of directories to run on all files in them
- "." to run in the current directory

formats are:
- numeric: "chr" followed by a number.                  ex chr4
- zeropad: "chr" followed by number padded to 2 places. ex chr04
- roman: "chr" followed by a Roman numeral.             ex chrIV
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description
 

def run():   
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-i', action='store', type='string', dest='in_format', default='autodetect',
                      help='Format input data is in. Default autodetect.')
    parser.add_option('-o', action='store', type='string', dest='out_format', default='numeric',
                      help='Format to output data in. Default numeric.')
    (options, args) = parser.parse_args()
        
    if options.in_format not in FORMATS+['autodetect']  or options.out_format not in FORMATS:
        parser.error('%s is not a valid method. Use -h option for a list of valid methods.' % options.method)
        
    if not args:
        parser.print_help()
        sys.exit(1)
        
    for path in args:
        if not os.path.exists(path):
            parser.error('Path %s does not exist.' % path)
        if os.path.isdir(path):
            for fname in os.listdir(path):
                fpath = os.path.join(path, fname)
                if os.path.isfile(fpath):
                    process_file(fpath, options.in_format, options.out_format)
        else:
            process_file(path, options.in_format, options.out_format)
    

if __name__ == '__main__':
    run()
    