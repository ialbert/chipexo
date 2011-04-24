""" Test harness for testing scripts """

from config import SCRIPTS
from optparse import OptionParser
import subprocess, shutil


def run_script(script, inpath, outpath):
    command = 'python %s.py data/%s data/temp.txt' % (script, inpath)
    p = subprocess.Popen(command, shell=True)
    p.wait() # Wait for script to complete
    expected_output = open('data/%s' % outpath, 'rt').read()
    actual_output = open('data/temp.txt', 'rt').read()
    if expected_output != actual_output:
        actual_path = 'data/actual-%s' % outpath
        shutil.copy('data/temp.txt', actual_path)
        print '%s FAILED. Actual output saved to %s' % (script, actual_path)
    else:
        print '%s passed.' % script


if __name__ == '__main__':
    parser = OptionParser(usage='%prog [options] [script]')
    (options, args) = parser.parse_args()
    if len(args) > 1:
        parser.error('Too many arguments')
    if len(args) == 1: # Run only the script provided
        script = args[0]
        if script in SCRIPTS:
            (infile, outfile) = SCRIPTS[script]
            run_script(script, infile, outfile)
        else:
            parser.error('No tests configured for %s' % script)
    else: # Run all scripts
        for script, (infile, outfile) in SCRIPTS.items():
            run_script(script, infile, outfile)