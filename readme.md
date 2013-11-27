
This is a standalone version of GeneTrack.
It was written by a Penn State freshman student Pindi Albert during an internship with the Pugh lab at Penn State.

It implements the peak calling algorithm as described
in the paper "GeneTrack--a genomic data processing and visualization framework" published in
Bioinformatics 2008, http://www.ncbi.nlm.nih.gov/pubmed/18388141 It does not however implement
the visualization interface that was described in the paper.

Usage:

    python genetrack/genetrack.py

This will print a help on usage.

The input files should be in `BED`, `GFF` or the internal `.idx` format.

Detailed usage:

    Usage: genetrack.py [options] input_paths

	input_paths may be:

		- a file to run on
		- "-" to run on standard input

	example usage:

		python genetrack.py -s 10 /path/to/a/file.txt
		python genetrack.py -s 5 -e 50 -

	Options:
	  -h, --help     show this help message and exit
	  -s SIGMA       Sigma to use when smoothing reads to call peaks. Default 5
	  -e EXCLUSION   Exclusion zone around each peak that prevents others from
					 being called. Default 20.
	  -u UP_WIDTH    Upstream width of called peaks. Default uses half exclusion
					 zone.
	  -d DOWN_WIDTH  Downstream width of called peaks. Default uses half exclusion
					 zone.
	  -F FILTER      Absolute read filter; outputs only peaks with larger peak
					 height. Default 3.
	  -c CHROMOSOME  Chromosome (ex chr11) to limit to. Default process all.
	  -k CHUNK_SIZE  Size, in millions of base pairs, to chunk each chromosome
					 into when processing. Each 1 million size uses approximately
					 20MB of memory. Default 10.
	  -o FORMAT      Output format for called peaks. Valid formats are gff and
					 txt. Default gff.
	  -b             Output bed graph tracks.
	  -v             Verbose mode: displays debug messages.


