BEGIN {
	FS = OFS = "\t"
	if (!minval) minval = 3
}

{
	if ($5 >= minval) {
		print $1, $3, $4, "peak", $5, $2
	}
}