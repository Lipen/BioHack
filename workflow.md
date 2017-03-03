1. Parse bed-file. Peaks -> nametuple<chr::str, start::int, end::int, feature::str, strand::('+' or '-')>

2. Make some Panyushev tea

3. Binarize data as follows:
	* Choose window size = median of peaks widths
	* Put a 1 in any window contains any peak, and 0 otherwise
