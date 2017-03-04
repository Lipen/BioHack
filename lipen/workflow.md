### Possible workflow:

1. Parse bed-file.
	Peaks ~> nametuple\<chr::str, start::int, end::int, feature::str, strand::('+' or '-')\>

	* TODO: Maybe reparse file for each chromosome, because all-chromosomes-data can not to fit into RAM

	* TODO: Maybe make robust available chromosomes getter...?

	> Done

2. Binarize data as follows:

	* Peak width: `width = end-start+1`

	For each chromosome do:

	* Choose window size `w = 2 * median(widths)`
	* Put a `1` in each window which contains any peak, and `0` otherwise
	* ~> `list [ {0,1} ]` of length `k=ceil(N/w)`

	* TODO: think about overlapping windows (slide 2*k window with 50% overlap with previous?)

	* TODO: maybe put `1/r` in each window, where `r` = number of windows where according peak lies (the cake is a lie)

	* TODO: do not count peaks with low height (need robust progressive threshold...)

	> Done

3. Countize (clusterize) as follows:
	* Pick number of windows: `k = 20` or more
	* Window size: `w = ceil(n / k)`

	For each chromosome do:

	* Slide a window and put value at each cluster: amount of binary_peaks inside window
	* ~> `list [ m_i ]` of length `k`

	> Done

4. Drink some coffee.

	> Done

5. Remove clusters which contain 0 in all cell-lines (including validation dataset) -- they are reduntant.
	* (But if there is >0 in validation set in according cluster, then maybe it is very important... or it's just an small error...)
	* TODO: maybe pick some threshold (such as value = 0..5) to remove cluster

6. Repeat for all cell-lines


### Sh~t. It doesnt work.

1. Split genome (each chromosome) into regions around (fully) every gene of given 20k.

2. For each region count number of good and bad peaks from input datasets -- Chip-seq histone, Chip-seq TF
	*
