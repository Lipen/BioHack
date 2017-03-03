1. Parse bed-file.
	Peaks ~> nametuple\<chr::str, start::int, end::int, feature::str, strand::('+' or '-')\>

2. Make some Panyushev tea

3. Binarize data as follows:
	* Peak width: `width = end-start+1`
	* Choose window size `k = median(widths)` (maybe twice as much)
	* Put a 1 in each window which contains any peak, and 0 otherwise
	* `~> list [ {0,1} ] of length N//k`

	* TODO: think about overlapping windows (slide 2*k window with 50% overlap with previous?)

4. Countize (clusterize) as follows:
	* Window size: large (around 1kk maybe) `w = 1_000_000`
	* Put a value at each cluster: amount of binary_peaks inside window
	* `~> list [ n_i ] of length M//w`

5. Remove clusters which contain 0 in all cell-lines (including validation dataset) -- they are reduntant.
	* (But if there is >0 in validation set in according cluster, then maybe it is very important... or it's just an small error...)
	* TODO: maybe pick some threshold (such as value = 0..5) to remove cluster
