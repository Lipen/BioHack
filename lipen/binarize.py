import csv
import time
import math
import pickle
from collections import namedtuple, defaultdict
import numpy as np

input_filename = '../data/dm6_rRNA.bed'
output_filename = 'C:/TMP/dm6_rRNA_bin_{chromosome}.pkl'
CHROMOLIST = ['chr4', 'chrX', 'chrY', 'chr2L', 'chr2R', 'chr3L', 'chr3R']
# input_filename = 'C:/TMP/H3K4me1.IMR90.mono.bed'
# output_filename = 'C:/TMP/H3K4me1.IMR90.mono_bin_{chromosome}.pkl'
# CHROMOLIST = ['chr11', 'chr2', 'chrM', 'chr13', 'chrX', 'chr6', 'chr15', 'chr16', 'chr3', 'chr12', 'chr20', 'chr7', 'chr17', 'chr1', 'chr19', 'chr9', 'chr4', 'chr5', 'chrY', 'chr18', 'chr21', 'chr14', 'chr22', 'chr8', 'chr10']

Item = namedtuple('Peak', ['chr', 'start', 'end', 'score'])
data = defaultdict(list)

print('[*] Parsing <{}>...'.format(input_filename))
time_start = time.time()
with open(input_filename) as f:
    reader = csv.reader(f, delimiter='\t')
    for chromo, start, end, stuff, score, strand in reader:
        if chromo in CHROMOLIST:
            data[chromo].append(Item(chromo, int(start), int(end), score))
print('[+] Parsed in {} seconds'.format(time.time() - time_start))

print('[*] Working...')
time_start = time.time()

for chromo in CHROMOLIST:
    meta = data[chromo]
    if len(meta):
        print('[*] Processing <{}> with {} peaks...'.format(chromo, len(meta)))
    else:
        print('[*] Skipping <{}> with no peaks'.format(chromo))
        continue

    # w = twice the median of peaks widths
    w = 2 * int(np.median([item.end - item.start + 1 for item in meta]))
    k = math.ceil(max(item.end for item in meta) / w)
    n = k * w
    print('  > Window width: <w = {}>'.format(w))
    print('  > Number of windows: <k = {}>'.format(k))
    print('  > Ceiled chromosome length: <n = {}>'.format(n))
    assert w > 40, "Maybe it would be tooooo large with such small window..."

    binarized = [False] * k

    for item in meta:
        # low/high indices of peaks covered by window
        low = item.start // w
        high = item.end // w
        for j in range(low, high + 1):
            try:
                binarized[j] = True
            except:
                print(low, high, item)
                raise

    out = output_filename.format(chromosome=chromo)
    print('[*] Dumping to <{}>...'.format(out))
    with open(out, 'wb') as f:
        pickle.dump(binarized, f, pickle.HIGHEST_PROTOCOL)

print('[+] Done in {} seconds!'.format(time.time() - time_start))
