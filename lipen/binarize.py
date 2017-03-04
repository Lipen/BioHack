import csv
import time
import math
import pickle
from collections import namedtuple, defaultdict
import numpy as np

input_filename = '../data/dm6_rRNA.bed'
output_filename = 'C:/TMP/dm6_rRNA_bin_{chromosome}.pkl'

Item = namedtuple('Peak', ['chr', 'start', 'end', 'stuff', 'score', 'strand'])
CHROMOLIST = ['chr4', 'chrX', 'chrY', 'chr2L', 'chr2R', 'chr3L', 'chr3R']

data = defaultdict(list)

print('[*] Parsing <{}>...'.format(input_filename))
with open(input_filename) as f:
    reader = csv.reader(f, delimiter='\t')
    for chromo, start, end, stuff, score, strand in reader:
        if chromo in CHROMOLIST:
            data[chromo].append(Item(chromo, int(start), int(end),
                                     stuff, score, strand))

print('[*] Working...')
time_start = time.time()

for chromo in CHROMOLIST:
    meta = data[chromo]
    print('[*] Processing <{}> with {} peaks'.format(chromo, len(meta)))

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
