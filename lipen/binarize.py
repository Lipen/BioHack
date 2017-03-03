import csv
import time
from collections import namedtuple
import numpy as np
import pickle

Item = namedtuple('Item', ['chr', 'start', 'end', 'feature', 'score', 'strand'])

input_filename = '../data/dm6_rRNA.bed'
output_filename = 'C:/TMP/dm6_rRNA.pkl'


print('[*] Parsing <{}>...'.format(input_filename))
with open(input_filename) as f:
    reader = csv.reader(f, delimiter='\t')
    data = [Item(chr, int(start), int(end), feature, score, strand)
            for chr, start, end, feature, score, strand in reader]
print('[+] len(data) = {}'.format(len(data)))


print('[*] Processing...')
time_start = time.time()
window = int(np.median([item.end - item.start + 1 for item in data]))

print('[*] Binarizing with <window = {}>...'.format(window))
assert window > 20
# binarized :: {str: [bool]} == {chromosome: [binary]}
binarized = {
    'chr4': [False] * (1500000 // window),
    # 'chrM': [False] * (20000 // window),  # Mitochondrial, do not count
    'chrX': [False] * (24000000 // window),
    'chrY': [False] * (4000000 // window),
    'chr2L': [False] * (24000000 // window),
    'chr2R': [False] * (26000000 // window),
    'chr3L': [False] * (29000000 // window),
    'chr3R': [False] * (33000000 // window)
}

# DEBUG
chroms = len(set(item.chr for item in data))
print('CHROMOSOMS: {}'.format(chroms))
# DEBUG END

for item in data:
    low = item.start // window
    high = item.end // window
    for j in range(low, high + 1):
        if item.chr in binarized:
            binarized[item.chr][j] = True

print('[+] Done processing in {} seconds'.format(time.time() - time_start))

print('[*] Dumping to <{}>...'.format(output_filename))
with open(output_filename, 'wb') as f:
    pickle.dump(binarized, f, pickle.HIGHEST_PROTOCOL)
