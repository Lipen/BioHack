import csv
import time
import pickle
import bisect
from collections import namedtuple, defaultdict
import numpy as np

# input_filename = '../data/dm6_rRNA.bed'
input_filename = 'C:/TMP/human.txt'
annotations_filename = '../data/annotation_{chromosome}.txt'
# output_filename = 'C:/TMP/dm6_rRNA_reb_{xy}_{chromosome}.pkl'
output_filename = 'C:/TMP/human_reb_{chromosome}.pkl'
# CHROMOLIST_ = ['chr4', 'chrX', 'chrY', 'chr2L', 'chr2R', 'chr3L', 'chr3R']
CHROMOLIST_ = ['chr5', 'chr7']
CHROMOLIST = [c + '+' for c in CHROMOLIST_] + [c + '-' for c in CHROMOLIST_]

Item = namedtuple('Peak', ['start', 'end', 'score'])
data = defaultdict(list)

print('[*] Parsing <{}>...'.format(input_filename))
time_start = time.time()
with open(input_filename) as f:
    reader = csv.reader(f, delimiter='\t')
    for chromo, start, end, stuff, score, strand in reader:
        if chromo in CHROMOLIST_:  # and int(score) > 100:
            data[chromo + strand].append(Item(int(start), int(end), int(score)))
print('[+] Parsed <{}> peaks in {} seconds'.format(sum(map(len, data.values())), time.time() - time_start))

print('[*] Working...')
time_start = time.time()

for chromo in CHROMOLIST:
    meta = data[chromo]
    if len(meta):
        print('[*] Processing <{}> with {} peaks...'.format(chromo, len(meta)))
    else:
        print('[*] Skipping <{}> with no peaks'.format(chromo))
        continue

    anno = annotations_filename.format(chromosome=chromo)
    print('[*] Reading annotation <{}>...'.format(anno))
    Gene = namedtuple('Gene', ['gene_id', 'left', 'right'])
    genes = []
    with open(anno) as f:
        for line in f.readlines():
            gene_id, left, right = line.rstrip().split('\t')
            genes.append(Gene(gene_id, int(left), int(right)))
    print('[+] Annotated genes: <{}>'.format(len(genes)))

    print('[*] Merging <{}> genes...'.format(len(genes)))
    genes.sort(key=lambda g: (g.left, g.right))
    merged = []  # :: (id, left, right)
    for higher in genes:
        if not merged:
            merged.append((higher.gene_id, higher.left, higher.right))
        else:
            i, a, b = merged[-1]
            if higher.left <= a:
                merged[-1] = (i, a, max(b, higher.right))
            else:
                merged.append((i, higher.left, higher.right))

    print('[*] Splitting <{}> "genes" into gene-regions...'.format(len(merged)))
    regions = []
    left = 0
    for i in range(len(merged) - 1):
        x, y = merged[i], merged[i + 1]
        right = (x[2] + y[1]) // 2
        r = (left, right)
        regions.append(r)
        left = right

    print('[*] Rebranding <{}> gene-regions...'.format(len(regions)))
    k = len(regions) + 1
    rebranded = np.zeros(2, k)
    rls = [r[0] for r in regions]

    for peak in meta:
        i = bisect.bisect_left(rls, peak.start)

        if 0 < peak.score < 2:
            rebranded[0, i] += 1
        else:
            rebranded[1, i] += 1

    out = output_filename.format(chromosome=chromo)
    print('[*] Dumping to <{}>...'.format(out))
    with open(out, 'wb') as f:
        pickle.dump(rebranded, f, pickle.HIGHEST_PROTOCOL)

print('[+] Done in {} seconds!'.format(time.time() - time_start))
