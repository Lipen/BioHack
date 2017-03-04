import csv
import time
import pickle
import bisect
from collections import namedtuple, defaultdict

input_filename = '../data/dm6_rRNA.bed'
annotations_filename = '../data/annotation_{chromosome}.txt'
output_filename = 'C:/TMP/dm6_rRNA_reb_{xy}_{chromosome}.pkl'
CHROMOLIST_ = ['chr4', 'chrX', 'chrY', 'chr2L', 'chr2R', 'chr3L', 'chr3R']
CHROMOLIST = [c + '+' for c in CHROMOLIST_] + [c + '-' for c in CHROMOLIST_]

Item = namedtuple('Peak', ['start', 'end', 'score'])
data = defaultdict(list)

print('[*] Parsing <{}>...'.format(input_filename))
time_start = time.time()
with open(input_filename) as f:
    reader = csv.reader(f, delimiter='\t')
    for chromo, start, end, stuff, score, strand in reader:
        if chromo in CHROMOLIST:  # and int(score) > 100:
            data[chromo].append(Item(int(start), int(end), int(score)))
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

    anno = annotations_filename.format(chromosome=chromo)
    print('[*] Reading annotation <{}>...'.format(anno))
    Gene = namedtuple('Gene', ['gene_id', 'left', 'right'])
    genes = []
    with open(anno) as f:
        gene_id, left, right = f.readline().rstrip().split('\t')
        genes.append(Gene(gene_id, int(left), int(right)))
    print('[+] Annotated genes: <{}>'.format(len(genes)))

    print('[*] Merging <{}> gene regions...'.format(len(genes)))
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
    print('[+] Merged into <{}> regions'.format(len(merged)))

    k = len(merged)
    rebranded_good = [0] * k
    rebranded_bad = [0] * k

    for peak in meta:
        i = bisect.bisect_left([item.left for item in merged], peak.start)

        if 0 < peak.score < 2:
            rebranded_bad[i] += 1
        else:
            rebranded_good[i] += 1

    out = output_filename.format(xy='good', chromosome=chromo)
    print('[*] Dumping to <{}>...'.format(out))
    with open(out, 'wb') as f:
        pickle.dump(rebranded_good, f, pickle.HIGHEST_PROTOCOL)
    out = output_filename.format(xy='bad', chromosome=chromo)
    print('[*] Dumping to <{}>...'.format(out))
    with open(out, 'wb') as f:
        pickle.dump(rebranded_bad, f, pickle.HIGHEST_PROTOCOL)
