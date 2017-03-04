import time
import math
import pickle

input_filename = 'C:/TMP/dm6_rRNA_bin_{chromosome}.pkl'
output_filename = 'C:/TMP/dm6_rRNA_clu_{chromosome}.pkl'
CHROMOLIST = ['chr4', 'chrX', 'chrY', 'chr2L', 'chr2R', 'chr3L', 'chr3R']
# input_filename = 'C:/TMP/H3K4me1.IMR90.mono_bin_{chromosome}.pkl'
# output_filename = 'C:/TMP/H3K4me1.IMR90.mono_clu_{chromosome}.pkl'
# CHROMOLIST = ['chr11', 'chr2', 'chrM', 'chr13', 'chrX', 'chr6', 'chr15', 'chr16', 'chr3', 'chr12', 'chr20', 'chr7', 'chr17', 'chr1', 'chr19', 'chr9', 'chr4', 'chr5', 'chrY', 'chr18', 'chr21', 'chr14', 'chr22', 'chr8', 'chr10']


def chunks(it, n):
    """Yield successive <n>-sized chunks from <it>."""
    for i in range(0, len(it), n):
        yield it[i:i + n]

print('[*] Working...')
time_start = time.time()

for chromo in CHROMOLIST:
    in_ = input_filename.format(chromosome=chromo)
    print('[*] Loading binarized dump <{}>...'.format(in_))
    try:
        with open(in_, 'rb') as f:
            binarized = pickle.load(f)
    except FileNotFoundError:
        print('[!] Okay, there is no such dump. Just continue...')
        continue

    n = len(binarized)
    k = 30  # number of clusters
    w = math.ceil(n / k)  # width of cluster
    clusters = []

    for chunk in chunks(binarized, w):
        clusters.append(sum(chunk))

    print('[+] <{}> clusters at <{}>:'.format(k, chromo))
    print('  > {}'.format(clusters))

    out = output_filename.format(chromosome=chromo)
    print('[*] Dumping to <{}>...'.format(out))
    with open(out, 'wb') as f:
        pickle.dump(clusters, f, pickle.HIGHEST_PROTOCOL)

print('[+] Done in {} seconds!'.format(time.time() - time_start))
