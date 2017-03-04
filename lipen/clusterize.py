import time
import math
import pickle


input_filename = 'C:/TMP/dm6_rRNA_bin_{chromosome}.pkl'
output_filename = 'C:/TMP/dm6_rRNA_clu_{chromosome}.pkl'

CHROMOLIST = ['chr4', 'chrX', 'chrY', 'chr2L', 'chr2R', 'chr3L', 'chr3R']


def chunks(it, n):
    """Yield successive <n>-sized chunks from <it>."""
    for i in range(0, len(it), n):
        yield it[i:i + n]

print('[*] Working...')
time_start = time.time()

for chromo in CHROMOLIST:
    in_ = input_filename.format(chromosome=chromo)
    print('[*] Loading dump <{}>...'.format(in_))
    with open(in_, 'rb') as f:
        binarized = pickle.load(f)

    n = len(binarized)
    k = 20  # number of clusters
    w = math.ceil(n / k)  # width of cluster
    clusters = []

    for chunk in chunks(binarized, w):
        clusters.append(sum(chunk))

    print('[+] Clusters at <{}>:'.format(chromo))
    print('  > {}'.format(clusters))

    out = output_filename.format(chromosome=chromo)
    print('[*] Dumping to <{}>...'.format(out))
    with open(out, 'wb') as f:
        pickle.dump(binarized, f, pickle.HIGHEST_PROTOCOL)

print('[+] Done in {} seconds!'.format(time.time() - time_start))
