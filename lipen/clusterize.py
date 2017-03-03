import pickle

input_filename = 'C:/TMP/dm6_rRNA.pkl'


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

print('[*] Loading dump <{}>...'.format(input_filename))
with open(input_filename, 'rb') as f:
    binarized = pickle.load(f)

print('[*] Processing with progressive window size...')
meta = {chrom: [] for chrom in binarized}

for chromosome, binary in binarized.items():
    n = len(binary)
    t = 20  # clusters amount for each chromosome
    w = n // t
    print('  > {}: w={}'.format(chromosome, w))
    for cluster in chunks(binary, w):
        meta[chromosome].append(sum(cluster))
        if meta[chromosome][-1] == 0 and len(meta[chromosome]) == 10:
            print('ZERO AT CLUSTER {}...'.format(cluster[:100]))

print('[+] Clusters:')
for chromosome, cluster in meta.items():
    print('  > {} :: {}'.format(chromosome, cluster))
