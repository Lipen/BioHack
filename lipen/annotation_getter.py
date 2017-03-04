from collections import defaultdict

input_filename = ''
output_filename = 'annotation_{chromosome}.txt'

data = defaultdict(list)

with open(input_filename) as f:
    while True:
        line = f.readline()
        if line[0] == '#':
            continue
        break
    for line in f.readlines():
        chromo, _, gene, left, right, _, strand, _, extra = line.split('\t')
        gene_id = extra.split(';')[0][9:-1]
        if chromo.startswith('chr') and gene == 'gene':
            data[chromo + strand].append((gene_id, int(left), int(right)))

for chromo, meta in data.items():
    out = output_filename.format(chromosome=chromo)
    print('[*] Writing to <{}>...'.format(out))
    with open(out, 'w') as f:
        f.write('\t'.join(map(str, meta)) + '\n')
