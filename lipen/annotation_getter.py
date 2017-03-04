input_filename = ''
output_filename = 'annotation_{chromosome}.txt'

data = dict()

with open(input_filename) as f:
    while True:
        line = f.readline()
        if line[0] == '#':
            continue
        break
    while line:
        chromo, _, gene, left, right, _, strand, _, extra = line.split('\t')
        gene_id = extra.split(';')[0][9:-1]
        data[chromo + strand] = (gene, int(left), int(right), gene_id)

for chromo, meta in data.items():
    out = output_filename.format(chromosome=chromo)
    print('[*] Writing to <{}>...'.format(out))
    with open(out, 'w') as f:
        f.write('\t'.join(map(str, [strand] + list(meta))) + '\n')
