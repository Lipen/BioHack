import os
import csv
import time
import glob
import pickle
from collections import namedtuple

input_filenames = '../data/*.tsv'
output_folder = 'C:/TMP/output'
annotations_filename = '../data/annotation_{chromosome}.txt'
CHROMOLIST = 'chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr2,chr20,chr21,chr22,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chrX'.split(',')


def parse(input_filename):
    print('[*] Parsing <{}>...'.format(input_filename))
    time_start = time.time()
    data = dict()

    with open(input_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for gene_id, log2_fold_change, _ in reader:
            try:
                log2_fold_change = float(log2_fold_change)
            except:
                log2_fold_change = 0
            data[gene_id] = log2_fold_change

    # print('[+] Parsed <{}> peaks in {} seconds'.format(sum(map(len, data.values())), time.time() - time_start))
    return data


def get_annotated_genes(chromo):
    anno = annotations_filename.format(chromosome=chromo)
    print('[*] Reading annotation <{}>...'.format(anno))
    Gene = namedtuple('Gene', ['gene_id', 'left', 'right'])
    genes = []

    with open(anno) as f:
        for line in f.readlines():
            gene_id, left, right = line.rstrip().split('\t')
            genes.append(Gene(gene_id, int(left), int(right)))

    print('[+] Annotated genes: <{}>'.format(len(genes)))
    return genes


genbase = {chromo: [x.gene_id for x in get_annotated_genes(chromo)] for chromo in CHROMOLIST}
for chromo in genbase:
    print('LEN({}) = {}'.format(chromo, len(genbase[chromo])))

for input_filename in glob.glob(input_filenames):
    data = parse(input_filename)
    for chromo in CHROMOLIST:
        meta = dict()
        for c, f in data.items():
            if c in genbase[chromo]:
                meta[c] = f
        out = 'DiffExpr_' + os.path.join(output_folder, os.path.splitext(os.path.split(input_filename)[1])[0].split('_')[-1]) + '_{}'.format(chromo) + '.pkl'
        with open(out, 'wb') as f:
            pickle.dump(meta, f)
