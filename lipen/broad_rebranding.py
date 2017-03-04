import os
import csv
import time
import glob
import pickle
import bisect
import argparse
from collections import namedtuple, defaultdict
import numpy as np


def parse(input_filename):
    print('[*] Parsing <{}>...'.format(input_filename))
    time_start = time.time()
    Peak = namedtuple('Peak', ['start', 'end', 'enrich'])
    data = defaultdict(list)

    with open(input_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for chromo, start, end, stuff, enrich, *_ in reader:
            if chromo in CHROMOLIST:
                data[chromo].append(Peak(int(start), int(end), float(enrich)))

    print('[+] Parsed <{}> peaks in {} seconds'.format(sum(map(len, data.values())), time.time() - time_start))
    return data


def work(data, input_filename):
    print('[*] Working...')
    time_start = time.time()

    for chromo in CHROMOLIST:
        process_chromosome(data, chromo, input_filename)

    print('[+] Done in {} seconds!'.format(time.time() - time_start))


def process_chromosome(data, chromo, input_filename):
    meta = data[chromo]
    if len(meta):
        print('[*] Processing <{}> with {} peaks...'.format(chromo, len(meta)))
    else:
        print('[*] Skipping <{}> with no peaks'.format(chromo))
        return

    genes = get_annotated_genes(chromo)
    genes.sort(key=lambda g: (g.left, g.right))

    lefts = [g.left - delta for g in genes]
    rebranded = rebrand(lefts, meta)  # {k x 3}
    blackhole = [(r, g.gene_id, g.left - delta, g.left + delta) for r, g in zip(rebranded, genes)]

    name = os.path.splitext(os.path.split(input_filename)[1])[0]
    out = os.path.join(output_folder, '{}_{chromosome}_d{delta}.pkl'.format(name, chromosome=chromo, delta=delta))
    print('[*] Dumping to <{}>...'.format(out))
    with open(out, 'wb') as f:
        pickle.dump(blackhole, f, pickle.HIGHEST_PROTOCOL)


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


def rebrand(lefts, meta):
    print('[*] Rebranding <{}> gene-delta(<{}>)-neighbours...'.format(len(lefts), delta))
    k = len(lefts) + 1
    rebranded = np.zeros((k, 3), dtype=int)
    meta_sorted = sorted(meta, key=lambda x: x.enrich)
    t = len(meta_sorted) / 10
    q10 = meta_sorted[int(t)]
    q90 = meta_sorted[int(t * 9)]

    for peak in meta:
        i = bisect.bisect_left(lefts, peak.start)

        if peak.enrich <= q10:
            rebranded[i, 0] += 1  # xi, bad
        elif peak.enrich <= q90:
            rebranded[i, 1] += 1  # yi, good
        else:
            rebranded[i, 2] += 1  # zi, idk

    return rebranded


def main():
    print('[*] Parsing args...')
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', default='../input/*.bed', help='Input folder and files in glob format')
    parser.add_argument('-output', default='../output')
    parser.add_argument('-annotation', default='../data/annotation_{chromosome}.txt')
    parser.add_argument('-chromolist', default='chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr2,chr20,chr21,chr22,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chrX,chrY')
    parser.add_argument('-delta', type=int, default=1000)  # bases

    global output_folder, annotations_filename, CHROMOLIST, delta
    args = parser.parse_args()
    input_glob = args.input
    output_folder = args.output
    annotations_filename = args.annotation
    CHROMOLIST = args.chromolist.split(',')
    delta = args.delta

    files = glob.glob(input_glob)
    if files:
        print('[*] Found input files:')
        for f in files:
            print('  > {}'.format(f))
    else:
        print('[!] There is no input files :c')
    for input_filename in files:
        data = parse(input_filename)
        work(data, input_filename)

if __name__ == '__main__':
    main()
