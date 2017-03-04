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
    Peak = namedtuple('Peak', ['start', 'end', 'score'])
    data = defaultdict(list)

    with open(input_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for chromo, start, end, stuff, score, strand in reader:
            if chromo + strand in CHROMOLIST:  # and int(score) > 100:
                data[chromo + strand].append(Peak(int(start), int(end), int(score)))

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

    # merged = merge(genes)
    # regions = get_gene_regions(merged)
    regions = get_gene_regions(genes)
    rebranded = rebrand(regions, meta)

    out = os.path.join(output_folder, '{}_reb_{}.pkl'.format(os.path.split(input_filename)[1], chromo))
    print('[*] Dumping to <{}>...'.format(out))
    with open(out, 'wb') as f:
        pickle.dump(rebranded, f, pickle.HIGHEST_PROTOCOL)

    # MAPPING
    out = os.path.join(output_folder, '{}_REGIONS_{}.pkl'.format(os.path.split(input_filename)[1], chromo))
    print('[*] Dumping to <{}>...'.format(out))
    with open(out, 'wb') as f:
        pickle.dump(regions, f, pickle.HIGHEST_PROTOCOL)


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


def merge(genes):
    print('[*] Merging <{}> genes...'.format(len(genes)))
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
    return merged


# def get_gene_regions(merged):
    # print('[*] Splitting <{}> "genes" into gene-regions...'.format(len(merged)))
    # regions = []
    # left = 0
    # for i in range(len(merged) - 1):
    #     x, y = merged[i], merged[i + 1]
    #     right = (x[2] + y[1]) // 2
    #     r = (left, right)
    #     regions.append(r)
    #     left = right
    # return regions
def get_gene_regions(genes):
    print('[*] Building regions around genes +-delta = <{}>'.format(delta))
    regions = [(g.left - delta, g.left + delta, g.gene_id) for g in genes]  # yes, BOTH from left edge!
    return regions


def rebrand(regions, meta):
    print('[*] Rebranding <{}> gene-regions...'.format(len(regions)))
    k = len(regions) + 1
    rebranded = np.zeros((k, 2), dtype=int)
    rls = [r[0] for r in regions]

    for peak in meta:
        i = bisect.bisect_left(rls, peak.start)

        if peak.score <= 3:
            rebranded[i, 0] += 1  # Bad
        else:
            rebranded[i, 1] += 1  # Good

    return rebranded


def main():
    print('[*] Parsing args...')
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', default='../input/*.bed', help='Input folder and files in glob format')
    parser.add_argument('-output', default='../output')
    parser.add_argument('-annotation', default='../data/annotation_{chromosome}.txt')
    parser.add_argument('-chromolist', default='chr1,chr2,chr3')
    parser.add_argument('-delta', default=1000)  # bases

    global input_folder, output_folder, annotations_filename, CHROMOLIST, delta
    args = parser.parse_args()
    input_folder = args.input
    output_folder = args.output
    annotations_filename = args.annotation
    CHROMOLIST_ = args.chromolist.split(',')
    CHROMOLIST = [c + '+' for c in CHROMOLIST_] + [c + '-' for c in CHROMOLIST_]
    print('CHROMOLIST: {}'.format(CHROMOLIST))
    delta = args.delta

    files = glob.glob(input_folder)
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
