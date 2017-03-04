import os
import time
import glob
import pickle
import argparse
from collections import defaultdict


def main():
    print('[*] Parsing args...')
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', default='../input/*.pkl', help='Input folder and files in glob format')
    parser.add_argument('-output', default='../output/{cl}_{chromosome}_d{delta}.pkl')

    global input_folder, output_folder, annotations_filename, CHROMOLIST, delta
    args = parser.parse_args()
    input_glob = args.input
    output_filename = args.output

    for chromo in 'chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr2,chr20,chr21,chr22,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chrX,chrY'.split(','):
        files = 'C:/TMP/output/IMR90_*_{chromosome}_d1000.pkl'.format(chromosome=chromo)
        blackholes = []
        for fn in glob.glob(files):
            print('FN = {}'.format(fn))
            with open(fn, 'rb') as f:
                blackhole = pickle.load(f)
                blackholes.append(blackhole)
        res = [[blackhole[j][0] for blackhole in blackholes] for j in range(len(blackholes[0]))]
        out = 'IMR90_{chromosome}_d1000.pkl'.format(chromosome=chromo)
        print('[*] Dumping back to <{}>...'.format(out))
        with open(out, 'wb') as f:
            pickle.dump(res, f)

    # data = defaultdict(lambda: defaultdict(list))

    # for input_filename in glob.glob(input_glob):
    #     cellline, chromo, delta = os.path.split(input_filename)[1].rsplit('_', 2)
    #     print(cellline, chromo, delta)
    #     with open(input_filename, 'rb') as f:
    #         blackhole = pickle.load(f)
    #     for rebranded, gene_id, left, right in blackhole:
    #         data[chromo][delta].append((cellline, rebranded, gene_id, left, right))

    # meta = []
    # warped = dict()
    # for chromo in data:
    #     for delta in data[chromo]:
    #         for cellline, rebranded, gene_id, left, right in data[chromo][delta]:
    #             warped[cellline].append(rebranded)


    # for chromo in data:
    #     for delta in data[chromo]:
    #         meta.append(())
    #             print(cellline, gene_id, left, right)

    #         # for cl in data[chromo][delta]:
    #         #     out = output_filename.format(cl=cl, chromosome=chromo, delta=delta)
    #         #     print('[*] Writing to <{}>...'.format(out))

    # # print('len(warped) = {}'.format(len(warped)))
    # # m = len(warped)  # experiments
    # # meta = [()]
    # # for cellline, (blackhole, chromo, delta) in warped.items():


if __name__ == '__main__':
    main()
