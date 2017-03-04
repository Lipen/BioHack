import glob
import pickle


def main():
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

if __name__ == '__main__':
    main()
