import os
import glob
import pickle
import numpy as np


def main():
    for delta, delta2 in {500: '500', 1000: '1k', 10000: '10k'}.items():
        for cl in ['IMR90', 'A549']:
            for chromo in 'chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr2,chr20,chr21,chr22,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chrX'.split(','):
                files = 'C:/TMP/output/{cl}*{chromosome}_d{delta}.pkl'.format(cl=cl, chromosome=chromo, delta=delta)
                chip_seqs = []
                for fn in glob.glob(files):
                    print('FN = {}'.format(fn))
                    with open(fn, 'rb') as f:
                        experiment = pickle.load(f)
                        chip_seqs.append(experiment)
                res = []
                for j in range(len(chip_seqs[0])):
                    res.append((np.array([chip_seqs[0][j][0],
                                          chip_seqs[1][j][0],
                                          chip_seqs[2][j][0],
                                          chip_seqs[3][j][0],
                                          chip_seqs[4][j][0]]),
                                chip_seqs[0][j][1],
                                chip_seqs[0][j][2],
                                chip_seqs[0][j][3]))
                out = 'C:/TMP/warped/{cl}_{delta2}/{cl}_{chromosome}_d{delta}.pkl'.format(cl=cl, chromosome=chromo, delta=delta, delta2=delta2)
                print('[*] Dumping back to <{}>...'.format(out))
                with open(out, 'wb') as f:
                    pickle.dump(res, f)

if __name__ == '__main__':
    main()
