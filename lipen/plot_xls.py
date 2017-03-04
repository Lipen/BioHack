import csv
import time
from collections import defaultdict
import matplotlib.pyplot as plt

input_filename = 'C:/TMP/input/IMR90_p65_TNFa_peaks.xls'


def parse(input_filename):
    print('[*] Parsing <{}>...'.format(input_filename))
    time_start = time.time()
    data = defaultdict(list)

    with open(input_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for _ in range(24):
            next(reader)
        for chromo, start, end, length, summit, tags, q, fold_enrichment, FDR in reader:
            data[chromo].append(float(fold_enrichment))

    print('[+] Parsed in {} seconds'.format(time.time() - time_start))
    return data

data = parse(input_filename)

plt.hist(data['chr1'], bins=100)
plt.show()
