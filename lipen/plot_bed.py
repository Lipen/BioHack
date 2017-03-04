import csv
import time
from collections import namedtuple, defaultdict
import matplotlib.pyplot as plt

input_filename = 'C:/TMP/input/H3K4me1.IMR90.mono.bed'
CHROMO = 'chr1+'


def parse(input_filename):
    print('[*] Parsing <{}>...'.format(input_filename))
    time_start = time.time()
    Peak = namedtuple('Peak', ['start', 'end', 'score'])
    data = defaultdict(list)

    with open(input_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for chromo, start, end, stuff, score, strand in reader:
            if chromo + strand == CHROMO:  # and int(score) > 100:
                data[chromo + strand].append(Peak(int(start), int(end), int(score)))

    print('[+] Parsed <{}> peaks in {} seconds'.format(sum(map(len, data.values())), time.time() - time_start))
    return data

data = parse(input_filename)

plt.hist(sorted([item.start for item in data[CHROMO]]), bins=100)
plt.show()
