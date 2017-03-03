import csv
from collections import namedtuple

Item = namedtuple('Peak', ['chr', 'start', 'end', ])

filename = '../data/dm6_rRNA.bed'

with open(filename) as f:
    reader = csv.reader(f, delimiter='\t')
    data = [Item(item) for item in reader]

print('LEN DATA: {}'.format(len(data)))
