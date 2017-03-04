import csv
import time
import pickle

input_filename = '../data/NFkB_diffExp_a549.tsv'
output_filename = 'C:/TMP/output/NFkB_diffExp_a549.pkl'


def parse(input_filename):
    print('[*] Parsing <{}>...'.format(input_filename))
    time_start = time.time()
    data = dict()

    with open(input_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for gene_id, log2_fold_change, padj in reader:
            try:
                log2_fold_change = float(log2_fold_change)
            except:
                log2_fold_change = 0
            try:
                padj = float(padj)
            except:
                padj = 0
            data[gene_id] = (log2_fold_change, padj)

    print('[+] Parsed <{}> peaks in {} seconds'.format(sum(map(len, data.values())), time.time() - time_start))
    return data

data = parse(input_filename)
with open(output_filename, 'wb') as f:
    pickle.dump(data, f)
