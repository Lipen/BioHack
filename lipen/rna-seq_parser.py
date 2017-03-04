import os
import csv
import time
import glob
import pickle

input_filenames = '../data/*.tsv'
output_folder = 'C:/TMP/output'


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


for input_filename in glob.glob(input_filenames):
    data = parse(input_filename)
    out = os.path.join(output_folder, os.path.splitext(os.path.split(input_filename)[1])[0]) + '.pkl'
    with open(out, 'wb') as f:
        pickle.dump(data, f)
