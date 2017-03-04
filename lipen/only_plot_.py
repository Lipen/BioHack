import time
import pickle
import matplotlib.pyplot as plt

input_filename = 'C:/TMP/dm6_rRNA_clu_{chromosome}.pkl'
CHROMOLIST = ['chr3L']

print('[*] Working...')
time_start = time.time()

for i, chromo in enumerate(CHROMOLIST):
    in_ = input_filename.format(chromosome=chromo)
    print('[*] Loading cluster dump <{}>...'.format(in_))
    try:
        with open(in_, 'rb') as f:
            clusters = pickle.load(f)
    except FileNotFoundError:
        print('[!] Okay, there is no such dump. Just continue...')
        continue

    plt.figure(figsize=(30, 12))
    plt.plot(clusters, '-', marker='o', markerfacecolor='red')
    # plt.plot(clusters)
    plt.grid(True)
    plt.title('{}'.format(chromo))
    plt.xlabel('# of cluster')
    plt.ylabel('Peaks in cluster')
    plt.tight_layout()
    plt.savefig('C:/TMP/clusters_{}.png'.format(chromo), dpi=300)

print('[+] Done in {} seconds!'.format(time.time() - time_start))
