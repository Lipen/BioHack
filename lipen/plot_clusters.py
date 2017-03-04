import time
import pickle
import matplotlib.pyplot as plt

input_filename = 'C:/TMP/dm6_rRNA_clu_{chromosome}.pkl'
CHROMOLIST = ['chr4', 'chrX', 'chrY', 'chr2L', 'chr2R', 'chr3L', 'chr3R']

print('[*] Working...')
time_start = time.time()
fig, axes = plt.subplots(2, 4, figsize=(16, 9))

for i, chromo in enumerate(CHROMOLIST):
    in_ = input_filename.format(chromosome=chromo)
    print('[*] Loading cluster dump <{}>...'.format(in_))
    try:
        with open(in_, 'rb') as f:
            clusters = pickle.load(f)
    except FileNotFoundError:
        print('[!] Okay, there is no such dump. Just continue...')
        continue

    plt.subplot(2, 4, i + 1)
    # plt.plot(clusters, '-', marker='o', markerfacecolor='red')
    plt.plot(clusters)
    plt.grid(True)
    plt.title('{}'.format(chromo))

print('[+] Done in {} seconds!'.format(time.time() - time_start))

# Dummy subplot for 'shared' labels
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.title('Peaks distributions')
plt.xlabel('# of cluster')
plt.ylabel('Peaks in cluster')
plt.savefig('C:/TMP/clusters.png', dpi=300)
# plt.show()
