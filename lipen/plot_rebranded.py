import time
import pickle
import matplotlib.pyplot as plt

input_filename = 'C:/TMP/output/H3K4me1.IMR90.mono.bed_reb_{chromosome}.pkl'
# input_filename = 'C:/TMP/output/H3K4me1.IMR90.mono.bed_REGIONS_{chromosome}.pkl'
annotations_filename = '../data/annotation_{chromosome}.txt'
CHROMOLIST = ['chr1+', 'chr2+', 'chr3+', 'chr1-', 'chr2-', 'chr3-']

print('[*] Working...')
time_start = time.time()
fig, axes = plt.subplots(2, 3, figsize=(16, 9))

for i, chromo in enumerate(CHROMOLIST):
    in_ = input_filename.format(chromosome=chromo)
    print('[*] Loading cluster dump <{}>...'.format(in_))
    try:
        with open(in_, 'rb') as f:
            rebranded = pickle.load(f)
    except FileNotFoundError:
        print('[!] Okay, there is no such dump. Just continue...')
        continue

    plt.subplot(2, 3, i + 1)
    # plt.hist(rebranded[:, 1], bins=100)
    plt.plot(rebranded[:, 1])
    # plt.hist([r[0] for r in rebranded], bins=100, range=(0, 100))
    # plt.plot([r[0] for r in rebranded])
    print(rebranded[:, 1])
    plt.grid(True)
    plt.title('{}'.format(chromo))

print('[+] Done in {} seconds!'.format(time.time() - time_start))

# Dummy subplot for 'shared' labels
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.title('MAYBE Peaks distributions')
plt.xlabel('MAYBE # of cluster')
plt.ylabel('MAYBE Peaks in cluster')
plt.savefig('C:/TMP/rebranded.png', dpi=300)
# plt.show()
