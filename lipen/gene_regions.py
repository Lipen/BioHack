# genes :: {(a, b)}
genes = [(1, 2), (5, 8), (10, 12), (18, 24), (23, 28)]
genes_sorted = sorted(genes)

print('[*] Merging...')
merged = []
i = 0
while i < len(genes) - 1:
    x = genes[i]
    y = genes[i + 1]
    if x[1] < y[0]:
        merged.append(x)
        merged.append(y)
        i += 1
    else:
        merged.append((x[0], y[1]))
    i += 1

print('[+] Merged:')
for g in merged:
    print('  > {}'.format(g))

print('[*] Splitting')
regions = []
left = 0
for i in range(len(merged) - 1):
    x = merged[i]
    y = merged[i + 1]
    right = (x[1] + y[0]) // 2
    r = (left, right)
    regions.append(r)
    left = right

print('[+] Regions:')
for r in regions:
    print('  > {}'.format(r))
