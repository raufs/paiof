import os
import sys

paiof_file = 'paiof_Results/Full_RBH_Listings.tsv'
genomes = ['GCA_000282715.1_ASM28271v1.faa', 'GCA_000384275.1_ASM38427v1.faa', 'GCA_001013905.1_sleC34.faa']
names = ['Amy', 'Micro', 'Strepto']

print('comparison\tahba_genes\taai')
for i, g1 in enumerate(genomes):
    n1 = names[i]
    for j, g2 in enumerate(genomes):
        n2 = names[j]
        if i <= j: continue
        comparison = n1 + '_vs_' + n2 
        cmd = ['python', 'createBoxplotDistributions.py', paiof_file, g1, g2, comparison]
        os.system(' '.join(cmd))
