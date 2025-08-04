import os
import sys

gca_to_genus = {}
for line in open('Common_Genera_Representatives.txt'):
    line = line.strip()
    ls = line.split('\t')
    gca_to_genus[ls[-1]] = ls[0]

# Genome A       Genes in A      Genome B        Genes in B      # orthologous genes     Mean AAI        Std AAI Orthologous fraction (OF)

gca_to_group = {}
for line in open('Listings/Multicellular1.txt'):
    line = line.strip()
    gca = '_'.join(line.split('_')[:2]) 
    gca_to_group[gca] = 'Multicellular1'

for line in open('Listings/Multicellular2.txt'):
    line = line.strip()
    gca = '_'.join(line.split('_')[:2])
    gca_to_group[gca] = 'Multicellular2'

for line in open('Listings/Simple.txt'):
    line = line.strip()
    gca = '_'.join(line.split('_')[:2])
    gca_to_group[gca] = 'Simple'

for line in open('Listings/Outgroups.txt'):
    line = line.strip()
    ls = line.split('\t')
    gca_to_group[gca] = 'Outgroup'

select = set([])
for line in open('Closer_to_Group2.txt'):
    line = line.strip()
    select.add(line)

print('Focal\tComparator\tComparator_Group\tAAI\tShared_Prop')
for i, line in enumerate(open('Pairwise_AAI_Estimates.tsv')):
    # Genome_1        Genome_2        AAI     Stdev_AAI       Genes_Shared    Genes_Shared_Normalized_by_Genome_1_Genes
    if i == 0: continue
    line = line.strip()
    ls = line.split('\t')
    gcaA = '_'.join(ls[0].split('_')[:2])
    gcaB = '_'.join(ls[1].split('_')[:2])
    if not gcaA in gca_to_group or not gcaB in gca_to_group: continue
    genusA = gca_to_genus[gcaA]
    genusB = gca_to_genus[gcaB]
    groupA = gca_to_group[gcaA]
    groupB = gca_to_group[gcaB]
    if groupA != 'Simple' and groupB != 'Simple': continue
    if groupA == 'Simple' and groupB.startswith('Multi'):
        if genusA in select: 
            aai = ls[2]
            shared_prop = ls[-1]
            print(genusA + '\t' + genusB + '\t' + groupB + '\t' + aai + '\t' + shared_prop)
