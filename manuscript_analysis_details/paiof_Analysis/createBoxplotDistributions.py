import os
import sys

paiof_full_tsv = sys.argv[1]
genome1_id = sys.argv[2]
genome2_id = sys.argv[3]
comparison_id = sys.argv[4]

ahba_genes = set(['GCA_000282715.1_ASM28271v1|BUPY_000639', 'GCA_000282715.1_ASM28271v1|BUPY_000640', 'GCA_000282715.1_ASM28271v1|BUPY_000641', 'GCA_000282715.1_ASM28271v1|BUPY_000642', 'GCA_000282715.1_ASM28271v1|BUPY_000643', 'GCA_000384275.1_ASM38427v1|AXGX_005050', 'GCA_000384275.1_ASM38427v1|AXGX_005051', 'GCA_000384275.1_ASM38427v1|AXGX_005071', 'GCA_000384275.1_ASM38427v1|AXGX_005052', 'GCA_000384275.1_ASM38427v1|AXGX_005053', 'BXKM_001028', 'BXKM_001027', 'BXKM_001026', 'BXKM_001025', 'BXKM_001024'])

with open(paiof_full_tsv) as opft:
    for line in opft:
        line = line.strip()
        ls = line.split('\t')
        if ls[0] == genome1_id and ls[1] == genome2_id:
            is_ahba_gene = False
            if ls[2] in ahba_genes and ls[3] in ahba_genes:
                is_ahba_gene = True
            print(comparison_id + '\t' + str(is_ahba_gene) + '\t' + ls[4])

