import os
import sys
from ete3 import Tree

gca_to_gen = {}
for line in open('Common_Genera_Representatives.txt'):
    line = line.strip()
    ls = line.split('\t')
    gca_to_gen[ls[-1]] = ls[0]

t = Tree(sys.argv[1], format=1)
for node in t.traverse('postorder'):
    if not node.is_leaf(): continue
    gca = '_'.join(node.name.replace("'", "").replace('"', '').split('_')[:2])
    node.name = gca_to_gen[gca]

R = t.get_midpoint_outgroup()
t.set_outgroup(R)

t.write(format=1, outfile=sys.argv[2])
