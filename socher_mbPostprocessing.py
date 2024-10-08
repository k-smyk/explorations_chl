import pandas as pd
from numpy import *
import subprocess,os
from ete3 import Tree


p1PW = pd.read_table('dialign+socher+original/socher_albanoRomancePW.run1.p',
                     sep='\t', skiprows=1)
p2PW = pd.read_table('dialign+socher+original/socher_albanoRomancePW.run2.p',
                     sep='\t', skiprows=1)
p1SW = pd.read_table('dialign+socher+original/socher_albanoRomanceSW.run1.p',
                     sep='\t', skiprows=1)
p2SW = pd.read_table('dialign+socher+original/socher_albanoRomanceSW.run2.p',
                     sep='\t', skiprows=1)
p1NW = pd.read_table('original+socher+original/socher_albanoRomance.run1.p', sep='\t', skiprows=1)
p2NW = pd.read_table('original+socher+original/socher_albanoRomance.run2.p', sep='\t', skiprows=1)

prPW = pd.concat([p1PW[p1PW.Gen>1000000],p2PW[p2PW.Gen>1000000]])
prSW = pd.concat([p1SW[p1SW.Gen>1000000],p2SW[p2SW.Gen>1000000]])
prNW = pd.concat([p1NW[p1NW.Gen>1000000],p2NW[p2NW.Gen>1000000]])

treesPW = []
treesSW = []
treesNW = []

with open('dialign+socher+original/socher_albanoRomancePW.run1.tre') as f:
    for i,ln in enumerate(f.readlines()):
        if i >100:
            t = Tree(ln.strip())
            treesPW.append(t)
with open('dialign+socher+original/socher_albanoRomancePW.run2.tre') as f:
    for i,ln in enumerate(f.readlines()):
        if i >100:
            t = Tree(ln.strip())
            treesPW.append(t)
with open('dialign+socher+original/socher_albanoRomanceSW.run1.tre') as f:
    for i,ln in enumerate(f.readlines()):
        if i >100:
            t = Tree(ln.strip())
            treesSW.append(t)
with open('dialign+socher+original/socher_albanoRomanceSW.run2.tre') as f:
    for i,ln in enumerate(f.readlines()):
        if i >100:
            t = Tree(ln.strip())
            treesSW.append(t)
with open('original+socher+original/socher_albanoRomance.run1.tre') as f:
    for i,ln in enumerate(f.readlines()):
        if i >100:
            t = Tree(ln.strip())
            treesNW.append(t)
with open('original+socher+original/socher_albanoRomance.run2.tre') as f:
    for i,ln in enumerate(f.readlines()):
        if i >100:
            t = Tree(ln.strip())
            treesNW.append(t)


taxaPW = array(treesPW[0].get_leaf_names())
taxaSW = array(treesSW[0].get_leaf_names())
taxaNW = array(treesNW[0].get_leaf_names())

romancePW = array([x for x in taxaPW if not 'ALBANIAN' in x])
romanceSW = array([x for x in taxaSW if not 'ALBANIAN' in x])
romanceNW = array([x for x in taxaNW if not 'ALBANIAN' in x])

for t in treesPW:
    t.set_outgroup('ALBANIAN')
    t.set_outgroup(t.get_common_ancestor([t&l for l in romancePW]))
for t in treesSW:
    t.set_outgroup('ALBANIAN')
    t.set_outgroup(t.get_common_ancestor([t&l for l in romanceSW]))
for t in treesNW:
    t.set_outgroup('ALBANIAN')
    t.set_outgroup(t.get_common_ancestor([t&l for l in romanceNW]))

with open('dialign+socher+original/socher_albanoRomancePW.posterior.tree', 'w') as f: f.write('')
with open('dialign+socher+original/socher_albanoRomancePW.posterior.tree', 'a') as f:
    for t in treesPW:
        f.write(t.write(format=1)+'\n')
with open('dialign+socher+original/socher_albanoRomanceSW.posterior.tree', 'w') as f: f.write('')
with open('dialign+socher+original/socher_albanoRomanceSW.posterior.tree', 'a') as f:
    for t in treesSW:
        f.write(t.write(format=1)+'\n')
with open('original+socher+original/socher_albanoRomance.posterior.tree', 'w') as f: f.write('')
with open('original+socher+original/socher_albanoRomance.posterior.tree', 'a') as f:
    for t in treesNW:
        f.write(t.write(format=1)+'\n')


for t in treesPW:
    t.prune([t&l for l in romancePW])
for t in treesSW:
    t.prune([t&l for l in romanceSW])
for t in treesNW:
    t.prune([t&l for l in romanceNW])


with open('dialign+socher+original/socher_romancePW.posterior.tree', 'w') as f: f.write('')
with open('dialign+socher+original/socher_romancePW.posterior.tree', 'a') as f:
    for t in treesPW:
        f.write(t.write(format=1)+'\n')
with open('sw+socher+original/socher_romanceSW.posterior.tree', 'w') as f: f.write('')
with open('sw+socher+original/socher_romanceSW.posterior.tree', 'a') as f:
    for t in treesSW:
        f.write(t.write(format=1)+'\n')
with open('original+socher+original/socher_romanceNW.posterior.tree', 'w') as f: f.write('')
with open('original+socher+original/socher_romanceNW.posterior.tree', 'a') as f:
    for t in treesNW:
        f.write(t.write(format=1)+'\n')

