from numpy import *
import pandas as pd
from subprocess import Popen
import os


ccPW = pd.read_csv('dialign+original_pipeline/albanoRomanceCCbinPW.csv', index_col=0, dtype='str')
ccSW = pd.read_csv('sw+original_pipeline/albanoRomanceCCbinSW.csv', index_col=0, dtype='str')

romancePW = array([x for x in ccPW.index if not 'ALBANIAN' in x])
romanceSW = array([x for x in ccSW.index if not 'ALBANIAN' in x])

ccPW = ccPW.loc[romancePW]
ccSW = ccSW.loc[romanceSW]

with open('dialign+original_pipeline/romanceCC_PW.tsv', 'w') as f:
    for i in ccPW.index:
        f.write(i+'\t')
        f.write('\t'.join(ccPW.loc[i].values)+'\n')

with open('sw+original_pipeline/romanceCC_SW.tsv', 'w') as f:
    for i in ccSW.index:
        f.write(i+'\t')
        f.write('\t'.join(ccSW.loc[i].values)+'\n')

paupCommandsPW = """#Nexus
Begin Paup;
set incr=auto;
gettrees file = romancePW.posterior.tree;
savetrees file = romancePW.posterior.nex.tree replace=yes brlen=user;
q;
End;
"""

paupCommandsSW = """#Nexus
Begin Paup;
set incr=auto;
gettrees file = romanceSW.posterior.tree;
savetrees file = romanceSW.posterior.nex.tree replace=yes brlen=user;
q;
End;
"""

with open('dialign+original_pipeline/convertRomancePosteriorPW.paup', 'w') as f:
    f.write(paupCommandsPW)

with open('sw+original_pipeline/convertRomancePosteriorSW.paup', 'w') as f:
    f.write(paupCommandsSW)

pPW = Popen('paup4 convertRomancePosteriorPW.paup>/dev/null',shell=True)
os.waitpid(pPW.pid,0)
pSW = Popen('paup4 convertRomancePosteriorSW.paup>/dev/null',shell=True)
os.waitpid(pSW.pid,0)

#p.169
btCommands = """1
1
seed 12345;
mlt 1;
ga 4;
run"""

with open('asrCC.bt','w') as f:
    f.write(btCommands)

pPW = Popen('BayesTraitsV4 romancePW.posterior.nex.tree romanceCC_PW.tsv < asrCC.bt>/dev/null',
          shell=True)
os.waitpid(pPW.pid,0)
pSW = Popen('BayesTraitsV4 romanceSW.posterior.nex.tree romanceCC_SW.tsv < asrCC.bt>/dev/null',
          shell=True)
os.waitpid(pSW.pid,0)


resultsPW = pd.read_table('romanceCC_PW.tsv.log.txt',sep='\t',
                        skiprows=33)
resultsSW = pd.read_table('romanceCC_SW.tsv.log.txt',sep='\t',
                        skiprows=33)

clPW = [x for x in resultsPW.columns if 'P(1)' in x]
clSW = [x for x in resultsSW.columns if 'P(1)' in x]

resultsPW = resultsPW[clPW]
resultsPW.columns = ccPW.columns
resultsSW = resultsSW[clSW]
resultsSW.columns = ccSW.columns

conceptsPW = unique([x.split(':')[0] for x in resultsPW.columns])
conceptsSW = unique([x.split(':')[0] for x in resultsSW.columns])

winnersPW = []
winnersSW = []
for c in conceptsPW:
    cChars = [x for x in resultsPW.columns if x.split(':')[0]==c]
    winnersPW.append(resultsPW.mean()[cChars].idxmax())
for c in conceptsSW:
    cChars = [x for x in resultsSW.columns if x.split(':')[0]==c]
    winnersSW.append(resultsSW.mean()[cChars].idxmax())

winnersPW = pd.Series(winnersPW,index=conceptsPW)
winnersSW = pd.Series(winnersSW,index=conceptsSW)

winnersPW.to_csv('asrCC_PW.csv', header=False)
winnersSW.to_csv('asrCC_SW.csv', header=False)
