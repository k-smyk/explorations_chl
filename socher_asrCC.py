from numpy import *
import pandas as pd
from subprocess import Popen
import os


ccPW = pd.read_csv('socher_albanoRomanceCC.bin.csv', index_col=0, dtype='str')
ccSW = ccPW
ccNW = ccPW

romancePW = array([x for x in ccPW.index if not 'ALBANIAN' in x])
romanceSW = array([x for x in ccSW.index if not 'ALBANIAN' in x])
romanceNW = array([x for x in ccNW.index if not 'ALBANIAN' in x])

ccPW = ccPW.loc[romancePW]
ccSW = ccSW.loc[romanceSW]
ccNW = ccNW.loc[romanceNW]

with open('dialign+socher+original/socher_romanceCC_PW.tsv', 'w') as f:
    for i in ccPW.index:
        f.write(i+'\t')
        f.write('\t'.join(ccPW.loc[i].values)+'\n')

with open('sw+socher+original/socher_romanceCC_SW.tsv', 'w') as f:
    for i in ccSW.index:
        f.write(i+'\t')
        f.write('\t'.join(ccSW.loc[i].values)+'\n')

with open('original+socher+original/socher_romanceCC.tsv', 'w') as f:
    for i in ccNW.index:
        f.write(i+'\t')
        f.write('\t'.join(ccNW.loc[i].values)+'\n')


paupCommandsPW = """#Nexus
Begin Paup;
set incr=auto;
gettrees file = socher_romancePW.posterior.tree;
savetrees file = socher_romancePW.posterior.nex.tree replace=yes brlen=user;
q;
End;
"""

paupCommandsSW = """#Nexus
Begin Paup;
set incr=auto;
gettrees file = socher_romanceSW.posterior.tree;
savetrees file = socher_romanceSW.posterior.nex.tree replace=yes brlen=user;
q;
End;
"""

paupCommandNW = """#Nexus
Begin Paup;
set incr=auto;
gettrees file = socher_romanceNW.posterior.tree;
savetrees file = socher_romanceNW.posterior.nex.tree replace=yes brlen=user;
q;
End;
"""

with open('dialign+socher+original/socher_convertRomancePosteriorPW.paup', 'w') as f:
    f.write(paupCommandsPW)

with open('sw+socher+original/socher_convertRomancePosteriorSW.paup', 'w') as f:
    f.write(paupCommandsSW)
with open('original+socher+original/socher_convertRomancePosteriorNW.paup', 'w') as f:
    f.write(paupCommandNW)

pPW = Popen('paup4 socher_convertRomancePosteriorPW.paup>/dev/null',shell=True)
os.waitpid(pPW.pid,0)
pSW = Popen('paup4 socher_convertRomancePosteriorSW.paup>/dev/null',shell=True)
os.waitpid(pSW.pid,0)
pNW = Popen('paup4 socher_convertRomancePosteriorNW.paup>/dev/null',shell=True)
os.waitpid(pNW.pid,0)

#p.169
btCommands = """1
1
seed 12345;
mlt 1;
ga 4;
run"""

with open('asrCC.bt','w') as f:
    f.write(btCommands)

pPW = Popen('BayesTraitsV4 socher_romancePW.posterior.nex.tree socher_romanceCC_PW.tsv < asrCC.bt>/dev/null',
          shell=True)
os.waitpid(pPW.pid,0)
pSW = Popen('BayesTraitsV4 socher_romanceSW.posterior.nex.tree socher_romanceCC_SW.tsv < asrCC.bt>/dev/null',
          shell=True)
os.waitpid(pSW.pid,0)
pNW = Popen('BayesTraitsV4 socher_romanceNW.posterior.nex.tree socher_romanceCC.tsv < asrCC.bt>/dev/null',
          shell=True)
os.waitpid(pNW.pid,0)


resultsPW = pd.read_table('socher_romanceCC_PW.tsv.log.txt',sep='\t',
                        skiprows=33)
resultsSW = pd.read_table('socher_romanceCC_SW.tsv.log.txt',sep='\t',
                        skiprows=33)
resultsNW = pd.read_table('socher_romanceCC.tsv.log.txt',sep='\t',
                        skiprows=33)

clPW = [x for x in resultsPW.columns if 'P(1)' in x]
clSW = [x for x in resultsSW.columns if 'P(1)' in x]
clNW = [x for x in resultsNW.columns if 'P(1)' in x]

resultsPW = resultsPW[clPW]
resultsPW.columns = ccPW.columns
resultsSW = resultsSW[clSW]
resultsSW.columns = ccSW.columns
resultsNW = resultsNW[clNW]
resultsNW.columns = ccNW.columns

conceptsPW = unique([x.split(':')[0] for x in resultsPW.columns])
conceptsSW = unique([x.split(':')[0] for x in resultsSW.columns])
conceptsNW = unique([x.split(':')[0] for x in resultsNW.columns])

winnersPW = []
winnersSW = []
winnersNW = []
for c in conceptsPW:
    cChars = [x for x in resultsPW.columns if x.split(':')[0]==c]
    winnersPW.append(resultsPW.mean()[cChars].idxmax())
for c in conceptsSW:
    cChars = [x for x in resultsSW.columns if x.split(':')[0]==c]
    winnersSW.append(resultsSW.mean()[cChars].idxmax())
for c in conceptsNW:
    cChars = [x for x in resultsNW.columns if x.split(':')[0]==c]
    winnersNW.append(resultsNW.mean()[cChars].idxmax())

winnersPW = pd.Series(winnersPW,index=conceptsPW)
winnersSW = pd.Series(winnersSW,index=conceptsSW)
winnersNW = pd.Series(winnersNW,index=conceptsNW)

winnersPW.to_csv('socher_asrCC_PW.csv', header=False)
winnersSW.to_csv('socher_asrCC_SW.csv', header=False)
winnersNW.to_csv('socher_asrCC_NW.csv', header=False)
