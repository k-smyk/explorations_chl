# -*- coding: utf-8 -*-
"""pmi_iwsa

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1XR80A5m-P9M-BTtWmq37hxowaNfMtxTP
"""

# -*- coding: utf-8 -*-

# !pip install random
# !pip install numpy
# !pip install pandas
# !pip install re
# !pip install Bio
# !pip install os
# !pip install multiprocessing
# !pip install Levenshtein

import random as pyrandom
pyrandom.seed(12345)
from numpy import *
import numpy as np
random.seed(12345)
import pandas as pd
import Levenshtein
import re
from Bio import pairwise2
from multiprocessing import Process,Manager
ncores = 50



# Function: nexCharOutput
# Description:
## This function takes a character array, a list of rownames
## and the name of the output nexus file as input
## and writes the character matrix into a nexus file.
## Missing entries are assumed to be coded as "-1"

def nexCharOutput(chMtx,names,outfile,datatype='STANDARD'):
    f = open(outfile,'w')
    f.write('#NEXUS\n\n')
    f.write('BEGIN DATA;\n')
    f.write('DIMENSIONS ntax='+str(len(chMtx))+' NCHAR='+str(len(chMtx.T))+';\n')
    f.write('FORMAT DATATYPE='+datatype+' GAP=? MISSING=- interleave=yes;\n')
    f.write('MATRIX\n\n')
    txLgth = max(map(len,names))
    for i in range(len(chMtx)):
        f.write(names[i].ljust(txLgth+2))
        for ch in chMtx[i]:
            if ch==-1: ch='-'
            else:
                ch = str(ch)
            f.write(ch)
        f.write('\n')
    f.write('\n;\n\nEND;\n')
    f.close()


data = pd.read_csv('dataset.tab',
                   index_col=0,na_filter=False,sep='\t')
data = data[data.wls_gen.isin(['ROMANCE','ALBANIAN'])]
data = data[data.index!='LATIN']

concepts100 = array(data.columns[9:])
nEntries = data[concepts100].apply(lambda x:sum(x!='')).sort_values()
concepts = nEntries.index[-40:]
data = data[concepts]
taxa = array(data.index)

def cleanASJP(word):
    """takes an ASJP string as argument
    and returns the string with all diacritics removed."""
    word = re.sub(r",|\%|\*|\"|\.~|\$(\d)|\s+", "", word)
    return word.replace('~', '')


dataWL = pd.DataFrame([(c,l,cleanASJP(w))
                       for c in concepts
                       for l in taxa
                       for w in data[c][l].split(',')
                       if data[c][l]!=''],
                      columns = ['concept','language','word'])

dataWL = dataWL.drop_duplicates(['concept','language'])

training = pd.DataFrame()
for c in concepts:
    cData = dataWL[dataWL.concept==c]
    cTraining = pd.DataFrame([cData.loc[[i,j]].word.values
                              for i in cData.index
                              for j in cData.index
                              if i<j])
    training = pd.concat([training, cTraining], ignore_index=True)

sounds = unique(concatenate(list(map(list,dataWL.word.values))))


def levalign(w):
    """takes a pair of strings as input
    and returns the Levenshtein alignment
    in column format. Gaps are removed."""
    x,y = w
    algn = zeros((0,2))
    if '0' in [x,y]: return algn
    e = Levenshtein.opcodes(x,y)
    for a in e:
        if a[0] in ['replace','equal']:
            x_a, x_e = a[1],a[2]
            y_a, y_e = a[3],a[4]
            ag = [list(x[x_a:x_e]), list(y[y_a:y_e])]
            algn = concatenate([algn,transpose(ag)])
    return algn


manager = Manager()
return_dict = manager.dict()
packages = array_split(training.values,ncores)


def doWork(i,pck):
    return_dict[i] = vstack([levalign(p) for p in pck])


jobs = []
for i,pck in enumerate(packages):
    p = Process(target=doWork,args=(i,pck))
    p.start()
    jobs.append(p)
for p in jobs:
    p.join()
alg0 = vstack(return_dict.values())


# count alignment frequencies
sFreqs = pd.crosstab(alg0[:,0],alg0[:,1])
sFreqs = sFreqs.reindex(sounds,fill_value=0).T
sFreqs = sFreqs.reindex(sounds,fill_value=0).T
# symmetrize
sFreqs = sFreqs.copy()+sFreqs.copy().T
# add-1 smoothing
sFreqs += 1
sProbs = sFreqs/sFreqs.sum().sum()

# extract relative sound frequencies
soundOccurrences = concatenate([list(w) for w in dataWL.word.values])
soundProbabilities = pd.value_counts(soundOccurrences,normalize=True)[sounds]
pmi0 = (log(sProbs).copy()-log(soundProbabilities)).T-log(soundProbabilities)


def sscore(a,b,pmiDict,gp1,gp2):
    """a,b: ASJP strings
    pmiDict: logodds dictionary
    gp1,gp2: gap penalties
    return PMI score of a/b
    """
    out = pairwise2.align.globalds(a,b,pmiDict,gp1,gp2)
    if len(out)==0: return nan
    return out[0][2]


def scoreNW(x,y,pmiDict,gp1,gp2):
    """x,y: sequences of ASJP strings, separated by '-'
    pmiDict: logodds dictionary
    gp1,g2: gap penalties
    returns maximal PMI score for the Cartesian product of x and y"""
    if '0' in [x,y]: return nan
    x1=x.split('-')
    y1=y.split('-')
    return max([sscore(xx,yy,pmiDict,gp1,gp2) for xx in x1 for yy in y1])


pmi0Dict = {(s1,s2):pmi0[s1][s2]
            for s1 in sounds for s2 in sounds}

gp1=-2.49302792222
gp2=-1.70573165621
th = 4.4451


'''
Start of Kateryna Smykovska's code
'''
def iwsa(a, b, w, IL):
    m, n = len(a), len(b)
    M = np.zeros((m + 1, n + 1))
    for i in range(1, m + 1):
        M[i][0] = M[i - 1][0] + w(a[i - 1], '-') * IL(a, b, i - 1, '-')
    for j in range(1, n + 1):
        M[0][j] = M[0][j - 1] + w('-', b[j - 1]) * IL(a, b, '-', j - 1)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            M[i][j] = min(
                M[i - 1][j - 1] + w(a[i - 1], b[j - 1]) * IL(a, b, i - 1, j - 1),
                M[i - 1][j] + w(a[i - 1], '-') * IL(a, b, i - 1, '-'),
                M[i][j - 1] + w('-', b[j - 1]) * IL(a, b, '-', j - 1)
            )
    return M[m][n]

def iwsa_score(x, y, pmiDict, gp1, gp2):
    if '0' in [x, y]:
        return nan
    x1 = x.split('-')
    y1 = y.split('-')
    return max([iwsa(xx, yy, lambda a, b: pmiDict.get((a, b), -inf), lambda a, b, i, j: 1) for xx in x1 for yy in y1])

def iwsa_align(w, pmiDict, gp1, gp2, th):
    x, y = w
    a = iwsa(x, y, lambda a, b: pmiDict.get((a, b), -inf), lambda a, b, i, j: 1)
    if a < th:
        return zeros((0, 2))
    algn = []
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i] != '-' and y[j] != '-':
                algn.append([x[i], y[j]])
    return array(algn)

def iwsa_align_star(crp, pmiDict, gp1, gp2, th):
    packages = array_split(crp, ncores)
    manager = Manager()
    return_dict = manager.dict()
    def doWork(i, pck):
        try:
            result = vstack([iwsa_align(w, pmiDict, gp1, gp2, th) for w in pck])
            return_dict[i] = result
        except Exception as e:
            print(f"Error in process {i}: {e}")
    jobs = []
    for i, pck in enumerate(packages):
        p = Process(target=doWork, args=(i, pck))
        p.start()
        jobs.append(p)
    for p in jobs:
        p.join()
    if return_dict:
        return vstack(return_dict.values())
    else:
        print("Warning: No results returned from child processes.")
        return np.array([])

'''
End of Kateryna Smykovska's code
'''

pmiDict = pmi0Dict.copy()
for i in range(10):
    alg = iwsa_align_star(training.values, pmiDict, gp1, gp2, th)
    sFreqs = pd.crosstab(alg[:, 0], alg[:, 1])
    sFreqs = sFreqs.reindex(sounds, fill_value=0).T
    sFreqs = sFreqs.reindex(sounds, fill_value=0).T
    sFreqs = sFreqs.copy() + sFreqs.copy().T
    sFreqs += 1
    sProbs = sFreqs / sFreqs.sum().sum()
    pmi = (log(sProbs).copy() - log(soundProbabilities)).T - log(soundProbabilities)
    pmiDict = {(s1, s2): pmi[s1][s2]
               for s1 in sounds for s2 in sounds}

pmi.to_csv('pmi-albanoRomance_IWSA.csv')
dataWL.to_csv('albanoRomanceASJP_IWSA.csv',index=False)

sc = pd.DataFrame(index=taxa)
for c in concepts:
    cData = dataWL[dataWL.concept==c]
    cTaxa = cData.language.unique()
    cWords = pd.Series([''.join(cData[cData.language==l].word.values)
                        for l in cTaxa],
                       index=cTaxa)
    cMtx = pd.DataFrame([[int(s in cWords[l]) for s in sounds]
                         for l in cTaxa],
                        index=cTaxa)
    cMtx = cMtx.reindex(taxa,fill_value='-')
    sc = pd.concat([sc,cMtx],axis=1)

nexCharOutput(sc.values,sc.index,'albanoRomanceSC_IWSA.nex','restriction')

def levalignFull(w):
    """takes a pair of strings as input
    and returns the Levenshtein alignment
    in column format."""
    x,y = w
    algn = zeros((0,2))
    e = Levenshtein.opcodes(x,y)
    for a in e:
        if a[0] in ['replace','equal']:
            x_a, x_e = a[1],a[2]
            y_a, y_e = a[3],a[4]
            ag = [list(x[x_a:x_e]), list(y[y_a:y_e])]
        elif a[0]=='delete':
            x_a, x_e = a[1],a[2]
            ag = [list(x[x_a:x_e]), ['-']*(x_e-x_a)]
        else:
            y_a, y_e = a[3],a[4]
            ag = [ ['-']*(y_e-y_a),list(y[y_a:y_e])]
        algn = concatenate([algn,transpose(ag)])
    return algn
