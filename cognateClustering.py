import random as pyrandom
pyrandom.seed(12345)
from numpy import *
random.seed(12345)
import pandas as pd
from Bio import pairwise2
from sklearn.metrics import precision_recall_curve
from sklearn.linear_model import LogisticRegression
import igraph
from scipy import stats
import tempfile
from ete3 import Tree
import subprocess
import os

gp1=-2.49302792222
gp2=-1.70573165621


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


data = pd.read_csv('albanoRomanceASJP.csv')
data['ID'] = range(len(data))

#pmi = pd.read_csv('pmi-albanoRomance.csv',index_col=0)
pmiPW = pd.read_csv('dialign+original_pipeline/pmi-albanoRomance-pw.csv', index_col=0)
pmiSW = pd.read_csv('sw+original_pipeline/pmi-albanoRomance-sw.csv', index_col=0)
#sounds = array(pmi.index)
soundsPW = array(pmiPW.index)
soundsSW = array(pmiSW.index)
# pmiDict = {(s1,s2):pmi[s1][s2]
#            for s1 in sounds for s2 in sounds}
pmiDictPW = {(s1,s2):pmiPW[s1][s2]
                for s1 in soundsPW for s2 in soundsPW}
pmiDictSW = {(s1,s2):pmiSW[s1][s2]
                for s1 in soundsSW for s2 in soundsSW}

taxa = data.language.unique()
lpairs = pd.DataFrame([(l1,l2)
                       for i,l1 in enumerate(taxa)
                       for j,l2 in enumerate(taxa)
                       if i<j])


wpairsPW = pd.DataFrame()
wpairsSW = pd.DataFrame()
for l1,l2 in lpairs.iloc[random.permutation(lpairs.index)].values:
    l1Data = data[data.language==l1]
    l2Data = data[data.language==l2]
    lpPairs = pd.DataFrame([list(l1Data[['concept','language','word','ID']].loc[i])+
                            list(l2Data[['concept','language','word','ID']].loc[j])
                            for i in l1Data.index
                            for j in l2Data.index],
                           columns=['concept1','language1','word1','ID1',
                                    'concept2','language2','word2','ID2'])
    wpairsPW = pd.concat([wpairsPW,lpPairs])
    wpairsSW = pd.concat([wpairsSW,lpPairs])

wpairsPW['target'] = array(wpairsPW.concept1==wpairsPW.concept2,int)
wpairsSW['target'] = array(wpairsSW.concept1==wpairsSW.concept2,int)

wpairsPW['PMI'] = [sscore(a,b,pmiDictPW,gp1,gp2)
                 for (a,b) in wpairsPW[['word1','word2']].values]
wpairsSW['PMI'] = [sscore(a,b,pmiDictSW,gp1,gp2)
                    for (a,b) in wpairsSW[['word1','word2']].values]

wpairsPW.to_csv('albanoRomancePW.wordpairs.csv',index=False)
wpairsSW.to_csv('albanoRomanceSW.wordpairs.csv',index=False)

lrPW = LogisticRegression()
lrSW = LogisticRegression()
lrPW.fit(c_[wpairsPW.PMI.values],wpairsPW.target.values)
lrSW.fit(c_[wpairsSW.PMI.values],wpairsSW.target.values)


synpairsPW = wpairsPW[wpairsPW.target==1][['concept1',
                                     'language1', 'language2',
                                     'word1','word2',
                                     'ID1','ID2','PMI']]

synpairsPW.columns = ['concept']+list(synpairsPW.columns[1:])
conceptsPW = data.concept.unique()
synpairsPW['prediction'] = lrPW.predict_proba(c_[synpairsPW.PMI.values])[:,1]

synpairsSW = wpairsSW[wpairsSW.target==1][['concept1',
                                     'language1', 'language2',
                                     'word1','word2',
                                     'ID1','ID2','PMI']]

synpairsSW.columns = ['concept']+list(synpairsSW.columns[1:])
conceptsSW = data.concept.unique()
synpairsSW['prediction'] = lrSW.predict_proba(c_[synpairsSW.PMI.values])[:,1]


ccDataPW = pd.DataFrame()
ccDataSW = pd.DataFrame()
th = .5
for c in conceptsPW:
    cData = data[data.concept==c].copy()
    cPairs = synpairsPW[synpairsPW.concept==c]
    cIDs = cData.ID.values
    simMtx = zeros((len(cIDs),len(cIDs)))
    id_index = pd.Index(cIDs)
    idx1 = id_index.get_indexer(cPairs.ID1.values)
    idx2 = id_index.get_indexer(cPairs.ID2.values)

    simMtx[idx1, idx2] = cPairs.prediction.values
    simMtx[idx2, idx1] = cPairs.prediction.values
    simMtx[simMtx<th]=0
    G = igraph.Graph.Weighted_Adjacency(list(simMtx))
    clusters = G.community_label_propagation(weights='weight')
    ccDict = {cIDs[x]:i for i,cl in enumerate(clusters)
              for x in cl}
    cData['cc'] = [c+':'+str(ccDict[i]) for i in cData.ID.values]
    ccDataPW = pd.concat([ccDataPW,cData])

for c in conceptsSW:
    cData = data[data.concept==c].copy()
    cPairs = synpairsSW[synpairsSW.concept==c]
    cIDs = cData.ID.values
    simMtx = zeros((len(cIDs),len(cIDs)))
    id_index = pd.Index(cIDs)
    idx1 = id_index.get_indexer(cPairs.ID1.values)
    idx2 = id_index.get_indexer(cPairs.ID2.values)

    simMtx[idx1, idx2] = cPairs.prediction.values
    simMtx[idx2, idx1] = cPairs.prediction.values
    simMtx[simMtx<th]=0
    G = igraph.Graph.Weighted_Adjacency(list(simMtx))
    clusters = G.community_label_propagation(weights='weight')
    ccDict = {cIDs[x]:i for i,cl in enumerate(clusters)
              for x in cl}
    cData['cc'] = [c+':'+str(ccDict[i]) for i in cData.ID.values]
    ccDataSW = pd.concat([ccDataSW,cData])


#taxa = ccData.language.unique()
taxaSW = ccDataSW.language.unique
taxaPW = ccDataPW.language.unique
ccMtxSW = pd.DataFrame(index=taxa)
ccMtxPW = pd.DataFrame(index=taxa)

for c in conceptsPW:
    cData = ccDataPW[ccDataPW.concept==c]
    cMtx = pd.crosstab(cData.language,cData.cc)
    cMtx[cMtx>1]=1
    cMtx = cMtx.reindex(taxa,fill_value='-')
    ccMtxPW = pd.concat([ccMtxPW,cMtx],axis=1)

for c in conceptsSW:
    cData = ccDataSW[ccDataSW.concept==c]
    cMtx = pd.crosstab(cData.language,cData.cc)
    cMtx[cMtx>1]=1
    cMtx = cMtx.reindex(taxa,fill_value='-')
    ccMtxSW = pd.concat([ccMtxSW,cMtx],axis=1)


ccMtxPW.to_csv('albanoRomanceCCbinPW.csv')
ccMtxSW.to_csv('albanoRomanceCCbinSW.csv')

nexCharOutput(ccMtxPW.values,ccMtxPW.index,'albanoRomanceCC_PW.nex')
nexCharOutput(ccMtxSW.values,ccMtxSW.index,'albanoRomanceCC_SW.nex')


ccDataPW.to_csv('albanoRomanceCC_PW.csv',index='False')
ccDataSW.to_csv('albanoRomanceCC_SW.csv',index='False')
