from numpy import *
import pandas as pd
from ete3 import Tree
from alignment import tCoffee
from subprocess import Popen
import os

import PyQt5 # for treestyle
from ete3 import TreeStyle


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


taxa_to_exclude = {'ALBANIAN', 'ALBANIAN_TOSK', 'ALBANIAN_GHEG'}
taxa_to_exclude_set = set(taxa_to_exclude)


def remove_taxa(tree, taxa_to_exclude=taxa_to_exclude_set):
    def recurse_remove(node):
        if node.is_leaf() and node.name in taxa_to_exclude:
            node.delete()
        else:
            for child in node.get_children():
                recurse_remove(child)
    recurse_remove(tree)
    return tree


guideTreePW = Tree('albanoRomancePW.mcc.nwk')
guideTreeSW = Tree('albanoRomanceSW.mcc.nwk')
guideTreePW_rt = remove_taxa(Tree('albanoRomancePW.mcc.nwk'))
guideTreeSW_rt = remove_taxa(Tree('albanoRomanceSW.mcc.nwk'))

ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True

guideTreePW.render('guideTreePW.png', tree_style=ts)
guideTreeSW.render('guideTreeSW.png', tree_style=ts)
guideTreePW_rt.render('guideTreePW_rt.png', tree_style=ts)
guideTreeSW_rt.render('guideTreeSW_rt.png', tree_style=ts)

asrCC_PW = pd.read_csv('dialign+original_pipeline/asrCC_PW.csv', header=None, index_col=0)[1]
asrCC_SW = pd.read_csv('sw+original_pipeline/asrCC_SW.csv', header=None, index_col=0)[1]

romancePW = array(guideTreePW_rt.get_leaf_names())
romanceSW = array(guideTreeSW_rt.get_leaf_names())

dataPW = pd.read_csv('dialign+original_pipeline/albanoRomanceCC_PW.csv', index_col=0)
dataSW = pd.read_csv('sw+original_pipeline/albanoRomanceCC_SW.csv', index_col=0)

dataPW = dataPW[(dataPW.language.isin(romancePW))&(dataPW.cc.isin(asrCC_PW.values))]
dataSW = dataSW[(dataSW.language.isin(romanceSW))&(dataSW.cc.isin(asrCC_SW.values))]


gp1=-2.49302792222
gp2=-1.70573165621

pmi = pd.read_csv('pmi-albanoRomance.csv',index_col=0)
sounds = array(pmi.index)

pmiDict = {(s1,s2):pmi[s1][s2]
           for s1 in sounds for s2 in sounds}


conceptsPW = asrCC_PW.index
conceptsSW = asrCC_SW.index


aBlocksPW, aBlocksSW = pd.DataFrame(), pd.DataFrame()


def aBlocks_binMtx_filler(concepts, data, asrCC, romance, aBlocks, guideTree):
    for c in concepts:
        cData = data[data.cc==asrCC[c]]
        if len(cData)>1:
            alg = tCoffee(cData.language.values,cData.word.values,
                          guideTree,pmiDict,gp1,gp2,sounds)
            cAlg = pd.DataFrame([list(x[1]) for x in alg],
                                index = [x[0] for x in alg])
            cAlg[cAlg=='-'] = '0'
        else:
            cAlg = pd.DataFrame(map(list,cData.word.values),index=cData.language.values)
        cAlg = cAlg.reindex(romance,fill_value='-')
        cAlg.columns = [c+':'+str(i) for i in cAlg.columns]
        aBlocks = pd.concat([aBlocks,cAlg],axis=1)

    binMtx = pd.DataFrame(index=romance)
    for i in aBlocks.columns:
        cl = aBlocks[i].values
        states = unique([x for x in cl if x != '-'])
        clMtx = pd.DataFrame([cl == s for s in states], dtype=int,
                             columns=romance,
                             index=[i + ':' + s for s in states]).T
        clMtx = clMtx.astype(object)
        clMtx[cl == '-'] = '-'
        binMtx = pd.concat([binMtx, clMtx], axis=1)
    return aBlocks, binMtx


aBlocksPW, binMtxPW = aBlocks_binMtx_filler(conceptsPW, dataPW, asrCC_PW, romancePW, aBlocksPW, guideTreePW_rt)
aBlocksSW, binMtxSW = aBlocks_binMtx_filler(conceptsSW, dataSW, asrCC_SW, romanceSW, aBlocksSW, guideTreeSW_rt)

#nexCharOutput(binMtxPW.values,binMtxPW.index,'romanceAlignmentsPW.nex')
nexCharOutput(binMtxPW.to_numpy(),binMtxPW.index,'romanceAlignmentsPW.nex')
#nexCharOutput(binMtxSW.values,binMtxSW.index,'romanceAlignmentsSW.nex')
nexCharOutput(binMtxSW.to_numpy(),binMtxSW.index,'romanceAlignmentsSW.nex')

binMtxPW.to_csv('romanceAlignmentsPW.tsv',sep='\t',header=False)
binMtxSW.to_csv('romanceAlignmentsSW.tsv',sep='\t',header=False)


with open('btAlignments.txt','w') as f:
    f.write('1\n1\n')
    f.write('seed 12345;\n')
    f.write('mlt 1;\n')
    f.write('run\n')


pPW = Popen('BayesTraitsV4 romancePW.posterior.nex.tree romanceAlignmentsPW.tsv < btAlignments.txt',
          shell=True)
os.waitpid(pPW.pid,0)
pSW = Popen('BayesTraitsV4 romanceSW.posterior.nex.tree romanceAlignmentsSW.tsv < btAlignments.txt',
          shell=True)
os.waitpid(pSW.pid,0)

resultsPW = pd.read_csv('romanceAlignmentsPW.tsv.log.txt',
                      skiprows=31,sep='\t')
resultsSW = pd.read_csv('romanceAlignmentsSW.tsv.log.txt',
                        skiprows=31,sep='\t')


def result_mean(results, binMtx):
    idc = [x for x in results.columns if 'P(1)' in x]
    print(idc)
    results = results[idc]
    results.columns = binMtx.columns
    return results.mean()


results_meanPW = result_mean(resultsPW, binMtxPW)
results_meanSW = result_mean(resultsSW, binMtxSW)


def recon_original(res, aBlocks, concepts):
    asr = []
    for x in aBlocks.columns:
        idc = [y for y in res.index if ':'.join(y.split(':')[:2])==x]
        asr.append(res[idc].sort_values().index[-1].split(':')[-1])
    asr = pd.Series(asr,index=aBlocks.columns)
    reconstruction = []
    for c in concepts:
        idc = [x for x in asr.index if x.split(':')[0]==c]
        reconstruction.append(''.join(asr[idc].values).replace('0',''))
    reconstruction = pd.Series(reconstruction,index=concepts)
    return reconstruction


reconstruction_original_PW = pd.DataFrame(recon_original(results_meanPW, aBlocksPW, conceptsPW))
reconstruction_original_SW = pd.DataFrame(recon_original(results_meanSW, aBlocksSW, conceptsSW))


asjp = pd.read_table('dataset.tab', index_col=0, sep='\t')


def recon_latin(reconstruction, asjp, concepts, filename):
    reconstruction['Latin'] = asjp.loc['LATIN'][concepts].values
    reconstruction.columns = ['reconstruction', 'Latin']
    reconstruction['concept'] = reconstruction.index
    reconstruction[['concept', 'Latin', 'reconstruction']].to_csv(f'{filename}.csv',
                                                                  index=False)
    return reconstruction


reconstruction_results_PW = recon_latin(reconstruction_original_PW, asjp, conceptsPW, 'reconstruction_results_PW')
reconstruction_results_SW = recon_latin(reconstruction_original_SW, asjp, conceptsSW, 'reconstruction_results_SW')
