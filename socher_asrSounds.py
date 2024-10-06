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


# guideTreePW = Tree('socher_albanoRomancePW.mcc.nwk')
# guideTreeSW = Tree('socher_albanoRomanceSW.mcc.nwk')
guideTreeNW = Tree('socher_albanoRomance.mcc.nwk')
# guideTreePW_rt = remove_taxa(Tree('socher_albanoRomancePW.mcc.nwk'))
# guideTreeSW_rt = remove_taxa(Tree('socher_albanoRomanceSW.mcc.nwk'))
guideTreeNW_rt = remove_taxa(Tree('socher_albanoRomance.mcc.nwk'))

ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True

# guideTreePW.render('socher_guideTreePW.png', tree_style=ts)
# guideTreeSW.render('socher_guideTreeSW.png', tree_style=ts)
guideTreeNW.render('socher_guideTreeNW.png', tree_style=ts)
# guideTreePW_rt.render('socher_guideTreePW_rt.png', tree_style=ts)
# guideTreeSW_rt.render('socher_guideTreeSW_rt.png', tree_style=ts)
guideTreeNW_rt.render('socher_guideTreeNW_rt.png', tree_style=ts)

# asrCC_PW = pd.read_csv('dialign+socher+original/socher_asrCC_PW.csv', header=None, index_col=0)[1]
# asrCC_SW = pd.read_csv('sw+socher+original/socher_asrCC_SW.csv', header=None, index_col=0)[1]
asrCC_NW = pd.read_csv('original+socher+original/socher_asrCC_NW.csv', header=None, index_col=0)[1]

# romancePW = array(guideTreePW_rt.get_leaf_names())
# romanceSW = array(guideTreeSW_rt.get_leaf_names())
romanceNW = array(guideTreeNW_rt.get_leaf_names())

# dataPW = pd.read_csv('socher_albanoRomanceCC.csv', index_col=0)
# dataSW = pd.read_csv('socher_albanoRomanceCC.csv', index_col=0)
dataNW = pd.read_csv('socher_albanoRomanceCC.csv', index_col=0)

# dataPW = dataPW[(dataPW.Language.isin(romancePW))&(dataPW.Cluster.isin(asrCC_PW.values))]
# dataSW = dataSW[(dataSW.Language.isin(romanceSW))&(dataSW.Cluster.isin(asrCC_SW.values))]
dataNW = dataNW[(dataNW.Language.isin(romanceNW))&(dataNW.Cluster.isin(asrCC_NW.values))]


gp1=-2.49302792222
gp2=-1.70573165621

pmi = pd.read_csv('pmi-albanoRomance.csv',index_col=0)
sounds = array(pmi.index)

pmiDict = {(s1,s2):pmi[s1][s2]
           for s1 in sounds for s2 in sounds}


# conceptsPW = asrCC_PW.index
# conceptsSW = asrCC_SW.index
conceptsNW = asrCC_NW.index


aBlocksPW, aBlocksSW, aBlocksNW = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()


def aBlocks_binMtx_filler(concepts, data, asrCC, romance, aBlocks, guideTree):
    for c in concepts:
        cData = data[data.Cluster==asrCC[c]]
        if len(cData)>1:
            alg = tCoffee(cData.Language.values,cData.Word.values,
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


# aBlocksPW, binMtxPW = aBlocks_binMtx_filler(conceptsPW, dataPW, asrCC_PW, romancePW, aBlocksPW, guideTreePW_rt)
# aBlocksSW, binMtxSW = aBlocks_binMtx_filler(conceptsSW, dataSW, asrCC_SW, romanceSW, aBlocksSW, guideTreeSW_rt)
aBlocksNW, binMtxNW = aBlocks_binMtx_filler(conceptsNW, dataNW, asrCC_NW, romanceNW, aBlocksNW, guideTreeNW_rt)

# nexCharOutput(binMtxPW.to_numpy(),binMtxPW.index,'socher_romanceAlignmentsPW.nex')
# nexCharOutput(binMtxSW.to_numpy(),binMtxSW.index,'socher_romanceAlignmentsSW.nex')
nexCharOutput(binMtxNW.to_numpy(),binMtxNW.index,'socher_romanceAlignmentsNW.nex')

# binMtxPW.to_csv('socher_romanceAlignmentsPW.tsv',sep='\t',header=False)
# binMtxSW.to_csv('socher_romanceAlignmentsSW.tsv',sep='\t',header=False)
binMtxNW.to_csv('socher_romanceAlignmentsNW.tsv',sep='\t',header=False)


with open('btAlignments.txt','w') as f:
    f.write('1\n1\n')
    f.write('seed 12345;\n')
    f.write('mlt 1;\n')
    f.write('run\n')


# pPW = Popen('BayesTraitsV4 socher_romancePW.posterior.nex.tree socher_romanceAlignmentsPW.tsv < btAlignments.txt',
#           shell=True)
# os.waitpid(pPW.pid,0)
# pSW = Popen('BayesTraitsV4 socher_romanceSW.posterior.nex.tree socher_romanceAlignmentsSW.tsv < btAlignments.txt',
#           shell=True)
# os.waitpid(pSW.pid,0)
pNW = Popen('BayesTraitsV4 socher_romanceNW.posterior.nex.tree socher_romanceAlignmentsNW.tsv < btAlignments.txt',
            shell=True)
os.waitpid(pNW.pid,0)

# resultsPW = pd.read_csv('socher_romanceAlignmentsPW.tsv.log.txt',
#                       skiprows=31,sep='\t')
# resultsSW = pd.read_csv('socher_romanceAlignmentsSW.tsv.log.txt',
#                         skiprows=31,sep='\t')
resultsNW = pd.read_csv('socher_romanceAlignmentsNW.tsv.log.txt',
                        skiprows=31,sep='\t')


def result_mean(results, binMtx):
    idc = [x for x in results.columns if 'P(1)' in x]
    print(idc)
    results = results[idc]
    results.columns = binMtx.columns
    return results.mean()


# results_meanPW = result_mean(resultsPW, binMtxPW)
# results_meanSW = result_mean(resultsSW, binMtxSW)
results_meanNW = result_mean(resultsNW, binMtxNW)


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


# reconstruction_original_PW = pd.DataFrame(recon_original(results_meanPW, aBlocksPW, conceptsPW))
# reconstruction_original_SW = pd.DataFrame(recon_original(results_meanSW, aBlocksSW, conceptsSW))
reconstruction_original_NW = pd.DataFrame(recon_original(results_meanNW, aBlocksNW, conceptsNW))


asjp = pd.read_table('dataset.tab', index_col=0, sep='\t')


def recon_latin(reconstruction, asjp, concepts, filename):
    reconstruction['Latin'] = asjp.loc['LATIN'][concepts].values
    reconstruction.columns = ['reconstruction', 'Latin']
    reconstruction['concept'] = reconstruction.index
    reconstruction[['concept', 'Latin', 'reconstruction']].to_csv(f'{filename}.csv',
                                                                  index=False)
    return reconstruction


# reconstruction_results_PW = recon_latin(reconstruction_original_PW, asjp, conceptsPW, 'socher_reconstruction_results_PW')
# reconstruction_results_SW = recon_latin(reconstruction_original_SW, asjp, conceptsSW, 'socher_reconstruction_results_SW')
reconstruction_results_NW = recon_latin(reconstruction_original_NW, asjp, conceptsNW, 'socher_reconstruction_results_NW')
