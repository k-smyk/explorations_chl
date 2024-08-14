# Results are saved in the pairwise_alignment folder, along with the ipynb used to execute.
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
from lingpy import *
import lingpy as lingpy
import ast
ncores = 50


def cleanASJP(word):
    """
    Adapted and shortened from Jaeger's cleaning in alignment.py.
    """
    word = re.sub(r"[,\%\"~]", "", word)
    word = re.sub(r"(.)(.)(.)\$", r"\2", word)
    return word.replace('~', '')

def convert_NW_PW(alignment_str, keep_score=False):
    """
    Cleaning up for NW alignment.
    """
    alignment_str = alignment_str.strip()
    try:
        parsed_input = ast.literal_eval(alignment_str)
    except (ValueError, SyntaxError) as e:
        raise ValueError("Wrong alignment format.") from e

    lists1, lists2, score = parsed_input
    seqA = ''.join(lists1)
    seqB = ''.join(lists2)

    if keep_score:
        return seqA, seqB, score

    return seqA, seqB


# print(convert_NW_PW("(['p', 'o', 'L', '-', '-', '-'], ['-', '-', '-', 'b', 'i', 'e'], -6)", keep_score=True))
# print(convert_NW_PW("(['p', 'o', 'L', '-'], ['-', '-', 'Z', 'u'], -1.44181953287)"))
#
# print("This is PW: \n")
# print(convert_NW_PW("(['p', 'o', 'L', '-', '-', '-'], ['-', '-', '-', 'b', 'i', 'e'], 0.0)"))
# print(convert_NW_PW("(['f', 'e', 'j', '3', '-', '-', '-'], ['-', '-', '-', '-', 'Z', 'u', 'L'], 0.0)"))


def convert_SW(alignment, keep_score=False):
    GAP_SYMBOL = '-'
    if isinstance(alignment, str):
        try:
            parsed_input = ast.literal_eval(alignment)
        except (ValueError, SyntaxError) as e:
            raise ValueError("Wrong alignment format.") from e
    elif isinstance(alignment, tuple):
        parsed_input = alignment
    else:
        raise ValueError("Unsupported alignment format.")

    lists1, lists2, score = parsed_input
    max_cols = max(len(lists1), len(lists2))
    output_lists1 = []
    output_lists2 = []

    for i in range(max_cols):
        seqA = lists1[i] if i < len(lists1) else []
        seqB = lists2[i] if i < len(lists2) else []
        lenA = len(seqA)
        lenB = len(seqB)
        max_len = max(lenA, lenB)
        padded_seqA = ''.join(seqA) + GAP_SYMBOL * (max_len - lenA)
        padded_seqB = ''.join(seqB) + GAP_SYMBOL * (max_len - lenB)
        output_lists1.append(padded_seqA)
        output_lists2.append(padded_seqB)
    final_seqA = ''.join(output_lists1)
    final_seqB = ''.join(output_lists2)

    if keep_score:
        return final_seqA, final_seqB, score

    return final_seqA, final_seqB



# Function: nexCharOutput
# Description:
## This function takes a character array, a list of rownames
## and the name of the output nexus file as input
## and writes the character matrix into a nexus file.
## Missing entries are assumed to be coded as "-1"

def nexCharOutput(chMtx, names, outfile, datatype='STANDARD'):
    f = open(outfile, 'w')
    f.write('#NEXUS\n\n')
    f.write('BEGIN DATA;\n')
    f.write('DIMENSIONS ntax=' + str(len(chMtx)) + ' NCHAR=' + str(len(chMtx.T)) + ';\n')
    f.write('FORMAT DATATYPE=' + datatype + ' GAP=? MISSING=- interleave=yes;\n')
    f.write('MATRIX\n\n')
    txLgth = max(map(len, names))
    for i in range(len(chMtx)):
        f.write(names[i].ljust(txLgth + 2))
        for ch in chMtx[i]:
            if ch == -1:
                ch = '-'
            else:
                ch = str(ch)
            f.write(ch)
        f.write('\n')
    f.write('\n;\n\nEND;\n')
    f.close()


data = pd.read_csv('dataset.tab',
                   index_col=0, na_filter=False, sep='\t')
data = data[data.wls_gen.isin(['ROMANCE', 'ALBANIAN'])]
data = data[data.index != 'LATIN']
concepts100 = array(data.columns[9:])
nEntries = data[concepts100].apply(lambda x: sum(x != '')).sort_values()
concepts = nEntries.index[-40:]
data = data[concepts]
taxa = array(data.index)

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


def recalculate_pmi(alignment, gp1, gp2):
    """Recalculate the PMI score based on gap penalties."""
    print(alignment)
    score = alignment[2]
    gap_count = 0
    seqA, seqB = alignment[0], alignment[1]
    for i in range(len(seqA)):
        if seqA[i] == '-' or seqB[i] == '-':
            gap_count += 1
            if gap_count == 1:
                score -= 1
                score += gp1
            else:
                score -= 1
                score += gp2
    return score


def sw_align_lingpy(w, scorer=False, gp1=None, gp2=None, keep_scores=False):
    x, y = w
    sw_aln = lingpy.align.pairwise.sw_align(x, y, scorer=scorer)
    sw_aln = convert_SW(sw_aln, keep_score=keep_scores)
    if gp1 and gp2:
        sw_aln = (sw_aln[0], sw_aln[1], recalculate_pmi(sw_aln, gp1, gp2))
    return sw_aln


def nw_align_lingpy(w, scorer=False, gp1=None, gp2=None, keep_scores=False):
    x, y = w
    nw_aln = lingpy.align.pairwise.nw_align(x, y, scorer=scorer)
    nw_aln = convert_NW_PW(nw_aln, keep_score=keep_scores)
    if gp1 and gp2:
        nw_aln = (nw_aln[0], nw_aln[1], recalculate_pmi(nw_aln, gp1, gp2))
    return nw_aln


def pw_align_lingpy(w, scorer=False, gp1=None, gp2=None, keep_scores=False):
    x, y = w
    pw_aln = lingpy.align.pairwise.pw_align(x, y, mode="dialign")
    pw_aln = convert_NW_PW(pw_aln, keep_score=keep_scores)
    if gp1 and gp2:
        sw_aln = (pw_aln[0], pw_aln[1], recalculate_pmi(pw_aln, gp1, gp2))
    return pw_aln



manager = Manager()
return_dict = manager.dict()
packages = array_split(training.values,ncores)


def doWork(i,pck):
    results = []
    for p in pck:
        alignment = nw_align_lingpy(p)
        if alignment:
            results.append(alignment)
    print(f'Package {i}: Number of successful alignments - {len(results)}') # Print number of successful alignments
    if results:
        return_dict[i] = vstack(results)
        print(f'Package {i}: Results added to return_dict') # Confirm results added to dictionary
    else:
        print(f"Warning: No successful alignments for package {i}")


jobs = []
for i,pck in enumerate(packages):
    p = Process(target=doWork,args=(i,pck))
    p.start()
    jobs.append(p)
for p in jobs:
    p.join()
if return_dict:
  alg0 = vstack(return_dict.values())
else:
  print("Error: the dict is empty")

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
pmi0 = (np.log(sProbs).copy()-np.log(soundProbabilities.to_numpy().reshape(-1, 1))).T-np.log(soundProbabilities.to_numpy().reshape(-1, 1))

def sscore(a, b, pmiDict, gp1, gp2):
    """a,b: ASJP strings
    pmiDict: logodds dictionary
    gp1,gp2: gap penalties
    return PMI score of a/b
    """
    out = pairwise2.align.globalds(a, b, pmiDict, gp1, gp2)
    if len(out) == 0: return nan
    return out[0][2]


def scoreNW(x, y, pmiDict, gp1, gp2):
    """x,y: sequences of ASJP strings, separated by '-'
    pmiDict: logodds dictionary
    gp1,g2: gap penalties
    returns maximal PMI score for the Cartesian product of x and y"""
    if '0' in [x, y]: return nan
    x1 = x.split('-')
    y1 = y.split('-')
    return max([sscore(xx, yy, pmiDict, gp1, gp2) for xx in x1 for yy in y1])


pmi0Dict = {(s1, s2): pmi0[s1][s2]
            for s1 in sounds for s2 in sounds}

def nw(x, y, pmiDict, gp1, gp2):
    """wrapper for Bio.pairwise2.align.globalds"""
    return pairwise2.align.globalds(x, y, pmiDict, gp1, gp2)


# def nwalign(w, pmiDict, gp1, gp2, th=-Inf):
# 	"""w: pair of ASJP strings
# 	pmiDict: dictionary of logodds (=PMI scores)
# 	gp1,gp2: gap penalties (non-positive)
# 	th: threshold; all pairs with a PMI-score <th will be ignored
# 	returns: array of pairwise alignment, with gaps removed"""
# 	x, y = w
# 	a = nw(x, y, pmiDict, gp1, gp2)
# 	if len(a) == 0: return zeros((0, 2))
# 	algn = []
# 	if a[0][2] < th:
# 		return zeros((0, 2))
# 	for aa in a:
# 		l = len(aa[0])
# 		aaa = [[aa[0][i], aa[1][i]] for i in range(l)]
# 		algn += [x for x in aaa if not '-' in x]
# 	return array(algn)

def align_no_gaps(w, pmiDict, gp1, gp2, th=-Inf, mode=None):
    """w: pair of ASJP strings
    pmiDict: dictionary of logodds (=PMI scores)
    gp1,gp2: gap penalties (non-positive)
    th: threshold; all pairs with a PMI-score <th will be ignored
    returns: array of pairwise alignment, with gaps removed"""
    if mode == 'nw':
        a = nw_align_lingpy(w, scorer=pmiDict, keep_scores=True)
    elif mode == 'sw':
        a = sw_align_lingpy(w, scorer=pmiDict, keep_scores=True)
    elif mode == 'pw':
        a = pw_align_lingpy(w, keep_scores=True)
    if len(a) == 0: return zeros((0, 2))
    algn = []
    # Check if the alignment score is below the threshold, handling potential NaN
    if a[2] < th if not isnan(a[2]) else False:
        return zeros((0, 2))
    # Extract aligned pairs without gaps
    for aa in a:
        l = len(aa[0])
        aaa = [[aa[0][i], aa[1][i]] for i in range(l)]
        algn += [x for x in aaa if not '-' in x]
    return array(algn)
    # if a[0][2] < th:
    #     return zeros((0, 2))
    # for aa in a:
    #     l = len(aa[0])
    #     aaa = [[aa[0][i], aa[1][i]] for i in range(l)]
    #     algn += [x for x in aaa if not '-' in x]
    # return array(algn)


def alignStar(crp, pmiDict, gp1, gp2, th, align_func):
    packages = array_split(crp, ncores)
    manager = Manager()
    return_dict = manager.dict()

    def doWork(i, pck):
        return_dict[i] = vstack([align_no_gaps(w, pmiDict, gp1, gp2, th, mode=align_func)
                                 for w in pck])

    jobs = []
    for i, pck in enumerate(packages):
        p = Process(target=doWork, args=(i, pck))
        p.start()
        jobs.append(p)
    for p in jobs:
        p.join()
    return vstack(return_dict.values())


def nwalignStar(crp, pmiDict, gp1, gp2, th):
    return alignStar(crp, pmiDict, gp1, gp2, th, 'nw')


def swalignStar(crp, pmiDict, gp1, gp2, th):
    return alignStar(crp, pmiDict, gp1, gp2, th, 'sw')


def pwalignStar(crp, pmiDict, gp1, gp2, th):
    return alignStar(crp, pmiDict, gp1, gp2, th, 'pw')


gp1 = -2.49302792222
gp2 = -1.70573165621
th = 4.4451

pmiDict = pmi0Dict.copy()
for i in range(10):
    nw_alg = nwalignStar(training.values, pmiDict, gp1, gp2, th)
    sw_alg = swalignStar(training.values, pmiDict, gp1, gp2, th)
    pw_alg = pwalignStar(training.values, pmiDict, gp1, gp2, th)
    for alg, filename in zip([nw_alg, sw_alg, pw_alg], ['nw', 'sw', 'pw']):
        sFreqs = pd.crosstab(alg[:, 0], alg[:, 1])
        sFreqs = sFreqs.reindex(sounds, fill_value=0).T
        sFreqs = sFreqs.reindex(sounds, fill_value=0).T
        sFreqs = sFreqs.copy() + sFreqs.copy().T
        sFreqs += 1
        sProbs = sFreqs / sFreqs.sum().sum()
        pmi = (log(sProbs).copy() - log(soundProbabilities)).T - log(soundProbabilities)
        pmiDict = {(s1, s2): pmi[s1][s2]
                   for s1 in sounds for s2 in sounds}
        pmi.to_csv(f'pmi-albanoRomance-{filename}.csv')

dataWL.to_csv('albanoRomanceASJP_new.csv', index=False)

sc = pd.DataFrame(index=taxa)
for c in concepts:
    cData = dataWL[dataWL.concept == c]
    cTaxa = cData.language.unique()
    cWords = pd.Series([''.join(cData[cData.language == l].word.values)
                        for l in cTaxa],
                       index=cTaxa)
    cMtx = pd.DataFrame([[int(s in cWords[l]) for s in sounds]
                         for l in cTaxa],
                        index=cTaxa)
    cMtx = cMtx.reindex(taxa, fill_value='-')
    sc = pd.concat([sc, cMtx], axis=1)

nexCharOutput(sc.values, sc.index, 'albanoRomanceSC_new.nex', 'restriction')
