from numpy import *
import pandas as pd
import Levenshtein
import re

from pairwise_alignment import iwsa


def cleanASJP(word):
    """takes an ASJP string as argument
    and returns the string with all diacritics removed."""
    word = re.sub(r",","-",word)
    word = re.sub(r"\%","",word)
    word = re.sub(r"\*","",word)
    word = re.sub(r"\"","",word)
    word = re.sub(r".~","",word)
    word = re.sub(r"(.)(.)(.)\$",r"\2",word)
    word = re.sub(r"\$","",word)
    word = re.sub(r"\s+","",word)
    return word.replace('~','')


asjp = pd.read_table('dataset.tab',index_col=0,
                     sep='\t',na_filter=False)

romance = array([x for x in asjp[asjp.wls_gen=='ROMANCE'].index
                 if x!='LATIN'])

reconstruction = pd.read_csv('reconstruction.csv',index_col=0)

concepts = array(reconstruction.index)

def ldn(a,b):
    return min([1.*Levenshtein.distance(x,y)/max(len(x),len(y))
                for x in a.split('-') for y in b.split('-')])


romanceCleaned = pd.DataFrame([[cleanASJP(x).split('-')[0] for x in y]
                               for y in asjp.loc[romance][concepts].values],
                              index=romance,
                              columns=concepts)

latinCleaned = pd.Series([cleanASJP(x) for x in reconstruction.Latin.values],
                         index=reconstruction.index)[concepts]

symbols = set(''.join(romanceCleaned.values.flatten()) + ''.join(latinCleaned.values.flatten()))
symbol_table = iwsa.PhoneticSymbolTable(symbols)
glo_corr_model = iwsa.CorrespondenceModel(symbol_table)
loc_corr_model = iwsa.CorrespondenceModel(symbol_table)
self_sim_model1 = iwsa.CorrespondenceModel(symbol_table)
self_sim_model2 = iwsa.CorrespondenceModel(symbol_table)
info_model1 = iwsa.CorrespondenceModel(symbol_table)
info_model2 = iwsa.CorrespondenceModel(symbol_table)

MATCH_SCORE = 2
MISMATCH_PENALTY = -1
GAP_PENALTY = -2


def initialize_model(model, symbols):
    for s1 in symbols:
        for s2 in symbols:
            if s1 == s2:
                model.set_score(s1, s2, MATCH_SCORE)
            else:
                model.set_score(s1, s2, MISMATCH_PENALTY)
    # gap penalties
    model.set_score(iwsa.PhoneticSymbolTable.EMPTY_SYMBOL, iwsa.PhoneticSymbolTable.EMPTY_SYMBOL, 0)
    for s in symbols:
        model.set_score(s, iwsa.PhoneticSymbolTable.EMPTY_SYMBOL, GAP_PENALTY)
        model.set_score(iwsa.PhoneticSymbolTable.EMPTY_SYMBOL, s, GAP_PENALTY)


# Initialize each correspondence model
initialize_model(glo_corr_model, symbols)
initialize_model(loc_corr_model, symbols)
initialize_model(self_sim_model1, symbols)
initialize_model(self_sim_model2, symbols)
initialize_model(info_model1, symbols)
initialize_model(info_model2, symbols)


def iwsa_distance(a, b):
    if not a or not b:  # check for empty strings
        return 1.0  # maximum distance if one of the strings is empty
    score, (a, b) = iwsa(a, b, glo_corr_model, loc_corr_model, self_sim_model1, self_sim_model2, info_model1, info_model2, -1, -1)
    return 1.0 - score / max(len(a), len(b))

print(iwsa_distance(['u', 'r', '-', 'a', 'k', '-', '-', '-', '-', 'a'],
['s', 'i', 'n', 'a', 'k', 's', 'i', 'n', 'a', 'k']))

reconEval = mean([iwsa_distance(x, y) for x, y in zip(reconstruction.reconstruction.values, latinCleaned.values)])

romanceEval = pd.Series(
    [mean([iwsa_distance(romanceCleaned.loc[l][c], latinCleaned[c]) for c in concepts]) for l in romance],
    index=romance
)
romanceEval.loc['Proto-Romance'] = reconEval

romanceEval.to_csv('romanceEvaluation.csv')
