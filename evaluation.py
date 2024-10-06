from numpy import *
import pandas as pd
import Levenshtein
import re

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

rec_socher_original = pd.read_csv('original+socher+original/socher_reconstruction_results_NW.csv',index_col=0)
rec_socher_dialign = pd.read_csv('dialign+socher+original/socher_reconstruction_results_PW.csv', index_col=0)
rec_socher_SW = pd.read_csv('sw+socher+original/socher_reconstruction_results_SW.csv', index_col=0)
rec_original_dialign = pd.read_csv('dialign+original_pipeline/reconstruction_results_PW.csv', index_col=0)
rec_original_sw = pd.read_csv('sw+original_pipeline/reconstruction_results_SW.csv', index_col=0)
rec_original = pd.read_csv('original_results/reconstruction.csv', index_col=0)

concepts_socher_original = array(rec_socher_original.index)
concepts_socher_dialign = array(rec_socher_dialign.index)
concepts_socher_SW = array(rec_socher_SW.index)
concepts_original_dialign = array(rec_original_dialign.index)
concepts_original_sw = array(rec_original_sw.index)
concepts_original = array(rec_original.index)


def evaluate_reconstruction(reconstruction, concepts, asjp, romance):
    def ldn(a, b):
        return min([1. * Levenshtein.distance(x, y) / max(len(x), len(y))
                    for x in a.split('-') for y in b.split('-')])

    romanceCleaned = pd.DataFrame([[cleanASJP(x).split('-')[0] for x in y]
                                   for y in asjp.loc[romance][concepts].values],
                                  index=romance,
                                  columns=concepts)

    latinCleaned = pd.Series([cleanASJP(x) for x in reconstruction.Latin.values],
                             index=reconstruction.index)[concepts]

    reconEval = mean([ldn(x, y) for x, y in zip(reconstruction.reconstruction.values,
                                                latinCleaned.values)])

    romanceEval = pd.Series([mean([ldn(romanceCleaned.loc[l][c], latinCleaned[c])
                                   for c in concepts])
                             for l in romance],
                            index=romance)
    romanceEval.loc['Proto-Romance'] = reconEval

    return romanceEval

reconstructions = {
    'socher_original': rec_socher_original,
    'socher_dialign': rec_socher_dialign,
    'socher_SW': rec_socher_SW,
    'original_dialign': rec_original_dialign,
    'original_sw': rec_original_sw,
    'original': rec_original
}

evaluation_results = pd.DataFrame({
    name: evaluate_reconstruction(reconstruction, globals()[f'concepts_{name}'], asjp, romance)
    for name, reconstruction in reconstructions.items()
})

evaluation_results.to_csv('romanceEvaluation_multiple.csv')


evaluation_results = pd.read_csv('romanceEvaluation_multiple.csv', index_col=0)
mean_scores = evaluation_results.mean()
best_pipeline = mean_scores.idxmax()
best_score = mean_scores.max()
print(f"The best-performing pipeline is: {best_pipeline} with a mean score of {best_score}")
# The best-performing pipeline is: original_dialign with a mean score of 0.6264740896358544
