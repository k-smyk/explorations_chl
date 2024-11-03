import pandas as pd
import re
from numpy import array, mean
import Levenshtein

def cleanASJP(word):
    word = re.sub(r"[,%*\"~]", "", word)
    word = re.sub(r"\s+", "", word)
    word = re.sub(r"(.)(.)(.)\$", r"\2", word)
    return word

asjp = pd.read_table('dataset.tab', index_col=0, sep='\t', na_filter=False)
romance = array([x for x in asjp[asjp.wls_gen=='ROMANCE'].index if x != 'LATIN'])

recon_files = {
    'socher_original': 'original+socher+original/socher_reconstruction_results_NW.csv',
    'socher_dialign': 'dialign+socher+original/socher_reconstruction_results_PW.csv',
    'socher_SW': 'sw+socher+original/socher_reconstruction_results_SW.csv',
    'original_dialign': 'dialign+original_pipeline/reconstruction_results_PW.csv',
    'original_sw': 'sw+original_pipeline/reconstruction_results_SW.csv',
    'original': 'original_results/reconstruction.csv',
    'weighted_PW': 'reconstruction_results_weighted_PW.csv',
    'weighted_SW': 'reconstruction_results_weighted_SW.csv'
}
reconstructions = {name: pd.read_csv(file, index_col=0) for name, file in recon_files.items()}
concepts = {name: array(reconstruction.index) for name, reconstruction in reconstructions.items()}

def ldn(a, b):
    return min([1.*Levenshtein.distance(x, y)/max(len(x), len(y))
                for x in a.split('-') for y in b.split('-')])

def evaluate_reconstruction(reconstruction, concepts, asjp, romance):
    romance_cleaned = pd.DataFrame([[cleanASJP(x).split('-')[0] for x in y]
                                   for y in asjp.loc[romance][concepts].values],
                                  index=romance, columns=concepts)
    latin_cleaned = pd.Series([cleanASJP(x) for x in reconstruction.Latin.values], index=reconstruction.index)[concepts]

    recon_eval = mean([ldn(x, y) for x, y in zip(reconstruction.reconstruction.values, latin_cleaned.values)])
    romance_eval = pd.Series([mean([ldn(romance_cleaned.loc[l][c], latin_cleaned[c]) for c in concepts]) for l in romance], index=romance)
    romance_eval.loc['Proto-Romance'] = recon_eval

    return romance_eval


evaluation_results = pd.DataFrame({name: evaluate_reconstruction(reconstruction, concepts[name], asjp, romance) for name, reconstruction in reconstructions.items()})
evaluation_results.to_csv('romanceEvaluation_multiple.csv')

proto_romance_scores = evaluation_results.loc['Proto-Romance']
best_pipeline = proto_romance_scores.idxmin()
best_score = proto_romance_scores.min()
worst_pipeline = proto_romance_scores.idxmax()
worst_score = proto_romance_scores.max()

print(f"The best-performing pipeline for Proto-Romance is: {best_pipeline} with a score of {best_score}")
print(f"The worst-performing pipeline for Proto-Romance is: {worst_pipeline} with a score of {worst_score}")

# one file with all pipelines
merged_df = reconstructions['socher_original'].merge(
    reconstructions['socher_dialign'],
    on=['concept', 'Latin'], suffixes=('_socher_original', '_socher_dialign'))
for name, df in reconstructions.items():
    if name not in ['socher_original', 'socher_dialign']:
        merged_df = merged_df.merge(df, on=['concept', 'Latin'], suffixes=('', f'_{name}'))
merged_df.to_csv('all_reconstruction_results.csv', index=False)
# merged_df = reconstructions['socher_original'].merge(
#     reconstructions['socher_dialign'],
#     on=['concept', 'Latin'], suffixes=('_original', '_dialign')
# )
#
# merged_df = merged_df.merge(reconstructions['socher_SW'], on=['concept', 'Latin'])
# merged_df.rename(columns={'reconstruction': 'reconstruction_SW'}, inplace=True)
#
# # Export differing reconstructions between pipelines
# diff_df = merged_df[(merged_df['reconstruction_original'] != merged_df['reconstruction_dialign']) |
#                     (merged_df['reconstruction_original'] != merged_df['reconstruction_SW']) |
#                     (merged_df['reconstruction_dialign'] != merged_df['reconstruction_SW'])]
# diff_df.to_csv('socher_reconstruction_differences.csv', index=False)
