import numpy as np
import pandas as pd

df = pd.read_csv('socher_albanoRomanceCC.csv')
pre_clusters = pd.Categorical(df['Cluster'], categories=df['Cluster'].unique(), ordered=True)
clusters = list(pre_clusters.categories)
languages = sorted(df['Language'].unique())

binary_matrix = np.zeros((len(languages), len(clusters)), dtype=int)

for i, language in enumerate(languages):
    language_clusters = df[df['Language'] == language]['Cluster'].tolist()
    for cluster in language_clusters:
        j = clusters.index(cluster)
        binary_matrix[i, j] = 1

bin_file_path = "socher_albanoRomanceCC.bin.csv"
# the bin file here doesn't have '-' for missing values, just a 0. it's because the socher cc I have doesn't handle them
binary_df = pd.DataFrame(binary_matrix, index=languages, columns=clusters)
binary_df.to_csv(bin_file_path, index=True)

nex_file_path = "socher_albanoRomanceCC.nex"


def create_nexus_from_matrix(matrix, names, output_file, datatype='STANDARD'):
    with open(output_file, 'w') as f:
        f.write('#NEXUS\n\n')
        f.write('BEGIN DATA;\n')
        f.write(f'DIMENSIONS ntax={len(matrix)} NCHAR={matrix.shape[1]};\n')
        f.write(f'FORMAT DATATYPE={datatype} GAP=? MISSING=- interleave=yes;\n')
        f.write('MATRIX\n\n')
        max_name_length = max(map(len, names))
        for i in range(len(matrix)):
            f.write(names[i].ljust(max_name_length + 2))
            for ch in matrix[i]:
                if ch == -1:
                    f.write('-')
                else:
                    f.write(str(ch))
            f.write('\n')
        f.write('\n;\n\nEND;\n')


create_nexus_from_matrix(binary_matrix, languages, nex_file_path)
