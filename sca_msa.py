import lingpy
import pandas as pd

# path = "albanoRomanceASJP.csv"
# df = pd.read_csv(path, delimiter=',', quotechar='"')
#
# df['ID'] = range(1, len(df)+1)
# df['ipa'] = df['word']
# df['cogid'] = range(1, len(df) + 1)
#
# df.to_csv('albanoRomanceASJP_with_ID.csv', index=False)
# path = "albanoRomanceASJP_with_ID.csv"
# wl = lingpy.basic.wordlist.get_wordlist(path, delimiter=',', quotechar='"', row="concept", col="language")
# alms = lingpy.align.sca.Alignments(wl)
# alms.align(model='asjp')
# alms.output('tsv', filename='albanoRomanceASJP_lingpy_msa')


# msa_file = 'albanoRomanceASJP_lingpy_msa.tsv'
# msa_data = pd.read_csv(msa_file, sep='\t', comment='#')
# msa_data = msa_data[['CONCEPT', 'DOCULECT', 'ALIGNMENT']]
# msa_data['BINARY'] = msa_data['ALIGNMENT'].replace({r"[a-zA-Z]+": "1", r"[^\d-]": "0", "-": "-"}, regex=True)
#
# binary_alignment = (
#     msa_data['ALIGNMENT']
#     .replace({r"[a-zA-Z]": "1", "-": "0"}, regex=True)
#     .apply(lambda x: ''.join(['1' if char.isalpha() else ('0' if char == '-' else char) for char in x]))
# )
#
# binary_df = msa_data.groupby('DOCULECT')['BINARY'].apply(lambda x: ''.join(x)).reset_index()
#
# for index, row in binary_df.iterrows():
#     binary_string = ''
#     for char in row['BINARY']:
#         if char == '1':
#             binary_string += '1'
#         elif char == '0':
#             binary_string += '0'
#         else:
#             binary_string += '-'
#     binary_df.at[index, 'BINARY'] = binary_string
#
# # Convert to NEXUS format
# output_lines = [
#     "#NEXUS\n\nBEGIN DATA;",
#     f"DIMENSIONS ntax={binary_df['DOCULECT'].nunique()} NCHAR={binary_df['BINARY'].str.len().max()};",
#     "FORMAT DATATYPE=STANDARD GAP=? MISSING=- interleave=yes;",
#     "MATRIX\n"
# ]
#
# for _, row in binary_df.iterrows():
#     output_lines.append(f"{row['DOCULECT']:<20} {row['BINARY']}")
#
# output_lines.append(";\nEND;")
#
# output_content = "\n".join(output_lines)
# print(output_content)
#
# with open('albanoRomanceASJP_lingpy_msa.nex', 'w') as f:
#     f.write(output_content)
