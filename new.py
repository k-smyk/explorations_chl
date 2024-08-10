import lingpy
import pandas as pd
import re
from lingpy import *
from Bio import pairwise2
import numpy as np
import ast


def cleanASJP(word):
    """
    Adapted and shortened from Jaeger's cleaning in alignment.py.
    """
    word = re.sub(r"[,\%\"~]", "", word)
    word = re.sub(r"(.)(.)(.)\$", r"\2", word)
    return word.replace('~', '')


data = pd.read_csv('albanoRomance.wordpairs.csv')
data = data.drop(columns=['PMI'])
data['word1'] = data['word1'].apply(cleanASJP)
data['word2'] = data['word2'].apply(cleanASJP)

print(data.head())

# pmi dict
pmi_df = pd.read_csv('pmi-albanoRomance.csv', index_col=0)
pmi_dict = {}
for i in range(pmi_df.shape[0]):
    for j in range(pmi_df.shape[1]):
        pmi_dict[(pmi_df.columns[j], pmi_df.index[i])] = pmi_df.iloc[i, j]



def align_word_pairs(df):
    alignments = []

    for index, row in df.iterrows():
        seqA = row['word1']
        seqB = row['word2']
        try:
            nw_aln = align.pairwise.nw_align(seqA, seqB, scorer=pmi_dict)
            sw_aln = align.pairwise.sw_align(seqA, seqB, scorer=pmi_dict)
            pw_aln = align.pairwise.pw_align(seqA, seqB, mode="dialign")
        except:
            continue

        alignments.append({
            'concept1': row['concept1'],
            'language1': row['language1'],
            'word1': row['word1'],
            'ID1': row['ID1'],
            'concept2': row['concept2'],
            'language2': row['language2'],
            'word2': row['word2'],
            'ID2': row['ID2'],
            'nw_alignment': nw_aln,
            'sw_alignment': sw_aln,
            'pw_alignment': pw_aln
        })

    return pd.DataFrame(alignments)


alignments_df = align_word_pairs(data)
alignments_df.to_csv('aligned_lingpy.csv', index=False)





# turchin only returns 0 or 1, not useful

def convert_NW_PW(alignment_str):
    """
    Cleaning up for NW alignment.
    """
    alignment_str = alignment_str.strip()
    parts = re.findall(r"\[([^\]]*)\]", alignment_str)

    if len(parts) < 2:
        print("Error: Unexpected format.")
        pass

    seqA_parts = [p.strip().replace("'", "").replace(" ", "") for p in parts[0].split(',')]
    seqB_parts = [p.strip().replace("'", "").replace(" ", "") for p in parts[1].split(',')]
    seqA = ''.join(seqA_parts)
    seqB = ''.join(seqB_parts)

    return seqA, seqB


print(convert_NW_PW("(['p', 'o', 'L', '-', '-'], ['-', '-', '-', 't', 'i'], -5)"))
print(convert_NW_PW("(['p', 'o', 'L', '-'], ['-', '-', 'Z', 'u'], -1.44181953287)"))

print("This is PW: \n")
print(convert_NW_PW("(['p', 'o', 'L', '-', '-', '-'], ['-', '-', '-', 'b', 'i', 'e'], 0.0)"))
print(convert_NW_PW("(['f', 'e', 'j', '3', '-', '-', '-'], ['-', '-', '-', '-', 'Z', 'u', 'L'], 0.0)"))


def convert_SW(alignment_str):
    GAP_SYMBOL = '-'
    try:
        parsed_input = ast.literal_eval(alignment_str)
    except (ValueError, SyntaxError) as e:
        raise ValueError("Input string does not match the expected format.") from e

    lists1, lists2, _ = parsed_input
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
    return final_seqA, final_seqB


print(convert_SW("((['p', 'o', 'L'], [], []), (['t', 'i'], [], []), 0)"))
print(convert_SW("((['f'], ['o'], ['L', '3']), (['a', 'i', 'g'], ['o'], []), 2.06525264777)"))


