import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


def evaluate_pmi_differences(original_file, alignment_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    df1 = pd.read_csv(original_file, header=0)
    alignment_files = [f for f in os.listdir(alignment_dir) if f.startswith('pmi-albanoRomance')]

    for alignment_file in alignment_files:
        df2 = pd.read_csv(os.path.join(alignment_dir, alignment_file), header=0)
        common_rows = df1['col_0'].isin(df2['col_0'])
        common_columns = df1.columns.intersection(df2.columns)
        df1_common = df1[common_rows].reset_index(drop=True).set_index('col_0')[common_columns[1:]]
        df2_common = df2[df2['col_0'].isin(df1['col_0'])].reset_index(drop=True).set_index('col_0')[common_columns[1:]]
        df2_common = df2_common.loc[df1_common.index]

        plt.figure(figsize=(12, 8))
        plt.subplot(2, 2, 1)
        pmi_differences = df1_common - df2_common
        sns.heatmap(pmi_differences, cmap='coolwarm', center=0)
        plt.title(f'Heatmap of PMI Differences for {alignment_file}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'pmi_difference_{alignment_file}.png'), bbox_inches='tight')
        plt.close()

        avg_pmi_original = df1_common.mean().mean()
        avg_pmi_new = df2_common.mean().mean()
        avg_pmi_difference = (df2_common - df1_common).mean().mean()
        print(f"File: {alignment_file}")
        print(f"Average PMI Score (Original): {avg_pmi_original:.4f}")
        print(f"Average PMI Score {alignment_file}: {avg_pmi_new:.4f}")
        print(f"Average Difference in PMI Scores ({alignment_file} - Original): {avg_pmi_difference:.4f}")


evaluate_pmi_differences('pmi-albanoRomance.csv', 'pairwise_alignment', 'pairwise_alignment/evaluation_output')

"""
File: pmi-albanoRomance-pw.csv
Average PMI Score (Original): -2.3855
Average PMI Score pmi-albanoRomance-pw.csv: -0.4156
Average Difference in PMI Scores (pmi-albanoRomance-pw.csv - Original): 1.9699

File: pmi-albanoRomance-sw.csv
Average PMI Score (Original): -2.3855
Average PMI Score pmi-albanoRomance-sw.csv: -1.9573
Average Difference in PMI Scores (pmi-albanoRomance-sw.csv - Original): 0.4282
"""
