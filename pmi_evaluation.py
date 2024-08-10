import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df1 = pd.read_csv('pmi-albanoRomance.csv', header=0)
df2 = pd.read_csv('pairwise_alignment/unfinished_iwsa/pmi-albanoRomance_IWSA.csv', header=0)
common_rows = df1['col_0'].isin(df2['col_0'])
common_columns = df1.columns.intersection(df2.columns)
df1_common = df1[common_rows].reset_index(drop=True).set_index('col_0')[common_columns[1:]]
df2_common = df2[df2['col_0'].isin(df1['col_0'])].reset_index(drop=True).set_index('col_0')[common_columns[1:]]
df2_common = df2_common.loc[df1_common.index]

plt.figure(figsize=(12, 8))
plt.subplot(2, 2, 1)
pmi_differences = df1_common - df2_common
sns.heatmap(pmi_differences, cmap='coolwarm', center=0)
plt.title('Heatmap of PMI Differences')
plt.tight_layout()
plt.savefig('pairwise_alignment/pmi_difference.png', bbox_inches='tight')

avg_pmi_original = df1_common.mean().mean()
avg_pmi_iwsa = df2_common.mean().mean()
print(f"Average PMI Score (Original): {avg_pmi_original:.4f}")
print(f"Average PMI Score (IWSA): {avg_pmi_iwsa:.4f}")
avg_pmi_difference = (df2_common - df1_common).mean().mean()
print(f"Average Difference in PMI Scores (IWSA - Original): {avg_pmi_difference:.4f}")

# Average PMI Score (Original): -2.3855
# Average PMI Score (IWSA): 1.2793
# Average Difference in PMI Scores (IWSA - Original): 3.6648

# IWSA doesn't work well on the small datasets - maybe that's why.

