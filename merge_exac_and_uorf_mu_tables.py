#%%
import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

import matplotlib
matplotlib.use('module://backend_interagg')

import matplotlib.pyplot as plt

#%%

df1 = pd.read_table("exac_2015/data/recurrence/fordist_1KG_mutation_rate_table.txt", sep=" ")
df1.loc[:, 'to'] = df1.loc[:, 'to'].apply(lambda s: s[1:2])
df1.set_index(['from', 'to'], inplace=True)
#%%

df1

#%%

df2 = pd.read_table("uORFs/code/SynonymousPropSingleton_byTrimer_NEW.txt")
df2.rename(columns={'Context': 'from', 'Alt': 'to'}, inplace=True)
df2.set_index(['from', 'to'], inplace=True)

#%%
df2


#%%

RC = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def reverse_complement(s):
    return "".join([RC[b] for b in s][::-1])

assert reverse_complement('A') == 'T'
assert reverse_complement('ACGT') == 'ACGT'
assert reverse_complement('GGA') == 'TCC'

#%%

len(df1)

len(df2)

#%%

df1_dict = {(r['from'], r['to']): r.mu_snp for _, r in df1.reset_index().iterrows()}

df2_dict = {(r['from'], r['to']): (r.mu_snp, r.Count_variants) for _, r in df2.reset_index().iterrows() if r.Methyl_bin == 0}
df2_dict_methyl_bin_1 = {(r['from'], r['to']): (r.mu_snp, r.Count_variants) for _, r in df2.reset_index().iterrows() if r.Methyl_bin == 1}
df2_dict_methyl_bin_2 = {(r['from'], r['to']): (r.mu_snp, r.Count_variants) for _, r in df2.reset_index().iterrows() if r.Methyl_bin == 2}

for key in df2_dict_methyl_bin_2.keys():
    total_variant_count = df2_dict[key][1] + df2_dict_methyl_bin_1[key][1] + df2_dict_methyl_bin_2[key][1]
    df2_dict[key] = (
        (df2_dict[key][0]*df2_dict[key][1] + df2_dict_methyl_bin_1[key][0]*df2_dict_methyl_bin_1[key][1] + df2_dict_methyl_bin_2[key][0]*df2_dict_methyl_bin_2[key][1])/total_variant_count,
        total_variant_count,
    )

df1_with_reverse_complement_scores = {}
for _, r in df1.reset_index().iterrows():
    rc_from = reverse_complement(r['from'])
    rc_to = reverse_complement(r['to'])
    if (rc_from, rc_to) in df1_with_reverse_complement_scores:
        continue

    is_cpg = (r['from'][1:] == 'CG' and r['to'] == 'T') or (r['from'][:2] == 'CG' and r['to'] == 'A')
    is_transition = (
        r['from'][1] == 'A' and r['to'] == 'G') or (
        r['from'][1] == 'G' and r['to'] == 'A') or (
        r['from'][1] == 'C' and r['to'] == 'T') or (
        r['from'][1] == 'T' and r['to'] == 'C')
    #is_tranversion = not is_transition

    exac_mu_snp_reverse_complement = df1_dict.get((rc_from, rc_to))
    assert exac_mu_snp_reverse_complement

    uorf_2_tuple = df2_dict_methyl_bin_2.get((r['from'], r['to']), (None, None))
    if not uorf_2_tuple[0]:
        uorf_2_tuple = df2_dict.get((r['from'], r['to']), (None, None))
    if not uorf_2_tuple[0]:
        uorf_2_tuple = df2_dict.get((rc_from, rc_to), (None, None))

    uorf_mu_snp = uorf_2_tuple[0]
    uorf_variant_count = uorf_2_tuple[1]

    assert uorf_mu_snp

    df1_with_reverse_complement_scores[(r['from'], r['to'])] = {
        'from': r['from'],
        'to': r['to'],
        'snp_type': "cpg" if is_cpg else ("transition" if is_transition else "transversion"),
        'exac_mu_snp': r.mu_snp,
        'exac_mu_snp_reverse_complement': exac_mu_snp_reverse_complement,
        'uorf_mu_snp': uorf_mu_snp,
        'uorf_variant_count': uorf_variant_count,
    }

# cpg
# (mutations$ref == 'C' & mutations$alt == 'T' & substr(mutations$context, 2, 3) == 'CG') |
# (mutations$ref == 'G' & mutations$alt == 'A' & substr(mutations$context, 1, 2) == 'CG')
# transition
# (mutations$ref == 'A' & mutations$alt == 'G') |
# (mutations$ref == 'G' & mutations$alt == 'A') |
# (mutations$ref == 'C' & mutations$alt == 'T') |
# (mutations$ref == 'T' & mutations$alt == 'C')

#%%
df3 = pd.DataFrame(df1_with_reverse_complement_scores.values())

df3.to_csv("exac_mu_snp_with_reverse_complements.tsv", index=False, sep="\t", header=True)


#%%

COLOR_MAP = {
    'cpg': '#2E9FFD',
    'transition': '#458A03',
    'transversion': '#EA4545',
}

#%%
p = df3.plot.scatter('exac_mu_snp', 'exac_mu_snp_reverse_complement', color=df3["snp_type"].map(COLOR_MAP))
p.set_title("ExAC 2015 - mu_snp vs. mu_snp of reverse-complement")
plt.show()

#%%
p = df3.plot.scatter('exac_mu_snp', 'uorf_mu_snp', color=df3["snp_type"].map(COLOR_MAP))
p.set_title("exac_mu_snp vs. uorf_mu_snp")
plt.show()



#%%
#%%
