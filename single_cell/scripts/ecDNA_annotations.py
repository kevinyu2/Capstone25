import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt


# ecDNA_df = pd.read_csv('../data/ecDNA/ecDNA_preds_cnaObj_HT268B1-Th1H3.all.tsv', sep = '\t')
gene = 'MIR4743'
sample = 'HT268B1-Th1H3'
boolean_ecDNA = True
fig_dir = '../data/ecDNA/figs'

# Read just the column
count_df = pd.read_csv(
    f'../data/raw/cnaObj_{sample}.tsv',
    sep='\t',
    usecols=[gene],
)
count_df.index = count_df.index.str.split('#').str[-1]

# Better visuals
count_dict = count_df[gene].to_dict()
for key in count_dict.keys() :
    if count_dict[key] > 20 :
        count_dict[key] = 20
print(np.mean(list(count_dict.values())))

# Make boolean
if boolean_ecDNA :
    count_dict_2 = {}
    for key in count_dict.keys() :
        if count_dict[key] > 6 :
            count_dict_2[key] = 1
        else :
            count_dict_2[key] = 0
    count_dict = count_dict_2

coords_folder = f'./{sample}_UMAPs/'

for file in glob.glob(f"{coords_folder}/UMAP_coords*") :
    print(file)

    type = file.split('_')[-1].split('.')[0]

    df = pd.read_csv(file)
    df["score"] = df["barcode"].map(count_dict)

    plt.figure(figsize=(6,6))
    plt.scatter(
        df["UMAP1"],
        df["UMAP2"],
        c=df["score"],
        cmap="viridis",
        s=4
    )

    plt.colorbar(label="score")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("UMAP score")
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/{sample}_{gene}_{type}_UMAP.png")

    