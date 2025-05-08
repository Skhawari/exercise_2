import pandas as pd
import os

input_files = snakemake.input
output_file = snakemake.output[0]

dfs = []
for file in input_files:
    sample = os.path.basename(file).split(".")[0]
    df = pd.read_csv(file, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["ref", "length", sample])
    df = df.set_index("ref")
    dfs.append(df[[sample]])

summary = pd.concat(dfs, axis=1)
summary.to_csv(output_file, sep="\t")
