import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--regions", type=str, dest="regions")
parser.add_argument("--denovo_counts", type=str, dest="denovo_counts")
parser.add_argument("--output", type=str, dest="output")
args = parser.parse_args()

# Reading files
df_regions = pd.read_table(args.regions, header=None)
df_denovo = pd.read_table(args.denovo_counts, header=None)

results = pd.DataFrame()

# Counting missense mutations in NDD patients in depleted regions
for iter_region in range(0, len(df_regions)):
    transcript = df_denovo.loc[df_denovo[0] == df_regions.iloc[iter_region, 1], ]
    # within the start and stop of an exon
    # Uncomment the line below for finding count per exon for top MAF

    region_mutations = transcript.loc[transcript[1].isin(range(int(df_regions.iloc[iter_region, 3]),
                                                               int(df_regions.iloc[iter_region, 4]+1)))]

    # Uncomment the line below for finding count per sig region from S
    """
    region_mutations = transcript.loc[transcript[1].isin(range(int(df_regions.iloc[iter_region, 2]),
                                                               int(df_regions.iloc[iter_region, 3] + 1)))]
    """
    sum_mutations = float(region_mutations[2].sum())
    # Uncomment the line below for finding count per exon for top MAF

    line = pd.DataFrame([[df_regions.iloc[iter_region, 0], df_regions.iloc[iter_region, 1],
                          df_regions.iloc[iter_region, 2], df_regions.iloc[iter_region, 3],
                          df_regions.iloc[iter_region, 4], df_regions.iloc[iter_region, 5],
                          df_regions.iloc[iter_region, 6], df_regions.iloc[iter_region, 7],
                          df_regions.iloc[iter_region, 8], df_regions.iloc[iter_region, 9],
                          df_regions.iloc[iter_region, 10], df_regions.iloc[iter_region, 11], sum_mutations]])
    """
    # Uncomment the line below for finding count per sig region from S
    line = pd.DataFrame([[df_regions.iloc[iter_region, 0], df_regions.iloc[iter_region, 1],
                          df_regions.iloc[iter_region, 2], df_regions.iloc[iter_region, 3],
                          df_regions.iloc[iter_region, 4], df_regions.iloc[iter_region, 5],
                          df_regions.iloc[iter_region, 6], sum_mutations]])
    """
    results = results.append(line)

# Saving files
results.to_csv(args.output, sep='\t', index=False, header=False)
