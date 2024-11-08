import pandas as pd

# Load the AMR and resFinder tables
amr_df = pd.read_csv('amrp_results_report-core_n287.tsv', index_col=0, sep="\t")
resfinder_df = pd.read_csv('res-and-pointFinder_results_report.tsv', index_col=0, sep="\t")

# Merge the tables on gene names (index)
merged_abundance = amr_df.add(resfinder_df, fill_value=0)

# Create a metadata DataFrame
metadata = pd.DataFrame(index=merged_abundance.index, columns=merged_abundance.columns)

for gene in merged_abundance.index:
    for sample in merged_abundance.columns:
        tags = []
        if gene in amr_df.index and sample in amr_df.columns and amr_df.at[gene, sample] > 0:
            tags.append('AMR')
        if gene in resfinder_df.index and sample in resfinder_df.columns and resfinder_df.at[gene, sample] > 0:
            tags.append('resFinder')
        
        metadata.at[gene, sample] = '|'.join(tags) if tags else 'none'

# Save the merged tables
merged_abundance.to_csv('amrc_resFinder_merged_abundance.tsv', sep="\t")
metadata.to_csv('amrc_resFinder_metadata.tsv', sep="\t")
