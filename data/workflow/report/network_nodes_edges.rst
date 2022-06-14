SparCC network of {{ snakemake.wildcards.groupcategory }} category, part of the group {{ snakemake.wildcards.group }}. Other parameters are as follows:

- Core level: {{ snakemake.config["threshold"] }}
- P-value cut off to consider: p> {{ snakemake.config["sparcc_pvalue"] }}
- Correlation cut off level: {{ snakemake.config["sparcc_corr"] }}

