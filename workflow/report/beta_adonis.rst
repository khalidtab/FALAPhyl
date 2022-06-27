ADONIS (also known as non-parametric multivariate analysis of variance NPMANOVA or PERMANOVA) of the {{ snakemake.wildcards.distance }} dissimilarity matrix. The variable {{ snakemake.wildcards.group }} was used for grouping the samples.

Null hypothesis: there are no differences in the presence/absence or relative magnitude of a set of variables among objects from different groups or treatments.

ADONIS is sensitive to dispersion of the samples that it is possible that two groups can have the same centroid but would still be significantly different because of the beta-dispersion of the samples. Therefore, it should be interpreted with the beta dispersion (permdisp).