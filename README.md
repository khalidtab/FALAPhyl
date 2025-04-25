![image](FALAPhyl.png)

# FALAPhyl: *F*orays into *A*utomating *L*aborious *A*nalyses of *Phyl*ogeny

This is a pipeline that fully automates some bioinformatic analysis using well-recognised packages such as PhyloSeq, ggplot and others. It is built on top of snakemake.

The rationale for this is to perform exploratory analysis of your feature data, and produce publication quality plots for beta diversity (Principal Coordinates Analysis - PCoA, Non-metric Multidimensional Scaling - NMDS), alpha diversity, the network connectivity Zi-Pi plot (Zi within-module connectivity, Pi among-module connectivity).

You need a features table in a biom file format, and a mapping file. If compositional beta diversity is specified (PhILR), a rooted tree is not needed as one is automatically generated using Ward Hierarchical clustering. Note that the PhILR script only uses the tree to identify which nodes to compare against other nodes, and not the distances between them.

## Example

If you are using the docker image, you will need to mount the OS's folder as a volume:

> docker run --rm -it -v ~/Path/To/Your/Folder:/data khalidtab/falaphyl:latest start

This will copy the following files to your path (if they are not there already):

- `FALAPhyl_environments.txt`: Contains the paths to the snakemake environments used by the pipeline. It is helpful if you want to run specific commands on the environment yourself.
- `FALAPhyl_metadata.txt`: template file to use
- `falaphyl.yaml`: configuration file

Before executing the code, ensure the following:

1. `filename.biom`, where filename` is the file name that is exactly the same as the other input files (mapping file)
2. Mapping file in the `filname..txt` format, where `filename` is the file name that is exactly the same as the other input files (`biom`). File is in the format typically used for Qiime2. That is, first column "#SampleID", contains sample IDs that match those in the biom file, then the other columns are the additional variables you're interested in analyzing. It might be best if you use the `FALAPhyl_metadata.txt` as your template, then modify it to be the same name as the `biom` file. Save it as `filename.txt`.
3. `falaphyl.yaml`: the file with the input parameters. The `filename` must be listed in there.

To run the pipeline. Run the following command

> snakemake alpha beta breakdown network diff subject_alpha subject_beta subject_diff --cores all --use-conda --keep-going --retries 5 --scheduler greedy --rerun-incomplete

## Command explanations

1. `alpha`: alpha diversity, and violin plots of the graphs
2. `beta`: dissimilarity matrices, ADONIS, ANOSIM, Beta dispersion, PCoA and NMDS
3. `breakdown`: Jaccard and Bray-Curtis broken down to their corresponding components (per the R package Betapart). ANOSIM on the two components. Probability density function of the components through permutation.
4. `network`: nodes and edges based on SparCC. Zi-Pi graphs to determine which nodes are most interesting.
5. `diff`: differential abundance using multiple methods, using the R package DAtest. Will do pairwise comparisons of selected variables. Will output differential abundance, power analysis, and the following graphs to help identify the best method: AUC, FDR, Power, and Scores.
6. `subject_alpha`: For paired/repeated measures. Alpha diversity differences but restricted to differences within a subject.
7. `subject_beta`: For paired/repeated measures. Same as the breakdown option above, but the comparison is restricted to differences within the subject.
8. `subject_diff`: For paired/repeated measures. Same as the diff option above, but the comparisons is restricted to differences within the subject.
