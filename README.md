
# 🧆 FALAPhyl: *F*orays into *A*utomating *L*aborious *A*nalysis of *Phyl*ogeny 
This is a pipeline that fully automates some bioinformatic analysis using well-recognised packages such as Qiime2, PhyloSeq, ggplot and others. It is built on top of snakemake.

The rationale for this is to perform exploratory analysis of your feature data, and produce publication quality plots for beta diversity (Principal Coordinates Analysis - PCoA, Non-metric Multidimensional Scaling - NMDS), alpha diversity, the network connectivity Zi-Pi plot (Zi within-module connectivity, Pi among-module connectivity).

You need a features table in a biom file format, and a mapping file. If compositional beta diversity is specified (PhILR), a rooted tree is optional. If no tree is provided, then one is automatically generated using Ward Hierarchical clustering.

## How do you run it?
If you are using the docker image, you will need to mount the OS's folder as a volume:

> docker run -ti -v ~/MyLocalFolder/:/data/ khalidtab/falaphyl:latest bash

The volume mounted above should include the following files and folders:

> /biom/ 

This is the folder with the biom files to be processed. Files must be in [filename].biom extension, where [filename] is the file name that is exactly the same as the other input files (biom, mapping file, and optional tree file)
> /map/ 
 
Optional folder with the tree file. Files must be in [filename].tre extension

> input.yaml

The file with the input parameters. This can be downloaded from this github repo, and it is also packaged in the docker image. If you are using the file in the docker image, it must be copied to the volume, and modified before usage. Any parameter that is not wanted should be commented with a #

To run the pipeline. Run the following command
> snakemake [filename].final --cores all --use-conda

[filename] must be the exact same name as the biom and the mapping files, with the .final at the end
