{{ snakemake.wildcards.sample }} –  differential abundance.

Details of the methods can be found in the official DAtest page (https://github.com/Russel88/DAtest/wiki/Methods). In summary:
- ds2: DESeq2: negative binomial. Uses manual geomeetric means and the results normalized as relative log-expression
- ds2x: DESeq2: negative binomial. The results are normalized as relative log-expression
- erq: EdgeR - Quasi likelihood. Negative binomial. Normalization as Trimmed Mean by M-value (TMM)
- erq2: EdgeR - Quasi likelihood2. Negatie binomial. Normalization as Relative Log Expression (RLE)
- lli: LIMMA log. Gaussian. Log transformation is done before normalization with Total sum scaling.
- lli2: LIMMA log #2. Gaussian. Log transformation is done after normalization with Total sum scaling.
- neb: GLM - Negative binomial. 
- poi: GLM - Poisson
- vli: LIMMA voom. Gaussian. Data normalized with Total sum scaling, and transformed with Voom
- zig: MetagenomeSeq ZIG. Zero-inflated Gaussian. Data normalized with CSS and transformed with Log	