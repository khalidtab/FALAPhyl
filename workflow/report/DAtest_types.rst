{{ snakemake.wildcards.sample }} –  differential abundance.

Details of the methods can be found in the official DAtest page (https://github.com/Russel88/DAtest/wiki/Methods). In summary:
- DESeq2 based tests
- ds2: DESeq2: negative binomial. Uses manual geometric means and the results normalized as relative log-expression
- ds2x: DESeq2: negative binomial. The results are normalized as relative log-expression
- EdgeR (ere,ere2,erq,erq2 – Exact test,Exact test2,Quasi likelihood,Quasi likelihood2)
	- erq: EdgeR - Quasi likelihood. Negative binomial. Normalization as Trimmed Mean by M-value (TMM)
	- erq2: EdgeR - Quasi likelihood2. Negatie binomial. Normalization as Relative Log Expression (RLE)
- LIMMA based tests(lia,lic,lli,lli2,vli = LIMMA(ALR,CLR,log,log #2,voom))
	- lli: LIMMA log. Gaussian. Log transformation is done before normalization with Total sum scaling.
	- lli2: LIMMA log #2. Gaussian. Log transformation is done after normalization with Total sum scaling.
- t-test based  tests
	- tta,ttc = t-test (ALR,CLR)
	- ttt,ltt,ltt2,ttr = Welch(t.test,t.test log,t.test log #2,t.test Rank-normalization)
- neb: GLM - Negative binomial. 
- poi: GLM - Poisson
- vli: LIMMA voom. Gaussian. Data normalized with Total sum scaling, and transformed with Voom
- zig: MetagenomeSeq ZIG. Zero-inflated Gaussian. Data normalized with CSS and transformed with Log	
- msf: MetagenomeSeq featuremodel
- per = Permutation test
- adx = ALDEx2
- abc = ANCOM-BC
- bay = baySeq
- neb,poi,qpo,znb,zpo = GLM (Negative binomial,Poisson,Quasi-poisson,ZI Negative Binomial,ZI Poisson)
- wil = Wilcoxon
- From Tweedieverse - CPLM: Compound poisson linear model
