FROM snakemake/snakemake:v7.32.4
WORKDIR ../../
ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN git clone https://github.com/khalidtab/FALAPhyl && mv ./FALAPhyl/workflow/build.snakefile ./build.snakefile && cp -r FALAPhyl/workflow ./workflow && cp FALAPhyl/data/input.yaml ./input.yaml && cp FALAPhyl/workflow/snakefile ./snakefile && cp FALAPhyl/LICENSE ./LICENSE && cp FALAPhyl/README.md ./README.md && rm -r FALAPhyl && snakemake snakefile.final --use-conda --cores all --conda-create-envs-only --snakefile  ./build.snakefile && rm ./build.snakefile && conda clean -a -y && rm ./snakefile
