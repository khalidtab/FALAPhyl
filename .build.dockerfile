FROM snakemake/snakemake:latest
WORKDIR ../../
RUN git clone https://github.com/khalidtab/FALAPhyl && mv ./FALAPhyl/workflow/build.snakefile ./build.snakefile &&cp -r FALAPhyl/workflow ./workflow &&cp FALAPhyl/data/input.yaml ./input.yaml && cp FALAPhyl/data/snakefile ./snakefile && cp FALAPhyl/LICENSE ./LICENSE && cp FALAPhyl/README.md ./README.md && rm -r FALAPhyl && snakemake snakefile.final --use-conda --cores all --conda-create-envs-only --snakefile  ./build.snakefile && rm ./build.snakefile && conda clean -a -y
