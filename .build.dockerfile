FROM snakemake/snakemake:v7.32.4
WORKDIR ../../

RUN git clone https://github.com/khalidtab/FALAPhyl && \
    mv ./FALAPhyl/workflow/build.snakefile ./build.snakefile && \
    cp -r FALAPhyl/workflow ./workflow && cp FALAPhyl/data/falaphyl.yaml ./falaphyl.yaml && cp FALAPhyl/FALAPhyl_metadata.txt ./FALAPhyl_metadata.txt && \
    cp FALAPhyl/workflow/snakefile ./snakefile && cp FALAPhyl/LICENSE ./LICENSE && \
    cp FALAPhyl/README.md ./README.md && \
    cp FALAPhyl/FALAPhyl_metadata.txt ./FALAPhyl_metadata.txt && \
    cp FALAPhyl/FALAPhyl.png ./FALAPhyl.png && \
    rm -r FALAPhyl && \
    snakemake snakefile.final --use-conda --cores all --conda-create-envs-only --snakefile  ./build.snakefile 2>&1 | tee -a FALAPhyl_environments.txt && \
    rm ./build.snakefile && conda clean -a -y && rm ./snakefile

COPY start.sh /usr/local/bin/start
RUN chmod +x /usr/local/bin/start
CMD ["/usr/local/bin/start"]

