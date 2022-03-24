rule all:
   input: "{sample}.final"
   shell: "touch {input}"

rule biom:
   conda:
      "./workflow/envs/biom.yaml"
   message: "Creating biom"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}biom"))
   shell:
      "touch {output}"

rule csvkit:
   conda:
      "./workflow/envs/csvkit.yaml"
   message: "Creating csvkit"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}csvkit"))
   shell:
      "touch {output}"

rule fastspar:
   conda:
      "./workflow/envs/fastspar.yaml"
   message: "Creating fastspar"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}fastspar"))
   shell:
      "touch {output}"


rule ggpubr:
   conda:
      "./workflow/envs/ggpubr.yaml"
   message: "Creating ggpubr"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}ggpubr"))
   shell:
      "touch {output}"


rule ggrepel:
   conda:
      "./workflow/envs/ggrepel.yaml"
   message: "Creating ggrepel"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}ggrepel"))
   shell:
      "touch {output}"

rule philr:
   conda:
      "./workflow/envs/philr.yaml"
   message: "Creating philr"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}philr"))
   shell:
      "touch {output}"

rule phyloseq:
   conda:
      "./workflow/envs/phyloseq_vegan_tidyverse.yaml"
   message: "Creating phyloseq"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}phyloseq"))
   shell:
      "touch {output}"

<<<<<<< Updated upstream
=======
rule betapart:
   conda:
      "./workflow/envs/betapart.yaml"
   message: "Creating betapart"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}betapart"))
   shell:
      "touch {output}"

>>>>>>> Stashed changes

      
rule results:
   version: "1.0"
   input:
      rules.betapart.output,
      rules.biom.output,
      rules.csvkit.output,
      rules.fastspar.output,
      rules.ggpubr.output,
      rules.ggrepel.output,
      rules.philr.output,
      rules.phyloseq.output
   message: "Cleaning up…"
   output:
      "{sample}.final"
   shell:
      "echo Done building environments"