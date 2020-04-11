rule alpha_div: # Calculates whether the two groups have different centroids, susceptible to the groups having different dispersions. Therefore, interpret along with PERDISP (also known as betadisper)
   version: "1.0"
   conda:
      "../../workflow/envs/phylotoast.yaml"
   input:
      biom="data/biom/{sample}.biom",
      map="data/map/{sample}.txt"
   output:
      sig=report("data/alpha_div/alpha_div_{sample}.txt"),
      svg=report(expand("data/plots/alpha_div_{{sample}}+{x}.svg",x=config["alpha"]))
   params: 
      alpha=expand("{alpha}",alpha=config["alpha"]),
      group=config["group"][0],
      color=config["color"][0]
   shell:
      "mkdir -p data/alpha_div data/scripts &&"
      "echo 'for x in {params.alpha}; do diversity.py -i {input.biom} -m {input.map} -c {params.group} -d $x --o data/plots/alpha_div_{wildcards.sample} --image_type svg --color_by {params.color} >> data/alpha_div/alpha_div_{wildcards.sample}.txt &&"
      "mv data/plots/alpha_div_{wildcards.sample}/$x.svg data/plots/alpha_div_{wildcards.sample}+$x.svg &&"
      "rm -rf data/plots/alpha_div_{wildcards.sample} "
      "; done' > data/scripts/alpha_div_{wildcards.sample}.sh |"
      "chmod +x data/scripts/alpha_div_{wildcards.sample}.sh |"
      "bash data/scripts/alpha_div_{wildcards.sample}.sh"
