rule alpha_div: # Calculates whether the two groups have different centroids, susceptible to the groups having different dispersions. Therefore, interpret along with PERDISP (also known as betadisper)
   version: "1.0"
   conda:
      "../../workflow/envs/phylotoast.yaml"
   input:
      biom="data/biom/{sample}.biom",
      map="data/map/{sample}.txt"
   output:
      sig=report(expand("data/alpha_div/alpha_div_{{sample}}+{x}.txt",x=config["alpha"])),
      svg=report(expand("data/plots/alpha_div_{{sample}}+{x}.svg",x=config["alpha"]))
   params: 
      alpha=expand("{alpha}",alpha=config["alpha"]),
      group=config["group"][0],
      color=config["color"][0]
   shell:
      "echo 'for x in {params.alpha}; do mkdir -p data/plots/alpha_div_{wildcards.sample}+$x && "
      "xvfb-run --auto-servernum diversity.py -i {input.biom} -m {input.map} -c {params.group} -d $x --o data/plots/alpha_div_{wildcards.sample}+$x --image_type svg --color_by {params.color} >> data/alpha_div/alpha_div_{wildcards.sample}+$x.txt && "
      "mv data/plots/alpha_div_{wildcards.sample}+$x/d.svg data/plots/alpha_div_{wildcards.sample}+$x.svg && "
      "rm -rf data/plots/alpha_div_{wildcards.sample}+$x "
      "; done' > tmp/alpha_div_{wildcards.sample}.sh | "
      "chmod +x tmp/alpha_div_{wildcards.sample}.sh | "
      "bash tmp/alpha_div_{wildcards.sample}.sh"
