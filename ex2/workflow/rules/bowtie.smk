import pandas as pd

configfile: "../config/config.yaml"


samples = pd.read_csv("../config/samples.tsv", sep="\t")
# Load samples table
SAMPLES = samples["sample"].tolist()
SAMPLE_READS = dict(zip(samples["sample"], zip(samples["fq1"], samples["fq2"])))

# List of expected index files
BT2_INDEX_FILES = expand("{ref}.{ext}", ref=config["reference"], ext=[
    "1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"
])

rule build_bowtie2_index:
    input:
        fasta = config["reference"]
    output:
        BT2_INDEX_FILES
    threads: 1
    conda:
        "../envs/bowtie.yaml"
    shell:
        """
        bowtie2-build {input.fasta} {input.fasta}
        """

rule bowtie2_map:
    input:
        reads = lambda wildcards: SAMPLE_READS[wildcards.sample],
        index = BT2_INDEX_FILES
    output:
        "../resources/tiny/{sample}.sam"
    threads: 4
    conda:
        "../envs/bowtie.yaml"
    shell:
        """
        bowtie2 {config[bowtie2_params][sensitivity]} \
                {config[bowtie2_params][seed_mismatch]} \
                -x {config[reference]} \
                -1 {input.reads[0]} -2 {input.reads[1]} \
                -S {output} \
                -p {threads}
        """

