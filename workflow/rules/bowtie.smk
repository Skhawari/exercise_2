rule bowtie_idx:
    input:
        ref = config["ref"]
    output:
        "resources/reference.fa.1.bt2",
        "resources/reference.fa.2.bt2",
        "resources/reference.fa.3.bt2",
        "resources/reference.fa.4.bt2",
        "resources/reference.fa.rev.1.bt2",
        "resources/reference.fa.rev.2.bt2"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build {input.ref} resources/reference.fa"


rule bowtie_map:
    input:
        r1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards: samples.at[wildcards.sample, 'fq2'],
        index = [
            "resources/reference.fa.1.bt2",
            "resources/reference.fa.2.bt2",
            "resources/reference.fa.3.bt2",
            "resources/reference.fa.4.bt2",
            "resources/reference.fa.rev.1.bt2",
            "resources/reference.fa.rev.2.bt2"
        ]

    output:
        "results/sam/{sample}.sam"
    threads: 4
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 {config[bowtie2_params]} \
        -x resources/reference.fa \
        -1 {input.r1} -2 {input.r2} \
        -S {output} -p {threads}
        """