rule sam_to_bam:
    input:
        lambda wildcards: f"data/tiny/{wildcards.sample}.sam"
    output:
        "results/bam/{sample}.bam"
    threads: 1
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -Sb {input} > {output}"

rule sort_bam:
    input:
        "results/bam/{sample}.bam"
    output:
        "results/bam_sorted/{sample}.sorted.bam"
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule index_bam:
    input:
        "results/bam_sorted/{sample}.sorted.bam"
    output:
        "results/bam_sorted/{sample}.sorted.bam.bai"
    threads: 1
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule idxstats:
    input:
        bam="results/bam_sorted/{sample}.sorted.bam",
        bai="results/bam_sorted/{sample}.sorted.bam.bai"
    output:
        "results/stats/{sample}.idxstats.txt"
    threads: 1
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools idxstats {input.bam} > {output}"

rule aggregate_idxstats:
    input:
        expand("results/stats/{sample}.idxstats.txt", sample=SAMPLES)
    output:
        "results/stats/aggregated_idxstats.tsv"
    run:
        import pandas as pd

        aggregated = None

        for infile in input:
            sample = os.path.basename(infile).replace(".idxstats.txt", "")
            df = pd.read_csv(infile, sep="\t", header=None,
                             names=["seqname", "seqlen", "mapped", "unmapped"])
            df = df[["seqname", "seqlen", "mapped"]]
            df = df.rename(columns={"mapped": sample})

            if aggregated is None:
                aggregated = df
            else:
                aggregated = pd.merge(aggregated, df[["seqname", sample]], on="seqname")

        # Move "seqlen" column to front if not already
        cols = ["seqname", "seqlen"] + [col for col in aggregated.columns if col not in ["seqname", "seqlen"]]
        aggregated = aggregated[cols]

        aggregated.to_csv(output[0], sep="\t", index=False)
