rule sam_to_bam:
  input:
    "results/sam/{sample}.sam"
  output:
    "results/bam/{sample}.bam"
  threads: 4
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
    "samtools sort -o {output} {input}"

rule index_bam:
  input:
    "results/bam_sorted/{sample}.sorted.bam"
  output:
    "results/bam_sorted/{sample}.sorted.bam.bai"
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
  conda:
    "../envs/samtools.yaml"
  shell:
    "samtools idxstats {input.bam} > {output}"
