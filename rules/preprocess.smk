# Build minimap2 index to reference genome
rule Minimap2Index:
    input:
        genome = ReferenceFasta
    output:
        index = "Results/Minimap2/%s.mmi" %(ReferenceFasta.split("/")[-1])
    #log:
    #    "logs/downsample/{sample}/{sample}_bbnorm.log"
    params:
        opts = config["minimap2_index_opts"],
    threads:
        config["threads"]
    resources:
        memory = 16,
        time = 1
    conda:
        "../envs/env.yaml"
    shell:
      "minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}"

# Build known splice junction bed file
rule SpliceJunctionIndex:
    input:
        GenomeGFF
    output:
        junc_bed = "ReferenceData/junctions.bed"
    shell:
        "paftools.js gff2bed {input} > {output.junc_bed}"

# Filter basecalled raw ONT reads
rule FilterReads:
    input:
        lambda wildcards: "RawData/Fastq/%s.fastq" %(SAMPLES[wildcards.sample])
    output:
        "FilteredData/{sample}.fastq",
    resources:
        memory = 32,
        time = 1
    params:
        config["min_mean_q"]
    shell:
        "filtlong --min_mean_q {params} {input} > {output}"
