# Build minimap2 index to reference genome
rule Minimap2Index:
    input:
        genome = ReferenceFasta
    output:
        index = "Results/Minimap2/%s.mmi" %(ReferenceFasta.split("/")[-1])
    params:
        opts = config["minimap2_index_opts"],
    #log:
    #    "logs/downsample/{sample}/{sample}_bbnorm.log"
    threads:
        config["threads"]
    resources:
        memory = 8,
        time = 1
    conda:
        "../envs/minimap2.yaml"
    shell:
      "minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}"

# Build known splice junction bed file
rule SpliceJunctionIndex:
    input:
        GenomeGFF
    output:
        junc_bed = "ReferenceData/junctions.bed"
    conda:
        "../envs/minimap2.yaml"
    shell:
        "paftools.js gff2bed {input} > {output.junc_bed}"

# Filter basecalled raw ONT reads

def FilterReadsInput(wildcards):
    if(os.path.exists("RawData/Fastq/%s.fastq" %(SAMPLES[wildcards.sample]))):
        return "RawData/Fastq/%s.fastq" %(SAMPLES[wildcards.sample])
    else:
        assert(os.path.exists("RawData/Fastq/%s.fastq.gz" %(SAMPLES[wildcards.sample])))
        return "RawData/Fastq/%s.fastq.gz" %(SAMPLES[wildcards.sample])

rule FilterReads:
    input:
        FilterReadsInput
    output:
        "FilteredData/{sample}.fastq",
    params:
        config["min_mean_q"]
    resources:
        memory = 8,
        time = 1
    conda:
        "../envs/filtlong.yaml"
    shell:
        "filtlong --min_mean_q {params} {input} > {output}"
