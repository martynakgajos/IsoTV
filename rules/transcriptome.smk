
# Build minimap2 index for transcriptome
rule Transcriptome2Index:
    input:
        genome = rules.PrepareCorrectedTranscriptomeFastaNonred.output.fasta,
    output:
        index = "Results/Minimap2/Transcriptome.mmi"
    threads: config["threads"]
    resources:
        memory = 16,
        time = 2
    shell:
      "minimap2 -t {threads} -I 1000G -d {output.index} {input.genome}"

# Map reads using minimap2
rule Minimap2Genome:
    input:
       index = rules.Minimap2Index.output.index,
       fastq = rules.FilterReads.output,
       use_junc = rules.SpliceJunctionIndex.output.junc_bed if config["minimap2_opts_junction"] else ""
    output:
       bam = "IGV/{sample}.genome.bam"
    resources:
        memory = 32,
        time = 4
    params:
        opts = config["minimap2_opts"]
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} --junc-bed {input.use_junc} {input.index} {input.fastq}\
    | samtools view -F 260 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """

# Map to transcriptome
rule Map2Transcriptome:
    input:
       index = rules.Transcriptome2Index.output.index,
       fastq = rules.FilterReads.output
    output:
       bam = "Results/Quantification/{sample}.bam",
       sbam = "Results/Quantification/{sample}.sorted.bam"
    threads: config["threads"]
    resources:
        memory = 32,
        time = 4
    params:
        msec = config["maximum_secondary"],
        psec = config["secondary_score_ratio"],
    priority: 10
    shell:"""
        minimap2 -ax map-ont -t {threads} -p {params.psec} -N {params.msec} {input.index} {input.fastq} | samtools view -Sb > {output.bam};
        samtools sort -@ {threads} {output.bam} -o {output.sbam};
        samtools index {output.sbam};
        """

rule ExtractPrimaryMapping:
    input:
        transcriptome = rules.Map2Transcriptome.output.sbam,
        genome = rules.Minimap2Genome.output.bam
    resources:
        memory = 4,
        time = 1
    output:
        "IGV/{sample}.transcriptome.bam"
    shell:"""
        samtools view -b -F 260 {input.transcriptome} > {output};
        samtools index {output}
        """
