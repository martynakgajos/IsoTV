# Pychopper
rule Pychopper:
    input:
        rules.FilterReads.output
    output:
        pdf = "Results/Pychopper/{sample}.pychopper_report.pdf",
        fastq = "Results/Pychopper/{sample}.pychop.fastq",
        stats = "Results/Pychopper/{sample}.pychop.stats",
        scores = "Results/Pychopper/{sample}.pychop.scores",
        unclass = "Results/Pychopper/{sample}.unclassified.fastq"
    resources:
        memory=32,
        time = 4
    params:
        opts = config["porechop_heu_stringency"]
    run:
        shell("cdna_classifier.py -b ReferenceData/cdna_barcodes.fas -r {output.pdf} -S {output.stats} -A {output.scores} -u {output.unclass} {input} {output.fastq}")

# Map reads using minimap2
rule Minimap2Pinfish:
    input:
       index = rules.Minimap2Index.output.index,
       fastq = expand("Results/Pychopper/{sample}.pychop.fastq", sample=SAMPLES) if (config["pychopper"]==True) else expand("FilteredData/{sample}.fastq", sample=SAMPLES),
       use_junc = rules.SpliceJunctionIndex.output.junc_bed if config["minimap2_opts_junction"] else ""
    output:
       bam = "Results/Minimap2/merged.mapping.bam"
    resources:
        memory = 32,
        time = 4
    params:
        opts = config["minimap2_opts"],
        min_mq = config["minimum_mapping_quality"],
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} --junc-bed {input.use_junc} {input.index} {input.fastq}\
    | samtools view -q {params.min_mq} -F 2304 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """
