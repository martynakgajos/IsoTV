
# Quantify transcripts
rule CountTranscripts:
    input:
        bam = rules.ExtractPrimaryMapping.output
    output:
        quant = "Results/Quantification/{sample}.counts"
    resources:
        memory = 4,
        time = 1
    shell:"""
        echo -e "counts\ttranscript" > {output.quant};
        samtools view {input.bam} | cut -f3 | uniq -c | grep -v "*" | sed -e 's/^[ \t]*//' | sed 's/ /\t/' >> {output.quant}
        """

# Concatenate transcript counts
rule ConcatenateCounts:
    input:
        annotation = GenomeGFF,
        transcriptome = rules.GffCompare.output.nanopore_gtf,
        counts = expand("Results/Quantification/{sample}.counts", sample=SAMPLES),
    resources:
        memory = 4,
        time = 1
    output:
        counts = "Results/Quantification/counts.txt"
    script:
        "../scripts/concatenate.py"

# Normalize isoform levels
rule CalculateDeSeq2Norm:
    input:
        counts = rules.ConcatenateCounts.output
    output:
        norm_counts = "Results/Quantification/counts_deseq2norm.txt"
    resources:
        memory = 4,
        time = 1
    script:
        "../scripts/deseq2norm.R"
