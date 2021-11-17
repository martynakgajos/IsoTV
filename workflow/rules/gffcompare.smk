rule GffCompare:
    input:
        reference = GenomeGFF,
        exptgff = rules.PinfishCollapsePolishedPartials.output.pol_gff_col
    output:
        nanopore_gtf = "Results/GffCompare/nanopore.combined.gtf"
    resources:
        memory = 16,
        time = 1
    conda:
        "../envs/gffcompare.yaml"
    shell:
        "gffcompare -r {input.reference} -R -A -K -o " + "Results/GffCompare/nanopore {input.exptgff}"

# Generate corrected transcriptome
rule PrepareCorrectedTranscriptomeFasta:
    input:
        genome = ReferenceFasta,
        gff = rules.PinfishCollapsePolishedPartials.output.pol_gff_col,
    output:
        fasta = "Results/Pinfish/corrected_transcriptome_polished_collapsed.fas",
    resources:
        memory = 8,
        time = 1
    conda:
        "../envs/gffcompare.yaml"
    shell:
        "gffread -g {input.genome} -w {output.fasta} {input.gff}"

rule PrepareCorrectedTranscriptomeFastaNonred:
    input:
        rules.PrepareCorrectedTranscriptomeFasta.output.fasta,
        rules.GffCompare.output.nanopore_gtf
    output:
        fasta = "Results/Pinfish/corrected_transcriptome_polished_collapsed_nonred.fas",
    resources:
        memory = 8,
        time = 1
    conda:
        "../envs/pythonpackages.yaml"
    script:
        "../scripts/nonredundant.py"

rule ReadsToTranscripts:
    input:
         rules.PrepareCorrectedTranscriptomeFastaNonred.output.fasta,
         rules.GffCompare.output.nanopore_gtf
    output:
        "Results/Pinfish/corrected_transcriptome_polished_collapsed_nonredundant_tcons.fas"
    resources:
        memory = 10,
        time = 1,
        tmpdir = 0
    priority: 8
    group: 
        "preprocess_sequences"
    conda:
        "../envs/pythonpackages.yaml"
    script:
        "../scripts/reads_to_transcripts.py"
