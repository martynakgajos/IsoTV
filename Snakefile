# Siddharth Annaldasula
# required packages: Bio (SeqIO), os, pandas, numpy,  subprocess, string, argparse, seaborn, matplotlib

import os
import pandas as pd
import shutil
from matplotlib.backends.backend_pdf import PdfPages

##### Logistics
workdir: config["outdir"]

ReferenceFasta = config["genome_fasta"]
GenomeGFF  = config["genome_annot"]
SAMPLES = config["samples"]


##### All
rule all:
  input:
    #"Results/Output/" + config["output_plots"]
    "Results/Quantification/counts_deseq2norm.txt"


# Package versions
rule dump_versions:
    output:
        ver = "versions.txt"
    conda:
        "environment.yaml"
    shell:
        "conda list > {output.ver}"

# include: "basecalling.smk"

include: "rules/preprocess.smk"

include: "rules/pychopper.smk"

include: "rules/pinfish.smk"

include: "rules/gffcompare.smk"

include: "rules/transcriptome.smk"

include: "rules/quantification.smk"

###### Isoform Transcript Processing

# Use correct file paths
CountsData =  rules.CalculateDeSeq2Norm.output if config["preprocess"] else config["counts_data"]
Transcripts =  rules.ReadsToTranscripts.output if config["preprocess"] else config["polished_reads"]
NanoporeGTF = rules.GffCompare.output.nanopore_gtf if config["preprocess"] else config["nanopore_gtf"]

# Filter genes that are not found in the analysis, or if their analysis has already been done
checkpoint gene_filter:
    input:
        NanoporeGTF if config["annotation"] else Transcripts
    output:
        directory("Results/Genes_temp")
    params:
        genes_final_dir = "Results/Genes"
    script:
        "scripts/filter_genes.py"

# Extract transcripts for each gene
rule gene_to_transcripts:
    input:
        NanoporeGTF,
        Transcripts,
        "Results/Genes_temp/{gene}"
    output:
        "Results/Genes_temp/{gene}/{gene}_seq.fa"
    params:
        gene = "{gene}"
    resources:
        memory = 10,
        time = 1,
        tmpdir = 0
    priority: 8
    threads: 1
    group: "preprocess_sequences"
    script:
        "scripts/gene_to_transcripts.py"

# Transcript isoform translation
checkpoint translation_to_protein:
    input:
        seq = rules.gene_to_transcripts.output
    output:
        transcripts = directory("Results/Genes_temp/{gene}/transcripts"),
        stats = "Results/Genes_temp/{gene}/{gene}_transcripts_stats.txt"
    params:
        gene = "{gene}"
    resources:
        memory = 10,
        time = 1,
        tmpdir = 0
    priority: 10
    threads: 1
    #group: "preprocess_sequences"
    script:
        "scripts/translation_to_protein.py"

###### Translated Isoform Feature Analysis

# Disorder region analysis
rule iupred2a_analysis:
    input:
        protein = "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa",
        iupred2a = config["iupred2a_path"]
    output:
        seq = "Results/Genes_temp/{gene}/transcripts/{transcript}_sequence.txt",
        idr = "Results/Genes_temp/{gene}/transcripts/{transcript}_func_iupred2a.txt"
    resources:
        memory = 10,
        time = 1,
        tmpdir = 0
    priority: 8
    threads: 1
    group: "sequence_analysis"
    shell:
        "tail +2 {input.protein} >> {output.seq}; {input.iupred2a} -a {output.seq} long >> {output.idr}"

# Domain analysis
rule interpro_scan:
    input:
        translation = "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa"
    output:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_func_domains.gff3"
    params:
        db = "Pfam",
        pfam = config["interproScan_path"],
    resources:
        memory = 16,
        time = 1,
        tmpdir = 0
    priority: 8
    threads: 1
    group: "sequence_analysis"
    shell:
        "sh  {params.pfam} -i {input.translation} -o {output} -f GFF3 -appl {params.db} -dra"

# Secondary structure analysis
rule brewery_analysis:
    input:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa"
    output:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa.ss3"
    params:
        brewery = config["brewery_path"],
        lock1 = config["lock1"],
        lock2 = config["lock2"],
    resources:
        memory = 8,
        time = 1,
        tmpdir = 16
        #memory = 72,
        #time = 1,
        #tmpdir = 0
    priority: 10
    threads: 16
    #group: "sequence_analysis"
    shell:
        "python3 {params.brewery} -i {input} --cpu 8 --noTA --noSA --noCD --fast"
#        'set -ve; export ORIGDIR="$(/bin/pwd)"; flock {params.lock1} bash -ve -c "scp {input} $MXQ_JOB_TMPDIR/"; cd $MXQ_JOB_TMPDIR/; python3 {params.brewery} -i *.fa --cpu 16 --noTA --noSA --noCD --fast; flock {params.lock2} bash -ve -c "cp *.fa.ss3 $ORIGDIR/{output}"'
#        "echo {params.brewery} > {output}"

# PTM analysis
rule functional_site_analysis:
    input:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa"
    output:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_func_sites.txt"
    params:
        ps_scan = config["prositeScan_path"],
        prosite_dat = config["prositeDat_path"],
        pf_scan = config["pfScan_path"],
    resources:
        memory = 10,
        time = 1,
        tmpdir = 0
    priority: 8
    threads: 1
    group: "sequence_analysis"
    shell:
        "perl {params.ps_scan} -d {params.prosite_dat} --pfscan {params.pf_scan} {input} > {output}"

# Combine functional feature analysis on individual isoforms
def functional_files():
    files = []
    if (config["aa"]): files.append("Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa")
    if (config["iupred2a"]): files.append(rules.iupred2a_analysis.output.idr)
    if (config["pfam"]): files.append(rules.interpro_scan.output)
    if (config["brewery"]): files.append(rules.brewery_analysis.output)
    if (config["pfScan"]): files.append(rules.functional_site_analysis.output)
    return files
rule individual_transcript_analysis:
    input:
        functional_files()
    output:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_analysis.txt"
    params:
        gene = "{gene}"
    resources:
        memory = 10,
        time = 1,
        tmpdir = 0
    priority: 5
    threads: 1
    script:
        "scripts/individual_transcript_analysis.py"

def aggregate_transcript_analysis_func(wildcards):
    checkpoint_output_gene = checkpoints.gene_filter.get(**wildcards).output[0]
    checkpoint_output = checkpoints.translation_to_protein.get(**wildcards).output[0]
    return expand("Results/Genes_temp/{gene}/transcripts/{transcript}_analysis.txt",
        gene = wildcards.gene,
        transcript=glob_wildcards(os.path.join(checkpoint_output,"{transcript}_protein.fa")).transcript)

rule aggregate_transcript_analysis:
    input:
        func = aggregate_transcript_analysis_func,
        stats = "Results/Genes_temp/{gene}/{gene}_transcripts_stats.txt",
        gene_temp = directory("Results/Genes_temp/{gene}")
    output:
        analysis = "Results/Genes/{gene}/{gene}_transcripts_filtered_analysis.txt",
        stats = "Results/Genes/{gene}/{gene}_transcripts_stats.txt"
    params:
        gene = "{gene}"
    resources:
        memory = 4,
        time = 1,
        tmpdir = 0
    priority: 5
    threads: 1
    shell:
        "cat {input.func} > {output.analysis}; cat {input.stats} > {output.stats}; rm -rf {input.gene_temp};"

rule protein_coding_potential_analysis:
    input:
        rules.aggregate_transcript_analysis.output.analysis,
        rules.aggregate_transcript_analysis.output.stats,
        CountsData,
        NanoporeGTF,
    output:
        "Results/Genes/{gene}/{gene}_functional_analysis.pdf"
    params:
        gene = "{gene}"
    resources:
        memory = 16,
        time = 1,
        tmpdir = 0
    priority: 8
    threads: 1
    script:
        "scripts/visualization.py"

def aggregate_protein_coding_potential_analysis(wildcards):
    checkpoint_output = checkpoints.gene_filter.get(**wildcards).output[0]
    return expand("Results/Genes/{gene}/{gene}_functional_analysis.pdf",
           gene=glob_wildcards(os.path.join(checkpoint_output,"{gene}.txt")).gene)

rule output_combine_files:
    input:
        aggregate = aggregate_protein_coding_potential_analysis,
        gene_temp = directory("Results/Genes_temp")
    output:
        plot = "Results/Output/" + config["output_plots"]
    resources:
        memory = 4,
        time = 1,
        tmpdir = 0
    priority: 5
    threads: 1
    shell:
        "rm -rf {input.gene_temp}; pdfunite {input.aggregate} {output.plot};"
