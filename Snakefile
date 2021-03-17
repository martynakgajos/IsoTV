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


### All
rule all:
  input:
    "Results/Output/" + config["output_plots"]

# Package versions
rule dump_versions:
    output:
        ver = "versions.txt"
    conda:
        "environment.yaml"
    shell:
        "conda list > {output.ver}"

###### Raw ONT Reads Processing

# Basecalling
rule Basecalling:
    input:
        "RawData/{sample}/{id}.fast5"
    output:
        "RawData/QC/{sample}/{id}.txt",
        "RawData/Fast5/{sample}/{id}/{id}.fast5",
        fastq = "RawData/Basecalled/{sample}/{id}.fastq"
    resources:
        memory = 10,
        time = 5,
        tmpdir = 4
    params:
        guppy = config["guppy"],
        flowcell = config["flowcell"],
        kit = config["kit"],
        work = config["outdir"]
    priority: 10
    group: "basecalling"
    threads: 8
	shell: """
        mkdir -p {params.work}/RawData/Fast5/{wildcards.sample}/{wildcards.id}

		ln -sf {params.work}/{input} {params.work}/RawData/Fast5/{wildcards.sample}/{wildcards.id}/{wildcards.id}.fast5

		{params.guppy} -i {params.work}/RawData/Fast5/{wildcards.sample}/{wildcards.id}/ -s {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id}/ \
		--flowcell {params.flowcell} --kit {params.kit} --num_callers 1 --cpu_threads_per_caller {threads} --fast5_out

		mv -f {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id}/workspace/{wildcards.id}.fast5 {params.work}/RawData/Fast5/{wildcards.sample}/{wildcards.id}/{wildcards.id}.fast5

		mv {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id}/sequencing_summary.txt {params.work}/RawData/QC/{wildcards.sample}/{wildcards.id}.txt

		find {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id} -name '*.fastq' -exec mv {{}} {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id}.fastq \;"""

rule aggregate_fast5:
    input:
        f5 = lambda wildcards: expand("RawData/Basecalled/{{sample}}/{id}.fastq", id = glob_wildcards("RawData/"+wildcards.sample+"/{id}.fast5")[0]),
    output:
        fq = "RawData/Fastq/{sample}"+".fastq"
    resources:
        memory = 50,
        time = 4,
        tmpdir = 4
    priority: 9
    run:
        with open(output.fq,'w') as fout:
            for fn in input.f5:
                with open(fn,'r') as fin:
                    shutil.copyfileobj(fin, fout)

		#f5 = lambda wildcards: expand("RawData/Basecalled/{sample}/{id}.fastq"),

# Build minimap2 index to reference genome
rule Minimap2Index:
    input:
        genome = ReferenceFasta
    output:
        index = "Results/Minimap2/"+ReferenceFasta.split("/")[-1]+".mmi"
    resources:
        memory = 32,
        time = 1
    params:
        opts = config["minimap2_index_opts"],
    threads:
        config["threads"]
    shell:
      "minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}"

# Filter basecalled raw ONT reads
rule FilterReads:
    input:
        lambda wildcards: "RawData/Fastq/"+SAMPLES[wildcards.sample]+".fastq"
    output:
        "FilteredData/{sample}.fastq",
    resources:
        memory = 32,
        time = 1
    params:
        config["min_mean_q"]
    shell:
        "filtlong --min_mean_q {params} {input} > {output}"

# Pychopper
rule Pychopper:
    input:
        "FilteredData/{sample}.fastq",
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

# Build known splice junction bed file
rule SpliceJunctionIndex:
    input:
        GenomeGFF
    output:
        junc_bed = "ReferenceData/junctions.bed"
    shell:
        "paftools.js gff2bed {input} > {output.junc_bed}"

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

# Convert BAM to GFF
rule PinfishRawBAM2GFF:
    input:
        bam = rules.Minimap2Pinfish.output.bam
    output:
        raw_gff = "Results/Pinfish/raw_transcripts.gff"
    resources:
        memory = 8,
        time = 1
    params:
        opts = config["spliced_bam2gff_opts"]
    threads: config["threads"]
    shell:
        "spliced_bam2gff {params.opts} -t {threads} -M {input.bam} > {output.raw_gff}"

# Cluster transcripts in GFF
rule PinfishClusterGFF:
    input:
        raw_gff = rules.PinfishRawBAM2GFF.output.raw_gff
    output:
        cls_gff = "Results/Pinfish/clustered_pol_transcripts.gff",
        cls_tab = "Results/Pinfish/cluster_memberships.tsv",
    resources:
        memory = 8,
        time = 1
    params:
        c = config["minimum_cluster_size"],
        d = config["exon_boundary_tolerance"],
        e = config["terminal_exon_boundary_tolerance"],
        min_iso_frac = config["minimum_isoform_percent"],
    threads: config["threads"]
    shell:
        "cluster_gff -p {params.min_iso_frac} -t {threads} -c {params.c} -d {params.d} -e {params.e} -a {output.cls_tab} {input.raw_gff} > {output.cls_gff}"

# Collapse clustered read artifacts
rule PinfishCollapseRawPartials:
    input:
        cls_gff = rules.PinfishClusterGFF.output.cls_gff
    output:
        cls_gff_col = "Results/Pinfish/clustered_transcripts_collapsed.gff"
    resources:
        memory = 8,
        time = 1
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"],
    shell:
       "collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.cls_gff} > {output.cls_gff_col}"

# Polish read clusters
rule PinfishPolishClusters:
    input:
        cls_gff = rules.PinfishClusterGFF.output.cls_gff,
        cls_tab = rules.PinfishClusterGFF.output.cls_tab,
        bam = rules.Minimap2Pinfish.output.bam
    output:
        pol_trs = "Results/Pinfish/polished_transcripts.fas"
    conda:
        "/project/owlmayerTemporary/Sid/nanopore-analysis/environment.yaml"
    resources:
        memory = 8,
        time = 1
    params:
        c = config["minimum_cluster_size"]
    threads: config["threads"]
    shell:
        "polish_clusters -t {threads} -a {input.cls_tab} -c {params.c} -o {output.pol_trs} {input.bam}"

# Map polished transcripts to genome
rule MinimapPolishedClusters:
    input:
       index = rules.Minimap2Index.output.index,
       fasta = rules.PinfishPolishClusters.output.pol_trs,
    output:
       pol_bam = "Results/Minimap2/polished_reads_aln_sorted.bam"
    resources:
        memory = 8,
        time = 1
    params:
        extra = config["minimap2_opts_polished"],
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} {params.extra} -ax splice {input.index} {input.fasta}\
    | samtools view -Sb -F 2304 | samtools sort -@ {threads} - -o {output.pol_bam};
    samtools index {output.pol_bam}
    """

# Convert BAM of polished transcripts to GFF
rule PinfishPolishedBAM2GFF:
    input:
        bam = rules.MinimapPolishedClusters.output.pol_bam
    output:
        pol_gff = "Results/Pinfish/polished_transcripts.gff"
    resources:
        memory = 8,
        time = 1
    params:
        extra = config["spliced_bam2gff_opts_pol"]
    threads: config["threads"]
    shell:
        "spliced_bam2gff {params.extra} -t {threads} -M {input.bam} > {output.pol_gff}"

# Collapse polished read artifacts
rule PinfishCollapsePolishedPartials:
    input:
        pol_gff = rules.PinfishPolishedBAM2GFF.output.pol_gff
    output:
        pol_gff_col = "Results/Pinfish/polished_transcripts_collapsed.gff"
    resources:
        memory = 8,
        time = 1
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"],
    shell:
        "collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.pol_gff} > {output.pol_gff_col}"

rule GffCompare:
    input:
        reference = GenomeGFF,
        exptgff = rules.PinfishCollapsePolishedPartials.output.pol_gff_col
    output:
        nanopore_gtf = "Results/GffCompare/nanopore.combined.gtf"
    resources:
        memory = 16,
        time = 1
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
    script:
        "scripts/nonredundant.py"

rule reads_to_transcripts:
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
    threads: 1
    group: "preprocess_sequences"
    script:
        "scripts/reads_to_transcripts.py"

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
       fastq = "FilteredData/{sample}.fastq",
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
       fastq = "FilteredData/{sample}.fastq"
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

# Quantify transcripts
rule CountTranscripts:
    input:
        bam = rules.ExtractPrimaryMapping.output
    resources:
        memory = 4,
        time = 1
    output:
        quant = "Results/Quantification/{sample}.counts"
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
        "scripts/concatenate.py"

# Normalize isoform levels
rule CalculateDeSeq2Norm:
    input:
        counts = rules.ConcatenateCounts.output
    resources:
        memory = 4,
        time = 1
    output:
        norm_counts = "Results/Quantification/counts_deseq2norm.txt"
    script:
        "scripts/deseq2norm.R"

###### Isoform Transcript Processing

# Use correct file paths
CountsData =  rules.CalculateDeSeq2Norm.output if config["preprocess"] else config["counts_data"]
Transcripts =  rules.reads_to_transcripts.output if config["preprocess"] else config["polished_reads"]
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
        java = config["java"],
        pfam = config["interproScan_path"],
    resources:
        memory = 16,
        time = 1,
        tmpdir = 0
    priority: 8
    threads: 1
    group: "sequence_analysis"
    shell:
        "source {params.java} ; sh  {params.pfam} -i {input.translation} -o {output} -f GFF3 -appl {params.db} -dra"

# Secondary structure analysis
rule porter_analysis:
    input:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa"
    output:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa.ss3"
    params:
        brewery = config["brewery_path"],
        lock1 = config["lock1"],
        lock2 = config["lock2"],
    resources:
        memory = 72,
        time = 1,
        tmpdir = 16
        #memory = 16,
        #time = 1,
        #tmpdir = 0
    priority: 10
    threads: 16
    #group: "sequence_analysis"
    shell:
        "python3 {params.brewery} -i {input} --cpu 16 --noTA --noSA --noCD --fast"
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

# Prion analysis
# rule prion_analysis:
#    input:
#        "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa"
#    output:
#        "Results/Genes_temp/{gene}/transcripts/{transcript}_func_prion.txt"
#    params:
#        plaac = config["plaac_path"],
#        gene = "{gene}"
#    resources:
#        memory = 10,
#        time = 1,
#        tmpdir = 0
#    priority: 8
#    threads: 1
#    group: "sequence_analysis"
#    shell:
#        "java -jar {params.plaac} -i {input} -p all > {output}"

# Combine functional feature analysis on individual isoforms
def functional_files():
    files = []
    if (config["aa"]): files.append("Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa")
    if (config["iupred2a"]): files.append(rules.iupred2a_analysis.output.idr)
    if (config["pfam"]): files.append(rules.interpro_scan.output)
    if (config["porter"]): files.append(rules.porter_analysis.output)
    if (config["pfScan"]): files.append(rules.functional_site_analysis.output)
#    if (config["prion"]): files.append(rules.prion_analysis.output)
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
        aggregate_transcript_analysis_func
    output:
        "Results/Genes_temp/{gene}/{gene}_transcripts_filtered_analysis.txt"
    params:
        gene = "{gene}"
    resources:
        memory = 4,
        time = 1,
        tmpdir = 0
    priority: 5
    threads: 1
    shell:
        "cat {input} > {output}"

rule protein_coding_potential_analysis:
    input:
        rules.aggregate_transcript_analysis.output,
        rules.translation_to_protein.output.stats,
        CountsData,
        NanoporeGTF,
    output:
        "Results/Genes_temp/{gene}/{gene}_functional_analysis.pdf"
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

rule copy_genes_temp:
    input:
        gene_plot = "Results/Genes_temp/{gene}/{gene}_functional_analysis.pdf"
    output:
        gene_folder = directory("Results/Genes/{gene}"),
        plot = "Results/Plots/{gene}_functional_analysis.pdf"
    resources:
        memory = 16,
        time = 1,
        tmpdir = 0
    priority: 8
    threads: 1
    shell:
        "mkdir {output.gene_folder}; cp {input.gene_plot} {output.plot}; cp -r {input} {output.gene_folder};"

def aggregate_protein_coding_potential_analysis(wildcards):
    checkpoint_output = checkpoints.gene_filter.get(**wildcards).output[0]
    return expand("Results/Plots/{gene}_functional_analysis.pdf",
           gene=glob_wildcards(os.path.join(checkpoint_output,"{gene}.txt")).gene)

rule output_combine_files:
    input:
        aggregate = aggregate_protein_coding_potential_analysis,
        gene_temp = directory("Results/Genes_temp"),
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
