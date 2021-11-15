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
