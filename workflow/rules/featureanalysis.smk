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
    resources:
        memory = 16,
        time = 1,
        tmpdir = 0
    priority: 8
    threads: 1
    group: "sequence_analysis"
    conda:
        "../envs/interproscan.yaml"
    shell:
        "interproscan.sh -i {input.translation} -o {output} -f GFF3 -appl {params.db} -dra"

# Secondary structure analysis
rule brewery_analysis:
    input:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa"
    output:
        "Results/Genes_temp/{gene}/transcripts/{transcript}_protein.fa.ss3"
    params:
        brewery = config["brewery_path"],
    resources:
        memory = 8,
        time = 1,
        tmpdir = 16
    priority: 10
    threads:
        config["threads"]
    conda:
        "../envs/brewery.yaml"
    shell:
        "python3.7 {params.brewery} v {input} --cpu {threads} --noTA --noSA --noCD --fast"

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
        "../scripts/individual_transcript_analysis.py"
