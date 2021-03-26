# Determine which genes to analyze in the input

import pandas as pd
import os
from Bio import SeqIO


# Determines if a gene already exists, and create a dummy file and folder if it does
def gene_exists_file_func(gene):
    gene_temp_dir = os.path.join(snakemake.output[0],gene)
    gene_final_dir = os.path.join(snakemake.params["genes_final_dir"],gene)

    gene_exists_file = os.path.join(snakemake.output[0],"%s.txt" %(gene))
    if (not os.path.isfile(gene_exists_file)):
        f = open(gene_exists_file, "w+")
        f.close()

    if ((not os.path.isdir(gene_final_dir)) and (not os.path.isdir(gene_temp_dir))):
        os.makedirs(gene_temp_dir)


gene_file = snakemake.config["gene_file"]
genes = list(set(pd.read_csv(gene_file, sep = "\t", names = ["gene_symbol"])["gene_symbol"]))
os.makedirs(snakemake.output[0])

# If the annotation file is given
if (snakemake.config["annotation"]):
    annotate_df = pd.read_csv(snakemake.input[0], sep = "\t", header = None, comment='#')
    annotate_df = annotate_df[annotate_df[2]  != "exon"]
    annotate_lines = list(annotate_df[8])

    for ann in range(len(annotate_lines)):
        if ("gene_name" in annotate_lines[ann]):
            gene = annotate_lines[ann].split('gene_name "')[1].split('"')[0]
            if (gene in genes):
                gene_exists_file_func(gene)
else: # If only the transcriptome file is given
    transcripts_filename = snakemake.input[0]
    transcripts = SeqIO.index(transcripts_filename, "fasta")
    try:
        for transcript in transcripts:
            gene = transcipts[transcript].id.split("|")[1]
            if (gene in genes):
                gene_exists_file_func(gene)
    except:
        for gene in genes:
            gene_exists_file_func(gene)
