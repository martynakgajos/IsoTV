#!/usr/bin/env python3

# Extract the transcripts that are linked to their respective gene

import pandas as pd
import numpy as np
from Bio import SeqIO

# Import transcript isoform sequences
transcripts_filename = snakemake.input[1]
transcripts = SeqIO.index(transcripts_filename, "fasta")

output = []
gene_name = snakemake.params[0]

# If the annotation file is given
if (snakemake.config["annotation"]):
    # Import the nanopore annotation file
    annotation_filename = snakemake.input[0]
    annotate_df = pd.read_csv(annotation_filename,sep = "\t", header = None, comment='#')
    #annotation = pd.read_csv(sep = "\t", skiprows = 5, names = ["chr","type","info"],usecols = [0,2,8])
    #annotate_df = pd.read_csv(annotation_filename,sep = "\t", skiprows = 5, names = ["chr","type","info"],usecols = [0,2,8], header = None)

    annotate_df = annotate_df[annotate_df[2]  != "exon"]
    annotate_type = list(annotate_df[2])
    annotate_lines = list(annotate_df[8])

    # Mapping gene name to transcriptID
    gene_tID = dict()

    for ann in range(len(annotate_lines)):
        if (("gene_name" in annotate_lines[ann]) and ("transcript_id" in annotate_lines[ann]) and ("transcript" == annotate_type[ann])):
            line = annotate_lines[ann].split(";")
            #tID = line[0].split(" ")[-1][1:-1]
            tID = annotate_lines[ann].split('transcript_id "')[1].split('"')[0]
            #gene = line[2].split(" ")[-1][1:-1]
            gene = annotate_lines[ann].split('gene_name "')[1].split('"')[0]
            if (gene not in gene_tID): gene_tID[gene] = [tID]
            else: gene_tID[gene].append(tID)

    # Extracting isoforms from related genes
    try:
        gene_tID[gene_name]
    except KeyError as error:
        print("The gene %s was not found in the annotation file." %error)
        raise

    for tID in gene_tID[gene_name]:
        output.append(">" + tID)
        output.append(str(transcripts[tID].seq))
# If only the transcriptome file is given
else:
    for transcript in transcripts:
        try:
            if (gene_name == transcipts[transcript].id.split("|")[1]):
                output.append(">" + transcipts[transcript].id.split("|")[0])
                output.append(str(transcipts[transcript].seq))
        except:
            output.append(">" + transcipts[transcript].id.split("|")[0])
            output.append(str(transcipts[transcript].seq))

    if (output == []):
        print("The gene %s was not found in the transcriptome file.")
        raise

output_filename = snakemake.output[0]
output_file = open(output_filename,"w+")
output_file.write("\n".join(output))
output_file.close()
