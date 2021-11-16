#!/usr/bin/env python3

# Select only non-redundant transcripts

import pandas as pd
from Bio import SeqIO

transcripts_filename = snakemake.input[0]
transcripts = SeqIO.index(transcripts_filename, "fasta")

ds=pd.read_csv(snakemake.input[1],delimiter='\t',usecols=[2,8],header=None, names=['type','info'])
ds=ds[ds["type"]=='transcript']
ds['transcript_id']=ds.apply(lambda x: x['info'].split('oId "')[1].split('";')[0],1)
non_redundant_tid = set(ds["transcript_id"])
output = []
for transcript in transcripts:
    if (transcript in non_redundant_tid):
        output.append(">" + transcript.strip())
        output.append(str(transcripts[transcript].seq).strip())
output_filename = snakemake.output[0]
output_file = open(output_filename, "w+")
output_file.writelines("\n".join(output))
output_file.close()
