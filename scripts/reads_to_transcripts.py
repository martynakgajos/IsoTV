#!/usr/bin/env python3

# Converts transcriptome file from reads (oID names) to transcripts name (TCONS)

from Bio import SeqIO
import pandas as pd

transcripts_filename = snakemake.input[0]
transcripts = SeqIO.index(transcripts_filename, "fasta")

annotation_filename = snakemake.input[1]

ds=pd.read_csv(annotation_filename,delimiter='\t',usecols=[2,8],header=None, names=['type','info'])
ds=ds[ds["type"]=='transcript']
ds['transcript']=ds.apply(lambda x: x['info'].split('oId "')[1].split('";')[0],1)
ds['transcript_id']=ds.apply(lambda x: x['info'].split('transcript_id "')[1].split('";')[0],1)
test = dict()
index = list(ds['transcript'])
refs = list(ds["transcript_id"])
for i in range(len(index)):
    test[index[i]] = refs[i]


output = []
for transcript in transcripts:
    try:
        output.append(">" + str(test[transcripts[transcript].id]))
        output.append(str(transcripts[transcript].seq))
    except:
        pass

output_filename = snakemake.output[0]
output_file = open(output_filename,"w+")
output_file.write("\n".join(output))
output_file.close()
