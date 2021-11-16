#!/usr/bin/env python3

# Concatentate and annotate isoforms and their expression levels 

import pandas as pd

# Concatenate datasets
l=[]
for fn in snakemake.input[2:]:
    df=pd.read_csv(fn,delimiter='\t',usecols=[0,1],index_col=[1])
    df=df.rename(columns={'counts':fn.split('.counts')[0].split('/')[-1]})
    l.append(df)
df=l[0]
for i in range(1,len(l)):
    df=df.join(l[i],how='outer')
df=df.fillna(0)

# Get connection to reference
ds=pd.read_csv(snakemake.input[1],delimiter='\t',usecols=[2,8],header=None, names=['type','info'])
ds=ds[ds["type"]=='transcript']
ds['transcript']=ds.apply(lambda x: x['info'].split('oId "')[1].split('";')[0],1)
ds['transcript_id']=ds.apply(lambda x: x['info'].split('transcript_id "')[1].split('";')[0],1)
ds['class_code']=ds.apply(lambda x: x['info'].split('class_code "')[1].split('";')[0],1)
ds['gene_name']=ds.apply(lambda x: '' if x['class_code']=='u' else x['info'].split('gene_name "')[1].split('";')[0],1)
ds['ref_transcript']=ds.apply(lambda x: '' if x['class_code']=='u' else x['info'].split('cmp_ref "')[1].split('";')[0],1)
ds=ds.set_index('transcript')
df=df.join(ds,how='left')

# Get gene names
d=[]
with open(snakemake.input[0],'r') as af:
    for linia in af:
        if linia[0]!='#':
            ln=linia.strip().split()
            if ln[2]=='gene':
                gene_name=linia.split('gene_name "')[1].split('"')[0]
                gene_id=linia.split('gene_id "')[1].split('"')[0]
                gene_type=linia.split('gene_type "')[1].split('"')[0]
                d.append([gene_name,gene_id,gene_type])
d=pd.DataFrame(data=d,columns=['gene_name','gene_id','gene_type'])

df=df.merge(d,how='left',on='gene_name')
df=df.drop(columns=['info','type'])
cols = list(df.columns)
cols_new = ["class_code"] + cols[:cols.index("class_code")] + cols[cols.index("class_code")+ 1:]
df=df.reindex(columns = cols_new)
df.index = df["transcript_id"]
df.sort_index().to_csv(snakemake.output[0],index=False)
