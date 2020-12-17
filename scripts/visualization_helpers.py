#!/usr/bin/env python3

import string
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# change font
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
plt.rc('font', size=14)

alph = string.ascii_lowercase
alph += "0123456789!@#$%^&*?;'{}-=+,.~_[]:<>`qwertyuiopasdfghjklzxcvbnm1234567890qwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnm"

class transcript:
    def __init__(self, tcons, ns_decay, aa, idr = [], dom = [], ss3 = [], plc = [], pss = []):
        self.tcons = tcons
        self.ns_decay = ns_decay
        self.aa = aa
        self.idr = idr
        self.dom = dom
        self.ss3 = ss3
        self.plc = plc
        self.pss = pss

class Exon:
    def __init__(self, tcons, old_id, start, stop, new_id = ""):
        self.tcons = tcons
        self.old_id = old_id
        self.new_id = new_id
        self.start = start
        self.stop = stop

    def __lt__(self, other):
         return self.start < other.start

    def add_new_id(self, id):
        self.new_id = id

    def __str__(self):
        return str(self.new_id) + " " + str(self.start) + " " + str(self.stop)

class Parser(object):
    def __init__(self, csv, gtf, samples):
        self.csv = csv
        self.gtf = gtf
        self.samples = samples

def preprocessArguments(args):
    if (args.gtf == None):
        annotation = None
    else:
        annotation = pd.read_csv(args.gtf, delimiter='\t', header=None, usecols=[0,2,3,4,6,8],names=['chrm','type','start','stop','strand','more'])
        annotation['transcript_id'] = annotation.apply(lambda x: x['more'].split('transcript_id "')[1].split('"')[0],1)
        annotation["transcript_id"] = annotation["transcript_id"].apply(lambda x: x.replace("_",""))
        annotation["exon_number"] = annotation.apply(lambda x: x["more"].split('exon_number "')[1].split('"')[0] if "exon_number" in x["more"] else 0,1)
        annotation = annotation.drop(columns='more')

    if (args.csv == None):
        data = None
    else:
        data = pd.read_csv(args.csv)

    samples = [sample for sample in args.samples]
    conditions = sorted(list(set([sample.split("_")[0] for sample in samples])))
    number_replicates={}
    for cond in conditions:
        number_replicates[cond]  = len([sample for sample in samples if sample.startswith(cond)])
    return annotation, data, samples, conditions, number_replicates

def calculateStatistics(df,conds,nreps):
    for cond in conds:
        df['mean'+cond]=df.filter(like=cond+'_').mean(1)
        df['stdn'+cond]=df.filter(like=cond+'_').std(1)/np.sqrt(nreps[cond])
    df=df.sort_index()
    return df

def chooseIsoforms2Plot(df,minTPM,minPct,maxIso):
    df['minimum']=df.filter(regex='^mean').min(axis=1)
    df['maximumPct']=df.filter(regex='^Pct').max(axis=1)
    df=df[df['maximumPct']>minPct]
    df['maximum']=df.filter(regex='^mean').max(axis=1)
    df=df[df['maximum']>minTPM]
    df['avg'] = df.filter(regex='^Pct').mean(axis=1)
    df=df.sort_values('avg',ascending=False)
    df=df.head(maxIso)
    return df

def plotTotal(conditions, df_gene, ax, colors, number_replicates, continuous):
    conditions_x = [i + 1 for i in range(len(conditions))]
    if (continuous):
        plt.errorbar(conditions_x, df_gene.filter(like = 'mean').iloc[0], yerr = df_gene.filter(like='stdn').iloc[0], color = '#FA3D36', linewidth = 3, label = "Total Expression", markersize=6)
    else:
        x = []
        boxes = []
        count = 0
        for i in range(len(conditions)):
            x += [i + 1] * number_replicates[conditions[i]]
            boxes.append(df_gene.filter(like = '_').values[0][count:count+number_replicates[conditions[i]]])
            count += number_replicates[conditions[i]]
        bp = plt.boxplot(boxes, patch_artist=True, showfliers=False, zorder=1)
        for patch in bp['boxes']:
            patch.set_facecolor("#F7F7F7")
        for median in bp['medians']:
            median.set(color ='black',linewidth = 1)

        plt.scatter(x + np.random.normal(0, 0.075, len(x)), df_gene.filter(like = '_'), color = "black", marker = '.', s = 200, alpha = 0.5, linewidths = 0, zorder=2)

    ax.set_ylim(bottom = -0.5)
    ax.set_xticks(conditions_x)
    ax.set_xticklabels(conditions)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.spines['left'].set_bounds(0,ax.get_yticks()[-2])
    ax.spines['bottom'].set_bounds(min(conditions_x),max(conditions_x))
    ax.tick_params(axis = "both", which = "major", labelsize = 14)
    plt.xlabel("Condition")
    plt.ylabel("Expression [TPM]")
    #plt.legend(loc = "upper center", bbox_to_anchor=(0.5,1), frameon=False)
    ax.set_title("Total Gene Expression", loc = "center", fontsize = 14, weight='bold')

def plotIsoformProfiles(conditions, df, ax, colors, continuous):
    x = [i + 1 for i in range(len(conditions))]
    last_bottom = [0] * len(conditions)
    for j in range(df.shape[0]):
        row = df.iloc[j]
        if (continuous):
            plt.errorbar(x + np.random.normal(0, 0.03, len(x)), row.filter(regex='^mean'),yerr=row.filter(like='stdn'),color=colors[j],linewidth=2, label = "")
        else:
            mean = row.filter(regex = '^mean')
            err = row.filter(like='stdn')
            jitter = 0.5 * (j + 1)/(df.shape[0] + 1) - 0.25
            plt.bar(x, mean, width = 0.5, bottom = last_bottom, color = colors[j])
            for i in x:
                plt.plot([i + jitter,i + jitter],[last_bottom[i - 1] + mean.iloc[i - 1] - err.iloc[i - 1],last_bottom[i - 1] + mean.iloc[i - 1] + err.iloc[i - 1]],linewidth = 2,color = "black")
            last_bottom += mean

    ax.set_ylim(bottom=-0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(conditions)
    ax.tick_params(axis = "both", which = "major", labelsize = 14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.spines['left'].set_bounds(0,ax.get_yticks()[-2])
    ax.spines['bottom'].set_bounds(min(x),max(x))
    plt.xlabel("Condition", fontsize = 14)
    plt.ylabel("Expression [TPM]", fontsize = 14)
    plt.legend(loc = "upper center", bbox_to_anchor=(0.5,1), frameon=False)
    ax.set_title("Individual Isoform Expression", loc = "center", fontsize = 14, weight='bold')

def plotIsoformStacked(conditions, df, ax, colors, continuous):
    x = [i + 1 for i in range(len(conditions))]
    if (continuous):
        plt.stackplot(x,df.filter(regex='^Pct').values,colors=colors[:df.shape[0]])
    else:
        x = [i + 1 for i in range(len(conditions))]
        last_bottom = [0] * len(conditions)
        for j in range(df.shape[0]):
            row = df.iloc[j]
            pct = row.filter(regex = '^Pct')
            #err = row.filter(like='stdn')
            err = np.divide(np.multiply(row.filter(like='stdn').values,row.filter(regex = '^Pct').values),row.filter(regex = '^mean').values)
            jitter = (j - df.shape[0]//2)/10
            plt.bar(x, pct, width = 0.5, bottom = last_bottom, color=colors[j])
            #for i in x:
            #    plt.plot([i + jitter,i + jitter],[last_bottom[i - 1] + pct.iloc[i - 1] - err[i - 1],last_bottom[i - 1] + pct.iloc[i - 1] + err[i - 1]],linewidth = 2,color = "black")
            last_bottom += pct

    ax.set_xticks(x)
    ax.set_xticklabels(conditions)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_bounds(0,100)
    ax.spines['bottom'].set_bounds(min(x),max(x))
    ax.axis(top = [-0.5,100])
    ax.set_yticks([i for i in range(0,101,20)])
    ax.tick_params(axis = "both", which = "major", labelsize = 14)
    plt.xlabel("Condition")
    plt.ylabel("Expression [%]")
    ax.set_title("Individual Isoform Usage", loc = "center", fontsize = 14, weight='bold')

def prepareAnnotation(annotation,transcripts_plot):
    cut = annotation[annotation['transcript_id'].isin(transcripts_plot)]
    strand = cut.iloc[0]['strand']
    if (strand == '+'):
        start = cut['start'].min()
        cut['plot_start'] = cut['start'] - start
        cut['plot_stop'] = cut['stop'] - start
    else:
        start = cut['stop'].max()
        cut['plot_start'] = (cut['start'] - start) * (-1)
        cut['plot_stop'] = (cut['stop'] - start) * (-1)
    return cut
