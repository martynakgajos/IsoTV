import copy
import string
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.transforms import Bbox
from matplotlib.collections import PatchCollection
import matplotlib.gridspec as gridspec
from visualization_helpers import *


gene = snakemake.params[0]
alph = string.ascii_uppercase
alph += "0123456789!@#$%^&*?;'{}-=+,.~_[]:<>`qwertyuiopasdfghjklzxcvbnm1234567890qwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnm"

#grey, orange, turquoise blue, blue green, yellow, royal blue, vermillion, red purple,
#pastel yellow, brown orange, burgundy, cerulean blue, dark purple, light purple, dark blue green
#pastel green, harvest brown, brown, sky blue, pink, dark pink
colors=["#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7","#F0D648",
        '#db6d00','#920000', '#6db6ff', '#b66dff', '#490092', '#ffff6d', '#b66dff', '#004949',
        '#24C324', '#924900', '#9a6324', 'ff6db6', 'CC79A7']
colors_ptm = ["#DB6D00","#920000","#B66DFF","#0072B2","#24C324"]
colors_ss3 = ['#FFC20A','#0C7BDC',"white"]
ss3_abbvs = ["H","E","C"]
colors_aa = ['#99BFFF', '#F56345', '#E0C830', '#8700F5', 'white']
aa_abbvs = {"A":"H","C":"S","D":"N","E":"N","F":"H","G":"S","H":"P","I":"H","K":"P","L":"H","M":"H","N":"L","P":"H","Q":"L","R":"P","S":"L","T":"L","V":"H","W":"H","Y":"H"}
aa_groups = ["P","N","L","S","H"]

# VISUALIZATION FUNCTIONS

# change font
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
plt.rc('font', size=14)

# Plotting Annotation
def plotAnnotation(annotation, transcripts_plot, colors, longest_length_protein, fig, gs, start_row_panel):
    curr_row_panel = start_row_panel

    exon_dict = dict()
    for tcons_curr in transcripts_plot:
        transcript_annotation = annotation[annotation['transcript_id'] == tcons_curr]

        transcript_exon_start = 0
        transcript_exon_end = -1
        if (transcript_annotation.shape[0] > 0):
            if (transcript_annotation.iloc[0]["strand"] == "-"):
                transcript_annotation = transcript_annotation.iloc[::-1]
                plot_start = copy.deepcopy(transcript_annotation["plot_start"])
                transcript_annotation["plot_start"] = transcript_annotation["plot_stop"]
                transcript_annotation["plot_stop"] = plot_start
            for idx,row in transcript_annotation.iterrows():
                if (row['type'] != 'transcript'):
                    exon_old_id = row["exon_number"]
                    start,stop = row['plot_start'],row['plot_stop']
                    transcript_exon_start = transcript_exon_end + 1
                    transcript_exon_end += abs(start - stop) + 1
                    exon = Exon(tcons_curr, exon_old_id, start, stop)

                    if (start in exon_dict):
                        exon_dict[start].append(exon)
                    else:
                        exon_dict[start] = [exon]

                    if (stop in exon_dict):
                        exon_dict[stop].append(exon)
                    else:
                        exon_dict[stop] = [exon]

    all_exons = []
    curr_exons = []
    count = 0
    master = dict()
    for key,value in sorted(exon_dict.items()):
        end_exons = []
        for exon in value:
            curr_exons.append(exon) if (exon not in curr_exons) else end_exons.append(exon)
        end_exons.sort()
        prev_start,prev_end = -1,-1
        for exon in end_exons:
            if (prev_start != exon.start or prev_end != exon.stop):
                count += 1
                exon.add_new_id(count)
                all_exons.append(exon)
            else:
                exon.add_new_id(count)
            prev_start,prev_end = exon.start,exon.stop
            if (exon.tcons not in master):
                master[exon.tcons] = dict()
            master[exon.tcons][exon.old_id] = exon.new_id

    count = -1
    previous_start,previous_stop = -1,-1
    exon_group = []
    exon_distances = []
    for exon in all_exons:
        if (exon.start > previous_stop):
            exon_group.append([exon])
            exon_distances.append([exon.start,exon.stop])
            previous_start,previous_stop = exon.start,exon.stop
            count += 1
        else:
            exon_group[count].append(exon)
            if (exon.stop > previous_stop):
                previous_stop = exon.stop
                exon_distances[count][1] = exon.stop
            if (exon.start < previous_start):
                previous_start = exon.start
                exon_distances[count][0] = exon.start

    exon_distances_original = copy.deepcopy(exon_distances)
    previous_stop = -1
    for exon_dist in exon_distances:
        start,stop = exon_dist[0],exon_dist[1]
        if ((start - previous_stop) > 200):
            exon_dist[0] = previous_stop + 200
            exon_dist[1] = exon_dist[0] + stop - start
        previous_stop = exon_dist[1]

    for exon_similar_group in range(len(exon_group)):
        for exon in exon_group[exon_similar_group]:
            start,stop = exon.start,exon.stop
            if (exon.start != exon_distances[exon_similar_group][0]):
                exon.start = exon_distances[exon_similar_group][0] + (exon.start - exon_distances_original[exon_similar_group][0])
                exon.stop = exon.start + (stop - start)

    for tcons_curr in transcripts_plot:
        transcript_annotation = annotation[annotation['transcript_id'] == tcons_curr]
        transcript_annotation["exon_number"] = transcript_annotation.apply(lambda x: master[tcons_curr][str(x["exon_number"])] if (x["type"] != "transcript") else 0, axis = 1)
        transcript_exons = transcript_annotation.iloc[1:]
        start = all_exons[transcript_exons.loc[transcript_exons['exon_number'].idxmin()]["exon_number"] - 1].start
        stop = all_exons[transcript_exons.loc[transcript_exons['exon_number'].idxmax()]["exon_number"] - 1].stop
        transcript_annotation["plot_start"] = transcript_annotation.apply(lambda x: all_exons[x["exon_number"] - 1].start if (x["type"] != "transcript") else start, axis = 1)
        transcript_annotation["plot_stop"] = transcript_annotation.apply(lambda x: all_exons[x["exon_number"] - 1].stop if (x["type"] != "transcript") else stop, axis = 1)
        annotation[annotation['transcript_id'] == tcons_curr] = transcript_annotation

    # Plot Transcript
    longest = max(annotation.loc[annotation['plot_start'].idxmax()]["plot_start"],annotation.loc[annotation['plot_stop'].idxmax()]["plot_stop"])
    for tcons_curr in transcripts_plot:
        transcript_color = colors[curr_row_panel - start_row_panel]
        ax = fig.add_subplot(gs[curr_row_panel, :])
        ax.set_xlim((-10, longest+1))
        ax.set_ylim((0, 1.5))
        ax.text(0,1,tcons_curr,horizontalalignment = "left",verticalalignment = "top", fontsize = 14, color = "black")

        transcript_annotation = annotation[annotation['transcript_id'] == tcons_curr]

        transcript_info = transcript_annotation.iloc[0]
        plt.plot([transcript_info["plot_start"],transcript_info["plot_stop"]],[0.5,0.5],color=transcript_color,linewidth=2)

        transcript_annotation = transcript_annotation.iloc[1:]
        if (transcript_annotation.iloc[0]["strand"] == "-"):
            transcript_annotation = transcript_annotation.iloc[::-1]
        prev_intron_start = -1
        for idx,row in transcript_annotation.iterrows():
            exon_start,exon_stop = row['plot_start'],row['plot_stop']
            exon = alph[int(row["exon_number"]) - 1]
            plt.plot([exon_start,exon_stop],[0.5,0.5],color = transcript_color,linewidth = 15)
            ax.annotate(exon,(exon_start/2 + exon_stop/2,0.49), fontsize = 14, color = "white", weight='bold', ha = "center", va = "center")

            if (prev_intron_start > 0 and exon_start - prev_intron_start >= 200):
                intron_shortening = [exon_distances[group - 1][1]/2 + exon_distances[group][0]/2 for group in range(len(exon_distances)) if (prev_intron_start <= exon_distances[group][0] and exon_start >= exon_distances[group][0])]
                for intron in intron_shortening:
                    plt.plot([intron - 10,intron + 10],[0.5,0.5],color = "white",linewidth = 3)
                    plt.plot([intron - 20,intron - 10],[0.3,0.7],color = transcript_color,linewidth = 2)
                    plt.plot([intron + 10,intron + 20],[0.3,0.7],color = transcript_color,linewidth = 2)
            prev_intron_start = exon_stop

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_yticks([])
        ax.set_xticks([])

        curr_row_panel += 1

# Plotting Functional Feature Analysis
def plotSequenceAnalysis(transcripts_plot, colors, longest_length_protein, fig, gs, start_row_panel, annotation = None):
    curr_row_panel = start_row_panel

    for isoform in range(len(transcripts_plot)):
        tcons_curr = transcripts[transcripts_plot[isoform]]

        transcript_color = colors[isoform]

        buffer_idr = 1.0
        buffer_domain = 0.4
        buffer_exon = 0.25
        buffer_aa = 0.3
        buffer_ss = 0.25
        buffer_phos = 0.075
        buffer_den = 0.5

        lg = fig.add_subplot(gs[curr_row_panel:curr_row_panel+total_features, 0])
        lg.set_xlim((-0.11,0.11))

        ax = fig.add_subplot(gs[curr_row_panel:curr_row_panel+total_features, 1:])
        ax.set_xlim((-0.5,longest_length_protein + 1))
        ax.patch.set_facecolor(transcript_color)
        ax.patch.set_alpha(0.15)

        lg.set_xlabel('AA Position', fontsize = 14)
        plt.plot((0,len(tcons_curr.aa)), (0,0),color = "black")
        bottom = 0
        exon_draw = []

        # IDR region
        if (snakemake.config["iupred2a"]):
            idr_df = pd.DataFrame(tcons_curr.idr[1:], columns = tcons_curr.idr[0])
            idr_df = idr_df.astype({'# POS': 'int32',"IUPRED2":"float"})
            idr_df["# POS"] = idr_df["# POS"] - 0.5
            plt.plot(idr_df["# POS"], bottom + idr_df["IUPRED2"], color = transcript_color)
            lg.annotate("Disordered",(0,bottom + buffer_idr * 0.65), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
            lg.annotate("Region",(0,bottom + buffer_idr * 0.35), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
            lg.plot((-0.1,0.1), (bottom + 1,bottom + 1),color = "black")
            bottom += 1.05
            exon_draw.append([bottom - 1.05,bottom])

        # Protein domain
        r = fig.canvas.get_renderer()
        if (snakemake.config["pfam"]):
            for i in tcons_curr.dom:
                start,stop,name = (int(i[0]) - 1),(int(i[1]) - 1),i[2].strip()
                name = name.replace(" domain","")
                plt.gca().add_patch(Rectangle((start, bottom), stop - start, buffer_domain, edgecolor = transcript_color, facecolor = transcript_color, alpha = 0.5))
                text = ax.annotate(name, (start/2.0 + stop/2.0, bottom + buffer_domain/2), fontsize = 14, color = "black", weight='bold', ha = "center", va = "center")
                bbtext =  Bbox(ax.transData.inverted().transform(text.get_window_extent(renderer = r)))
                if (bbtext.x0 < start and bbtext.x1 > stop):
                    text.set_visible(False)
                    name_split = name.split(" ")
                    name_top = " ".join(name_split[:int((len(name_split)/2))])
                    name_bottom = " ".join(name_split[int((len(name_split)/2)):])
                    ax.annotate(name_top, (start/2.0 + stop/2.0, bottom + buffer_domain * 3/4), fontsize = 14, color = "black", weight='bold', ha = "center", va = "center")
                    ax.annotate(name_bottom, (start/2.0 + stop/2.0, bottom + buffer_domain * 1/4), fontsize = 14, color = "black", weight='bold', ha = "center", va = "center")
            lg.annotate("Domain", (0, bottom + buffer_domain/2), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
            bottom += 0.5

        # Secondary structure prediction
        if (snakemake.config["brewery"]):
            ss3_df = pd.DataFrame(tcons_curr.ss3[1:], columns = tcons_curr.ss3[0])
            ss = list(ss3_df["SS"])
            for i in range(len(tcons_curr.aa)):
                plt.gca().add_patch(Rectangle((i,bottom),1,buffer_ss,edgecolor=colors_ss3[ss3_abbvs.index(ss[i])],facecolor=colors_ss3[ss3_abbvs.index(ss[i])]))
            lg.annotate("Structure",(0,bottom + buffer_ss/2), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
            bottom += 0.4

        # Amino acid sequence
        if (snakemake.config["aa"]):
            for i in range(len(tcons_curr.aa)):
                plt.gca().add_patch(Rectangle((i,bottom),1,buffer_aa, edgecolor=colors_aa[aa_groups.index(aa_abbvs[tcons_curr.aa[i]])],facecolor=colors_aa[aa_groups.index(aa_abbvs[tcons_curr.aa[i]])]))
            lg.annotate("Amino Acid",(0,bottom + buffer_aa/2), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
            bottom += 0.4

        # Phosphorlyation Site
        if (snakemake.config["pfScan"]):
            row_phos = 0
            sites_focus = ["AMIDATION","MYRISTYL","GLYCOSYLATION","HYDROXYL","PHOSPHO_SITE"]

            plt.plot([0,len(tcons_curr.aa)],[bottom + 0.1,bottom + 0.1], linewidth = 4, color = "#919191")
            for site_type in sites_focus:
                if (site_type in tcons_curr.pss):
                    for i in tcons_curr.pss[site_type]:
                        start,stop = (int(i[0]) - 1),(int(i[1]) - 1)
                        mid = start/2 + stop/2
                        radius_x = longest_length_protein/150 + 1
                        radius_y = -longest_length_protein/8000 + 0.22
                        patches = [Rectangle((mid - 0.5,bottom + 0.1),1,0.2,facecolor=colors_ptm[row_phos])]
                        patches.append(Ellipse((mid,bottom + 0.3),radius_x,radius_y,facecolor=colors_ptm[row_phos]))
                        p = PatchCollection(patches)
                        p.set_color(colors_ptm[row_phos])
                        ax.add_collection(p)
                row_phos += 1

            lg.annotate("PTM",(0,bottom + row_phos * buffer_phos/2), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
            bottom += 0.45

            # Phosphorlyation Site Density
            density = [0.0] * len(tcons_curr.aa)
            for site_type in tcons_curr.pss:
                for i in tcons_curr.pss[site_type]:
                    start,stop = (int(i[0]) - 1),(int(i[1]) - 1)
                    hamming_values = list(np.hamming(stop - start + 1 + 10))
                    for j in range(max(0,start - 5),min(stop + 1 + 5,len(tcons_curr.aa))):
                        density[j] += hamming_values[j - max(0,start - 5)]
            density = np.array(density)
            density = 0.5 * density/max(density)
            density_x = [i + 0.5 for i in range(len(tcons_curr.aa))]
            plt.plot(density_x, density + bottom, color =  transcript_color)
            lg.annotate("PTM Density",(0,bottom + buffer_den/2), fontsize = 13, color = "black", ha = "center", va = "center", weight='bold')
            bottom += 0.55
            exon_draw.append([bottom - 0.6,bottom])

        # Exons
        if (snakemake.config["annotation"]):
            transcript_annotation = annotation[annotation['transcript_id'] == transcripts_plot[isoform]]
            transcript_stats = transcript_stats_df[transcript_stats_df["transcript_id"] == transcripts_plot[isoform]]
            translate = False
            protein_start,protein_end = int(transcript_stats["start"]),int(transcript_stats["end"])
            transcript_exon_start = 0
            transcript_exon_end = -1
            protein_exon_start = 0
            protein_exon_end = 0
            if (transcript_annotation.shape[0] > 0):
                #if (transcript_annotation.iloc[0]["strand"] == "-"):
                #    transcript_annotation = transcript_annotation.iloc[::-1]
                for idx,row in transcript_annotation.iterrows():
                    if (row['type'] != 'transcript'):
                        start = row['plot_start']
                        stop = row['plot_stop']
                        transcript_exon_start = transcript_exon_end + 1
                        transcript_exon_end += abs(start - stop) + 1
                        exon = alph[int(row["exon_number"]) - 1]

                        if ((transcript_exon_start <= protein_end) and (protein_end <= transcript_exon_end)):
                            translate = False
                            protein_exon_start = protein_exon_end
                            protein_exon_end = (protein_end - protein_start)/3
                            plotExonLines(protein_exon_start,protein_exon_end,bottom,buffer_exon,exon_draw)
                            ax.annotate(exon,(protein_exon_start/2 + protein_exon_end/2 - 1,bottom + buffer_exon/2), fontsize = 14, color = "black", weight='bold', ha = "center", va = "center")
                            break

                        if (translate):
                            protein_exon_start = protein_exon_end
                            protein_exon_end += abs(start - stop)/3
                            plotExonLines(protein_exon_start,protein_exon_end,bottom,buffer_exon,exon_draw)
                            ax.annotate(exon,(protein_exon_start/2 + protein_exon_end/2 - 1,bottom + buffer_exon/2), fontsize = 14, color = "black", weight='bold', ha = "center", va = "center")
                            protein_exon_start = protein_exon_end

                        if ((transcript_exon_start <= protein_start) and (protein_start <= transcript_exon_end)):
                            translate = True
                            protein_exon_start = 0
                            protein_exon_end = abs(transcript_exon_end - protein_start)/3
                            plotExonLines(protein_exon_start,protein_exon_end,bottom,buffer_exon,exon_draw)
                            ax.annotate(exon,(protein_exon_start/2 + protein_exon_end/2 - 1,bottom + buffer_exon/2), fontsize = 14, color = "black", weight='bold', ha = "center", va = "center")
            lg.annotate("Exon",(0,bottom + buffer_exon/2), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
            bottom += 0.25

        plt.gca().add_patch(Rectangle((len(tcons_curr.aa), 0), longest_length_protein - len(tcons_curr.aa) + 1, bottom, edgecolor = "white", facecolor = "white"))

        plt.plot((0,len(tcons_curr.aa)), (bottom,bottom),color = "black")
        ax.set_title(tcons_curr.tcons + str("***" if tcons_curr.ns_decay else ""), loc = "left", fontsize = 14)
        # Finish Formatting
        lg.set_ylim((0, bottom))
        lg.spines['top'].set_visible(False)
        lg.spines['right'].set_visible(False)
        lg.spines['left'].set_visible(False)
        #lg.spines['bottom'].set_visible(False)
        lg.set_xticks([])
        lg.set_yticks([])

        ax.set_ylim((0, bottom))
        plt.plot((0,0),(0,bottom),color = "black")
        plt.plot((len(tcons_curr.aa),len(tcons_curr.aa)),(0,bottom),color = "black")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        #ax.spines['bottom'].set_visible(False)
        ax.set_yticks([])
        curr_row_panel += total_features + 1

# Plotting legend for the netire visualization
def plotLegend(plt,fig,gs,total_gs):
    ax = fig.add_subplot(gs[total_gs-2:total_gs, 0:5])
    ax.set_xlim((0, 13))
    ax.set_ylim((0, 2))
    mods_abbvs = ["Amidation","Myristoylation","Glycoslyation","Hydroxylation","Phosphorlyation"]
    for i in range(0,len(mods_abbvs)):
        plt.gca().add_patch(Rectangle((3,float(i*2)/len(mods_abbvs)),3,1.0/len(mods_abbvs),edgecolor=colors_ptm[i],facecolor=colors_ptm[i]))
        ax.text(7,float(i*2)/len(mods_abbvs) + 0.2,mods_abbvs[i],horizontalalignment = "left",verticalalignment = "top")
    ax.annotate("Post-translation Modification (PTM)",(6.75,2), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
    ax.axis("off")

    ax = fig.add_subplot(gs[total_gs-2:total_gs, 5:10])
    ax.set_xlim((0, 13))
    ax.set_ylim((0, 2))
    ss_abbvs = ["Alpha Helix","Beta Sheet","Other"]
    for i in range(0,len(ss_abbvs)):
        if (ss_abbvs[i] == "Other"):
            plt.gca().add_patch(Rectangle((3,float(i*2 + 4)/5),3,1.0/5,edgecolor="black",facecolor=colors_ss3[i]))
        else:
            plt.gca().add_patch(Rectangle((3,float(i*2 + 4)/5),3,1.0/5,edgecolor=colors_ss3[i],facecolor=colors_ss3[i]))
        ax.text(7,float(i*2 + 4)/5 + 0.25,ss_abbvs[i],horizontalalignment = "left",verticalalignment = "top")
    ax.annotate("Secondary Structure",(6.75,2), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
    ax.axis("off")

    ax = fig.add_subplot(gs[total_gs-2:total_gs, 10:15])
    ax.set_xlim((0, 13))
    ax.set_ylim((0, 2))
    aa_abbvs = ["Positively Charged","Negatively Charged","Polar","Special","Hydrophobic"]
    for i in range(0,len(aa_abbvs)):
        if (aa_abbvs[i] == "Hydrophobic"):
            plt.gca().add_patch(Rectangle((3,float(i*2)/len(aa_abbvs)),3,1.0/len(aa_abbvs),edgecolor="black",facecolor=colors_aa[i]))
        else:
            plt.gca().add_patch(Rectangle((3,float(i*2)/len(aa_abbvs)),3,1.0/len(aa_abbvs),edgecolor=colors_aa[i],facecolor=colors_aa[i]))
        ax.text(7,float(i*2)/len(aa_abbvs) + 0.2,aa_abbvs[i],horizontalalignment = "left",verticalalignment = "top")
    ax.annotate("Amino Acid",(6.75,2), fontsize = 14, color = "black", ha = "center", va = "center", weight='bold')
    ax.axis("off")

########

# Input/preprocessing
gene_functional_analysis_file = open(snakemake.input[0], "r")
lines = gene_functional_analysis_file.readlines()
gene_functional_analysis_file.close()
transcript_stats_df = pd.read_csv(snakemake.input[1])
transcript_stats_df["transcript_id"] = transcript_stats_df["transcript_id"].apply(lambda x: x.replace("_",""))

counts_data_filename = snakemake.input[2] if snakemake.config["quantification"] else None
annotation_filename = snakemake.input[3] if snakemake.config["annotation"] else None
args = Parser(counts_data_filename,annotation_filename,snakemake.config["samples"])
continuous = snakemake.config["continuous"]
minimumTPM = snakemake.config["minIsoTPM"]
maximumIso = snakemake.config["maxIsoNum"]
minimumPct = snakemake.config["minIsoPct"]
(annotation, data, samples, conditions, number_replicates) = preprocessArguments(args)

# Assign each transcript their features
transcripts,total_features = functionalFeatureParser(lines)

# Determine longest transcript
longest_length = 0
longest_tcon = ""
for ids in transcripts:
    if (len(transcripts[ids].aa) > longest_length):
        longest_length = len(transcripts[ids].aa)
        longest_tcon = ids

# Determine which transcripts to visualize and preprocess the quanitifcation data
total_gs = 1
transcripts_plot = sorted(list(transcripts.keys()))
if (snakemake.config["quantification"]):
    #df_temp = data[data["gene_name"] == gene]
    #print(transcripts_plot)
    #print(data["transcript_id"])
    data["transcript_id"] = data["transcript_id"].apply(lambda x: x.replace("_",""))
    df_temp = data[data["transcript_id"].isin(transcripts_plot)]
    data_gene = df_temp[samples].sum().to_frame().transpose()
    data_gene = calculateStatistics(data_gene,conditions,number_replicates)
    df_temp = calculateStatistics(df_temp,conditions,number_replicates)
    df_temp = (df_temp.filter(like='mean').div(data_gene.filter(like='mean').values[0],1)*100).add_prefix('Pct_').join(df_temp)
    df_temp = chooseIsoforms2Plot(df_temp,minimumTPM,minimumPct,maximumIso)
    df_temp["transcript_id"] = df_temp["transcript_id"].apply(lambda x: x.replace("_",""))
    transcripts_plot = list(pd.unique(df_temp["transcript_id"]))
    total_gs += 4

if (snakemake.config["annotation"]):
    transcripts_delete = []
    for transcript_id in transcripts_plot:
        transcript_annotation = annotation[annotation['transcript_id']==transcript_id]
        if ((transcript_annotation.shape[0] < 3)):
            transcripts_delete.append(transcript_id.replace("_",""))
    for transcript_id in transcripts_delete:
        transcripts_plot.remove(transcript_id)
    total_gs += len(transcripts_plot) + 1
total_gs += (total_features + 1) * len(transcripts_plot) + 2


# Visualization Code
fig = plt.figure(figsize = (21,total_gs))
gs = fig.add_gridspec(total_gs, 15)

fig.tight_layout(pad = 5.0)
fig.subplots_adjust(top = 0.95,bottom = 0.05,left = 0.05, right = 0.95)
fig.suptitle(gene,fontsize=30)

total_gs = 1
if (snakemake.config["quantification"]):
    #plot total expression
    axg = fig.add_subplot(gs[0:4, 0:4])
    plotTotal(conditions, data_gene, axg, colors, number_replicates, continuous)
    total_gs += 4

if (len(transcripts_plot) > 0):
    if (snakemake.config["quantification"]):
        df_temp = df_temp[df_temp["transcript_id"].isin(transcripts_plot)]

        #plot isoform expression
        axg = fig.add_subplot(gs[0:4, 5:10])
        plotIsoformProfiles(conditions, df_temp, axg, colors, continuous)

        #plot isoform expression percentage
        axg = fig.add_subplot(gs[0:4, 11:15])
        plotIsoformStacked(conditions, df_temp, axg, colors, continuous)

    if (snakemake.config["annotation"]):
        #prepare annotation
        annotation_cut = prepareAnnotation(annotation,transcripts_plot)

        #plot annotation
        plotAnnotation(annotation_cut, transcripts_plot, colors, longest_length, fig, gs, total_gs)
        total_gs += len(transcripts_plot) + 1

        plotSequenceAnalysis(transcripts_plot, colors, longest_length, fig, gs, total_gs, annotation_cut)
    else:
        plotSequenceAnalysis(transcripts_plot, colors, longest_length, fig, gs, total_gs)

    total_gs += (total_features + 1) * len(transcripts_plot) + 2
    plotLegend(plt,fig,gs,total_gs)

plt.subplots_adjust(wspace=0.05, hspace=0.05)
fig.savefig(snakemake.output[0])
