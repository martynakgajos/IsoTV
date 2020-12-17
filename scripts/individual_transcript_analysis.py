#!/usr/bin/env python3

# Individual transcript analysis

count = 0
output = [">" + snakemake.wildcards["transcript"]]

for filename in snakemake.input:
    file = open(filename,"r")
    lines = file.readlines()
    file.close()

    if("ss" in filename):
        output.append("#####SS Analysis")
        for position in lines:
            output.append(position.strip())

    elif ("iupred2a" in filename):
        output.append("#####IUPred2A Analysis")
        lines = lines[7:]
        for position in lines:
            output.append(position.strip())

    elif("domains" in filename):
        output.append("#####Domain Analysis")
        lines = lines[5:]
        for region in lines:
            if (not region.startswith("#")):
                region_split = region.strip().split("\t")
                if (region_split[1] == "Pfam"):
                    region_short = "\t".join([region_split[3],region_split[4],region_split[8].split(";")[3].split("=")[-1]])
                    output.append(region_short.strip())
            else:
                break

    elif("prion" in filename):
        output.append("#####Prion Analysis")
        for prionsite in lines:
            if (not prionsite.startswith("#")):
                prionsite_split = prionsite.strip().split("\t")
                prionsite_short = "\t".join([prionsite_split[2],prionsite_split[15]])
                output.append(prionsite_short)

    elif("sites" in filename):
        output.append("#####Prosite Analysis")
        modification = ""
        site_only = False
        for site in lines:
            if (site.startswith(">")):
                modification = site.strip().split(" ")[3]
                site_only = "site" in site
            else:
                site_split = site.strip().replace("  "," ").split(" ")
                if ((len(site_split) > 2) and site_only):
                    site_short = "\t".join([site_split[0],site_split[2],modification])
                    output.append(site_short)

    elif ("protein" in filename):
        output.append("#####AA Sequence")
        lines = lines[1:]
        output.append(lines[0].strip())

output.append("")

output_file = open(snakemake.output[0],"w+")
output_file.write("\n".join(output))
output_file.close()
