###### Raw ONT Reads Processing

# Basecalling
rule Basecalling:
    input:
        "RawData/{sample}/{id}.fast5"
    output:
        "RawData/QC/{sample}/{id}.txt",
        "RawData/Fast5/{sample}/{id}/{id}.fast5",
        fastq = "RawData/Basecalled/{sample}/{id}.fastq"
    resources:
        memory = 10,
        time = 5,
        tmpdir = 4
    params:
        guppy = config["guppy"],
        flowcell = config["flowcell"],
        kit = config["kit"],
        work = config["outdir"]
    priority: 10
    group: "basecalling"
    threads: 8
	shell: """
        mkdir -p {params.work}/RawData/Fast5/{wildcards.sample}/{wildcards.id}

		ln -sf {params.work}/{input} {params.work}/RawData/Fast5/{wildcards.sample}/{wildcards.id}/{wildcards.id}.fast5

		{params.guppy} -i {params.work}/RawData/Fast5/{wildcards.sample}/{wildcards.id}/ -s {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id}/ \
		--flowcell {params.flowcell} --kit {params.kit} --num_callers 1 --cpu_threads_per_caller {threads} --fast5_out

		mv -f {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id}/workspace/{wildcards.id}.fast5 {params.work}/RawData/Fast5/{wildcards.sample}/{wildcards.id}/{wildcards.id}.fast5

		mv {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id}/sequencing_summary.txt {params.work}/RawData/QC/{wildcards.sample}/{wildcards.id}.txt

		find {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id} -name '*.fastq' -exec mv {{}} {params.work}/RawData/Basecalled/{wildcards.sample}/{wildcards.id}.fastq \;"""

rule aggregate_fast5:
    input:
        f5 = lambda wildcards: expand("RawData/Basecalled/{{sample}}/{id}.fastq", id = glob_wildcards("RawData/"+wildcards.sample+"/{id}.fast5")[0]),
    output:
        fq = "RawData/Fastq/{sample}"+".fastq"
    resources:
        memory = 50,
        time = 4,
        tmpdir = 4
    priority: 9
    run:
        with open(output.fq,'w') as fout:
            for fn in input.f5:
                with open(fn,'r') as fin:
                    shutil.copyfileobj(fin, fout)

		#f5 = lambda wildcards: expand("RawData/Basecalled/{sample}/{id}.fastq"),
