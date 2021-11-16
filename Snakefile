import datetime
import time
import os
import os.path
import json

curr_date = datetime.datetime.today().strftime("%y%m%d")


REF_DIR = "refs/NC_045512.2"
REF_SEQ = f"{REF_DIR}/reference.fasta"
INPUT_DIR = "input"
GISAID_DIR = f"{INPUT_DIR}/gisaid"

pipelines_free = {
    "A01": "seqs/{dir}/msa/kl-ma_0.01/snp-sites/catalog.cleaned.vcf.gz",
    "B01": "seqs/{dir}/msa/nkl-ma_0.01/msa2vcf/catalog.norm.filt_maxlen.vcf.gz",
}

rule all:
    input:
        expand(
            "vcf/{curr_date}-SARS-CoV-2-{source}-pipeline_{pip}/{file}",
            curr_date = [curr_date],
            source = ["NCBIVirusS-filters_50_100"],
            pip = pipelines_free.keys(),
            file = [
                "catalog.vcf.gz",
                "reference.fasta",
                "GCF_009858895.2_ASM985889v3_genomic.gff",
                "info.json",
                "status.json"
            ]
        ),


rule copy_from_ref:
    input:
        "refs/NC_045512.2/{file}"
    output:
        "vcf/{vcfdir}/{file}"
    shell: "cp -a {input} {output}"

rule download_ncbi_virus_seqs:
    output:
        archive = temp("seqs/{curr_date}-SARS-CoV-2-NCBIVirus/archive.zip")
    log:
        out = "seqs/{curr_date}-SARS-CoV-2-NCBIVirus/dataset.log",
        out_bz2 = "seqs/{curr_date}-SARS-CoV-2-NCBIVirus/dataset.log.bz2",
    shadow: "shallow"
    shell: """
    rm -f archive.zip &> {log.out}
    ./software/datasets download virus genome taxon SARS2 --host human --complete-only --exclude-cds --exclude-gpff --exclude-pdb --exclude-protein --filename archive.zip &>> {log.out}
    mv archive.zip {output.archive}
    bzip2 -9k {log.out}
    """

rule prepare_ncbi_virus_seqs:
    input:
        archive = "seqs/{curr_date}-SARS-CoV-2-NCBIVirus/archive.zip"
    output:
        fa = "seqs/{curr_date}-SARS-CoV-2-NCBIVirus/sequences.fa.gz",
        ids = "seqs/{curr_date}-SARS-CoV-2-NCBIVirus/accessions.txt.bz2",
    log:
        out = "seqs/{curr_date}-SARS-CoV-2-NCBIVirus/dataset-pre.log",
    shadow: "shallow"
    conda: "envs/pigz.yml"
    threads: 8
    shell: """
    rm -rf ncbi_dataset/ &> {log.out}
    unzip {input.archive} ncbi_dataset/data/genomic.fna &>> {log.out}
    pigz -9k -p {threads} ncbi_dataset/data/genomic.fna
    mv ncbi_dataset/data/genomic.fna.gz {output.fa}
    grep '^>' ncbi_dataset/data/genomic.fna | cut -d ' ' -f 1 | sed 's/>//' | bzip2 -9 > {output.ids}
    """

rule pangolin_ncbi_virus_seqs:
    input:
        seqs = f"seqs/{curr_date}-SARS-CoV-2-NCBIVirus/sequences.fa",
    output:
        lineages = f"seqs/{curr_date}-SARS-CoV-2-NCBIVirus/lineage_report.csv",
        lineages_bz2 = f"seqs/{curr_date}-SARS-CoV-2-NCBIVirus/lineage_report.csv.bz2",
    params:
        outdir = f"seqs/{curr_date}-SARS-CoV-2-NCBIVirus",
        min_len = 29000,
        max_ambig = 0.01,
    log:
        out = f"seqs/{curr_date}-SARS-CoV-2-NCBIVirus/pangolin.log",
        out_bz2 = f"seqs/{curr_date}-SARS-CoV-2-NCBIVirus/pangolin.log.bz2",
    conda: "envs/pangolin.yml"
    threads: 8
    shell: """
    pangolin --update &> {log.out}
    pangolin -v &>> {log.out}
    pangolin -pv &>> {log.out}
    pangolin --no-temp --alignment --outdir {params.outdir} --outfile lineage_report.csv --threads {threads} --max-ambig {params.max_ambig} --min-length {params.min_len} {input.seqs} &>> {log.out}
    bzip2 -9k {output.lineages} {log.out}
    """

rule subsample_ncbi_virus_seqs:
    input:
        seqs = "seqs/{curr_date}-SARS-CoV-2-NCBIVirus/sequences.fa.gz",
        lineages = "seqs/{curr_date}-SARS-CoV-2-NCBIVirus/lineage_report.csv",
    output:
        report = "seqs/{curr_date}-SARS-CoV-2-NCBIVirusS-filters_{minout}_{maxout}/selected_lineage_report.csv",
        selected = "seqs/{curr_date}-SARS-CoV-2-NCBIVirusS-filters_{minout}_{maxout}/selected_seqs.txt",
        report_bz2 = "seqs/{curr_date}-SARS-CoV-2-NCBIVirusS-filters_{minout}_{maxout}/selected_lineage_report.csv.bz2",
        selected_bz2 = "seqs/{curr_date}-SARS-CoV-2-NCBIVirusS-filters_{minout}_{maxout}/selected_seqs.txt.bz2",
        seqs = "seqs/{curr_date}-SARS-CoV-2-NCBIVirusS-filters_{minout}_{maxout}/sequences.fa.gz",
    params:
        seqs = "seqs/{curr_date}-SARS-CoV-2-NCBIVirusS-filters_{minout}_{maxout}/sequences.fa",
    log:
        out = "seqs/{curr_date}-SARS-CoV-2-NCBIVirusS-filters_{minout}_{maxout}/subsample.log",
    conda:
        "envs/tidyverse.yml"
    shell: """
    ./software/subsample-free-seqs.R {input.lineages} {output.report} {output.selected} {wildcards.minout} {wildcards.maxout} &> {log.out}
    ./software/faSomeRecords <( zcat {input.seqs} | tr " ," "__" ) {output.selected} {params.seqs} &>> {log.out}
    gzip -9 {params.seqs}
    bzip2 -9k {output.report} {output.selected}
    """

rule index_fasta:
    input:  "{dir}/{file}.fasta"
    output: "{dir}/{file}.fasta.fai"
    conda: "envs/samtools.yml"
    shell: "samtools faidx {input}"

rule gunzip_sequences:
    input:
        "{dir}/sequences.fa.gz",
    output:
        "{dir}/sequences.fa",
    shell:
        "gunzip -k {input}"

rule gunzip_msa:
    input:
        "{dir}/output.msa.gz",
    output:
        temp("{dir}/output.msa"),
    shell:
        "gunzip -k {input}"

rule gunzip_vcf:
    input:
        "{dir}/{catalog}.vcf.gz",
    output:
        temp("{dir}/{catalog}.vcf"),
    shell:
        "gunzip -k {input}"

rule mafft_keeplen:
    input:
        seqs = "seqs/{dir}/sequences.fa",
        ref = REF_SEQ,
    output:
        msa = "seqs/{dir}/msa/kl-ma_{ma}/output.msa.gz",
    log:
        out = "seqs/{dir}/msa/kl-ma_{ma}/mafft.log",
        time = "seqs/{dir}/msa/kl-ma_{ma}/mafft.time",
    threads:
        8
    conda:
        "envs/mafft.yml"
    shell:
        "/usr/bin/time -o {log.time} -v mafft --thread {threads} --6merpair --maxambiguous {wildcards.ma} --keeplength --addfragments {input.seqs} {input.ref} 2> {log.out} | gzip -9 > {output.msa}"

rule mafft_nokeeplen:
    input:
        seqs = "seqs/{dir}/sequences.fa",
        ref = REF_SEQ,
    output:
        msa = "seqs/{dir}/msa/nkl-ma_{ma}/output.msa.gz",
    log:
        out = "seqs/{dir}/msa/nkl-ma_{ma}/mafft.log",
        time = "seqs/{dir}/msa/nkl-ma_{ma}/mafft.time",
    threads:
        8
    conda:
        "envs/mafft.yml"
    shell:
        "/usr/bin/time -o {log.time} -v mafft --thread {threads} --6merpair --maxambiguous {wildcards.ma} --addfragments {input.seqs} {input.ref} 2> {log.out} | gzip -9 > {output.msa}"

rule fill_msa:
    input:
        fill_msa = "software/fill_msa",
        msa = "seqs/{dir}/output.msa.gz",
    output:
        filled = "seqs/{dir}/output.msa.filled.gz",
    log:
        out = "seqs/{dir}/fill_msa.log",
        time = "seqs/{dir}/fill_msa.time",
    threads:
        2 # 1 for the program, 1 for gzip
    shell:
        "/usr/bin/time -o {log.time} -v ./{input.fill_msa} {input.msa} 2> {log.out} | gzip -9 > {output.filled}"

rule snpsites:
    input:
        filled = "seqs/{dir}/output.msa.filled.gz",
    output:
        vcf = "seqs/{dir}/snp-sites/catalog.vcf.gz",
    log:
        out = "seqs/{dir}/snp-sites/snp-sites.log",
        time = "seqs/{dir}/snp-sites/snp-sites.time",
    threads:
        2 # 1 for the program, 1 for gzip
    conda:
        "envs/snpsites.yml"
    shell:
        "/usr/bin/time -o {log.time} -v snp-sites -cv {input.filled} 2> {log.out} | bgzip > {output.vcf}"

rule clean_vcf:
    input:
        vcf = "seqs/{dir}/snp-sites/catalog.vcf",
    output:
        vcf = "seqs/{dir}/snp-sites/catalog.cleaned.vcf.gz",
    log:
        out = "seqs/{dir}/snp-sites/clean.log",
        time = "seqs/{dir}/snp-sites/clean.time",
    threads:
        2 # 1 for the program, 1 for gzip
    conda:
        "envs/pysam.yml"
    shell:
        "/usr/bin/time -o {log.time} -v ./software/format_vcf.py clean {input.vcf} 2> {log.out} | cut -f 1-9,11- | bgzip > {output.vcf}"


def translate_pipeline(wildcards):
    return {
        "refidx": REF_SEQ + ".fai",
        "vcf": (pipelines_free[wildcards.pip]).format(**wildcards),
    }

rule filter_vcf:
    input:
        unpack(translate_pipeline)
    output:
        vcf = "vcf/{dir}-pipeline_{pip}-filters_minsample_{minsample}/catalog.vcf.gz",
    log:
        out = "vcf/{dir}-pipeline_{pip}-filters_minsample_{minsample}/filter.log",
        time = "vcf/{dir}-pipeline_{pip}-filters_minsample_{minsample}/filter.time",
    threads:
        2 # 1 for the program, 1 for gzip
    conda:
        "envs/pysam.yml"
    shell:
        "/usr/bin/time -o {log.time} -v ./software/format_vcf.py freq {input.vcf} {input.refidx} {wildcards.minsample} 2> {log.out} | bgzip > {output.vcf}"

rule no_filter_vcf:
    input:
        unpack(translate_pipeline)
    output:
        vcf = "vcf/{dir}-pipeline_{pip}/catalog.vcf.gz",
    log:
        out = "vcf/{dir}-pipeline_{pip}/filter.log",
        time = "vcf/{dir}-pipeline_{pip}/filter.time",
    threads:
        2 # 1 for the program, 1 for gzip
    conda:
        "envs/pysam.yml"
    shell:
        "/usr/bin/time -o {log.time} -v ./software/format_vcf.py freq {input.vcf} {input.refidx} 0 2> {log.out} | bgzip > {output.vcf}"

rule no_filter_vcf_info_json:
    input:
        unpack(translate_pipeline)
    output:
        json = "vcf/{dir}-pipeline_{pip}/info.json",
    run:
        vcfdir = f"{wildcards.dir}-pipeline_{wildcards.pip}"
        info = {
	        "filename": f"/jobs/vcf/{vcfdir}/catalog.vcf.gz",
	        "id": vcfdir,
	        "description": f"SARS-CoV-2 precomputed VCF and ref downloaded from NCBI on {curr_date}.",
	        "submission_time": int(time.time()),
	        "alias": vcfdir,
	        "reference": f"/jobs/vcf/{vcfdir}/reference.fasta",
	        "gtf": f"/jobs/vcf/{vcfdir}/GCF_009858895.2_ASM985889v3_genomic.gff",
	        "internal_ref": {
		        "snpEff": {
			        "id": "NC_045512.2"
		        },
		        "reference": {
			        "file": "NC_045512.2/reference.fasta"
		        },
		        "id": "NC_045512.2",
		        "annotation": {
			        "file": "NC_045512.2/GCF_009858895.2_ASM985889v3_genomic.gff"
		        },
		        "alias": "SARS-CoV-2, Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1",
		        "pangolin": True
	        }
        }
        with open(output.json, "w") as f:
            json.dump(info, f)

rule status_json:
    input:
        unpack(translate_pipeline)
    output:
        json = "vcf/{dir}-pipeline_{pip}/status.json",
    run:
        vcfdir = f"{wildcards.dir}-pipeline_{wildcards.pip}"
        status = {
	        "status": "Precomputed",
	        "last_time": datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S"),
	        "output": {
		        "vcf": f"/jobs/vcf/{vcfdir}/catalog.vcf.gz",
		        "reference": f"/jobs/vcf/{vcfdir}/reference.fasta"
	        }
        }
        with open(output.json, "w") as f:
            json.dump(status, f)

rule msa2vcf:
    input:
        msa = "seqs/{dir}/output.msa",
        ref = REF_SEQ,
    output:
        file_list = "seqs/{dir}/msa2vcf/base-files.list",
        vcf_dir = directory("seqs/{dir}/msa2vcf/base_vcfs"),
    log:
        out = "seqs/{dir}/msa2vcf/msa2vcf.log",
        time = "seqs/{dir}/msa2vcf/msa2vcf.time",
    threads:
        16 # 1 for the program, 1 for gzip
    conda:
        "envs/msa2vcf.yml"
    shell: """
    rm -rf {output.vcf_dir}
    mkdir -p {output.vcf_dir}
    /usr/bin/time -o {log.time} -v software/msa2vcf.py {input.ref} {input.msa} {output.vcf_dir} {threads} 2> {log.out}
    ls {output.vcf_dir} | while read f; do echo "{output.vcf_dir}/$f"; done > {output.file_list}
    """

rule index_vcf:
    input:
        "seqs/{dir}/msa2vcf/base-files.list",
    output:
        file_list = "seqs/{dir}/msa2vcf/indexed-files.list",
        vcf_dir = directory("seqs/{dir}/msa2vcf/indexed_vcfs"),
    conda:
        "envs/bcftools.yml"
    shell: """
    rm -rf {output.vcf_dir}
    mkdir -p {output.vcf_dir}
    while read f; do
        bgzip -c $f > {output.vcf_dir}/$(basename $f)
    done < {input}
    ls {output.vcf_dir} | while read f; do echo "{output.vcf_dir}/$f"; done > {output.file_list}
    while read f; do
        bcftools index $f
    done < {output.file_list}
    """

checkpoint compute_chunks:
    input:
        "seqs/{dir}/msa2vcf/indexed-files.list",
    output:
        chunks = directory("seqs/{dir}/msa2vcf/chunks"),
    shell: """
    rm -rf {output.chunks}
    mkdir -p {output.chunks}
    split -l 1000 -a 4 --additional-suffix .txt {input} {output.chunks}/chunk-
    """

rule merge_chunk_vcf:
    input:
        chunk="seqs/{dir}/msa2vcf/chunks/chunk-{idx}.txt",
    output:
        vcf="seqs/{dir}/msa2vcf/chunk_vcfs/chunk-{idx}.vcf.gz",
        idx="seqs/{dir}/msa2vcf/chunk_vcfs/chunk-{idx}.vcf.gz.csi",
    log:
        out="seqs/{dir}/msa2vcf/merge-chunk-{idx}.out",
        time="seqs/{dir}/msa2vcf/merge-chunk-{idx}.time",
    conda:
        "envs/bcftools.yml"
    shell: """
        /usr/bin/time -o {log.time} -v bcftools merge --threads {threads} -m all -0 -l {input} -O z -o {output.vcf} 2> {log.out}
        bcftools index {output.vcf} 2>> {log.out}
    """

rule bound_vcf:
    input:
        vcf="{dir}/catalog{suffix}.vcf.gz",
    output:
        vcf="{dir}/catalog{suffix}.bound.vcf.gz",
        isx="{dir}/catalog{suffix}.bound.vcf.gz.csi",
    log:
        out="{dir}/bound.{suffix}.out",
    conda:
        "envs/bcftools.yml"
    shell: """
        bcftools filter -i "POS>260 & POS<29680" {input.vcf} -Oz -o {output.vcf} 2> {log.out}
        bcftools index {output.vcf} 2>> {log.out}
    """

def merge_all_input(wildcards):
    msa2vcf_out = checkpoints.compute_chunks.get(**wildcards).output.chunks
    return expand(
        "seqs/{dir}/msa2vcf/chunk_vcfs/chunk-{sid}.vcf.gz",
        sid=glob_wildcards(os.path.join(msa2vcf_out, "chunk-{sid}.txt")).sid,
        **wildcards
    )

rule merge_all:
    input:
        merge_all_input,
    output:
        vcf="seqs/{dir}/msa2vcf/catalog.vcf.gz",
        idx="seqs/{dir}/msa2vcf/catalog.vcf.gz.csi",
    log:
        out="seqs/{dir}/msa2vcf/merge-all.out",
        time="seqs/{dir}/msa2vcf/merge-all.time",
    conda:
        "envs/bcftools.yml"
    threads: 4
    shell: """
        /usr/bin/time -o {log.time} -v bcftools merge --threads {threads} -m all -0 -O z -o {output.vcf} {input} 2> {log.out}
        bcftools index {output.vcf} 2>> {log.out}
    """

rule norm_merge_any:
    input:
        vcf="seqs/{dir}/msa2vcf/catalog.vcf.gz",
        idx="seqs/{dir}/msa2vcf/catalog.vcf.gz.csi",
        ref = REF_SEQ,
    output:
        vcf="seqs/{dir}/msa2vcf/catalog.norm.vcf.gz",
        idx="seqs/{dir}/msa2vcf/catalog.norm.vcf.gz.csi",
    log:
        out="seqs/{dir}/msa2vcf/norm.out",
        time="seqs/{dir}/msa2vcf/norm.time",
    conda:
        "envs/bcftools.yml"
    threads: 2
    shell: """
        /usr/bin/time -o {log.time} -v bcftools norm --threads {threads} -c w -d all -f {input.ref} -m +any -O z -o {output.vcf} {input.vcf} 2> {log.out}
        bcftools index {output.vcf} 2>> {log.out}
    """

rule filt_maxlen:
    input:
        vcf = "seqs/{dir}/msa2vcf/catalog.norm.vcf.gz",
        idx = "seqs/{dir}/msa2vcf/catalog.norm.vcf.gz.csi",
    output:
        vcf = "seqs/{dir}/msa2vcf/catalog.norm.filt_maxlen.vcf.gz",
        idx = "seqs/{dir}/msa2vcf/catalog.norm.filt_maxlen.vcf.gz.csi",
    log:
        out = "seqs/{dir}/msa2vcf/norm.filt_maxlen.out",
        time = "seqs/{dir}/msa2vcf/norm.filt_maxlen.time",
    params:
        maxlen = 10,
    conda:
        "envs/bcftools.yml"
    threads: 2
    shell: """
        /usr/bin/time -o {log.time} -v bcftools view --threads $(({threads} -1)) -a -e "STRLEN(REF)>{params.maxlen} || MAX(STRLEN(ALT[*]))>{params.maxlen}" -o {output.vcf} -O z {input.vcf} 2> {log.out}
        bcftools index {output.vcf} 2>> {log.out}
    """
