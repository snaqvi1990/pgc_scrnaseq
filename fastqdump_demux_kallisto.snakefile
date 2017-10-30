#########
# Pipeline for mapping scRNA-seq data from Li et al, Cell Stem Cell 2017.
# Uses: snakemake v3.4.2, kallisto v0.43, umis package https://github.com/vals/umis, fastq-dump 2.8.1
# SAMPLE_PATH variable is the path to a list of the SRA run ids from Li et al 2017
# KSTO_IDX is the path to a kallisto index for the human transcriptome
# WRKDIR is the path to a folder where all the temporary processing of kallisto-mapped pseudobam files will happen
# FQDIR is the path to a folder where raw fastq files will be downloaded
# if multithreading on an LSF cluster is not available, set all "threads" variables to 1
# The output is, for each run ID provided, a matrix of gene by cellular barcode counts (RUNID_CELLBARCODE).
# These matrices need to be combined prior to processing by Seurat and SCDE
#########

SAMPLE_PATH='/lab/solexa_page/snaqvi/pgc_scrnaseq/runids.all.txt'
KSTO_IDX = '/lab/page_human_data/gtex/index/hg38.gencode_v24.basic_ccds_nopar.transcripts.ksto'
WRKDIR = '/lab/solexa_page/snaqvi/pgc_scrnaseq/'
GENEMAP = '/lab/solexa_page/snaqvi/pgc_scrnaseq/gencode.v24.annotation.basic_ccds_nopar.genesymbol_tx.txt'
FQDIR = '/lab/page_downloads/li2017/'

def get_inputs():
    f = open(SAMPLE_PATH, "r")

    inputs = []
    for line in f:
        sample_id = line.rstrip()
        if (sample_id) == '' : continue
        inputs.append(sample_id)

    f.close()
    #print(inputs)
    return inputs

rule all:
    input: expand((WRKDIR+"{sample_id}.gene.cb.counts.txt"),sample_id=get_inputs())

rule fastqdump:
    input: SAMPLE_PATH
    params: sampid = "{sample_id}", fqdir = FQDIR
    output: fwd = FQDIR+"{sample_id}_1.fastq.gz", rev = FQDIR+"{sample_id}_2.fastq.gz"
    shell:
        """
        cd {params.fqdir}
        fastq-dump --split-3 -gzip {params.sampid}
        """

rule fastqtransform:
    input: fwd = FQDIR+"{sample_id}_1.fastq.gz", rev = FQDIR+"{sample_id}_2.fastq.gz"
    params: transf = WRKDIR+"transform.json"
    output: fwd = temp(WRKDIR+"{sample_id}_1.transform.fastq.gz"), rev = temp(WRKDIR+"{sample_id}_2.transform.fastq.gz")
    threads: 8
    shell:
        """umis fastqtransform --cores {threads} --fastq1out {output.fwd} --fastq2out {output.rev} {params.transf} {input.fwd} {input.rev}"""

rule cb_filterfwd:
    input: fwd = WRKDIR+"{sample_id}_1.transform.fastq.gz"
    params: cb = WRKDIR+"cell_barcode_only.txt"
    output: temp(WRKDIR+"{sample_id}_1.transform.filtered.fastq.gz")
    threads: 8
    shell:
        """umis cb_filter --bc1 {params.cb} --cores {threads} --nedit 0 {input.fwd} | gzip -f > {output}"""

rule cb_filterrev:
    input: rev = WRKDIR+"{sample_id}_2.transform.fastq.gz"
    params: cb = WRKDIR+"cell_barcode_only.txt"
    output: temp(WRKDIR+"{sample_id}_2.transform.filtered.fastq.gz")
    threads: 8
    shell:
        """umis cb_filter --bc1 {params.cb} --cores {threads} --nedit 0 {input.rev} | gzip -f > {output}"""

rule kallisto:
    input: fwd = WRKDIR+"{sample_id}_1.transform.filtered.fastq.gz",rev = WRKDIR+"{sample_id}_2.transform.filtered.fastq.gz"
    params: idx = KSTO_IDX
    output: temp(WRKDIR+"{sample_id}.transform.filtered.bam")
    shell:
        """kallisto quant -i {params.idx} -o temp --pseudobam {input.fwd} {input.rev} | samtools view -F 16 -Sb - > {output}"""

rule tagcount:
    input: bam = WRKDIR+"{sample_id}.transform.filtered.bam"
    params: genemap = GENEMAP
    output: WRKDIR+"{sample_id}.gene.cb.counts.txt"
    shell:
        """umis tagcount --genemap {params.genemap} {input.bam} {output}"""
