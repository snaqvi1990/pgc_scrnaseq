# pgc_scrnaseq
Code for analyzing single-cell RNA-seq data from Li et al, Cell Stem Cell 2017

fastqdump_demux_kallisto.snakefile is a snakemake file that downloads and pseudomaps raw data from Li et al, Cell Stem Cell 2017. It writes a gene x cell matrix file for each run (corresponding to one dissected/sorted embryo) separately, however; a combine matrix is needed for further downstream analyses. It is provided at: http://pagelab.wi.mit.edu/page/papers/Nicholls_et_al_2017/alldat.noembryo.countmat.txt


