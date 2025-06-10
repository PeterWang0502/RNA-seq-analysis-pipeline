# RNA-seq-analysis-pipeline

This is a bioinformatics pipeline used for analyzing RNA sequencing (RNA-seq) data. The pipeline performs quality controls (QC), trimming, alignment, and quantifying transcript counts. Following differential gene expression analysis are provided in R scripts.
1. Read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([Cutadapt](https://cutadapt.readthedocs.io/en/stable/))
3. Alignment ([STAR](https://github.com/alexdobin/STAR))
4. Transcript quantification ([RSEM](https://github.com/deweylab/RSEM))
5. Extensive QC ([Qualimap](http://qualimap.conesalab.org/))
6. Differential expression analysis:
     1. [DeSeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
     2. [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

# Pipeline setup and run

The pipeline is run with Snakemake and Singularity. Setting up a virtual environment to run this pipeline with mamba is recommended.
```
mamba env create --name snakemake-tutorial --file environment.yaml
```
Singularity containers can be downloaded from here: https://depot.galaxyproject.org/singularity/f26fb5d65f9c

The Snakemake pipeline can also be run without singularity if it's not supported in the cluster. In that case, you will need to install all necessary softwares (FastQC, Cutadapt, STAR, etc) and compile them manually. Remove singularity parts from each rule in `Snakefile` and add software paths to `config.yaml` and cooresponding labels to `Snakefile`.

# DE analysis

DeSeq2 and edgeR can be applied to output of the pipeline, and differentially expressed genes are nominated. Genes with |log2FC| > 1.5 and adjusted p-value (FDR) < 0.05 are selected as differentially expressed genes.
