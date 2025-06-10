#!/bin/bash
echo $1
echo $2

sample=$1
dir=$2

### use specific adapter sequence for Cutadapt (Illumina usually uses ACTGTCTCTTATACACATCT), or you can use a universal adapter sequence AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

cutadapt -a ACTGTCTCTTATACACATCT -q 20 -m 50 -o $dir/cutadapt_RNAseq/"$sample"_trimmed_R1.fastq.gz -p $dir/cutadapt_RNAseq/"$sample"_trimmed_R2.fastq.gz \
	$dir/RNAseq__fastq/"$sample"_R1.fastq.gz $dir/RNAseq__fastq/"$sample"_R2.fastq.gz > $dir/cutadapt_RNAseq/$sample.qc.txt \
	2> ~/logs/cutadapt/$sample.log

### STAR requires generating genome indices before use. similarly, RSEM requires preparing reference with rsem-prepare-reference (details in RSEM GitHub: https://github.com/deweylab/RSEM)

#STAR --runThreadN 8 --runMode genomeGenerate --genomeFastaFiles ~/resources/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbOverhang 100 \
#       --sjdbGTFfile ~/resources/Homo_sapiens.GRCh38.83.gtf --genomeDir ~/resources/star

STAR --sjdbOverhang 100 --twopassMode Basic --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMunmapped Within \
	--outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts TranscriptomeSAM --limitBAMsortRAM 20000000000 --sjdbGTFfile ~/resources/Homo_sapiens.GRCh38.83.gtf \
	--runThreadN 16 --genomeDir ~/resources/star --readFilesIn $dir/cutadapt_RNAseq/"$sample"_trimmed_R1.fastq.gz $dir/cutadapt_RNAseq/"$sample"_trimmed_R2.fastq.gz --readFilesCommand zcat \
	--outFileNamePrefix $dir/star_RNAseq/$sample/ --outStd Log > ~/logs/star/$sample.log

### if need to merge bam files for technical replicates

#samtools merge $dir/star_RNAseq/"$sample"_1+2/"$sample"_merged.out.bam $dir/star_RNAseq/"$sample"_1/Aligned.toTranscriptome.out.bam $dir/star_RNAseq/"$sample"_2/Aligned.toTranscriptome.out.bam

mkdir /$dir/rsem_RNAseq/$sample

rsem-calculate-expression --num-threads 8 --no-bam-output --paired-end --bam $dir/star_RNAseq/$sample/Aligned.toTranscriptome.out.bam ~/resources/rsem/human_ensembl \
	$dir/rsem_RNAseq/$sample/$sample.Quant.rsem > ~/logs/rsem/$sample.log
