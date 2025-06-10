{\rtf1\ansi\ansicpg1252\cocoartf2639
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;\f1\fnil\fcharset0 Menlo-Bold;}
{\colortbl;\red255\green255\blue255;\red46\green174\blue187;\red0\green0\blue0;\red47\green180\blue29;
\red180\green36\blue25;\red200\green20\blue201;\red159\green160\blue28;\red64\green11\blue217;}
{\*\expandedcolortbl;;\cssrgb\c20199\c73241\c78251;\csgray\c0;\cssrgb\c20241\c73898\c14950;
\cssrgb\c76411\c21697\c12527;\cssrgb\c83397\c23074\c82666;\cssrgb\c68469\c68012\c14211;\cssrgb\c32309\c18666\c88229;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f0\fs22 \cf2 \CocoaLigature0 #!/bin/bash\cf3 \
sample\cf4 =
\f1\b \cf5 $1
\f0\b0 \cf3 \
dir\cf4 =
\f1\b \cf5 $2
\f0\b0 \cf3 \
\
\
cutadapt
\f1\b \cf6  -a
\f0\b0 \cf3  ACTGTCTCTTATACACATCT
\f1\b \cf6  -q
\f0\b0 \cf3  20
\f1\b \cf6  -m
\f0\b0 \cf3  50
\f1\b \cf6  -o
\f0\b0 \cf3  
\f1\b \cf5 $dir
\f0\b0 \cf3 /cutadapt_RNAseq/
\f1\b \cf7 "$sample"
\f0\b0 \cf3 _trimmed_R1.fastq.gz
\f1\b \cf6  -p
\f0\b0 \cf3  
\f1\b \cf5 $dir
\f0\b0 \cf3 /cutadapt_RNAseq/
\f1\b \cf7 "$sample"
\f0\b0 \cf3 _trimmed_R2.fastq.gz \cf4 \\\cf3 \
        
\f1\b \cf5 $dir
\f0\b0 \cf3 /RNAseq_fastq/
\f1\b \cf7 "$sample"
\f0\b0 \cf3 _R1.fastq.gz 
\f1\b \cf5 $dir
\f0\b0 \cf3 /RNAseq_fastq/
\f1\b \cf7 "$sample"
\f0\b0 \cf3 _R2.fastq.gz \cf4 >\cf3  
\f1\b \cf5 $dir
\f0\b0 \cf3 /cutadapt_RNAseq/
\f1\b \cf5 $sample
\f0\b0 \cf3 .qc.txt \cf4 \\\cf3 \
        2\cf4 >\cf3  ~/logs2/cutadapt/
\f1\b \cf5 $sample
\f0\b0 \cf3 .log\
\
STAR
\f1\b \cf6  --sjdbOverhang
\f0\b0 \cf3  100
\f1\b \cf6  --twopassMode
\f0\b0 \cf3  Basic
\f1\b \cf6  --outFilterType
\f0\b0 \cf3  BySJout
\f1\b \cf6  --outFilterMultimapNmax
\f0\b0 \cf3  20
\f1\b \cf6  --alignSJoverhangMin
\f0\b0 \cf3  8
\f1\b \cf6  --alignSJDBoverhangMin
\f0\b0 \cf3  1 \cf4 \\\cf3 \

\f1\b \cf6         --outFilterMismatchNmax
\f0\b0 \cf3  999
\f1\b \cf6  --outFilterMismatchNoverReadLmax
\f0\b0 \cf3  0.04
\f1\b \cf6  --alignIntronMin
\f0\b0 \cf3  20
\f1\b \cf6  --alignIntronMax
\f0\b0 \cf3  1000000
\f1\b \cf6  --alignMatesGapMax
\f0\b0 \cf3  1000000
\f1\b \cf6  --outSAMunmapped
\f0\b0 \cf3  Within \cf4 \\\cf3 \

\f1\b \cf6         --outSAMattributes
\f0\b0 \cf3  NH HI AS NM MD
\f1\b \cf6  --outSAMtype
\f0\b0 \cf3  BAM SortedByCoordinate
\f1\b \cf6  --quantMode
\f0\b0 \cf3  GeneCounts TranscriptomeSAM
\f1\b \cf6  --limitBAMsortRAM
\f0\b0 \cf3  20000000000
\f1\b \cf6  --sjdbGTFfile
\f0\b0 \cf3  ~/resources/Homo_sapiens.GRCh38.83.gtf \cf4 \\\cf3 \

\f1\b \cf6         --runThreadN
\f0\b0 \cf3  16
\f1\b \cf6  --genomeDir
\f0\b0 \cf3  ~/resources/star
\f1\b \cf6  --readFilesIn
\f0\b0 \cf3  
\f1\b \cf5 $dir
\f0\b0 \cf3 /cutadapt_RNAseq/
\f1\b \cf7 "$sample"
\f0\b0 \cf3 _trimmed_R1.fastq.gz 
\f1\b \cf5 $dir
\f0\b0 \cf3 /cutadapt_RNAseq/
\f1\b \cf7 "$sample"
\f0\b0 \cf3 _trimmed_R2.fastq.gz
\f1\b \cf6  --readFilesCommand
\f0\b0 \cf3  zcat \cf4 \\\cf3 \

\f1\b \cf6         --outFileNamePrefix
\f0\b0 \cf3  
\f1\b \cf5 $dir
\f0\b0 \cf3 /star_RNAseq/
\f1\b \cf5 $sample
\f0\b0 \cf3 /
\f1\b \cf6  --outStd
\f0\b0 \cf3  Log \cf4 >\cf3  ~/logs2/star/
\f1\b \cf5 $sample
\f0\b0 \cf3 .log\
\

\f1\b \cf8 mkdir
\f0\b0 \cf3  /
\f1\b \cf5 $dir
\f0\b0 \cf3 /rsem_RNAseq/
\f1\b \cf5 $sample
\f0\b0 \cf3 \
\
rsem-calculate-expression
\f1\b \cf6  --num-threads
\f0\b0 \cf3  8
\f1\b \cf6  --no-bam-output --paired-end --bam
\f0\b0 \cf3  
\f1\b \cf5 $dir
\f0\b0 \cf3 /star_RNAseq/
\f1\b \cf5 $sample
\f0\b0 \cf3 /Aligned.toTranscriptome.out.bam /home/bow012/resources/rsem/human_ensembl \cf4 \\\cf3 \
        
\f1\b \cf5 $dir
\f0\b0 \cf3 /rsem_RNAseq/
\f1\b \cf5 $sample
\f0\b0 \cf3 /
\f1\b \cf5 $sample
\f0\b0 \cf3 .Quant.rsem \cf4 >\cf3  /home/bow012/logs2/rsem/
\f1\b \cf5 $sample
\f0\b0 \cf3 .log}