---
title: 差异分析 RSEM
date: 2018-10-12
toc: true
description: 呦呦呦
categories: Bioinformatics
tags:
  - 软件和包
  - 差异分析
---

RSEM利用的是transcripts而非genome。我们有两种方式来构建RSEM转录参考，其一是利用参考基因组来构建；另外一种方式是从许多转录本中构建。

<!--more-->

## 1. 下载对应的基因组和注释文件
```
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
ftp://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/Mus_musculus.GRCm38.82.chr.gtf.gz
```
## 2. 构建索引

```bash
# 使用bowtie2索引
gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz
gunzip ref/Mus_musculus.GRCm38.82.chr.gtf.gz
# 从基因组构建
rsem-prepare-reference --gtf Mus_musculus.GRCm38.82.chr.gtf \
					   --bowtie2 --bowtie2-path bowtie2 \
					   Mus_musculus.GRCm38.dna.toplevel.fa mouse_ref

# 从转录本构建
rsem-prepare-reference --transcript-to-gene-map mouse_ref_mapping.txt \
					   --bowtie2 --bowtie2-path bowtie2 \
					   ref/mouse_ref.fa ref/mouse_ref
```

## 3. 表达定量

```bash
rsem-calculate-expression -p 8 --paired-end \
						--bowtie2 --bowtie2-path bowtie2 \
						--estimate-rspd \
						--append-names \
						--output-genome-bam \
						read1.fastq read2.fastq \
						ref/mouse_ref path/to/output/prefix
# output => prefix.genes.results, prefix.genome.bam, prefix.genome.sorted.bam,
#			prefix.genome.sorted.bam.bai, prefix.isoforms.results, prefix.stat,
#			prefix.transcript.bam, prefix.transcript.sorted.bam, prefix.transcript.sorted.bam.bai

# --output-genome-bam	只在由基因组构建索引时有效
# --bam 				当输入不是fastq文件而是bam文件时,指定该参数
# --calc-ci 			指定了会计算CQV值, coefficient of quartile variation，这将回答该基因的测序深度是否足够用于差异检测
# --no-bam-output 		不输出bam文件
# --single-cell-prior 	使用针对单细胞的预设条件运行
```

## 4. 结果输出

```r
# in R
data = read.table('prefix.gene.results', header=T, stringsAsFactors=F)
idx = order(data[,'TPM'], decreasing=T)
data[idx[1:10], c('gene_id', 'expected_count', 'TPM')]
```

```bash
# in bash
rsem-plot-model /path/to/output/prefix prefix.diagnostic.pdf
rsem-plot-transcript-wiggles --gene-list --show-unique \
							prefix gene_ids.txt aim_gene_transcript_wiggle.pdf
# gene_ids.txt: 含有gene identifier of aim_gene

rsem-bam2wig prefix.genome.sorted.bam prefix.wig prefix
```

## 5. 生成表达矩阵

```bash
rsem-generate-data-matrix sample1.genes.results sample2.genes.results \
						  sample3.genes.results ... \
						  sampleN.genes.results > gene_matrix.txt
rsem-run-ebseq gene_matrix.txt NumberOfgroup1,numberofgroup2 gene_matrix.ebseq.results
rsem-control-fdr gene_matrix.ebseq.result 0.05 gene_matrix.de.txt
```

## 6. 检测转录表达差异(differentially expressed isoforms)

```bash
rsem-generate-ngvector mouse_ref.transcripts.fa mouse_ref
# mouse_ref.transcripts.fa是mouse的转录组参考序列
# output=> mouse_ref.ngvec
rsem-generate-data-matrix sample1.isoforms.results \
						  sample2.isoforms.results sample3.isoforms.results \
						  ...
						  sampleN.isoforms.results > isoform_matrix.txt
rsem-run-ebseq --ngvector mouse_ref.ngvec isoform_matrix.txt control_sample_num,treate_sample_num \
						  isoform_ebseq.results
rsem-control-fdr isoform_ebseq.results 0.05 isoform.de.txt
```

## 7. 检测测序深度是否足够

```bash
# 1. 在rsem-calculate-expression中指定--calc-ci
# 2. 如果基因的cov值高于0.05表明测序深度不够,需要多少？
rsem-simulate-reads ../ref/mouse_ref prefix.stat/prefix.model \
					prefix.isoforms.results 0.36 20000000 prefix_sim_2M \
					--seed 0
# ../ref/mouse_ref指定参考文件的位置
# 0.36是预设的背景噪音比率
# output => prefix_sim_2M_1.fastq, prefix_sim_2M_2.fastq

rsem-calculate-expression -p 8 --paired-end \
						  --bowtie2 --bowtie2-path bowtie2 \
						  --estimate-rspd \
						  --append-names \
						  --no-bam-output \
						  --calc-ci \
						  prefix_sim_2M_1.fastq prefix_sim_2M_2.fastq \
						  ../ref/mouse_ref prefix_sim_2M
# 查看output的六个文件中prefix_sim_2m.genes.results,如果基因的cov低于0.05,表明20000000的测序深度足够.
```
