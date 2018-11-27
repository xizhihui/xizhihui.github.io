---
title: 比对 STAR
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 软件和包
---

STAR 的比对速率要比 bowtie 快那么一丢丢。

<!--more-->

## 1. 生成索引

```bash
STAR --runThreadN numberOfThreads \
	--runMode genomeGenerata \
	--genomeDir /path/to/output/index \
	--genomeFastaFiles /path/to/ref.fasta1 /path/to/ref.fasta2... \
	--sjdbGTFfile /path/to/annotations.gtf \
	--sjdbOverhang ReadLength-1

# 使用的参考基因组：一般来说，不应该包含patches和可变单倍体。可以使用的比如：
# EMSEMBL：ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# NCBI：ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# annotation: 来源要和参考基因组匹配,比如注释和基因组都是来源ENSEMBL的,或者都是来自UCSC的;千万不要混用。这是因为它们对染色体的命名不一样。
# 如果你手动更改成一样的也可以使用。


# --sjdbGTFfile	<anno.gtf>		选项默认只处理GTF中exon的行;
# --sddbGTFfeatureExon exon 	默认值为exon,使只处理exon的行.
# --sjdbGTFtagExonParent 		如果使用的是GTF3格式文件,需要指定本选项。
# --sjdbGTFtagExonParentTranscript <transcript_id> 	如果指定,将只会处理能比对到指定transcript(transcript_id)的exon
# --sjdbFileChrStartEnd <sjdbFile.txt> 指定可变剪接的注释文件
# sjdbOverhang: 对于2x100b的Illumina双端测序来说, 是100-1=99.如果length是变化的，就用max(length)-1
# --genomeSAindexNbases <int> 	对于较小的基因组,需要用此选项指定N的值.计算方式: min(14, log2(genomeLength)/2-1). 1Mb的基因组一般为7.
# --genomeChrBinNbits <int>		对于较大的基因组,需要用此选项指定N的值,计算方式: min(18, log2(genomeLength/NumberofReference)).对于3GB含100000个染色体/scaffolds的基因组来说,取15.
```

## 2. 进行比对

```bash
STAR --runThreadN numberofthreads \
	--genomeDir index_path \
	--readFilesIn read1.fastq read2.fastq
	--outFileNamePrefix path/prefix

# --readFilesCommand <command>		如果是压缩文件,可以提供解压缩命令.zcat, gzip -c, bzip2 -c.依据压缩文件的压缩形式而不同.
# --outSAMstrandField introMotif 	指定值为introMotif将会为unstranded RNA-seq数据生成可用于cufflinks和cuffdiff的文件(含XS strand attribute)
# --outFilterIntroMotifs RemoveNoncanonical 	后续进行cufflinks的话,推荐设定此选项.
# --outSAMtype BAM Unsorted 		输出未排序的aligned.out.bam,可以直接用于HTSeq,不用进行name sorting
# --outSAMtype BAM SortedByCoordinate 输出根据坐标排序的aligned.sortedByCoord.out.bam,与samtools sort命令的输出一致
# --outSAMtype BAM Unsorted SortedByCoordinate	两种情况均输出.
# --quantMode TranscriptomeSAM		将会输出翻译成转录本坐标的bam文件,aligned.toTranscriptome.out.bam.
									这个可以用于转录本定量软件的输入.比如RSEM, eXpress.
```

## 3. 2-pass mapping with re-generated genome  
在初次比对后,为了提高比对到新的剪接位点的reads数,推荐再进行一次比对.新的比对将需要提供上一次比对产生的剪接位点文件.

```bash
# 1. 把SJ.out.tab文件合并,过滤掉在线粒体染色体上的junction、non-canonical junction。
# 1. 如果使用了注释的话,这里只需考虑新的剪接位点.

# 2. 把1.中得到的合并文件通过--sjdbFileChrStartEnd指定输入, --sjdbGTFile指定注释,重新生成索引
# 3. 使用新的索引再次进行比对
```

## 4. 2-pass mapping without re-generated genome  
要进行此项,在首次进行比对时,就要设置--twopass1readsN选项,值为所有reads数.此外，你还得设置--sjdbOverhang参数.
