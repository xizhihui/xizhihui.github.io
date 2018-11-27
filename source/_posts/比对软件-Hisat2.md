---
title: 比对软件 Hisat2
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 软件和包
---

hisat2是快速灵敏的比对软件,可用于全基因组测序，转录组测序，外显子测序的数据比对.基于GCSA（bwt的拓展），我们设计了graph FM index用于比对。hisat2的比对结果是sam格式文件，你可以使用samtools，GATK等软件进行后续的分析.

<!--more-->

## 下载与安装


```bash
conda install hisat2
```

## [Getting started](https://ccb.jhu.edu/software/hisat2/manual.shtml)

### 1. 构建索引

如果你用到了--snp, --ss, --exon等参数,对于人类基因组的大小来说,hisat2需要200G RAM。这几个参数慎用。[预构建的索引](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data).

```bash
# 输出8个以index为前缀的.n.ht2的文件,n=1-8
hisat2-build ref.fa index
```

### 2. 进行比对

```bash
# single-end
hisat2 -f -x index -U reads.fa -S align.sam > align.log
# pair-end
hisat2 -f -x index -1 read_1.fa -2 read_2.fa -S align.sam > align.log
```

### 3. [使用samtools进行下游分析](http://samtools.sourceforge.net/mpileup.shtml)

```bash
# 转为bam
samtools view -bS align.sam > align.bam
# 对bam进行排序
samtools sort align.bam -o align.sorted.bam
# variant calling
samtools mpileup -uf ref.fa align.sorted.bam | bcftools view -bvcg - > align.raw.bcf
bcftools view align.raw.bcf
```

## options
### hisat2：进行比对

```bash
hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]

# 选项
-x <hisat2-idx> 			索引的前缀
-1 <ml> 					以逗号分隔的多个pair-end read1
-2 <ml>						以逗号分隔的多个pair-end read2
-U <r>						以逗号分隔的多个single-end reads
--sra-acc <number>			指定sra accession number
-S <hit>					指定输出文件

# 输入可选参数
-q 							说明输入文件格式是fastq文件
--qseq 						说明输入文件格式是QSEQ文件
-f 							说明输入文件是fasta文件
-r 							说明输入文件是每行一个序列，除此无他
-c 							reads是以序列的形式用逗号分隔输入的
-s,--skip <int>				reads的前int个跳过
-u <int>					只比对int个reads
-5,--trim5 <int>			比对前截去5‘端的int个base
-3,--trim3 <int>			比对前截去3’端的int个base
--phred33					说明fastq的碱基质量体系phred33
--phred64 					说明fastq的碱基质量体系是phred64
--solexa-quals				说明质量体系是solexa,并且转换成phred

# 比对可选参数
--n-ceil <func>				指定每条reads允许N碱基的个数的函数;超过将被丢弃
--ignore-quals				在对错配罚分时,考虑该位置的碱基质量
--nofw/--norc 				指定后,在pair-end无法比对时,不会试图去比对参考序列的forward链(nofw)和reverse链(norc)

# 打分选项
--mp MX,MN 					指定错配时的最大和最小罚分
--sp MX,MN 					指定发生soft clipped的碱基的最大和最小罚分
--no-softclip 				禁用softclip
--np <int>					指定N碱基的罚分
--rdg m1,m2					指定比对时read gap open(m1)和gap extend(m2)的罚分
--rfg m1,m2 				指定比对时reference gap open(m1)和gap extend(m2)的罚分
--score-min <func>			指定比对得分的函数, 当超过计算所得分数时才算一个成功比对

```

### hisat2-index: 构建索引

```bash
hisat2-build [options]* <reference_in> <ht2_base>

# 选项
-f 				指定ref的文件格式
-c 				以逗号分隔的形式指定有多个ref 序列，而非fasta文件的列表
--large-index	指定构建大型索引即便参考序列短于 4 billion bp
-a 				禁止自动使用参数（--bmax，--dv，--packed）
-r 				不构建3.ht2,4.ht2
-3 				只构建3.ht2,4.ht2
-p 				使用多少个线程
--snp <path>	提供snps的列表
--haplotype <path>	提供haplotypes的文件
--ss 			提供间接位点(splice sites)的文件.与--exon联用
--exon <path>	提供外显子文件,与--ss联用
--cutoff <int>  只对参考序列int个base构建索引,丢弃其余部分
-q 				静默运行
```

### hisat2-inspect:从index提取源参考序列

```bash
hisat2-inspect [options]* <ht2_base>

# 选项
-a <int>				输出时每隔60base就换行
-n 						输出参考序列名,每行一个
-s 						打印进行构建索引的参数
--snp 					打印snps
--ss 					打印splice site
--exon 					打印exon
```
