---
title: 生信基础 unique比对的获取
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 生信基础
---




## 1.Sam文件各标签含义（tophat/hisat2)

+ NH:i:<N>: N=1时 为unique。常用于tophat/hisat2产生的sam文件unique read筛选。
+ CC:Z: 当为‘=’为map到同一条基因上，一般在map基因组时由于内含子存在而容易出现，他只代表两种不同的方式，计数时应记为1。此处一般为其他基因的名字。CP:i 和HI：i标签为map到第i条基因及起始位置。
+ YT:Z:S 代表的含义与bowtie产生的sam也不同。具体还未知！其他标签AS，XN,XM,XO,XG,NM,MD等如下图可以看出都相同。

<!--more-->

> 对于tophat/hisat2比对产生的sam文件我们可以直接筛选NH标签。


```bash
grep‘NH:i:1’ out.sam >unique.sam
```

## 2.Sam文件各标签含义（bowtie2)
+ AS:i:<N>Alignmentscore.可以为负的，在local下可以为正的。 只有当Align≥1 time才出现
+ XS:i:<N>Alignmentscorefor second-best alignment. 当Align>1 time出现
+ YS:i:<N>Alignmentscorefor opposite mate in the paired-end alignment. 当该read是双末端测序中的另一条时出现
+ XN:i:<N>Thenumber of ambiguous bases in the reference covering this alignment.（推测是指不知道错配发生在哪个位置，推测是针对于**和缺失，待查证）
+ XM:i:s错配碱基的数目
+ XO:i:<N>Thenumberof gap opens(针对于比对中的**和缺失)
+ XG:i:<N>Thenumberof gap extensions(针对于比对中的**和缺失延伸数目)
+ NM:i:<N>Theeditdistance。（edits:**/缺失/替换数目)
+ YF:Z:s该reads被过滤掉的原因。可能为LN(错配数太多，待查证)、NS(read中包含N或者．)、SC(match bonus低于设定的阈值)、QC(failing quality control，待证)
+ YT:Z:s值为UU表示不是pair中一部分、CP是pair且可以完美匹配、DP是pair但不能很好的匹配、UP是pair但是无法比对到参考序列上。
+ MD:Z:s比对上的错配碱基的字符串表示。

> 由于bowtie2产生的sam文件并没有NH标签，所以提取uniqueread可能比较麻烦。首先提取“AS”标签表示能比对上的read（>=1 time），然后利用grep反正则表达式过滤掉XS标签得到我们需要的unique read。


```bash
grep “AS:” aligned.sam | grep –v “XS:” >unique_alignments.sam
```

对于双端测序用bowtie2比对筛选unique concordant pair时则需要在上一步的基础上增加如下命令：


```bash
grep ‘YT:Z:CP’ unique.sam>pair-end_unique.sam
```

## 3.Bwa获取unique


```bash
samtools view bwa.bam | grep "XT:A:U"
```
