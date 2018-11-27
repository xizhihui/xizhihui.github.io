---
title: 比对 Bowtie2
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 软件和包
---

bowtie2是个超快的、内存占用少的序列比对工具，善于比对相对较长的基因组。bowtie2有gapped、pair-end和local比对模式，可以多线程进行。它是许多pipeline的首个步骤，例如变异检测，CHIP-seq，RNA-seq，BS-seq等等。
bowtie2不像常规目的的比对工具如MUMmer，Blast等。它在大的参考基因组的比对上表现更好，因为它针对当前各个测序平台的测序reads进行过优化。如果你的目的是比对很大的两个序列，比如基因组之间的比对，你应考虑使用MUMmer。如果你的目的是比对相对较短的序列如大肠杆菌的基因组，用bowtie2可以大大减少你的时间。

<!--more-->

## bowtie1和bowtie2的区别

+ 对于长于50bp的序列,bowtie2速度更快也更灵敏，内存占用更少；但是Bowtie1在短于50bp的比对上有时会更快或更灵敏
+ bowtie2通过gap penalty支持gapped alignment
+ bowtie2支持局部比对(soft clipped, not end-to-end)
+ bowtie1的read length上限是1000bp，而bowtie2则没有上限
+ bowtie2支持参考基因组的跨Ns比对
+ 对于pair-end reads,如果比对不上,bowtie2尝试将单个reads进行比对
+ bowtie2不支持colorspace reads的比对

## Getting started

### 1. 对参考基因组建立索引

```bash
# prefix是输出的索引文件前缀, .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2
bowtie2-build path/to/ref.fa prefix
```

### 2. 进行比对

```bash
# 对于单端测序(single end)
bowtie2 -x prefix -U Single_end.fq -S align.sam > align.log
# 对于双端测序(paire-end)
bowtie2 -x prefix -1 read1.fq -2 read2.fq -S align.sam > align.log
# 如果要使用局部比对模式
bowtie2 --local -x prefix -U sigle_end.fq -S align.sam > align.log
```

### 3. 从index里获取原始ref.fa 

```bash
bowtie2-inspect prefix
```
## 构建好的index下载
总卷 | 分卷
---|---
[H. sapiens, UCSC hg18](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg18.zip) | [part1](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg18.1.zip), [part2](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg18.2.zip), [part3](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg18.3.zip)
[H. sapiens, UCSC hg19](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip) | [part1](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.1.zip), [part2](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.2.zip), [part3](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.3.zip)
[H. sapiens, NCBI GRCh38](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz) | none
[M. musculus, UCSC mm10](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip) | [part1](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.1.zip), [part2](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.2.zip), [part3](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.3.zip)
[M. musculus, UCSC mm9](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm9.zip) | [part1](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm9.1.zip), [part2](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm9.2.zip), [part3](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm9.3.zip)
[R. norvegicus, UCSC rn4](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/rn4.zip) | [part1](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/rn4.1.zip), [part2](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/rn4.2.zip), [part3](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/rn4.3.zip)


## bowtie2相关术语解读
### 1. 比对模式

bowtie2里面有两种比对模式：end-to-end, local。前者在比对时要read的首尾比对上,中间没有比对上的作为gap进行罚分.而local比对模式在比对时，优先保证read内部而非两端的比对，就比对效果上看，local模式下的read两端没有比对上，类似切了两端，所以又称为soft clipped。

### 2. scores: higher = more similar

比对得分指示着read与参考序列的相似性，得分越高相似性越高。比对得分的控制可以通过指定各个情况的得分实现。如下：

```bash
--ma:	match bonus,比对上的得分
--mp:   mismatch penalty,比对不上的罚分
--np:	read/ref序列中有N的罚分
--rdg:	affine read gap penalty,比对时read出现gap的罚分;注意gap与mismatch的不同
--rfg:	affine reference gap penalty,比对时ref上出现gap的罚分
--score-min:	指定的最低得分,高于该得分的比对才算比对成功
```

### 3. mapping quality: higher = more unique

比对质量指示着read比对到参考序列上的唯一性，得分越高，越是唯一比对。如果一个read有多种比对情况,那么某个比对情况score越高，我们就说这个比对越唯一。

### 4. mixed mode

如果pair-end的reads比对不上，bowtie2会尝试使用单个read进行比对，这个情况称为mixed mode。这样的结果就是得到的比对率可能会比不使用该模式的比对率高辣么一点点。你可以指定--no-mixed参数禁用。

### 5. reporting

bowtie2有三种结果模式,默认模式是搜索多个比对情况,然后输出最佳比对。-a模式会搜索并输出所有比对情况。-k模式则是搜索1至多个比对情况，并输出。如果多个比对情况得分相同的话，bowtie2会随机选一个。
