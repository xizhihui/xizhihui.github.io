---
title: 可变剪切 rMATS和rmats2sashimiplot
date: 2018-10-28
toc: true
categories: Bioinformatics
tags:
  - 软件和包
---

## 安装

`environment: py2.7.x`

+ rMATs

```
# in linux
wget http://rnaseq-mats.sourceforge.net/rMATS.4.0.2.tgz

# check version to use, in python
import sys
print sys.maxunicode
# 1114111: rMATS-turbo-xxx-UCS4
# 65535: rMATS-turbo-xxx-UCS2

"path/to/rMATs-turbo-xxx-UCS{x}/rMATS.py"
```

<!--more-->

+ rmats2sashimiplot

```
wget https://files.pythonhosted.org/packages/f7/15/cd6fa8dc70de55be94f4a55edd388e5ba99c77eca8fd73e556e66e9d8110/rmats2sashimiplot-2.0.2.tar.gz
tar -zxvf rmats2sashimiplot-2.0.2.tar.gz & cd rmats2sashimiplot-2.0.2
python setup.py install
```

## 使用

+ rMATs

```
# input as fastq, samples in each group list in s1.txt/s2.txt separated with comma
python rmats.py --s1 s1.txt --s2 s2.txt --gtf gtfFile --bi STARindexFolder -od outDir -t readType -readLength readLength [options]*
# input as sorted.bam, samples in each group list in s1.txt/s2.txt separated with comma
python rmats.py --b1 b1.txt --b2 b2.txt --gtf gtfFile --od outDir -t readType --nthread nthread --readLength readLength --tstat tstat [options]*
```

+ rmats2sashimiplot

```
# input as sam
$rmats2sashimiplot --s1 s1_rep1.sam[,s1_rep2.sam]* --s2 s2.rep1.sam[,s2.rep2.sam]* -t eventType -e eventsFile --l1 SampleLabel1 --l2 SampleLabel2 --exon_s exonScale --intron_s intronScale -o outDir

# input as bam
$rmats2sashimiplot --b1 s1_rep1.bam[,s1_rep2.bam]* --b2 s2.rep1.bam[,s2.rep2.bam]* -c coordinate:annotaionFile --l1 SampleLabel1 --l2 SampleLabel2 --exon_s exonScale --intron_s intronScale -o outDir

# a group file provided
$rmats2sashimiplot --b1 s1_rep1.bam[,s1_rep2.bam]* --b2 s2.rep1.bam[,s2.rep2.bam]* -c coordinate:annotaionFile --l1 SampleLabel1 --l2 SampleLabel2 --exon_s exonScale --intron_s intronScale -o outDir --group-info gf.gf
```

## rMATS 结果
### 识别的剪切事件类型

rMATS analyzes skipped exon (SE), alternative 5' splice site (A5SS), alternative 3' splice site (A3SS), mutually exclusive exons (MXE), and retained intron (RI) events. 

![5种ASE](https://s1.ax1x.com/2018/10/28/ic5Jp9.jpg)
![5种ASE](https://s1.ax1x.com/2018/10/28/ic5BkD.png)

### 结果文件释义
+ AS_Event.MATS.JC.txt, 只使用 Junction counts (counts of reads that span splicing junctions) 检测到的 AS_Event 结果
+ JS.raw.input.AS_Event.txt, 进行 AS_Event.MATS.JC 分析的源数据
+ AS_Event.MATS.JCEC.txt, 同时使用 junction counts 和 reads on target 检测到的 AS_Event 结果
+ JCEC.raw.input.AS_Event.txt，进行 AS_Event.MATS.JCEC.txt 分析的源数据
+ fromGTF.AS_Event.txt, 来源于 GTF 和 RNA 的所有可能可变剪切事件

### 结果文件各列释义
+ IJC_SAMPLE_*: inclusion junction counts
+ SJC_SAMPLE_*: skipping junction counts
+ IC_SAMPLE_*: inclusion counts
+ SC_SAMPLE_*: skipping counts
+ IncFormLen: length of inclusion form, used for normalization
+ SkipFormLen: length of skipping form, used for normalization
+ IncLevel1: inclusion level for SAMPLE_1 replicates (comma separated) calculated from normalized counts
+ IncLevel2: inclusion level for SAMPLE_2 replicates (comma separated) calculated from normalized counts
+ IncLevelDifference: average(IncLevel1) - average(IncLevel2)
+ P-Value: Significance of splicing difference between two sample groups. (Only available if statistical model is on)
+ FDR: False Discovery Rate calculated from p-value. (Only available if statistical model is on)
+ upstreamES/upstreamEE: upstreamExonStart, upstreamExonEnd
+ downstreamES/downstreamEE: like above

![IC and SC](https://s1.ax1x.com/2018/10/28/ic5am6.png)
![upstreamES/EE](https://s1.ax1x.com/2018/10/28/ic5t61.png)

### 相关计算

`IncLevel1 = (IC_SAMPLE_1 / IncFormLen) / (IC_SAMPLE_1 / IncFormLen + SC_SAMPLE_1 / SkipFormLen)`

`IncFormLen = 2 * (Junction_length - read_length + 1)`

`SkipFormLen = Junction_length - read_length + 1`

`Junction_length = 2 * (read_length - anchor), anchor = 8 bp(default)`

### 结果解读：
1. JunctionCountOnly

在对剪切事件发生的差异分析过程中，rMATs 采用了 2 种定量方式，即 JunctionCountOnly 和 ReadOnTargetAndJunctionCounts. 对于前者，它把那些会全部比对到 alternatively spliced exon(ASE) 的 reads 进行剔除，而后者不会。也就是说，Junction counts 是 reads 覆盖 (span) 到 splicing site 的 reads count，read on target 就是被 JunctionCountOnly 剔除的 reads，也就是 这些 reads 整体只比对到 ASE 区域的某个部分。

2. PSI

PSI, 全称为 Percent-splice-in，可以针对 isoform，exon，ASE 进行计算。对于 ASE 来说，PSI = splice_in / (splice_in + splice_out), 在 rMATS 里，splice_in 和 splice_out 是支持 splice_in 和 splice_out 发生的 reads 数目，可以同基因表达的 reads count 作类比。

## rmats2sashimiplot

### 结果展示

![可视化结果](https://raw.githubusercontent.com/Xinglab/rmats2sashimiplot/master/img/plotwitheventgf.png)

### 结果解读

+ Y-轴释义

![Y-轴释义](https://raw.githubusercontent.com/Xinglab/rmats2sashimiplot/master/img/RPKM.png)

## 参考

+ [rMATS](http://rnaseq-mats.sourceforge.net/user_guide.htm)
+ [rmats2sashimiplot](https://github.com/Xinglab/rmats2sashimiplot)
+ [biostars](https://www.biostars.org/p/256949/)
