---
title: 扩增子测序 mothur使用示例
date: 2018-10-15
toc: true
description: mothur, 搞那么长的命名后缀，让你想骂你 “mothur”
categories: Bioinformatics
tags:
  - 软件和包
  - 扩增子
---



## 1 初始说明

+ 测序数据类型

Illumina Miseq paired-end reads

+ 实验设计

断奶后365天(dpw 365)的小鼠排泄物，比较初始10天(dpw 10)和中间10天(dpw140-150)的排泄物的微生物组的稳定性(肠道微生物组的变化情况)。为了简化操作，只用到一只小鼠的十个时间点(前5后5)的数据。这里还有模拟了由21种细菌组成的菌群的全基因组测序数据。先用小鼠的排泄物测序数据学习分析微生物群落，然后用模拟的菌落判断分析的错误率和它在其他分析中的作用。

+ 关于软件

mothur既提供交互模式(像python)，也提供命令行模式；后者可以进行批量操作。

<!--more-->

```bash
# 解析模式一
mothur # 进入解析模式
make.files(....)
quit() # 退出解析模式

# 命令行模式
mothur '#make.config(file=stability.files, processors=8)'
# 批量操作,把多个命令写入cmd.batch文件里
mothur cmd.batch
```

## 2 进入mothur交互模式

```bash
mothur
set.dir(tempdefault=MiSeq_SOP)
set.dir(inputdir=MiSeq_SOP)
help(query_cmd)
```

## 3 stability.files

这个文件记录了每个样品对应的read1和read2文件。

```bash
make.file(inputdir=MiSeq_SOP, type=fastq, prefix=metadata)
```

## 4 序列合并和质控：减少测序和PCR错误

使用make.configs合并序列,合并逻辑是先对两个reads比对，找到overlap；根据overla的base quality score，如果是base-gap，那么base的质量值要大于25才取base否则认为该位置没有base；如果是base-base，base的质量值要大于6才取base；否则置为N。


```bash
# 1. 合并序列
make.contigs(file=metadata.files, processor=8)
		# MiSeq_SOP/metadata.trim.contigs.[fasta,qual]
		# MiSeq_SOP/metadata.scrap.contigs.[fasta, qual]
		# MiSeq_SOP/metadata.contigs.[report, groups]
# 2. 统计合并后的序列
summary.seqs(fasta=metadata.trim.contigs.fasta)
		# MiSeq_SOP/metadata.trim.contigs.summary

# 3. 筛选合并后的序列, 排除含有N和长度大于275的序列
screen.seqs(fasta=metadata.trim.contigs.fasta, 
        summary=metadata.trim.contigs.summary, 
        group=metadata.contigs.groups, 
        maxambig=0, maxlength=275)
		# MiSeq_SOP/metadata.contigs.[pick, good].groups
		# MiSeq_SOP/metadata.trim.contigs.[good.summary, good.fasta, bad.accnos]
```

## 5 序列的去冗余


```bash
# 1.去测序冗余
unique.seqs(fasta=metadata.trim.contigs.good.fasta)
		# MiSeq_SOP/metadata.trim.contigs.good.[names, unique.fasta]
# 2.生成count矩阵
count.seqs(name=metadata.trim.contigs.good.names, 
        group=metadata.contigs.good.groups)
		# MiSeq_SOP/metadata.trim.contigs.good.count_table

# 3.对count矩阵进行统计
summary.seqs(fasta=metadata.trim.contigs.good.unique.fasta, 
        count=metadata.trim.contigs.good.count_table)
		# MiSeq_SOP/metadata.trim.contigs.good.unique.summary

# 4.去非扩增区域冗余; 通过设置参考数据库中你的PCR扩增区域来去除非扩增区域的序列，以减少比对量
#   你可以使用oligos来指定primer文件进行去除，keepprimer
pcr.seqs(fasta=../silva.bacteria/silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=8)
		# silva.bacteria/silva.bacteria.pcr.fasta
rename.file(input=../silva.bacteria/silva.bacteria.pcr.fasta, new=silva.v4.fasta)
		# MiSeq_SOP/silva.v4.fasta

# 5. 对silva.v4.fasta进行统计
summar.seqs(fasta=silva.v4.fasta)
		# MiSeq_SOP/silva.v4.summary
```

## 6 序列的预聚类

要进行预聚类，我们需要先进行序列比对，获得相似序列。

### 6.1 序列比对

```bash
# 1.比对并统计
align.seqs(fasta=metadata.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)
		# MiSeq_SOP/metadata.trim.contigs.good.unique.align.['', report]
summary.seqs(fasta=metadata.trim.contigs.good.unique.align,
        count=metadata.trim.contigs.good.count_table)
		# MiSeq_SOP/metadata.trim.contigs.good.unique.summary

# 2. 根据比对统计结果，大部分序列比对在1968-11550, polymer最大为8；再进行一次筛选序列
screen.seqs(fasta=metadata.trim.contigs.good.unique.align, count=metadata.trim.contigs.good.count_table, summary=metadata.trim.contigs.good.unique.summary, start=1968, end=11550, maxhomop=8)
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.summary
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.align
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.accnos
		# MiSeq_SOP/metadata.trim.contigs.good.good.count_table
summary.seqs(fasta=current, count=current)
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.summary

# 3. 把比对结果中未必对上部分序列和表示gap的-去掉
filter.seqs(fasta=metadata.trim.contigs.good.unique.good.align, vertical=T, trump=.)
		# Length of filtered alignment: 376 
		# Number of columns removed: 13049
		# Length of the original alignment: 13425
		# Number of sequences used to construct filter: 16299
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.fasta
		# MiSeq_SOP/metadata.filter

# 4. 去冗余, 在上一步剪切序列时可能产生
unique.seqs(fasta=metadata.trim.contigs.good.unique.good.filter.fasta, count=metadata.trim.contigs.good.good.count_table)
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.count_table
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.fasta
```

### 6.2 预聚类

这里使用pre.cluster进行预聚类。它首先根据group对序列分组，然后依据序列丰度进行排序，最后把差异在2个核苷酸及以内的序列进行合并。一般是100bp取1个核苷酸的差异,这里的序列长度是250左右，所以为2。

```bash
pre.cluster(fasta=metadata.trim.contigs.good.unique.good.filter.unique.fasta, count=metadata.trim.contigs.good.unique.good.filter.count_table, diffs=2)
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.precluster.[SampleID].map
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.precluster.[fasta,count_table]
```

## 7 去除嵌合体 chimeras

chimera.vsearch将会把数据按样品分组，检测嵌合体的存在。我们倾向于使用高丰度序列作为reference。


```bash
# 1. 检测嵌合体序列；输出的count_table依据去除了
chimera.vsearch(fasta=metadata.trim.contigs.good.unique.good.filter.unique.precluster.fasta,
        count=metadata.trim.contigs.good.unique.good.filter.unique.precluster.count_table, 
        dereplicate=T)
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.[pick.count_table, chimeras,accnos]
# 2. 去除嵌合体序列
remove.seqs(fasta=metadata.trim.contigs.good.unique.good.filter.unique.precluster.fasta,
        accnos=metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
		# MiSeq_SOP/metadata.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta
# 3. 统计
summary.seqs(fasta=current, count=current)
	# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.summary
```

## 8 进行聚类（分类）


```bash
# 1. 分类
classify.seqs(fasta=metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta,
        count=metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table,
        reference=../trainset9_032012.pds.fasta, 
        taxonomy=trainset9_032012.pds.tax, cutoff=80)
	# MiSeq_SOP/metadata.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.[taxonomy, tax.summary]
	
# 2. 去除污染; 这里的unknown是在上一步分类时无法根据训练集和准确进行准确分类的一类
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table

# 3. 统计
summary.tax(taxonomy=current, count=current)
		# MiSeq_SOP/metadata.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary
```

## 9 检测错误率

如果你在测序的时候有一些模拟群落的数据的话，你就可以检测一下你的数据的错误率。


```bash
# 1. 提取出mock的分组
get.groups(count=current, fasta=current, groups=Mock)
		# MiSeq_SOP/metadata.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table
		# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta

# 2. 检测错误率
seq.error(fasta=current, count=current, reference=HMP_MOCK.v35.fasta, aligned=F)
		# Overall error rate:	6.5108e-05
	# MiSeq_SOP/metadata.trim.contigs.good.unique.good.filter.unique.
	# precluster.pick.pick.pick.error.[summary,seq,chimera,seq.forward, seq.reverse, count, matrix, ref]
```

## 10 生成Mock OTUs


```bash
# 1. 计算距离
dist.seqs(fasta=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta,
    cutoff=0.03)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist

# 2.进行聚类
cluster(column=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist,
    count=metadata.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.[list, steps, sensspec]

# 3. 生成shared文件
make.shared(list=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.list,
    count=metadata.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, 
    label=0.03)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared

# 4. 稀释
rarefaction.single(shared=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared)
		# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.groups.rarefaction
```

## 11 准备分析

接下来我们有2件事要做。第一件就是把序列与OTU对应起来，第二件是把序列和系统发育表型对应起来。
当前我们需要把mock数据从我们的数据集里面去掉。


```bash
remove.groups(fasta=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,
    count=metadata.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,
    taxonomy=metadata.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, 
    groups=Mock)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy
```

## 12 生成OTUs

对于较小的数据集，可以使用dist.seqs和cluster来进行，就像前面Mock OTUs一样。 对于大的数据集，使用cluster.split。在该方法中，会先根据分类信息对序列进行分箱，然后在各个箱内进行聚类。对于指定的分类水平，可以通过指定taxlevel=level来实现。


```bash
# 1. 生成OTUs
# dist.seqs + cluster
dist.seqs(fasta=curent, count=current)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist

cluster(column=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist, count=current)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.[list,steps,sensspec]

# cluster.split
cluster.split(fasta=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta,
    count=metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table,
    taxonomy=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy, 
    splitmethod=classify, cutoff=0.03, taxlevel=4)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.[list,steps,sensspec]
```


```bash
# 2. 生成shared文件(每个样品为行，每个OTU为列，OTU数目矩阵)
makd.shared(list=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.list,
    count=metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, 
    label=0.03)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
```

然后使用classify.otu给每个OTU同对应的物种分类联系起来。


```bash
# 3. otu被观察到的次数（OTU丰度）和对应的物种分类信息
classify.otu(list=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.list,
    count=metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, 
    label=0.03)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.
	# pick.pick.pick.opti_mcc.0.03.cons.[taxonomy, tax.summary]
```

## 13 生成phylotypes 系统发育表型信息

你可以使用phylotype来对序列按照他们的物种分类分箱成系统发育表型。


```bash
# 1. 生成系统发育表型表
phylotype(taxnomy=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.taxonomy)	# current
		# metadata.trim.contigs.good.unique.good.filter.unique.
		# precluster.pick.pds.wang.pick.pick.tx.[sabund, rabund, list]

# 2. 生成shared文件，这里生成level6，genus的;这里设置label是从1-7对应着genus-kingdom
make.shared(list=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.list,
    count=metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table,
    label=1)
	# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.shared
	# label   Group   numOtus Otu01   Otu02   Otu03   Otu04
	# 1       F3D0    62      1663    2626    216     57
	# 1       F3D1    62      1959    1249    246     35

# 3. 生成OTU丰度、对应的物种分类信息
classify.otu(list=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.list,
    count=metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, 
    label=1)
	# metadata.trim.contigs.good.unique.good.filter.unique.
	# precluster.pick.pds.wang.pick.pick.tx.1.cons.[taxonomy,tax.summary]
	# OTU     Size    Taxonomy
	# Otu01   22183   Bacteria(100);...;Lachnospiraceae_unclassified(100);
	# Otu02   54221   Bacteria(100);...;"Porphyromonadaceae"_unclassified(100);
```

## 14 生成系统发育树

我们可以使用clearcut来生成系统发育树，它既可以以比对后的序列作为输入，也可以以序列距离来作为输入。这里以序列距离作为输入。如果以比对后的序列作为输入，你需要指定序列是DNA还是蛋白。


```bash
# 1. 计算距离; fasta和计算OTUs的是一致的
dist.seqs(fasta=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, 
        output=lt, processors=8)
		# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.phylip.dist

# 2. 生成树
clearcut(phylip=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.phylip.dist)	# current
		# metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.phylip.tre
```

## 15 准备分析: OTU-based

下面我们要进行基于OTU的分析，基于phylotype的分析也是一样的操作. 为了简化，我们要把current的各个相关的默认参数指定为OTU的。


```bash
# 这里这样设置后，...0.03.cons.taxonomy就成了metadata.taxonomy
rename.files(taxonomy=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy,
    shared=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared)
```

如果要知道每个样品有多少条序列，可以这么做. 然后对样品进行稀释抽样

```bash
# 1. 对序列进行分组计数；得到样品的最小总序列数是2390
count.groups(shared=metadata.opti_mcc.shared)
		# metadata.opti_mcc.count.summary

# 2. 进行稀释抽样，抽样大小是2390
sub.sample(shared=metadata.opti_mcc.shared, size=2390)
		# metadata.opti_mcc.0.03.subsample.shared
```

## 16 基于OTU的分析: alpha diversity

+ 稀释曲线


```bash
# 生成的稀释表，需要你自己去画图
# calc指定计算哪个alpha diversity index, 有多个的话，sobs-chao-ace
#  ace, bootstrap, chao, coverage, default, heip, invsimpson, jack, npshannon, nseqs, shannon, shannoneven, shannonrange, simpson, simpsoneven, smithwilson, sobs
rarefection.single(shared=metadata.opti_mcc.shared, calc=sobs, freq=100)
		# metadata.opti_mcc.groups.rarefaction
```

+ 生成一个含序列数目、sample converage, observed OTUs, inverse simpson diversity estimation的表

基于这里的ave-std.summary，你还可以进行repeated-measures ANOVA进行分组的组间差异比较


```bash
summary.single(shared=metadata.opti_mcc.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=2390)
		# metadata.opti_mcc.groups.ave-std.summary
		# metadata.opti_mcc.groups.summary
```

## 17 基于OTU的分析: beta diversity

现在我们要比较样品的关系和组成了。

+ OTU 相对丰度热图


```bash
# 你可以直接使用metadata.opti_mcc.0.03.subsample.shared在R里面画热图
heatmap.bin(shared=metadata.opti_mcc.0.03.subsample.shared, scale=log2, numotu=50)
	# metadata.opti_mcc.0.03.subsample.0.03.heatmap.bin.svg
```

+ 计算样品的关系和组成相似性


```bash
dist.shared(shared=metadata.opti_mcc.shared, calc=thetayc-jclass, subsample=2390)
		# metadata.opti_mcc.thetayc.0.03.lt.[dist, ave.dist,std.dist]
		# metadata.opti_mcc.jclass.0.03.lt.[dist, ave.dist,std.dist]

heatmap.sim(calc=jcalss-thetayc)
#heatmap.sim(phylip=metadata.opti_mcc.jclass.0.03.lt.ave.dist)
		# metadata.opti_mcc.thetayc.0.03.lt.ave.heatmap.sim.svg
		# metadata.opti_mcc.jclass.0.03.lt.ave.heatmap.sim.svg
```

+ venn图

当组别在5以内时，利用veen才方便。


```bash
# 利用生成的.sharedotus到R里面画图
venn(shared=metadata.opti_mcc.0.03.subsample.shared, groups=F3D0-F3D1-F3D2-F3D3)
		# metadata.opti_mcc.0.03.subsample.0.03.sharedsobs.F3D0-F3D1-F3D2-F3D3.svg
		# metadata.opti_mcc.0.03.subsample.0.03.sharedsobs.F3D0-F3D1-F3D2-F3D3.sharedotus
```

+ 系统发育树

这里利用tree.shared来画系统发育树，以可视化样品之间的相似性。


```bash
# 这里选择jclass
tree.shared(phylip=metadata.opti_mcc.jclass.0.03.lt.ave.dist)
		# metadata.opti_mcc.jclass.0.03.lt.ave.tre
```

+ 基于系统发育树的差异比较

parsimony, unifrac.unweighted, unifrac.weighted.
我们需要构建一个design文件，以指示那个样品属于哪个组；然后进行计算


```bash
parsimony(tree=metadata.opti_mcc.jclass.0.03.lt.ave.tre, group=mouse.time.desin, groups=all)
		# Tree#   Groups  ParsScore       ParsSig
		# 1       Early-Late      1       <0.001
		# metadata.opti_mcc.jclass.0.03.lt.ave.tre.parsimony
		# metadata.opti_mcc.jclass.0.03.lt.ave.tre.psummary
```

+ PCoA和nmds分析/可视化


```bash
pcoa(phylip=metadata.opti_mcc.jclass.0.03.lt.ave.dist)
		# metadata.opti_mcc.jclass.0.03.lt.ave.pcoa.[axes, loadings]
		# metadata.opti_mcc.jclass.0.03.lt.ave.pcoa.loadings

nmds(phylip=metadata.opti_mcc.jclass.0.03.lt.ave.dist)
		# Number of dimensions:   2
		# Lowest stress : 0.142764
		# R-squared for configuration:    0.918996
		# metadata.opti_mcc.thetayc.0.03.lt.ave.nmds.[iters, stress, axes]
```

+ AMOVA分析/HOMOVA分析


```bash
amova(phylp=metadata.opti_mcc.thetayc.0.03.lt.ave.dist, design=mouse.time.design)
		# metadata.opti_mcc.thetayc.0.03.lt.ave.amova

homova(phylip=metadata.opti_mcc.thetayc.0.03.lt.ave.dist, design=mouse.time.desin)
		# metadata.opti_mcc.thetayc.0.03.lt.ave.homova
```

+ 核心OTU分析：导致差异的OTU


```bash
# 生成的corr.axes文件可用biplot图可视化
corr.axes(axes=metadata.opti_mcc.jclass.0.03.lt.ave.pcoa.axes, 
        shared=metadata.opti_mcc.0.03.subsample.shared, 
        method=spearman, numaxes=3)
		# metadata.opti_mcc.0.03.subsample.spearman.corr.axes

corr.axes(axes=metadata.opti_mcc.jclass.0.03.lt.ave.pcoa.axes, 
        metadata=mouse.dpw.metadata, 
        method=spearman, numaxes=3)
		# mouse.dpw.spearman.corr.axes
```

+ 是否可以将数据分成不同的群落类型？


```bash
# 得到的表Laplace最小值对应k=2。所以可以分为2个类型
get.communitytype(shared=metadata.opti_mcc.0.03.subsample.shared)
		# metadata.opti_mcc.0.03.subsample.0.03.dmm.mix.fit
		# metadata.opti_mcc.0.03.subsample.0.03.dmm.[1-5].mix.posterior
		# metadata.opti_mcc.0.03.subsample.0.03.dmm.[1-5].mix.relabund
		# metadata.opti_mcc.0.03.subsample.0.03.dmm.mix.design，指明各个样品所属类型
		# metadata.opti_mcc.0.03.subsample.0.03.dmm.mix.parameters
		# metadata.opti_mcc.0.03.subsample.0.03.dmm.mix.summary，指明主要是哪个OTU导致的分类型
```

## 18 群体水平的分析

+ metastats


```bash
metastats(shared=metadata.opti_mcc.0.03.subsample.shared, 
        design=mouse.time.design)
		# metadata.opti_mcc.0.03.subsample.0.03.Late_Early.metastats
```

+ lefse


```bash
lefse(shared=metadata.opti_mcc.0.03.subsample.shared, design=mouse.time.design)
```

## 19 基于系统发育的分析 phylogeny-based analysis

### 19.1 alpha diversity


```bash
rename.file(tree=metadata.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.phylip.tre)
rename.file(count=metadata.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table)
phylo.diversity(tree=metadata.tre, count=metadata.count_table, rarefy=T)
		# metadata.1.phylodiv.summary
		# phylodiv.rarefaction
```

### 19.2 beta diversity


```bash
# 生成.ave.dist是距离矩阵，可以像OTU-based analysis那样进行分析那些PCOA之类的。
unifrac.unweighted(tree=metadata.tre, count=metadata.count_table, distance=lt, processors=2, random=F, subsample=2390)
unifrac.weighted(tree=metadata.tre, count=metadata.count_table, distance=lt, processors=2, random=F, subsample=2390)
		# metadata.uwsummary
		# metadata.1.unweighted.ave.dist
		# metadata.1.unweighted.std.dist
		# metadata.tre1.unweighted.phylip.dist
```
