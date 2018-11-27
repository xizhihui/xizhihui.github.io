---
title: 差异分析 Ballgown
date: 2018-10-12
toc: true
description: 呦呦呦
categories: Bioinformatics
tags:
  - 软件和包
  - 差异分析
---

Ballgown是一款灵活的用于RNA-seq数据差异分析的软件，除了差异分析，他还可以进行转录本的组织、可视化和分析表达程度。

<!--more-->

## preprocessing
在正式使用之前,你应该做完以下几步：
+ RNA-seq的reads应该比对到参考基因组上了
+ 转录本也组装完成，或者你有一个参考转录本
+ 用于分析的表达矩阵应该是ballgown的可读格式

推荐的两个pipeline：
+ pipeline1: TopHat2 -> stringtie -> ballgown
+ pipeline2: TopHat2 -> Cufflinks -> Tablemaker -> ballgown

Tablemaker:     https://github.com/leekgroup/tablemaker

```bash
# pipeline1
tophat2 -G ref.gff -o output_path -p 6 ref_index reads.fastq
stringtie -B -G ref.gff -p 6 accepted_hist.bam -o stringtie.gff

# pipeline2
tophat2 -G ref.gff -o output_path -p 6 ref_index reads.fastq
cufflinks -g ref.gff -o output_path accepted_hits.bam
# tablemaker调用cufflinks
```

## 安装


```r
source('http://bioconductor.org/biocLite.R')
biocLite('ballgown')
```

## Ballgown需求的格式

https://github.com/alyssafrazee/ballgown#ballgown-readable-expression-output


```bash
# 文件目录
extdata/
	sample01/
		e2t.ctab
		e_data.ctab
		i2t.ctab
		i_data.ctab
		t_data.ctab

# e_data.ctab: exon-level expression measurements, one row per exon
# e_id 	chr strand	start	end rcount	ucount	mrcount	cov	cov_sd	mcov	mcov_sd

# i_data.ctab: intro-(ie. junction-)level expression measurements, one row per intron
# e_id 	chr	strand	start	end	rcount	ucount	mrcount

# t_data.ctab: transcript-level expression measurements, one row per transcript
# t_id	chr	strand	start	end	t_name	num_exons	length	gene_id	gene_name	cov 	FPKM

# e2t.ctab:	which exons belong to which transcripts
# e_id t_id

# i2t.ctab: which introns belong to which transcripts
# i_data	t_data
```

## 使用

### 1. 加载数据

```r
library(ballgown)

# 1.寻找数据目录
data_directory = system.file('extdata', package='ballgown')		# 找到Ballgown所在的目录下extdata目录
data_directory

# 1.1 番外：自己构造目录文件
samples_path = data.frame(c('sample01_path', 'sample02_path', ...))

# 2. 构造Ballgown对象
bg = ballgown(dataDir=data_directory, samplePattern='sample', meas='all')
# meas的参数是那些.ctab文件不重复列的列名
# bg = ballgown(samples=samples_path, samplePattern='sample', meas='all')
bg

# 如果数据量太大
# in load.R
library(ballgown)
data_directory = system.file('extdata', package='ballgown')
bg = ballgown(dataDir=data_directory, samplePattern='sample', meas='all')
save(bg, file='bg.rda')
# 然后在每次使用时,运行以下命令,然后通过load()加载bg.rda
R CMD BATCH load.R
```

### 2. 处理组装数据
一个Ballgown对象有6个分支, structure, expr, indexes, dirs, mergeDate, meas.


```r
# structure十分依赖GenomicRanges包。它指定了组装的转录本的基因位置，外显子内含子转录本的关系。
structure(bg)$exon
structure(bg)$intron
structure(bg)$trans

# expr包含了表达数据
# *expr(bg, meas_name)
transcript_fpkm = texpr(bg, 'FPKM')
transcript_cov = texpr(bg, 'cov')
whole_tx_table = texpr(bg, 'all')
exon_mcov = eexpr(bg, 'mcov')
junction_rcount = iexpr(bg)
whole_intron_table = iexpr(bg, 'all')
gene_expression = gexpr(bg)

# indexes
indexes(bg)
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=10))
exon_transcript_table = indexes(bg)$e2t
transcript_gene_table = indexes(bg)$t2g
phenotype_table = pData(bg)
head(transcript_gene_table)

# 其他
head(bg@dirs)
bg@mergeDate
bg@meas
```

### 3.数据特征可视化

+ 可视化组装的转录本结构

```r
plotTranscripts(gene='XLOC_000454',
				gown=bg, samples='sample12',
				meas='FPKM',
				colorby='transcript',
				main='transcripts from gene XLOC_000454: sample12, FPKM')
# 多个样本同时画图
plotTranscripts('XLOC_00054', bg,
				samples=c('sample01', 'sample02'),
				meas='FPKM',
				colorby='transcript')
```

+ 进行分组比较表达量

```r
plotMeans('XLOC_00054', bg,
		   groupvar='group',
		   meas='FPKM',
		   colorby='transcript')
```

### 4. 差异分析

ballgown默认使用parametric F-test comparing nested linear model进行统计.它同时使用两个模型,一个模型是包含感兴趣协变量,如case/control status, 如time;另一个不包含协变量.显著性的p-value表明包含协变量的模型得到的结果要比不包含的模型结果更为显著，也就说明具有差异。而q-value<0.05表明FDR被控制在5%。stattest函数可自动处理两组比较和多组比较和时间系列比较。

+ 组组比较，group分组


```r
stat_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group')
head(stat_resuts)
```

+ 时间系列比较


```r
# 如果你的时间点较少，比如低于5，建议将其作为分组变量进行组组比较
pData(bg) = data.frame(pData(bg), time=rep(1:10,2))		# dummy time corvariate
timecourse_results = stattest(bg, feature='transcript', meas='FPKM', covariate='time', timecourse=T)
```

+ 调整混杂因素


```r
# 你可以通过pData调整任何混杂因素,也可以像下面一样添加
group_adj_timecourse_results = stattest(bg, feature='transcript',
										meas='FPKM', covariate='time', 
										timecourse=T, adjustvars='group')
```

+ 自定义比较模型 

```r
# sex,age,group,time vs. group,time
	# 模拟数据
set.seed(43)
sex = sample(c('M', 'F'), size=nrow(pData(bg)), replace=T)
age = sample(21:52, size=nrow(pData(bg)), replace=T)
	# 创建design matrices
mod = model.matrix(~ sex + age + pData(bg)$group + pData(bg)$time)
mod0 = model.matrix(~ pData(bg)$group + pData(bg)$time)
	# 进行差异分析
adjusted_results = stattest(bg, feature='transcript', meas='FPKM', mod0=mod0, mod=mod)
head(adjusted_results)
```

+ 输出到其他差异统计软件如limma,limma voom, DESeq, DEXSeq, EdgeR, etc.


```r
# 见2.处理组装数据
*expr(bg, meas_name)
```

+ 简单地转录聚类
如此做是为了避免一个基因有多个相似的转录用作组装，会导致表达差异检测不准确。可以进行一下聚类看看。聚类使用的聚类是Jaccard distance;你也可以使用k-means clustering或者hierarchical clustering。

```r
clusterTranscripts(gene='XLOC_000454', gown=bg, k=2, method='kmeans')
# 也可以可视化
plotLatentTranscripts(gene='XLOC_00054', gown=bg, k=2, method='kmeans', returncluster=F)
```

+ 针对基因进行聚合表达量

```r
agg = collapseTranscripts(gene='XLOC_00054', gown=bg, k=2, method='kmeans')
stattest(gowntable=agg$tab, pData=pData(bg), feature='transcript_cluster',
		covariate='group', libadjust=F)
```
