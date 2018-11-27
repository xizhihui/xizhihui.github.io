---
title: 差异分析 DESeq2
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 生信基础
---

> 与DESeq类似的包有：edgeR、limma、DSS、EBSeq、baySeq

## 一、核心逻辑代码

```r
# countData: 以样品为列名,基因为行名的表达矩阵
# colData的形式: 以样品作为行名,列是样品对应的分组类型
# 生成count matrix
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ batch + condition)
# 生成DESeq数据集
dds <- DESeq(dds)
# 进行比较，获得结果
res <- results(dds, contrast=c('condition', 'treat', 'ctrl'))
resultsNames(dds)
res <- lfcShrink(dds, coef=2)

#  DESeqDataSetFromTximport:   由Salmon, Saifish, kallisto生成
#  DESeqDataSetFromHTSeq:      由htseq-counts生成
#  DESeqDataSet:               由RangedSummarizedExpriment生成
```

<!--more-->

DESeq2为count数据提供两类变换方法，使不同均值的方差趋于稳定，rlog和vst，这两个函数可以用于处理含有色散平均趋势的负二项数据（类如RNA-seq）。由于rlog计算量很大，与vst的效果相近，那么，在数据集小于30的时候，使用rlog；大数据集使用vst可以加快速度。但做这样的变换不是用于差异分析，而是用于pca，WGCNA，clustering等与聚类相关的分析才用得到。

## 二、标准流程 standard workflow

### 2.1 input data

+ DESeq的input data是不需要预先进行标准化的，因为软件包内部会自己根据文库大小等进行标准化。  
+ input data有四类，分别对应四个读取函数，见核心逻辑代码部分。
+ DESeqDataSet是包中存取read counts和统计中间值的对象，通常简写为dds。其参数必需要有一个design formula，以被后续估计dispersion和log2 fold changes使用。

### 2.2 DESeqDataSetFromTximport：.txi文件
> 这类数据是由Salmon、Saifish、kallisto、RSEM等软件产生，数据应该是estimated gene counts。  
    使用这些软件的好处在于：  
        1. 它们会修正样品间gene length的可能改变。  
        2. 它们使用的速度更快、占用CPU内存更少（相比于基于比对的软件来讲）。  
        3. 它们可以避免丢弃多匹配片段（基因同源序列），提高灵敏度。

### 2.3 使用DESeqDataSetFromMatrix: count matrix文件
> 这类数据来源于Rsubread包的featureCounts函数。使用DESeqDataSetFromMatrix，需要提供counts matrix、数据框格式的样品信息、design formula。

### 2.4 input: htseq-count
> 具体的使用可见前面，以pasilla来源数据为例

### 2.5 input: SummarizedExperiment

```r
# 使用airway包来源的数据进行示例
library('airway')
data('airway')
se <- airway

library('DESeq2')
ddsSE <- DESeqDataSet(se, design = ~cell + dex)
```

### 2.6 预过滤

```r
# 为了减少内存占用，先把count数很低的数据剔除
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```
### 2.7 指定factor levels
> R 将基于字母顺序默认参考水平，但实际通常是根据对照组作为参考水平。因此有必要时要设置


```r
dds$condition <- factor(dds$condition, levels=c('untreated', 'treated'))
# 也可以直接指定
dds$condition <- relevel(dds$condition, ref='untreated')
# 如果有时对dds取子集时，导致某些水平不含数据，那么这个水平就可以丢弃
dds$condition <- droplevels(dds$condition)
```

### 2.8 合并技术性重复：collapseReplicates函数

## 三、差异表达分析

> 函数DESeq用来进行差异表达分析，所有计算都整合在该函数里面；生成一个结果对象。你需要使用results函数对该结果对象取值，才能获取到包含了log2FoldChange、p-value、p-adj-value等的结果。


```r
dds <- DESeq(dds)
res <- results(dds)
```
> 如果你想修改coeffficient的数量的话，可以这样设置：


```r
res <- results(dds, contrast='2')
# 或者
resFC <- lfcShrink(dds, coef=2)
```

> 如果有太多样本参与比较，可以并行计算


```r
library('BiocParallel')
register(MulticoreParam(4))
dds <- DESeq(dds)
res <- results(dds)
resFc <- lfcShrink(dds, coef=2)
```
> 按照p-value排序结果、取摘要值、统计个数、调整p-value值


```r
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res0.05 <- results(dds, alpha=0.05)
```

## 四、探索和导出结果
### MA-plot
收缩估计(shrinkage estimates)通过在回归时缩小对结果影响较小的系数值来达到改善结果的目的.有多种收缩方法可供选择,如apeglm, ashr, normal.前二者需要安装同名bioconductor包,后者是DESeq2自带的.

```r
plotMA(res, ylim=c(-2,2))
# 使用lfcShrink来改变shrinkage estimators
resApe <- lfcShrink(dds, coef=2, type='apeglm')
resAsh <- lfcShrink(dds, coef=2, type='ashr')
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)   #改变yx轴的值范围
plotMA(resLFC, xlim=xlim, ylim=ylim, main="normal")
plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```

### plot counts
> 计算一个基因的reads counts在各个处理组的情况


```r
plotCounts(dds, gene=which.min(res$padj), intgroup='condition')
# 输出数据框以供ggplot使用
d <- plotCounts(dds, gene=which.min(res$padj),intgroup='condition', returnData=TRUE)
library('ggplot2')
ggplot(d, aes(x=condition,y=count)) +
    geom_point(position=position_jitter(w=0.1,h=0)) +
    scale_y_log10(breaks=c(25,100,400))
```

### 更多有关结果的信息

> p-value是NA？原因有三：某行(row)所有样品都是0 counts；某行的样品值极端大或者小；自动独立过滤引起的低矫正后counts，p-adj会被设为NA。


```r
mcols(res)$description
```
### 富可视化、生成报告

> 可以使用Reporting Tools、regionReport、Glimma、pcaExplore输出结果


```r
write.csv(as.data.frame(resOrdered), file='results.csv')
resSig <- subset(resOrdered, padj < 0.1)
write.csv(as.data.frame(resSig), file='sig_results.csv')
```

### 多因子设计

> 重点在于colData的生成和类似formula(~type+condition)的生成


```r
## DataFrame with 7 rows and 3 columns
##            condition        type sizeFactor
##             <factor>    <factor>  <numeric>
## treated1     treated single-read  1.6355014
## treated2     treated  paired-end  0.7612159
## treated3     treated  paired-end  0.8326603
## untreated1 untreated single-read  1.1383376
## untreated2 untreated single-read  1.7935406
## untreated3 untreated  paired-end  0.6494828
## untreated4 untreated  paired-end  0.7516005
colData(dds) # 存起备用 
ddsMF <- dds 
levels(ddsMF$type)  # 生成新的level 
# re-run DESeq
design(ddsMF) <- formula(~ type + condition) 
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
# 还可以这样干
resMFType <- results(ddsMF, contrast=c('type', 'single', 'paired'))
```

## 五、通过聚类可视化表征数据质量

### 表达矩阵的热图


```r
# ntd,vsd,rld分别是norm,vst,rlog变换后的数据

library('pheatmap')
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c('condition', 'type')])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
```

### 样品距离热图

```r
# 这类图直观地显示了样品间的相似度和差异性
sampleDists <- dist(t(assay(vsd)))

library('RColorBrewer')
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep='-')
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, 'Blues')) )(255)
pheatmap(sampleDistMatrix,
        clustering_distance_rows=sampleDists,
        clustering_distance_cols=sampleDists,
        col=colors)
```

### 样品的主成分分析

```r
vsd = vst(dds)
plotPCA(vsd, intgroup=c('condition', 'type'))   # 这里condition,type都是分组信息
# returnData可以返回主成分分析数据,进而使用其他方法画图
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```


## FAQs

### 分析成对的样品

> 使用多因子设计即可。

### 含有多个组，共同分析还是分开分析？

> 共同分析，两两比较时，使用result函数的contrast指定要比较的两个

### 可以分析无重复（without replicate）的数据吗？

> 可以，但是只推荐用于探索性的分析，不适用于比较分析。
