---
title: 差异分析 limma
date: 2018-10-12
toc: true
categories: Bioinformatics
tags:
  - 软件和包
  - 差异分析
---

limma最开始是用于芯片数据分析的,不过现在也支持RNA-seq等数据的差异分析，但是需要通过voom函数进行校正表达矩阵。

<!--more-->

## 分析芯片数据

在limma芯片数据分析的过程中,有需要注意的一个地方,那就是分组矩阵和差异比较矩阵的问题,使用差异比较矩阵时的代码是有不一样的,结果仍旧相同,当然你要指定好makeContrasts的参数咯.

使用topTable时,coef是用来指定提取特定比较的结果.例如你在design的时候,1-2的比较,coef可以设置为'1-2'.如果不是用差异比较矩阵,coef=1就是分组矩阵的column1-column2,以此类推

### 不使用差异比较矩阵

```r
library(limma)
groups = sort(rep(1:3, 4))
names(groups) = level(factor(groups))
design = model.matrix(~factor(groups))
fit.lm = lmfit(expr, design)
fit.ebs = eBayes(fit.lm)
topTable(fit.ebs, adjust='BH', coef=2, lfc=1, p.value=0.05, numbers=30000)
# 输出所有的差异分析结果
results = desideTests(fit.ebs)
```

### 使用差异比较矩阵

```r
groups = sort(rep(1:3, 4))
design = model.matrix(~ 0 + factor(groups))
colnames(design) = c('1','2','3')
fit.lm = lmfit(expr, design)
contrast = makeContrasts('1-2', '1-3', '2-3', levels=design)
fit.cts = contrasts.fit(fit.lm, contrast)
fit.ebs = eBayes(fit.cts)
# 选取top的基因查看
topTable(fit.ebs, adjust='BH', coef=2, lfc=1, p.value=0.05, numbers=30000)
# 输出所有的差异分析结果
results = desideTests(fit.ebs)
```

## 分析RNA-seq数据

用于差异分析的RNA-seq数据必需是 **raw_counts** 数据.


```r
suppressMessages(library(limma))
# 构建分组矩阵
expr_matrix = your_expr_matrix
groups = your_groups
design = model.matrix(~factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = colnames(expr_matrix)

# 进行校正
v = voom(expr_matrix, design, normalize='quantile')

# 进行差异分析
fit.lm = lmFit(v, design)
fit.ebs = eBayes(fit.lm)
DE.genes = topTable(fit.ebs, coef=2, n=Inf, lfc=1, p.value=0.05)
```

## 结果说明
+ topTable用于提取top基因和对应的contrast
+ logFC是log fold change
+ aveExpr是average log2-expression level for gene in all array
+ pvalue/adj.p.value是p-value和校正后的p-value
+ t 是 moderated t-statistic
+ B 假设B=1.5,则exp(1.5)=4.48,因此基因出现差异表达的可能性为4.48/(1+4.48)=0.8
