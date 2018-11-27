---
title: 差异分析 edgeR
date: 2018-10-12
toc: true
categories: Bioinformatics
tags:
  - 软件和包
  - 差异分析
---

下方的代码把edgeR的三种差异分析整合成一个函数, 调用时直接指定参数即可.

+ classcial
+ glm: likelihood ratio test/ quasi-likelihood F-test
+ + quasi-likelihood(qlf): 推荐用于差异分析，因为他对错误率限制较好。
+ + likelihood(lrt)：对与单细胞RNA-测序和没有重复的数据较好

<!--more-->


## 调用代码

```r
suppressPackageStartupMessages(library(edgeR))
setwd('~/practice/180716_edgeR/')

rawdata = read.table('rawdata.txt')
head(rawdata)
rawdata = rawdata[-(1:5),]

groups = grepl('01A', colnames(rawdata))
groups = ifelse(groups, 'tumor', 'normal')
table(groups)

edgeR_DGE(exprSet = rawdata, group=groups, type='classical')
```

## 定义代码

```r
edgeR_DGE <- function(exprSet, group, type, cpm=c(100,4)) {
    # 生成DEGList
    DGE = DGEList(counts=exprSet, group=group)
    DGE.old = DGE
    # 生成design
    design = model.matrix(~group)
    # cpm过滤
    keep = rowSums(edgeR::cpm(DGE) > cpm[1]) >= cpm[2]
    DGE = DGE[keep,]
    # 进行校正
    DGE = calcNormFactors(DGE)
    # 检测离群值和关系
    png('plotMDS.png')
    plotMDS(DGE, method='bcv', col=as.numeric(DGE$samples$group))
    legendCol = unique(as.numeric(DGE$samples$group))
    legendGroup = unique(as.character(DGE$samples$group))
    legend("bottomright", legendGroup, col=legendCol, pch=20)
    dev.off()
    if (type == 'classical') {
        # 计算离散度dispersion
        d = estimateCommonDisp(DGE)
        d = estimateTagwiseDisp(d)
        test = exactTest(d)
    } else {
        # 计算离散度dispersion
        d = estimateGLMCommonDisp(DGE)
        d = estimateGLMTrendedDisp(d)
        d = estimateGLMTagwiseDisp(d)
        if (type == 'qlf') {
            fit = glmQLFit(d, design)
            test = glmQLFTest(fit, coef=2)
        } else if (type == 'lrt') {
            fit = glmFit(d, design)
            test = glmLRT(fit, coef=2)
        }
    }
    png('plotBCV.png')
    plotBCV(d)
    dev.off()
    png('plotSmear.png')
    de = decideTestsDGE(test, adjust.method="BH", p.value = 0.05)
    tags = rownames(d)[as.logical(de)]
    plotSmear(test, de.tags=tags)
    abline(h=c(-4,4), col='blue')
    dev.off()
    finalDGE = topTags(test, n=nrow(exprSet))
    finalDGE = as.data.frame(finalDGE)
    write.csv(file='DGE_edgeR.txt', finalDGE, quote = F)
}
```
