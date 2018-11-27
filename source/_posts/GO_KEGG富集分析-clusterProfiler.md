---
title: GO_KEGG富集分析 clusterProfiler
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 软件和包
---

首先呢,要详细了解的话,需要看这篇文献(Ten Years of Pathway Analysis: Current Approaches and Outstanding Challenges),他把基本的信号通路分析方法进行了总结.

<!--more-->

## 1.对Ensembl ID进行转换,得到对应的基因名和EntrezID

不同的物种要选择不同的数据库,不然选择不上.

+ 斑马鱼: [org.Dr.eg.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html)
+ 拟南芥: [org.At.tair.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)
+ 小鼠: [org.Mm.eg.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html)
+ 大鼠: [org.Rn.eg.db](https://www.bioconductor.org/packages/release/data/annotation/html/org.Rn.eg.db.html)
+ 人类: [org.Hs.eg.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)


```r
suppressMessages(library(org.Hs.eg.db))
keytypes(org.Hs.eg.db)

# ensids的ID类型与select函数的keytype要一致.如果是ENTREZID,就不用执行这一步.
ensids = c('400212', '238240', '238204')    
cols <- c("SYMBOL", "GENENAME", 'ENTREZID')
gene = select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="ENTREZID")
```

## 2. 进行注释

同样的,OrgDb的参数数据库, 不同的物种要选择不同的数据库.

enrichKEGG的organism的参数要符合 http://www.genome.jp/kegg/catalog/org_list.html 所列.

### GO & KEGG

```r
library(clusterProfiler)

# 1.cell category
ego_CC <- enrichGO(gene = gene$ENTREZID, OrgDb= org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH",
                   minGSSize = 1, pvalueCutoff = 0.01, qvalueCutoff = 0.01, readable = TRUE)
# 2.biological progress
ego_BP <- enrichGO(gene = gene$ENTREZID, OrgDb= org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", minGSSize = 1,
                   pvalueCutoff = 0.01, qvalueCutoff = 0.01, readable = TRUE)
# 3.molecular function
ego_MF <- enrichGO(gene = gene$ENTREZID, OrgDb= org.Hs.eg.db, ont = "MF",
                   pAdjustMethod = "BH", minGSSize = 1, pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01, readable = TRUE)

# 4.KEGG
kk <- enrichKEGG(gene = gene$ENTREZID, organism ="human", pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01, minGSSize = 1,
                 #readable = TRUE, 
                 use_internal_data =FALSE)
```

### GSEA富集分析

```r
gse <- gseGO(gene = gene$ENTREZID, ont = "BP", 
            OrgDb = org.Hs.eg.db, keyType = "ENTREZID", exponent = 1,
            nPerm = 1000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
            pAdjustMethod = "BH", verbose = TRUE, seed = FALSE, by = "fgsea")
```

## 3.进行可视化

```r
results = c(ego_cc, ego_BP, ego_MF, kk)
names = c('GO enrich: CC', 'GO enrich: BP', 'GO_enrich_MF', 'KEGG enrich')
for (i in 1:length(results)) {
    barplot(results[i], showCategory=20,title=names[i])
    dotplot(results[i],title=names[i])
}

# GSEA plot
gseaplot(gse, geneSetID="GO:0004871")
```
