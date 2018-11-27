---
title: 加权共表达基因网络分析 WGCNA
date: 2018-10-12
toc: true
description: 呦呦呦
categories: Bioinformatics
tags:
  - 软件和包
  - 富集分析
---

在拿到高通量测序数据以后(特别是基因表达数据),通常要分析基因与表型之间的相关性,以探究基因对表型所起关键的调节功能；而加权基因共表达网络分析(Weighted Gene Coexpression Network Analysis)就是其中比较实用的一种分析方法。那么,什么是加权基因共表达网络分析呢？要了解这个,我们需要对以下前提有所了解。

+ 功能相关的基因,其表达水平/表达模式也基本上是相似的
+ 表达水平高度相关的基因具有潜在的共有调控机制或参与相似的生物学过程
+ 如果得到某个表达模式中涉及到许多基因,其中有些基因是已知的,基于上面两点,那么就可以通过该已知基因的功能推测未知基因的功能

根据上面三点,我们就要得到有许多基因涉及的表达模式,而这个可以通过对基因的表达数据进行聚类得到.讲到这里,WGCNA实质上也是一种聚类方法.而在WGCNA之前,现在也在使用的另一种基因与表型的分析方法,那就是基因共表达网络分析.而WGCNA相较于后者来讲,有何优势?

<!--more-->

## 一、总述
### 1.1 相较于GCNA,WGCNA的优势
1. 对Pearson coefficient的绝对值取β次幂构建相关性  
这样做的优势在于,如果两个基因相关性较高的话,取幂后相关性更高;相关性较低,取幂后相关性更低.如此,就更加凸显两个基因的相关性.
2. 相关性阈值的软设定  
使用相关性软阈值设定的方法,进一步涵盖潜在的显著相关性基因;而非一刀切,因为你无法判定0.89与0.90之间,前者就相关性不显著,后者就显著相关.如果一刀切的话,你的结果则依据这个相关性阈值的设定不同,而显著不同
3. 以TOM进一步表征基因的相关性  
在基因两两相关而构建的网络（图）中，通过判断两个基因与其他基因相关性的相似程度来确定这两个基因的实际相关性（topological overlap），以此聚类得出的分组更符合前面所述的前提。

### 1.2 相较于差异分析,WGCNA优势

+ 从更高层次理解基因表达与表型之间的因果关系  
相较于简单比较不同处理中基因表达的差异，WGCNA立足于基因模块（多个基因）的共有效应，更加符合实际情况，因为生命活动就是多中因素共同调节的结果。

### 1.3 基本术语
术语 | 解析
---|---
共表达网络(co-expression network) | 无向的、加权的表达网络。无符号的a(i,j) =|cor(xi, xj)|^β,有符号的a(i,j) = |(1 + cor(xi, xj))/2|^β
模块(module) | 高度相关的基因集
连接性(connectivity) | 网络中某个基因与其他基因连接值的总和
模块内部连接性(intramodular connectivity) | 模块中某个基因与同模块其他基因的连接值总和
模块特征(module eigengene) | 给定模块的第一个主成分
特征显著性(eigengene significance) | 模块与表型之间的相关系数
模块关系(基于模块特征的连接性) | 给定模块的模块特征与基因表达的相关性
关键基因(hub gene) | 给定模块中具有较大连接性的基因
基因显著性(gene significance) | 基因表达的显著性(-log(p-value))
模块显著性(module significant) | 给定模块中所有基因显著性的绝对值的平均值

### 1.4 基本过程

1. 计算基因的相关性(Pearson correlation的β次幂，connectivity, K)
2. 依据connectivity计算两两基因的topological overlap（relative connectivity）
3. 根据拓扑重叠值构建基因表达网络
4. 对此进行层次聚类，确定共表达模块（表达模式）
5. 选取感兴趣的模块进行后续分析：模块的功能富集、模块与性状的相关性、模块与样本间的相关系数
6. 挖掘模块的关键信息：寻找模块核心基因（hub gene)、预测基因功能

## 二、分析过程


## 2.1 导入数据、预处理

这里主要有两个目的，一个是检测数据集是否都满足WGCNA的要求，二是排除离群值，后者需要你根据数据产生的聚类树，选择合适的cutHeight来去除离群值。


```r
expr = read.table('./clean_data/expr_matrix.txt')
gsg = goodSamplesGenes(expr)
if (!gsg$allOK) expr = expr[gsg$goodSamples, gsg$goodGenes]
sampleTree = hclust(as.dist(expr), method='average')

# [ Visualization 1 ], 根据画出的图确定cutHeight,这里假设是120
pdf(file='visual_01_sampleTree_with_outlier.pdf')
plot(sampleTree,
	main='Sample Clustering to detect outliers',
	xlab='', sub='')
abline(h=120, col='red')
dev.off()

cutHeight = 120
sample_cut = cutreeStatic(sampleTree, cutHeight=cutHeight)
expr = expr[sample_cut == 1, ]
sampleTreeCut = hclust(as.dist(expr), method='average')

# [ Visualization 2], 可视化样本聚类树和表型数据
traits = read.table(file='./clean_data/traits.txt')
pdf(file='visual_02_sampleCutTree_with_traits.pdf')
colors = numbers2colors(traits)
plotDendroAndColors(sampleTreeCut, colors,
					groupLabels=colnames(traits),
                    main='Sample Clustering and traits')
dev.off()

save(sampleTree, sampleTreeCut, expr, file='01_sampleTree_cut_outlier.RData')
```

## 2.2 获取软阈值

这里的软阈值选取power significance 大于0.8或0.85,0.9的power值。一般来说都在7-10之间。如果数据较好，可以直接通过"beta=soft$EstiamteSoft"获取，否则通过可视化获取。如果数据展示出的scale-free network不太好，那就多选取一些数据，而非仅仅基于有差异表达的基因表达矩阵。


```r
powers = 1:20
powerSigR2 = 0.8
soft = pickSoftThreshold(expr, powerVectors=powers, verbose=3)
beta = soft$EstimateSoft

# [ visualizaiton 3], 可视化软阈值的选取
pdf(file='visual_03_soft_pick_threshold.pdf')
par(mfrow=c(1,2))
plot(soft$fitIndices$Power,  -sign(soft$fitIndices$slope) * soft$fitIndices$SFT.R.sq, 
     xlab='Soft Threshold (power)',
     ylab='Scale Free Topology Model Fit, signed R^2', type='n',
     main=paste('Scale Independce'))
text(soft$fitIndices$Power, -sign(soft$fitIndices$slope) * soft$fitIndices$SFT.R.sq,
     labels=soft$fitIndices$Power, cex=0.9, col='red')
abline(h=0.8, col='blue')
plot(soft$fitIndices$Power, soft$fitIndices$mean.k.,
     xlab='Soft Threshold (power)',
     ylab='Mean connectivity', type='n',
     main=paste('Mean connectivity'))
text(soft$fitIndices$Power, soft$fitIndices$mean.k.,
     labels=soft$fitIndices$Power, cex=0.9, col='red')
dev.off()

save(soft, beta, file='02_soft_threshold.RData')
```

## 2.3 构建网络

构建网络有一步法和多步法，由于后续分析和可视化的需要，这里使用多步法。在我们通过拓扑网络计算节点距离进行聚类后，还可以根据聚类后的合并模块特征靠近的或者较小的模块。

### 2.3.3 计算拓扑网络和内部连接性，基于此构建聚类树

```r
softConnect = softConnectivity(expr, power=beta)-1
adjacency = adjacency(expr, power=beta)
TOM = TOMsimilarity(adjacency) # 一步到位 TOM = TOMsimilarityFromExpr(expr,beta)
dissTOM = 1 - TOM # dissTOM = TOMdist(adjacency)
geneTree = hclust(as.dist(dissTOM), method='average')
# [visualization 4], 可视化拓扑图
	pdf(file='visual_04_scale_free_topology_plot.pdf')
	scaleFreePlot(softConnect, truncated=F,
		main=paste('Scale Free Topology, power=', beta,sep=''))
	dev.off()
# [visualization 5], 可视化拓扑网络聚类树
	pdf(file='visual_05_raw_gene_tree_on_dissTOM.pdf')
	plot(geneTree, xlab='', sub='',
		    main='Gene Clustering on TOM-based dissimilarity',
		    labels=F, hang=0.04)
	dev.off()
# [visualization 5.1], 可视化模块与模块之间的拓扑相异性，揭示基于拓扑重叠进行聚类的效果
	pdf(file='visual_05_dissTOM_between_module.pdf')
	TOMplot(dissTOM, geneTree)
	dev.off()
save(softConnect, adjacency, TOM, geneTree, file='03_conn_adj_TOM.RData')
rm(c(soft, softConnect))
```

### 2.3.2 把拓扑网络聚类树进行动态修剪，产生初始模块

```r
minModuleSize=30
dissTOM_cut = cutreeDynamic(dendro=geneTree, distM=dissTOM,
                            deepSplit=2, pamRespectsDendro=F,
                            minClusterSize=minModuleSize)
dynamicColors = labels2colors(dissTOM_cut)
# [visualization 6] 可视化拓扑重叠相异性的多维尺度图(mds-plot)
	pdf(file='visual_06_MDS_plot_on_dissTOM_cut.pdf')
	dimension_num = 2
	cmd = cmdscale(as.dist(dissTOM_cut),dimensional_num)
	plot(cmd, col=as.character(dynamicColors), main='MDS plot',
		xlab='Scaling Dimension 1', ylab='Scaling Dimension 2')
	dev.off()
```
### 2.3.3 计算模块特征，合并特征相似或靠近的模块，产生最终模块

```r
MEs = moduleEigengenes(expr, dynamicColors)$eigengenes
dissMEs = 1 - cor(MEs)
dissMEs_tree = hclust(as.dist(dissMEs), method='average')
MEDissThreshold = 0.25
mergedMEs = mergeCloseModules(expr, dynamicColors,
							cutHeight=MEDissThreshold, verbose=3)
# [visualization 6], 可视化拓扑网络聚类树的动态修剪和合并产生的模块比较
	pdf(file='visual_06_raw_tree_with_dynamicCut_merge.pdf')
	mergedColors = mergedMEs$colors
	plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
							c('Dynamic Tree cut', 'Dynamic Merged'),
							dendroLabels=F, hang=0.03,
	                    	addGuide=T, guideHang=0.05,
	                    	main='Gene Dendrogram with Dynamic cut & Merge')
	dev.off()
# [visualization 7], 可视化模块特征聚类图
	pdf(file='visual_07_module_eigengenes_cluster.pdf')
	plot(dissMEs_tree, main='Module Eigengenes Clustering')
	dev.off()
save(MEs, dissMEs_tree, mergedMEs, file='04_rawMEs_dissMEs_tree_mergedMEs.RData')
rm(c(MEs, dissMEs_tree))
```

## 2.4 关联分析

从“2.1-2.3”是我们进行WGCNA的分析过程，下面是一些具体的关联分析。所谓关联分析就是把我们得到的最终的聚类模块同其他数据如表型数据计算相关性，看有哪些模块与这些数据相关，那么模块中的基因就与这些数据相关。下面是一些示例。


```r
# 导入相关数据
traits = read.table('./clean_data/traits.txt')
# expr 如果不存在就load('01_sampleTree_cut_outlier.RData')
# mergedMEs如果不存在，就load('04_rawMEs_dissMEs_tree_mergedMEs.RData')

# 重新计算最终的模块特征, 并据此排序
mergedColors = mergedMEs$colors
MEs = moduleEigengenes(expr, mergedColors)$eigengenes
MEs = orderMEs(MEs)
nGenes = ncol(expr)
nSamples = nrow(expr)
modNames = substring(names(MEs), 3)
moduleColors = mergedColors
```

### 2.4.1 模块与模块相关性


```r
# 方法一
# plotMEpairs(MEs)
# 方法二
moduleModuleCor = cor(MEs, MEs, use='p')
moduleModulePvalue = corPvalueStudent(moduleModuleCor, nSamples)
moduleModuleText = paste0(signif(moduleModuleCor,2), '\n(',
						signif(moduleModulePvalue,2), ')')
dim(moduleModuleText) = dim(moduleModuleCor)
# [visualization 8] 模型与模块之间的相关性
pdf(file='visual_09_module_module_relationship.pdf')
labeledHeatmap(Matrix=moduleModuleCor, xLabels=names(MEs),
				yLabels=names(MEs), ySymbols=names(MEs),
				xSymbols=names(MEs), colorLabels=F,
				colors=greenWhiteRed(50), textMatrix=moduleModuleText,
				setStdMargins=F, zlim=c(-1,1),
				main=paste('Relationships between Module Eigengenes'))
dev.off()
MMdata = cbind(moduleModuleCor, moduleModulePvalue)
colnames(MMdata) = c(paste0(colnames(moduleModuleCor), '.Cor'),
					paste0(colnames(moduleModulePvalue), '.Pvalue'))
write.csv(MMdata, file='WGCNA_02_Modue_Module_relationship.csv')
```
### 2.4.2 模块与表型之间的相关性


```r
moduleTraitCor = cor(MEs, traits, use='p')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitText = paste0(signif(moduleTraitCor,2), '\n(',
						signif(moduleTraitPvalue,2), ')')
dim(moduleTraitText) = dim(moduleTraitCor)
# [visualization 9] 模型与表型之间的相关性
pdf(file='visual_08_module_traits_relationship.pdf')
labeledHeatmap(Matrix=moduleTraitCor, xLabels=names(traits),
				yLabels=names(MEs), ySymbols=names(MEs),
				colorLabels=F, colors=greenWhiteRed(50),
				textMatrix=moduleTraitText, setStdMargins=F,
				zlim=c(-1,1),
				main=paste('Relationships between module eigengenes and traits '))
dev.off()
MTdata = cbind(moduleTraitCor, moduleTraitPvalue)
colnames(MTdata) = c(paste0(colnames(moduleTraitCor), '.Cor'),
					paste0(colnames(moduleTraitPvalue), '.Pvalue'))
write.csv(MTdata, file='WGCNA_01_Modue_Traits_relationship.csv')
```

### 2.4.3 模块与基因表达量相关性

通过检测模块与表型之间的相关性，我们可以得到高度相关的模块与表型。使用这个模块和表型，获取要分析的基因。假定我们的模块是“blue”和表型是“weight”。


```r
geneModuleCor = cor(expr, MEs, use='p') 
# modNames = substring(names(MEs), 3)
geneModulePvalue = corPvalueStudent(geneModuleCor, nSamples)
names(geneModuleCor) = paste(modNames, '.Cor', sep='')
names(geneModulePvalue) = paste(modNames, '.Cor.pvalue', sep='')
write.csv(cbind(geneModuleCor, geneModulePvalue),
	file='WGCNA_03_gene_module_relationship.csv')

geneTraitCor = cor(expr, traits, use='p')
geneTraitPvalue = corPvalueStudent(geneTraitCor, nSamples)
names(geneTraitCor) = paste(colnames(traits), '.Cor', sep='')
names(geneTraitPvalue) = paste(colnames(traits), '.Cor.pvalue', sep='')
write.csv(cbind(geneTraitCor, geneTraitPvalue),
	file='WGCNA_04_gene_trait_relationship.csv')

module_aim = 'blue'
trait_aim = 'weight'
module_column = match(module_aim, modNames)
trait_column = match(trait_aim, colnames(traits))
geneModuleCorAim = geneModuleCor[moduleColors == module_aim, module_column]
geneTraitCorAim = geneTraitCor[moduleColors == module_aim, trait_column]
# [visualization 10] 模块与基因表达量的相关性
pdf(file=paste('visual_10_gene_significance_MM_in_', module_aim, '.pdf', sep=''))
# 这里应该也可以用labledHeatmap
verboseScatterplot(abs(geneModuleCorAim),
					abs(geneTraitCorAim),
					xlab=paste('Module Membership in', module_aim, ' module'),
					ylab=paste('Gene significance for', trait_aim),
					main=paste('Module membership (', module_aim, ') vs Gene significance (', trait_aim, ')'),
					col = module_aim)
dev.off()
```

### 2.4.4 模块与模块内总体基因显著值

“we define a gene significance variable as minus log10 of the univarite Cox regression pvalue for predicting survival on the basis of the gene epxression info“, 有文献定义了个gene significance, 然后可视化各个模块内该gene significance的分布情况。这说gen significance的分布情况。这说明gene significance是由我们自己定义的。我们这里使用每个基因的平均表达量作为其基因显著特征


```r
exprPerGene = colMeans(expr)
# [visualizaton 11] 模块与总体基因表达量柱状图
pdf(file='visual_11_total_geneSig_in_module.pdf')
plotModuleSignificance(exprPerGene, moduleColors, main="Module significance")
dev.off()
```

### 2.4.5 各个模块内的基因表达模式图


```r
pdf(file=paste('visual_12_gene_expression_pattern_in_all_modules.pdf', sep=''), w=12,h=9)
for(moduleName in modNames) {
    column = match(moduleName, modNames)
    sizeGrWindow(12,9)
    par(mfrow=c(2,1))
    par(mar=c(1,3,6,3))
    plotMat(t(expr[, moduleColors == moduleName]),
            nrgcols=30, clabels=rownames(traits));
    par(mar=c(2,3,2,1))
    barplot(MEs[, column], xlab = moduleName, ylab='Gene Expression pattern', col=moduleName)
}
dev.off()
```

### 2.4.6 模块内部基因连接性与基因显著的相关性

这里的基因显著定义还是与2.4.4是一样的。


```r
exprPerGene = colMeans(expr)
ConnectivityMeasure = intramodularConnectivity(adjacency, moduleColors)
pdf(file=paste('visual_13_intramodularConnect_gene_significance.pdf', sep=''))
for(moduleName in modNames) {
	verboseScatterplot(ConnectivityMeasure$kWithin[moduleName == moduleColors],
						exprPerGene[moduleName == moduleColors],
						col=moduleColors[moduleColors == moduleName],
						main=paste('module ', moduleName, sep=''),
						ylab='Gene significance',
						xlab='Intramodular K')
}
dev.off()
```

## 2.5 导出数据


```r
load('01_sampleTree_cut_outlier.RData')
load('03_conn_adj_TOM.RData')
load('04_rawMEs_dissMEs_tree_mergedMEs.RData')

# softConnect, adjacency, TOM, geneTree,MEs, dissMEs_tree, mergedMEs
probes = colnames(expr)
# 想要导出的模块
modules = c('green', 'black')
inModule = is.finite(match(mergedColors, modules))
modTOM = TOM[inModule, inModule]
modProbes = probes[inModule]
dimnames(modTOM) = list(modProbes, modProbes)
IMConn = softConnectivity(expr[, modProbes]);

nTop = 30;
top = (rank(-IMConn) <= nTop)

cyt = exportNetworkToCytoscape(modTOM[top, top],
        edgeFile=paste('CytoscapeInput-edges-', paste(modules, collapse='-'), '.txt', sep=''),
        nodeFile=paste('CytoscapeInput-nodes-', paste(modules, collapse='-'), '.txt', sep=''),
        weighted=T,
        threshold = 0.02,
        nodeNames = modProbes,
        nodeAttr=mergedColors[inModule])

vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )
```





































