---
title: 芯片专题 AffyMetrix
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 生信基础
---

## 一、芯片基础
### 1.1 芯片平台和包
平台 | 包
---|---
illumina | beadarray, lumi
affymatrix | affy,simpleaffy, oligo
agilent | 没有找到

<!--more-->

### 1.2 bioconductor包内置芯片数据

> http://www.bio-info-trainee.com/1399.html

### 1.3 GEO数据库

+ GEO platform: GPL
+ GEO sample:	GSM
+ GEO series:	GSE
+ GEO dataset:	GDS

处理GEO数据有两条思路,第一条速录是利用GEOquery来下载。另一条则是自己下载原始的芯片数据，然后转换成表达矩阵.当然，有些比较常用的芯片则是有对应的包,我们只需要把对应包进行安装就有相应的芯片数据了。

#### 1.3.1 GEOquery

getGEO(GEOID, destdir):
+ GDS858: 	下载数据集
+ GPL96：	下载芯片设计信息,probe-geneID对应表(soft文件)
+ GSE1009:	下载表达矩阵,在下载该数据时GPL96.soft也会自动下载,但无法读取.得用GPL96作为ID使用getGEO函数.这时候的getGPL干脆设为False

getGEOfile:		 下载GEO soft文件
getGEOSuppFiles: 下载原始数据

```
# 使用GDS858的ID下载
GDS = getGEO('GDS858', destdir='.')
expr_matrix = Table(GDS)	#得到表达矩阵
Des_info = Meta(GDS)		#得到描述信息(metadata)
eSet = GDS2eSet(GDS, do.log2=T)	#得到expressionSet,做log2变换

# 使用GSE1009的ID下载
# 返回的是含expressionSet对象的list
genes = geneNames(GSE[[1]])		# 得到gene名, geneNames是biobase包里面的
samples = sampleNames(GSE[[1]])	# 得到样品名
pdata = pData(GSE[[1]])			# 得到描述信息(metadata)
expr_matrix = exprs(GSE[[1]])	# 得到表达矩阵

# GPL16699的ID下载
Des_info = Meta(GPL)		# 得到描述信息
Annotation = Table(GPL)		# 得到芯片注释信息(geneID-probeID)

# 下载原始数据信息
rawdata = getGEOSuppFiles(GSE1009)

```
#### 1.3.2 GEO页面

![P3GY5j.png](https://s1.ax1x.com/2018/07/19/P3GY5j.png)

### 1.4 原始数据的预处理
1. 背景信号处理
由于图像处理软件对芯片划格后,取杂交点周围区域内每个像素的吸光度均值作为背景信号,但是芯片的不同区域像素不同则均值不同,各个区域扣除的背景则多少不一.所以要进行背景信号处理. 可以利用芯片最低信号轻度的点作为背景信号,或者把芯片非杂交点的背景吸光度均值作为整个芯片的背景信号.

2. 芯片数据清理
经背景校正后的芯片数据会出现异常值,通常的处理方法就是去除这些异常值,去除的方法有:标准值法, 奇异值法, 变异系数法, 前景值法(front value <200), 中位数法等等.Affymetrix的芯片则时直接将负值修正为一个固定值.

紧接着,对于缺失值的处理,方法有多种.一个时删除, 如果某行或某列的缺失值多于某个阈值,则删除该行/列. 另一种是补全.补全值的选择有使用0值补全,使用均值或中值补全,还可以时K邻近法进行补全.

3. 归一化
由于芯片数据值变化很大, 呈偏态,标准差大.但是后续的差异分析基本都是假定数据的正态分布,所以可以使用对数转换使数据趋于正态分布, 减少标准差.然后是消除技术误差/偏差.这时候使用的方法有平均数标准化, 中位数标准化等方法.此时就得到了表达矩阵.

4. 差异分析
芯片数据的差异分析有三种方法：
+ 倍数分析法(fold change)
+ 参数检验(t-test)
+ 非参数检验(经验贝叶斯,SAM法)

常用的软件有:limma、DESeq2、edgeR、GFOLD等。


## 二、使用Affy包分析芯片质量

```r
source('http://biconductor.org/biocLite.R')
biocLite(c('affy', 'simpleaffy'))
```
### 2.1 读取CEL文件,查看基本信息

```r
suppressPackageStartupMessages({
    library(affy))
    library(simpleaffy))

###################读取文件####################
filenames = list.files('./GSE30668', pattern='*.CEL')
the.data = ReadAffy(filenames=paste0('./GSE30668/', filenames, sep=''))

##############芯片基本信息####################
old.names = sampleNames(the.data) # 芯片名
sampleNames(the.data) = gsub('.CEL', '', filenames)

# 完美匹配探针,错配探针; perfect match probes / mismatch probes
pm.data = pm(the.data)
mm.data = mm(the.data)
head(pm.data)
head(mm.data)

pdata = pData(the.data)
probes = geneNames(the.data)
```

### 2.2 查看芯片质量
+ 灰度图、灰度值

```r
nmatrix = length(filenames)

#################### 画灰度图,并统计灰度值 ######################
par(mfrow=c(ceiling(nmatrix/2), 2))
par(mar=c(0.2,0.2,0.2,0.2))
# 设置调色板为灰度
pallette.gray = c(grep(gray(0:10/10), times=seq(1, 41,by=4)))

for (i in 1:nmatrix) {
	image(the.data[, i], col=pallette.gray)
	# 如果芯片图像有版斑块就有可能时坏片
}
par(mfrow=c(1,1))
par(mar=c(4,4,3,0.5))
par(cex=0.7)
if (nmatrix > 40) par(cex=0.5)
colors = rainbow(nmatrix*1.2)
boxplot(the.data, col=colors, xlab='Sample', ylab='Log intensity')
```
+ 曲线图，MA-plot，RNA降解分析

```r
##################### 画曲线图 ######################
par(mar=c(4,4,3,0.5))
hist(the.data, lty=1:3, col=colors)
legend('topright', legend=sampleNames(the.data),
		lty=1:3, col=colors, 
		box.col='white', xpd=T,
		cex=1.5)
box()

####################### MA-plot #####################
# 作图得到IQR,如果IQR较大的话芯片可能有问题
par(mfrow=c(ceiling(nmatrix/2), 2))
par(mar=c(3,3,2,0.5))
par(tcl=0.2)
par(mgp=c(2,0.5,0))
MAplot(the.data, cex=0.8)

####################### RNA降解分析 #################
# 理想情况下每个样品的线(分段)是平行的
par(mfrow=c(1,1))
par(mar=c(4,4,3,0.5))
RNAdeg = AffyRNAdeg(the.data)
summaryAffyRNAdeg(RNAdeg)
plotAffyRNAdeg(RNAdeg, cols=colors)
legend('topleft', legend=sampleNames(the.data),
		lty=1, col=colors,
		box.col='white', xpd=T)
box()
```
+ 背景信号值、样品缩放因子、阳性信号值（表达基因占比，内参占比）

```r
suppressPackageStartupMessages(library(simpleaffy))

################# 查看背景信号值 ####################
# 如果背景信号值太大,说明芯片结果并不好(太大是多大?)
qc.data = qc(the.data)
avbg.data = as.data.frame(sort(avbg(qc.data)))
max(avbg.data)

################# 样品scale factor###################
# 如果scale factor的最大值与最小值的比值超过3倍的话
# 芯片结果也不好
sfs.data = sort(sfs(qc.data))
ratio = max(sfs.data) / min(sfs.data)
ratio

################# 有表达的基因占比多少 ###############
# 表达基因占比太小表明芯片结果有问题
per.data = as.data.frame(percent.present(qc.data))
min(per.data)

################# 内参基因的表达占比 #################
rat.data = ratios(qc.data)
head(rat.data)
```

## 三、使用Affy包处理原始芯片数据

从原始芯片数据开始处理需要经过四个步骤：背景校正、归一化(标准化)、计算表达量、得到表达矩阵。你也可以通过使用mas5或rma或gcrma函数直接从读入的原始芯片数据(见2.1, the.data)计算表达量(见3.3),最后得到表达矩阵.

### 3.1.背景校正

```r
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(simpleaffy))
load('GSE30668_eSet.RData')
################ rma, mas 两种校正方法###########################
# mas校正后pm和mm的值都被重新计算了
# rma校正仅针对pm,mm的值不会变化
rma.data = bg.correct(the.data, method='rma')
mas.data = bg.correct(the.data, method='mas')
crt.data = list(rma.data, mas.data)
```

### 3.2.归一化处理
芯片的归一化方法有线性缩放法、非线性缩放法、分位数法、Cyclic loess法、Contrasts法。


```r
#这里取rma背景校正的结果进行后续的归一化.
data = crt.data[[1]]
############## 线性缩放 ############################
# 选择某个芯片作为参考,其他芯片与其分别做线性回归
# 得到的回归线对剩余其他芯片做缩放(同一缩放倍数)

ln.data = normalize(data, method='constant')
head(pm(ln.data)/pm(data))
head(mm(ln.data)/mm(data))
############# 非线性缩放 ###########################
# 仅使用部分芯片结果做非线性拟合,然后进行缩放(不同缩放倍数)
invar.data = normalize(data, method='invariantest')
head(pm(invar.data)/pm(data))
head(mm(invar.data)/mm(data))
############# 分位数(quantile) #####################
# 前提是假设每个芯片信号的经验分布函数一样
# 任两个芯片的QQ-plot会形成45度对角线
qt.data = normalize(data, method='quantiles')
head(pm(qt.data)/pm(data))
head(mm(qt.data)/mm(data))
############# Cyclic loess #########################
loe.data = normalize(data, method='loess')
############# contrasts ############################
cont.data = normalize(data, method='contrasts')
```

### 3.3 计算表达量

+ 使用affy::computeExprSet计算表达矩阵,它需要指定统计方法和PM校正方式
```
################# 显示PM校正方式种类 #######################
pmcorrect.methods()
################# 显示统计方法种类 ##########################
generateExprSet.methods()
################# 计算表达矩阵 ######################
eSet.pmo.liw = computeExprSet(ln.data, pmcorrect.method='pmonly', summary.method='liwong')
```
+ 使用mas5或rma或gcrma函数一步到位
```
eSet.mas5 = mas5(the.data)
eSet.rma = rma(the.data)
suppressPackageStartupMessages(library(gcrma))
eset.gcrma = gcrma(the.data)
```
### 3.4 得到表达矩阵

```r
############### 进行rma或mas5处理 ###############
#eSet.mas5 = mas5(the.data)
#eSet.rma = rma(the.data)
############### 获取表达矩阵 ################
# 注意，rma法计算表达量
expr.rma.lg2 = exprs(eSet.rma)
expr.rma.nlg = 2^expr.rma.lg2
expr.mas5.nlg = exprs(eSet.mas5)
expr.mas5.lg2 = log2(expr.mas5.nlg)

############### 筛选有表达的基因 ################
# mas5calls: P为present, A为absent, M为marginal
eSet.calls = mas5calls(the.data)
expr.calls = exprs(eSet.calls)
probes = apply(expr.calls, 1, function(x) any(x=='P'))
# 以rma为例
expr.chose = expr.rma.lg2[probes, ]
save(expr.chose, file='GEO30668_expr_rma_chose.RData')
```


### 3.5 差异分析和注释
差异分析可以使用limma直接对表达矩阵进行。注释的话，需要通过GEOquery得到的GPL文件中获取探针ID对应的基因ID。然后可以进行其他分析如富集分析、基因关联分析。
