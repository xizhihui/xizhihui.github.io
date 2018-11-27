---
title: ChIP-seq专题 MACS2_ChIPSeeker_deeptools
date: 2018-10-14
toc: true
description: 搞了个基于 snakemake 的 ChIP 分析流程，详见 github.com/xizhiui/bioinformatic_learning/pipelines/ChIP-seq
categories: Bioinformatics
tags: 软件和包
---

ChIP-seq是使用抗体捕获富集DNA片段和高通量测序技术来获得某些marker与DNA的结合位点的一项综合技术。ChIP是染色质免疫共沉淀, 通过特异抗体将DNA结合蛋白免疫沉淀, 用于捕获蛋白质的DNA靶点, 比如转录因子啊, 组蛋白修饰啊. 它主要分为以下四步：cross-linking、sonication、IP、Sequencing。在DNA与蛋白交联以后, 通过超声的方式随机打断染色体, 在利用抗体将目的交联物筛选出来, 再反交联获取DNA,最后上机测序.获取到测序数据后,典型的分析流程如图.

![ChIP分析流程](http://mmbiz.qpic.cn/mmbiz_jpg/MPBFtnFrw4m2jWmOk9ic1ZkZ8VrQBrr8oCMwfEZDZG75OY3ic9qMqHhibpiafgZVfJaCSp7gqI87ETibVltjV4GgXkQ/640?wx_fmt=jpeg&tp=webp&wxfrom=5&wx_lazy=1)

<!--more-->

文字版: rawdata -> QC -> mapping -> peak calling(binding site) -> visualization/annotation/MotifAnalysis.
在peak calling完成之后,我们首先会可视化看一下数据，接下来就会想知道这些peak都和什么样的基因有关. 

## 一、Overview

### 1.1 抗体的种类 
+ 转录因子或抑制因子(transcription factors/repressors)如CTCT
+ 组蛋白或组蛋白修饰(histones/histone modifications)如H3/H3K4me3
+ DNA修饰(DNA modfications)如DNA甲基化
+ 染色质重塑蛋白(chromatin remodelling proteins)如BMI1
+ 其他转录机制相关蛋白(transcription machinery)如Pol2

### 1.2 实验结果简要展示

![ChIP结果简要展示](https://s1.ax1x.com/2018/10/14/iUEKwq.png)

### 1.3 富集的类型

经过抗体富集的DNA可以是单个结合点，也可以是较宽的结合区域，当然也可能是任何地方。

![type of enrichment](https://s1.ax1x.com/2018/10/14/iUElkV.png)

### 1.4. 表征结果的形式

ChIP测的是富集, 而且是相对富集, 就是指A/B两个区域的富集信号的富集差别.如果要得到绝对值的话,那就要进行归一化或者校正.影响富集的因素有,比如mark结合的起始位点, 能结合的位点种类或数目,富集时的信号强度.

## 二、step1: 数据预处理

### 2.1 ChIP文库
+ potential technical problems
+ + adapter contimination
+ + PCR duplication
+ potential biological problems
+ + lack of enrichment
+ + other selection bias

### 2.2 质控

+ QC of reads: fastqc, multiqc
+ QC of alignment: multiqc
+ Filtering of alignment with MAPQ: samtools
+ Deduplication: picard
+ 基于read density去除outliers


```bash
# 用MAPQ进行过滤其值低于20的
samtools view -q 20 -b -o filtered.bam input.bam
```

#### what about deduplication?为何不直接去掉所有的重复？

重复可以是由PCR引起的，也可能是由enrich引起的同一序列存在多个的情况，一般来说，重复使得富集更为明显，对识别富集区域(可能的结合区域/位点)有一定帮助。基于此，我们不能一股脑地去掉所有重复。

## 三、step2：peak calling

### 3.1 Peak callers workflow
![Peak callers workflow](https://s1.ax1x.com/2018/10/14/iUVn3D.png)

+ 优化初始数据：校正forwar/reverse peak offset, deduplication
+ 构建模型：peak = Observed + Model
+ 滑动窗口: 窗口大小=fragmentSize/2, 保留counts数大于Model的窗口
+ 校正：若合并后的窗口counts数大于合并后model的counts，则合并相邻窗口生成总的候选peak set

### 3.2 Peak calling 失败的原因
Peak calling是富集程度、背景值、序列总数的综合结果。
+ 用于Observed的data是进行Model的data的子集(比如处理无效果，无处理的两组)
+ Observed data有peak的地方，Model data连数据都没有(With no input the region around the peak is used to model the background)

### 3.3 结果可靠吗？

+ 多数ChIP enrichment不是链特异性的，你应该在两条链上都能看到富集
+ 重复之间的结果应该是一致的

### 3.4 下游分析

+ composition and motif analysis: 前者分析peak所在DNA组成，后者分析序列潜在生物学功能
+ GO：peak所在基因的GO富集

## 四、step3：探索peak data

### 4.1 可视化peak在序列上的分布情况

+ Is there any enrichment?
+ What is the size / patterning of enrichment?
+ How well are my controls behaving?
+ What is the best way to quantitate this data?
+ Are there any technical artefacts?
+ Are my peaks narrow or broad
+ Do peak positions obviously correspond to existing features(like TSS)?

#### IGV查看enrichment

> 请注意，因为你需要使用IGV可视化查看enrichment的分布情况，这里IGV的input bam文件是sorted bam文件。


```bash
samtools sort -O bam -o $name.sorted.bam $name.bam
```

### 4.2 检查对照

如果有使用IgG或Mock IP，不同的片段化DNA的方法的话，检查它们各自的peak中，peak coverage是均匀的吗？peak的分布模式是一致的吗？如果不是，那么：
+ Low coverage：Repetitive unmappable regions, Holes in the assembly
+ High coverage: Mismapped reads from outside the assembly
+ Biases: DNA stability, GC content, Segmental Duplication

### 4.3 compare samples

+ 可视化比较peak分布
+ scatterplot比较input/chip，input/input，chip/chip
+ correlation metrix, correlation tree, pca-plot, tSNE plot
+ Cumulative Distribution Plot, Q-Q plot

### 4.4 把enrichment/peak同序列特征相联系：trend plots
可视化peak在features(Gene body,promoter,CpG islands, Tss)的分布情况，查看是否符合一般性质(比如甲基化修饰应该在CpG islands富集), 然后分析可能的feature.

+ 查看总的enrich情况，看是否有peak在某个feature
![overall average](https://s1.ax1x.com/2018/10/14/iUZSat.png)

+ 查看不同sample在单个feature的enrichment
![enrichment of single feature](https://s1.ax1x.com/2018/10/14/iUZAMQ.md.png)




## 五、使用MACS2进行peak calling


```bash
# 第一个方法,使用python的pip
# 1.先用conda info --envs查看当前的环境, 我的是有base(python3.7)和env_name(python2.7)
# 2.如果没有python2.7的环境,可以创建一个名为py27的环境
conda create -n py27 python=2.7
source activate py27
pip install MACS2
source deactivate py27
```


```bash
# 1. 进入python2.7环境
source activate env_name
macs2 callpeak -c control.bam -t treat.bam \
		-m 10 30 -p 1e-5 \
		-f BAM -g mm -n treat
    # -c <bam>		进行peak calling的对照
    # -t <bam>		进行peak calling的目标
    # -m 10 30		建立双峰模型的参数
    # -p <pvalue>	显著性阈值,1e-5
    # -f BAM 		输入的文件类型
    # -g <organism> 进行peak calling的物种, hg(人), mm(小鼠), rn(rat)
    # -n <string>   输出文件的前缀
nohup macs2 callpeak -c ../02alignment/IgGold.sorted.bam -t ../02alignment/suz12.sorted.bam -m 10 30 -p 1e-5 -f BAM -g mm -n suz12 2>suz12.masc2.log &

# 每个样本共输出四个文件:
    # 1. name_peaks.xls         存放peak信息,格式与bed一致;坐标从1开始,而bed从0开始.
    # 2. name_model.r           peak模型,可直接进行R作图.
    # 3. name_peaks.narrowpeak  peak矩阵,与.broadpeak相似.可导入R等进行数据分析.BED6+4 format file.
    # 4. name_summits.bed       每个peak的peak summits(极值点的位置),用此寻找结合位点的motif.
# 查看peaks数目
ls *.bed | xargs -i wc -l {} # 6514 suz12_summits.bed
```

## 六、使用ChIPseeker进行可视化
### 6.0 BED 文件格式

 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12
---|---|---|---|---|---|---|---|---|---|---|---
chrom | start | end | name | score | strand | thickStart | thickEnd | itemRgb | blockCount | blockSizes | blockStarts

+ 在这12列中,前3列是必需的字段.
+ score是在genome browser中上色的分值(0-1000),就是peak峰的峰高(summit height of fragment pileup)
+ strand值链的正负性, 但是在ChIP中我们无法区分链的正负性
+ thickStart/End是在genome browser画矩形的起点和中点
+ itemRgb是上色的颜色
+ block是子元件,比如外显子内含子,5'utr之类的.

### 6.1 包和文件

TxDb数据库依据你的ChIP实验动物的物种下载对应物种的数据库.也可以自己通过makeTxDbFromBiomart/makeTxDbFromUCSC包来准备TxDb类型的数据库.hg38(TxDb.Hsapiens.UCSC.hg38.knownGene)、hg19(TxDb.Hsapiens.UCSC.hg19.knownGene)、mm10(TxDb.Mmusculus.UCSC.mm10.knownGene)、mm9(TxDb.Mmusculus.UCSC.mm9.knownGene)


```r
library(ChIPseeker)							# for all
library(TxDb.Hsapiens.UCSC.hg19.knownGene)	# for 6.3
library(org.Hs.eg.db)						# for 6.3,6.5
library(clusterProfiler)					# for 6.5,6.6
library(ReactomePA)							# for 6.5,6.6
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
peak = readPeakFile(filename)
```

### 6.2 可视化peak在基因组上的富集区域(peaks coverage plot)

除了peak文件, GRangeList也支持作为出入,来比较不同bed文件的peak情况。


```r
covplot(peak, weightCol='V5')		# V5是peakfile的第五列(scores)
# 指定染色体区域
covplot(peak, weightCol='V5', chrs=c('chr10', 'chr12'), xlim=c(4.5e7, 5e7))
```

### 6.3 结合到TSS区域的peaks分析

#### 6.3.1 peaks热图分析

为了分析结合到TSS区域的peaks, 我们需要准备好表示TSS regions的文件, 然后再进行计算.注意,我们这里指定了TSS上下游3000bp的区域, 你也可以指定其他区域(getBioRegion和getTagMatrix联合使用)

```r
promoter = getPromoter(TxDb=txdb, upstream=3000, downstream=3000, by='gene') # by=gene/transcript
# getBioRegion(TxDb = NULL, upstream = 1000, downstream = 1000, by = "gene") # by=gene/transcript/exon/intron
tagMatrix = getTagMatrix(peak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color='red')
# 一步法热图
peakHeatmap(filename, TxDb=txdb, upstream=3000, downstream=3000, color='red')
```

#### 6.3.2 peaks的平均read count frequency

```r
plotAvgProf(tagMatrix, xlim=c(-3000,3000), 
			xlab="Genomic Region (5'->3')", ylab="Read Count Frequency",
			conf=0.95, resample=1000)	# 添加置信区间
# 由文件名一步到位
plotAvgProf2(filename, TxDb=txdb, upstream=3000, 
			downstream=3000, xlab="Genomic Region (5'->3')", 
			ylab="Read Count Frequency")
```

### 6.4 注释及其可视化

+ 获取peak附近最近的基因(6.4.2)
+ 注释peak所在的基因组区域(6.4.1)

就像拿到fastq进行比对,你需要提供全基因组的参考序列, peak的注释同样需要参考信息, 也就是注释信息, 这些信息包含基因的起始和结束位置, 在基因的哪个区域是内含子, 哪里是外显子等信息.ChIP支持所有的有基因位置注释信息的物种, 这些注释我们要存储在TxDb对象里面.  

注释有好几种种, 一种是genomic annotation(annotation列),另一种是nearest gene annotation(其他咧).第一种genomic annotation的注释peak的位置,在基因的什么地方,比如UTR,外显子内含子之类的, 这个位置可能就是调控的根本, 可变剪切的调控就属于这类.而nearest gene annotation的注释是最近的基因, 是离peak的距离最近的TSS, 这个TSS所在的基因,这个关注的是启动子区域, 这个基因最有可能被调控. 然后由于一个基因有多个TSS, 多个转录本, 所以有一列transcriptId的列. 还有第三种注释,那就是注释peak区域上下游某个范围, 看这个范围的基因都有些啥.


```r
# peakAnno是个csAnno实例, 可以用as.GRanges把它转化成GRanges实例.
# as.data.frame可以把csAnno转换成data.frame,以供写入文件
# annoDb可选,用以添加各种ID的列. 物种要对应正确
peakAnno = annotatePeak(filename, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb='org.Hs.eg.db')
```
#### 6.4.1 peaks分布区域的饼图, 柱形图


```r
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
# vennpie(peakAnno)
upsetplot(peakAnno, vennpie=T)
```

#### 6.4.2 peak(binding site)与最近基因的TSS距离

```r
plotDistToTSS(peakAnno, title='Distribution of transcription factor-binding loci\n relative to TSS')
```

### 6.5 进行功能富集分析

一旦我们掌握peak区域(TF-binding site)最近的基因信息后,就可以对这些基因进行富集分析啦, 看看它们有什么功能咯.这就有好多方向.Go、KEGG、DO(基因和人类疾病,DOSE)、Reactome（基因与pathways和reactions,ReactomePA).为了进行这些富集分析,就可能要用到seq2gene函数了.它可以以多对多的形式把基因组区域同基因联系起来.


```r
pathway1 = ReactomePA::enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1)
gene = seq2gene(peak, tssRegion=c(-1000, 1000), flankDistance=3000, TxDb=txdb)
pathway2 = ReactomePA::enrichPathway(gene)
head(pathway2, 2)
dotplot(pathway2)
```

### 6.6 比较不同的peak data set (多个peak结果一起分析咯)

#### 6.6.1 结合到TSS region的peak分析: 频率和热图

plotAvgProf,tagHeatmap都支持tagMatrixList的输入, plotAvgProf2和peakHeatmap也接受filenameList


```r
peaks = sapply(filenames, readPeakFile)
tagMatrixList = lapply(peaks, getTagMatrix, window=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000,3000))		# average profiles of ChIP peaks among different experiments
plotAvfProf(tagMatrixList, xlim=c(-3000,3000), conf=0.95, resample=500, facet='row')	# 置信区间和按行作图
tagHeatmap(tagMatrixList, xlim=c(-3000,3000), color=NULL)	# heatmap of peaks among different experiments
```

#### 6.6.2 peak注释比较

多个annotatePeak输出的结果整合成列表, 也可以输入到plotAnnoBar和plotDistToTSS里面去.


```r
peakAnnoList = lapply(filenames, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=F)
plotAnnoBar(peakAnnoList)	# Genomic Annotation among different ChIPseq data
plotDistToTSS(peakAnnoList)	# Distribution of Binding Sites among different ChIPseq data
```

#### 6.6.3 功能富集比较


```r
genes = lapply(peakAnnoList, function(x) as.data.frame(x)$geneId)
names(genes) = sub('_', '\n', names(genes))
compKEGG = clusterProfiler::compareCluster(geneCluster='genes', fun='enrichKEGG',
							# One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
							pvalueCutoff=0.05,
							pAdjustMethod='BH')
dotplot(compKEGG, showCategory=15, title='KEGG Pathway Enrichment Analysis')
```

#### 6.6.4 peak和注释后基因的重叠

这个可以用于实验重复/样品重复peak结果的一致性. 也可以看不同样品peak结果的交叉点.


```r
genes = lapply(peakAnnoList, function(x) as.data.frame(x)$geneId)
vennplot(genes)
```

### 6.7 ChIP-seq peak结果重叠的统计分析
如果两个不同的ChIP实验的peak结果有很大重叠的话, 表明这两个实验中用于IP的蛋白极有可能是分子伴侣(cooperative in regulation). 所以,通过对重叠进行统计分析可以确认这个事实的显著性.  
当然,这里面有一些理论前提.我们把基因组位置随机打断的话, 在把这些片段进行peak calling, 我们会得到一个peak的概率分布. 基于此, 我们可以统计ChIP-seq peak calling的显著性.


```r
p = GRanges(seqnames=c('chr1', 'chr3'), ranges=IRanges(start=c(1,100), end=c(50,130)))
shuffle(p, TxDb=txdb)
# 基于基因组坐标计算的overlap显著性
#（overlap significant of ChIP experiments based on the genome coordinations)
enrichPeakOverlap(queryPeak=files[[5]], targetPeak=unlist(files[1:4]),
				  TxDb=txdb, pAdjustMethod='BH',
				  nShuffle=50, chainFile=NULL,
				  verbose=F)

# 基于最近TSS的基因注释计算overlap显著性
# overlap sig. calculated based on their nearest gene annotation
enrichAnnoOverlap(queryPeak = filename, targetPeak = filename(s),
				  TxDb=txdb, pAdjustMethod='BH',
				  chainFile=NULL, distanceToTSS_cutoff=NULL)
```

### 6.8 GEO中ChIP-seq数据挖掘

我们已经可以计算不同ChIP-seq实验peak数据的重叠显著性了.而在GEO上有相当多的数据,那么我们可以拿自己的数据与里面的peak数据去找重叠,如果找到了新的蛋白-蛋白复合转录体,又是个新发现了.
这个包里面涵盖了17000个GEO bed文件, 可以通过getGEOspecies来获取这些文件的概要.


```r
getGEOspecies()	# 基于物种
getGEOgenomeVersion() # 基于基因组版本
# 获取详细信息
hg19 = getGEOInfo(genome='hg19', simplify=T)
# 下载数据
downloadGEObedFiles(genome='hg19', destDir='hg19')	# 根据基因组版本
downloadGSMbedFiles(gsm_vector, destDir='hg19')	    # 根据gsm accession来下载
```


## 七、使用deeptools进行可视化

deeptools可以用于处理比对结果的多项质控, 并基于比对文件生成校正后的bed或bigwig格式的覆盖度文件(normalized coverage file),这就允许多个处理组的比较了. 当然,利用这些文件, 你还可以进行多种可视化.

![deeptools in ChIP-seq](http://deeptools.readthedocs.io/en/latest/_images/start_workflow.png)

+ 多个Bam文件的相关性: multiBamSummary, plotCorrelation, multiBigWigSummary
+ 比对后read的覆盖度: plotCoverage(--ignoreDuplicates is useful), bamPEFragmentSize(for PE)
+ GC-bias: computeGCBias, correctGCBias
+ 估测ChIP信号强度:

### 7.1 质控和数据处理

1. 从bam文件生成bigwig
[bamCoverage](http://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html)

```bash
bamCoverage -b reads.bam -o coverage.bw \
			--binSize 10 \
			--normalizeUsing RPGC \
			--extendReads \
# -b,--bam 			输入的bam文件
# -o 				输出文件名
# -of 				输出格式, bigwig, bedgraph
```

2. 检测多个测序结果之间的可重复性

```bash
multiBamSummary, plotCorrelation
```

3. 检测和校正 GC bias

```bash
computeGCBias
```

4. 以特定文件校正ChIP的覆盖率度

```bash
bamCompare with ChIP = treatment, input=control
```

5. 检测不同ChIP实验的信号强度

```bash
plotFingerprint
```

6. 得到TSS(转录起始位点)富集基因的覆盖度热图

```bash
computeMatrix, plotHeatmap
```

7. 比较基因在常染色体与性染色体之间的信号差别

```bash
1. 过滤出欲比较的基因
2. computeMatrix
3. plotProfile
4. plotting the summary plots for multiple samples
```

### 7.2 函数的使用

#### 7.2.1 computeMatrix
对基因组上的任何一个区域进行打分(peak 数),生成供plotHeatMap和plotProfiles用的中间文件。它的输入的是多个打分文件如bigwig文件,bed文件。当然，这个函数还可以根据输入的bw/bed文件对基因组上的区域进行筛选.  
它有两种模式,参考点模式(reference-point)和大区域模式(scale-regions).

bw文件可以通过bamCoverage/bamCompare函数得到.

+ computeMatrix scale-regions模式

```bash
computeMatrix scale-regions -S <bw/bed files> -R <bed file> -b 1000

# -R,--regionsFileName 			你想可视化的区域的文件名,以BED/GTF格式.
#								多个文件的话,可以使用#进行分组.分组的文件将画在一起.
# -S,--scoreFileName 			打分文件, 以bw,bed格式.看这里,就需要打分文件啦.bigwig.
# -o,-out,--outFileName 		输出.gz文件以供plotHeatMap和plotProfile使用
# --outFileNameMatrix 			输出peak calling矩阵,可以进行差异分析
# --outFileSortedRegions 		在你设置一个筛选条件后输出的筛选后文件名
# -a,--downstream				转录起始位点下游的bp数
# -b,--upstream 				转录起始位点上游的bp数
# 具体筛选参数设定见handbook

computeMatrix scale-regions \
                -R genes_chr19_firstHalf.bed genes_chr19_secondHalf.bed \ # separate multiple files with spaces
                -S testFiles/log2ratio_*.bw  \ or use the wild card approach
                -b 3000 -a 3000 \
                --regionBodyLength 5000 \
                --skipZeros -o matrix2_multipleBW_l2r_twoGroups_scaled.gz \
                --outFileNameMatrix matrix2_multipleBW_l2r_twoGroups_scaled.tab \
                --outFileSortedRegions regions2_multipleBW_l2r_twoGroups_genes.bed
```

+ computeMatrix reference-point模式


```bash
computedMatrix reference-point -S <bw/bed file(s)> -R <bed file> -a 3000 -b 3000

# 参数可见scale-regions模式
# --referencePoint 		指定参考点,可选的参考点为TSS,TES,center.TSS(region start) TES(region end), center(region center)
						不管你用的啥参考点, 可视化时都是用TSS作为默认标签。

computeMatrix reference-point --referencePoint TSS \ # alternatives: TES, center
       -b 3000 -a 10000 \ # define the region you are interested in
       -R testFiles/genes.bed \
       -S testFiles/log2ratio_H3K4Me3_chr19.bw  \
       --skipZeros \
       -o matrix1_H3K4me3_l2r_TSS.gz \ # to be used with plotHeatmap and plotProfile
       --outFileSortedRegions regions1_H3K4me3_l2r_genes.bed
```

#### 7.2.2 multiBamSummary
multiBamSummary计算多个bam文件的基因组区域的read覆盖率,这个可以通过指定‘bins’模式来实现。你也可以通过‘BED-file’模式来计算指定区域的read覆盖度。它的标准输出是压缩后的numpy array(.npz). 可以使用plotCorrelation函数直接计算和可视化成对相关性。你还可以用plotPCA函数来对整个numpy文件进行主成分分析。我们并不推荐你只输入单个bw文件, 如果你只是想生成bedGraph文件的话, 倒是可以使用,但你要指定--outRawCounts选项.

+ multiBamSummary bins

```bash
multiBamSummary bins --bamfiles file1.bam file2.bam ... -o results.npz

# --bamfiles 			你输入的bam文件,用空格隔开
# -o 					输出文件名, read覆盖度矩阵. .npz文件
# -l,--labels			指定输入文件的标签,应该是在覆盖度矩阵里作为列名
# --smartLabels 		自动指定标签, 使用文件名(basename filename)
# -r,--region 			指定计算覆盖度的区域.形式为-region chr10, -region chr10:456700:891000
# --ignoreDuplicates 	忽略重复的reads
# --outRawCounts 		把counts per region以\t分隔的形式写入文件
# -e,--extendReads 		当read的长度不足,允许read延伸到与fragment一致的长度.
						不推荐在含有spliced-read的数据里使用, 如RNA-seq。
# --minMappingQuality	最低的比对质量分数,低于此的reads不进行计算
# --centerReads 		指定reads要覆盖到fragment长度的中心
# 其他设定fragment长短的选项见handbook
```

+ multiBamSummary BED-file 

```bash
multiBamSummary BED-file --BED selection.bed --bamfiles file1.bam file2.bam -o results.npz

# -b,--bamfiles 			索引后的bam文件,以空格分隔
# -o 						生成覆盖度的矩阵文件名
# --BED 					指定计算覆盖度的区域文件
```
