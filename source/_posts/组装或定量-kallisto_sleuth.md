---
title: 组装或定量 kallisto_sleuth
date: 2018-10-12
toc: true
categories: Bioinformatics
tags:
  - 软件和包
  - 定量
---

## 免比对的定量：kallisto
[kallisto](https://pachterlab.github.io/kallisto/manual.html)是一个align-free的测序结果定量工具。速度快得很.他可以在10min内完成index的构建，然后花3min完成30 million的human reads的定量。它不仅快速而且准确。

Pseudoalignment保存了用于定量的关键信息。An important feature of kallisto is that it outputs bootstraps along with the estimates of transcript abundances

<!--more-->

### 流程概览

```bash
# 建立一个索引
kallisto index -i transcripts.idx transcripts.fasta.gz
# 进行定量
kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq.gz reads_2.fastq.gz
# single-end
# kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 reads.fastq.gz

```

### 0. 安装

```bash
conda install kallisto
```

### 1. kallisto index: 建立索引


```bash
kallisto index [arguments] ref.fa[.gz]

# 必需参数
-i, --index=STRING          存放索引的目录名字

# 可选参数
-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: 31)
    --make-unique           Replace repeated target names with unique names
```

### 2. kallisto quant: 进行定量
#### 注意事项
+ 每次的输入文件只能是一个样本的文件, 输入多个文件仅仅是因为这个样本的数据被分成了多个文件.
+ 对于single-end,你必需通过-l参数指定fragment length


```bash
kallisto quant [arguments] reads.fastq

示例: pair-end / single-end
kallisto quant -i index -o output pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq
kallisto quant -i index -o output --single -l 200 -s 20 file1.fastq.gz file2.fastq.gz file3.fastq.gz

# 必需参数
-i, --index=STRING            索引所在目录
-o, --output-dir=STRING       输出目录。默认输出三个文件,abuncance.h5, abundance.tsv, run_info.json

Optional arguments:
    --bias                    Perform sequence based bias correction
-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)
    --seed=INT                Seed for the bootstrap sampling (default: 42)
    --plaintext               指定输出纯文本文件,而非HDF5格式
    --fusion                  搜索融合位点?(Search for fusions for Pizzly)
    --single                  输入的文件是single-end reads
    --single-overhang         对不再转录本上的read也进行定量
    --fr-stranded             Strand specific reads, first read forward
    --rf-stranded             Strand specific reads, first read reverse
-l, --fragment-length=DOUBLE  Estimated average fragment length
-s, --sd=DOUBLE               Estimated standard deviation of fragment length
                              (default: -l, -s values are estimated from paired
                               end data, but are required when using --single)
-t, --threads=INT             Number of threads to use (default: 1)
    --pseudobam               Save pseudoalignments to transcriptome to BAM file
    --genomebam               Project pseudoalignments to genome sorted BAM file
-g, --gtf                     GTF file for transcriptome information
                              (required for --genomebam)
-c, --chromosomes             Tab separated file with chrosome names and lengths
                              (optional for --genomebam, but recommended)
```

### 3. kallisto preudo: 伪比对

这个功能是进行伪比对这一步，主要是用于单细胞RNA-seq的。它的命令在形式上与目的上与kallisto quant相似，但是它不会利用EM-算法来计算reads abundance.而且它还有个选项指定多个细胞到一个批文件里面.


```bash
kallisto pseudo [arguments] FASTQ-files

# 必需参数
-i, --index=STRING            Filename for the kallisto index to be used for
                              pseudoalignment
-o, --output-dir=STRING       Directory to write output to

# 可选参数
-u  --umi                     First file in pair is a UMI file
-b  --batch=FILE              Process files listed in FILE
    --single                  Quantify single-end reads
-l, --fragment-length=DOUBLE  Estimated average fragment length
-s, --sd=DOUBLE               Estimated standard deviation of fragment length
                              (default: -l, -s values are estimated from paired
                               end data, but are required when using --single)
-t, --threads=INT             Number of threads to use (default: 1)
```


```bash
# batch file的形式, single-end只有一个fastq,2变成1
{	
	#id file1 file 2
	cell1 cell1_1.fastq.gz cell1_1.fastq.gz
	cell2 cell2_1.fastq.gz cell2_1.fastq.gz
	cell3 cell3_1.fastq.gz cell3_1.fastq.gz
}
# 对于--umi指定后的形式
{
	#id umi-file file-1
	cell1 cell_1.umi cell_1.fastq.gz
	cell2 cell_2.umi cell_2.fastq.gz
	cell3 cell_3.umi cell_3.fastq.gz
}
```

### 4. kallisto h5dump: 把HDF5文件转为纯文本格式

```bash
kallisto h5dump -o output_path abundance.h5
```

## 利用sleuth进行差异分析

sleuth可以与kallisto进行无缝对接，进行表达定量的后续差异分析；它在RStudio上工作，让你能交互式地进行数据探索。

### 安装

```r
# 通过bioconductor安装
source('http://bioconductor.org/biocLite.R')
biocLite('devtools')
biocLite('pachterlab/sleuth')
# 通过conda安装
# conda install --channel bioconda r-sleuth
```

### 特性
+ 既适用于transcript-level的分析，也适合gene-level的分析
+ 与kallisto兼容，让你的分析流程更快
+ 使用bootstraps来检测和校正实验的技术性变异
+ 用于数据分析的可交互应用

### [get started with RStudio](https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html)

#### 1. 构建实验设计表格

```r
# 不打印messages
suppressMessages({
	library('sleuth')
})
options(stringsAsFactors=F)

# 找到kallisto处理结果
sample_id <- dir(file.path('..', 'results'))
kal_dirs <- file.path('..', 'results', sample_id, 'kallisto')

# 找到对于的实验设计文件
s2c <- read.table(file.path('..', 'metadata', 'hiseq_info.txt'), header=T)
s2c <- dplyr::select(s2c, sample=run_accession, condition)
# 将path作为列添加上去，否则会报错
s2c <- dplyr::mutate(s2c, path=kal_dirs)

##      sample condition                          path
## 1 SRR493366  scramble ../results/SRR493366/kallisto
## 2 SRR493367  scramble ../results/SRR493367/kallisto
## 3 SRR493368  scramble ../results/SRR493368/kallisto
## 4 SRR493369   HOXA1KD ../results/SRR493369/kallisto
## 5 SRR493370   HOXA1KD ../results/SRR493370/kallisto
## 6 SRR493371   HOXA1KD ../results/SRR493371/kallisto
```

#### 2. 构建‘sleuth’对象  
sleuth对象不仅包含了实验信息，也包含了用于差异分析的模型和结果。

```r
# 1. 添加kallisto处理好的数据
so <- sleuth_prep(s2c, extra_bootstrap_summary=T)

# 2. 确定sleuth应答错误衡量(full)模型的参数, condition是前面构建的实验设计表格里面的
# 详细讲就是~condition作为formula进行线性回归，然后对样本数据进行平滑
so <- sleuth_fit(so, ~condition, 'full')

# 3. 确定sleuth归并模型的参数;前提是假定各个条件下的表达量是一致的
so <- sleuth_fit(so, ~1, 'reduced')

# 4. 使用似然率确定表达差异
so <- sleuth_lrt(so, 'reduced', 'full')

# 查看使用了哪些回归模型
models(so)

# 5. 进行差异分析
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all=F)
sleuth_significant <- dplyr::filter(sleuth_table, qval<=0.05)
head(sleuth_significant, 20)
```

#### 3. 基因名注释

```r
# 注意,这里使用的ENSEMBL human transcriptome, 所以使用biomaRt
source('http://bioconductor.org/biocLite.R')
biocLite('biomaRt')
```

```r
mart <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL',
						dataset='hsapiens_gene_ensembl',
						host='ensembl.org')
t2g <- biomaRt::getBM(attributes=c('ensembl_transcript_id',
									'ensembl_gene_id',
									'extenal_gene_name'),
					  mart=mart)
t2g <- dplyr::rename(t2g, 
					 target_id=ensembl_transcript_id,
					 ens_gene=ensembl_gene_id,
					 ext_gene=external_gene_name)

# 添加到sleuth对象上
so <- sleuth_prep(s2c, target_mapping=t2g)

# 再进行一次计算
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_rlt(so, 'reduced', 'full')
```

#### 4. 可视化结果
+ 对结果的查看和交互可以通过产生sleuth live site实现

```r
sleuth_live(so)
```

+ 箱图

```r
plot_bootstrap(so, 'ENST00000264734', units='est_counts', color_by='condition')
```

+ PCA plot
以实验条件的pca条件

```r
plot_pca(so, color_by='condition', text_labels=T)
```

+ 样品counts分布图

```r
plot_group_density(so, use_filtered=T, units='est_counts', trans='log',grouping=setdiff(colnames(so$sample_to_corvariates), 'sample'),offset=1)
```
+ MA-plot

```r
plot_ma(so, test='condition')
```

+ 均值方差分析

```r
plot_mean_var(so)
```

+ 热图

```r
library(pheatmap)
plot_transcript_heatmap(so, transcripts=sleuth_table$target_id[1:20])
```

+ 火山图

```r
plot_volcano(so, test='condition')
```
