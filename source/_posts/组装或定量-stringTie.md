---
title: 组装或定量 stringTie
date: 2018-10-12
toc: true
categories: Bioinformatics
tags:
  - 软件和包
  - 定量
  - 组装
---

stringTie 是用于 RNA-seq 的转录本组装和定量软件

<!--more-->

## 基本使用

输入文件是BAM格式的比对结果文件，该文件必需经过排序，排序的方式基因组位置。这些文件可以是来源于Tophat比对的结果文件accepted_hits.bam，也可以是hisat2的结果文件经过转换和排序的文件(使用samtools)。


```bash
stringtie <aligned_reads.bam> [options]*

# 选项
-o [<path/>]<out.gtf> 			设置输出文件;组装的结果将写入out.gtf
-p <int>						指定threads
-G <ref_ann.gff>				指定参考注释文件
--rf 							说明是fr-firststrand的文库
--fr 							说明是fr-secondstrand的文库
-l <label>						设定输出转录本的前缀,默认是STRG
-m <int> 						最短的预测转录本长度
-A <gene_abund.tab>				将基因表达量输出到gene_abund.tab文件(tab delimited format)
-C <cov_refs.gtf>				输出reads能完全覆盖的参考序列区域的所有转录本
-a <int>						没有剪切后reads比对上的junction最短长度
-j <float>						可比对到junction的剪切reads的最少个数
-t 								禁止在组装转录本的两端进行trim
-c <float> 						设置进行转录本预测的最小read覆盖度
-g <int>						locus的最小间隔值,如果reads比对区域距离小于该值,则会被合并
-B 								允许生成Ballgown的输入文件
-b <path> 						指定-B生成文件的存储路径
-e 								指定只进行比对文件的定量和输出能比对到参考转录本的组装转录本
-M <0.0-1.0>					设定多位置比对reads在给定位置序列的最大占比
-x <seqid_list>					忽略这里设置的参考序列区域,可以是逗号分隔的染色体. -x 'chrM,chrX,chrY'
--merge							指定使用转录合并模式,在本模式中，将使用GTF/GFF文件作为输入，然后把这些转录合并成非冗余的转录集合。
								主要是在多样本的RNA-seq数据结果中使用.
								本模式下有效的选项:-G, -o, -c,-m,-F,-T,-f,-i,-l
```


> 使用hisat2结果作为输入时, 首先在比对时要指定--dta选项,其次需要进行排序


```bash
samtools view -Su alns.sam | samtools sort - alns.sorted
```

## 输出文件

主要的输出:

+ .GTF文件, 组装的转录本

```bash
#seqname source      feature     start   end     score   strand  frame attributes
#chrX    StringTie   transcript  281394  303355  1000    +       .     gene_id "ERR188044.1"; transcript_id "ERR188044.1.1"; reference_id "NM_018390"; ref_gene_id "NM_018390";ref_gene_name "PLCXD1"; cov "101.256691"; FPKM "530.078918"; TPM "705.667908";
#chrX    StringTie   exon        281394  281684  1000    +       .     gene_id "ERR188044.1"; transcript_id "ERR188044.1.1"; exon_number "1"; reference_id "NM_018390";ref_gene_id "NM_018390"; ref_gene_name "PLCXD1"; cov "116.270836";
```

+ .tab文件, 基因表达量

```bash
#Gene ID     Gene Name   Reference   Strand  Start   End     Coverage    FPKM        TPM
#NM_000451   SHOX        chrX        +       624344  646823  0.000000    0.000000    0.000000
#NM_006883   SHOX        chrX        +       624344  659411  0.000000    0.000000    0.000000
```

+ .GTF文件, 能与参考注释匹配的完全覆盖转录
同组装的转录本

+ 用于Ballgown进行下游差异分析的文件
总共有5个文件, e2t.ctab, e_data.ctab, i2t.ctab, i_data.ctab, t_data.ctab

+ 在merge mode, merged_gtf文件
如果stringtie在merge mode下运行的话，会将输入的一系列GTF/GFF文件合并组装成无冗余的转录本。输出的文件只是包含转录本，没有其他的数据如coverage，FPKM，TPM。stringtie可以用这个新的转录本重新计算表达量，但是你得设置-e参数。如果你要寻找新的转录本或者寻找转录本的来源，你可以使用gffcompare软件来实现。

## 结合Ballgown进行差异分析
### 完全的差异分析
![差异分析](https://ccb.jhu.edu/software/stringtie/DE_pipeline.png)

推荐的流程:
1. 对每个样本, 使用hisat2 --dta选项与参考基因组进行比对。  
在比对时，我们强烈推荐加入参考注释信息。这可以通过hisat2 --ss/--exon实现,也可以通过hisat2的--know-splicesite-infile实现。值得注意的是，你一定要对输出的比对文件进行排序和转换成BAM文件。

2. 对得到的每个样本的比对文件，用strintie进行组装。  
我们推荐使用-G参数提供参考注释(如果有的话.)

3. 用生成组装文件，使用stringtie --merge生成无冗余的转录本  
同样的，如果有参考注释，也推荐使用

4. 对每个样本的比对文件, 使用stringtie -B/-b, -G和-e选项计算表达量和生成Ballgown所需的table文件。  
-G选项指定第三步生成无冗余转录本，这是唯一没有使用参考注释的情况。虽然说这里的-e选项不是必需的，但是使用它可以得到一个更准确的结果。

5. 使用Ballgown进行差异分析

### 简化的差异分析
![简化的差异分析](https://ccb.jhu.edu/software/stringtie/DE_pipeline_refonly.png)

---

## stringtie结合DESeq2和edgeR的差异分析

 1. 使用stringtie -e得到定量结果
 2. 使用[“prepDE.py”](https://ccb.jhu.edu/software/stringtie/dl/prepDE.py)提取定量结果生成counts矩阵


```python
# prepDE.py [options]
# -i <input>	 	指定gtfs文件的路径和样品ID的txt文件
# -g G 			指定gene count matrix的输出路径
# -t T 			指定transcripts count matrix的输出路径
# -l length 		指定平均read长度
# -p pattern		选择样品子路径的正则表达式
# -c 			是否合并重复的gene,虽然它们有不同的gene ID
# -s string 		指定stringtie添加geneID前缀的字符串,默认为MSTRG
# -k Key 		如果-c指定,则指定本脚本添加到gneIDs的前缀, 默认为prepG
# --legend=Legend 如果-c指定,转录比对到geneIDs的表头文件的输出路径,默认为legend.csv

#prepDE.py -i sample_lst.txt

 ## sample_lst.txt
 ## ERR188021 <PATH_TO_ERR188021.gtf>
 ## ERR188023 <PATH_TO_ERR188023.gtf>
 ## ERR188024 <PATH_TO_ERR188024.gtf>
```

 3. 使用DESeq2进行差异分析


```r
countData <- as.matrix(read.csv('gene_count_matrix.csv', row.names='gene_id'))
colData <- read.csv(PHENO_DATA, sep='\t', row.names=1) 	# 这里是你的表型文件
all(rownames(colData) %in% colnames(countData))			# 检测是否一致
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

dds <- DESeqDataSetFromMatrix(countData = countData,
							  colData = colData,
							  design=~CHOOSE_FEATURE)
dds <- DESeq(dds)
res <- results(dds)
```

---

## 组装super-reads
stringtie还可以把短的transcripts组装成长的contigs，我们把这个叫做super-reads。这一步可以省略，但我们还是推荐使用这一步，因为他可以提高转录本组装的正确性。不过你需要安装MaSuRCA genome assembler包。

1.super-reads的生成

```bash
superreads.pl reads_1.fastq reads_2.fastq <masurca_path> [options]*

# 选项
-t <num_threads>
-j <jf_size> 			MasuRCA需要运行Jellyfish,这里设置后者使用的hash size
-s <step>				打印运行成功的步骤
-r <paired_read_prefix> 设置paired read的前缀
-f <fragment_size>		设置mean library insert length
-d <stdev>				设置insert length的标准偏差
-l <name> 				设置组装的super-reads文件名
-u <prefix>				设置为组装的reads的文件名
```

2. super-reads比对到参考基因组上
可以使用你喜欢的比对软件，如tophat，hisat，bowtie等。比如：

```bash
tophat [options]* <genome_index_base> PE_reads_1.notAssembled.fq.gz,super_reads.fq PE_reads_2.notAssembled.fq.gz
```
