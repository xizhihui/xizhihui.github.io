---
title: 比对 subread
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 软件和包
---

subread是个套件,里面有subread aligner, subjunc aligner, featureCounts, exactSNP.

subread aligner可以用于DNA-seq和RNA-seq.当用于RNA-seq时,subread只适用于差异分析;对于检测基因组变异如可变剪接之类的,需要reads的完全比对,这时候可以使用subjunc进行比对.在比对RNA-seq数据时,subread不会取检测exon-exon junctions的存在,只会把exon-spanning eads的最大可比对区域作为比对结果.但是,如果只是进行差异分析的话,subread的结果足以进行.subread的比对上reads可能会比subjunc多.

<!--more-->

## 安装

```bash
# 下载源文件
https://sourceforge.net/projects/subread/files/
tar zxvf subread-1.6.2.tar.gz
make -f Makefile.Linux

# in R
biocLite('Rsubread')
```

## subread aligner

```bash
subread-buildindex -o my_index ref.fa  			#允许单染色体文件
# -B 					把index进行分割多个文件;会覆盖-M的设置
# -c 					构建color-space索引
# -f <int> 				去除高度重复序列的重复阈值,高于该阈值的重复序列将被去除
# -F 					构建一个完全索引,小鼠的完全索引达14G
# -M <int> 				可以用int MB的RAM
# -o <string>			index的basename

subread-align -t 1 -i my_index -r reads.fastq -o subread_results.bam
# -t <int>				进行比对的数据类型, 0是RNA-seq, 1是DNA-seq
# --multiMapping 		允许多比对
# -B <int> 				允许的多比对数指定为int
# -T <int>				指定用int个threads
# -I <int> 				检测的indel长度最长为int bp
# -d <int> 				minFragLength指定为50
# -D 600				maxFragLength指定为600
# -r fastq文件			输入的reads
# -R read2.fastq 		如果是paired-end,由此指定read2
# -a <string>			指定注释文件
# -A <string> 			指定参考基因组和注释之间的chr名对应关系的文件
#						第一列是注释里面的染色体名,第二列是参考基因组对应的染色体名.不需要列名.
```

```bash
# -b 					当比对color space文件时,输出的比对结果保持正常格式,而非color-space
# -F <string> 			指定注释文件的格式, 'GTF', 'GFF', 'SAF'(在Rsubread里默认这个)
# -m <int> 				指定一致性阈值,当比对一致性超过该阈值则认为比对上
# -n <int> 				允许的最大错配数
# -p <int> 				当是paired-end时,两个reads的比对一致性最低值,应该小于-m的指定值
# -P 3/6/33/64			指定使用的phred质量值,3指phred+33; 在Rsubread中,33值phred+33.

# 其他有关genomic variance的设置见handbook.

# 比对microRNA-seq reads
# 注释下载: http://www.mirbase.org/
subread-buildindex -F -B -o mm10_full_index mm10.fa
subread-align -t 1 -i mm10_full_index \
			  -n 35 -m 4 -M 3 -T 10 -I 0 \
			  --multiMapping -B 10 \
			  -r miRNA_reads.fastq -o results.sam
```


## subjunc

```bash
subjunc -i my_index -r rnaseq-reads1.txt -R rnaseq-read2.txt -o subjunc_result
```

## featureCounts定量
```
featureCounts -T 5 -a annotation.gtf -t exon -g gene_id \
			  -o counts.txt mapping_results_SE.sam
			  # 可以有多个bam文件
# paired-end
featureCounts -p -a annotation.gtf -t exon -g gene_id \
			  # -P -d 50 -D 600		# 指定了fragment长度
			  # -B 	# 不考虑fragment,但两个read都要比对上,才计数
			  # -C  # 排除嵌合(chimeric fragments)
			  -o counts.txt mapping_results_PE.bam


```

## SNP calling

```bash
exactSNP [options] -i input -g reference genome -o output

# -a <file> 			指定vcf格式的snp注释文件
# -b 					指定输入文件是bam文件
# -f <float> 			指定含SNP位点区域的错配碱基的最小区域
# -g <file> 			指定参考基因组文件, fasta格式
# -i <file>				指定输入文件,SAM或BAM,如果是BAM,指定-b选项
# -n <int> 				指定出现错配碱基的最小数目
# -Q <int> 				指定在50x测序深度是的SNPcalling的q-value的cutoff值
# -r <int> 				判定为SNp时需要的最少比对上reads数目
# -s <int> 				指定作为SNP位点的碱基质量值阈值
# -t <int> 				read两端切除的碱基数
# -T <int> 				Threads
# -x <int> 				判断为SNP的位点的最大测序深度
```

## 在R里面

```r
library(Rsubread)
buildindex(basename='my_index', reference='genome.fa')
align(index='my_index', type='dna',
	  readfile1='reads.txt.gz',
	  output_file='rsubread.bam',
	  nthreads=5,
	  indels=16,
	  unique=F, nBestLocations=3)
# pairedd-end
align(index='my_index',
	  readfile1='reads1.fq.fz', readfile2='reads2.fq.gz',
	  type='dna',
	  output_file='rsubread.bam',
	  minFragLength=50, maxFragLength=600)

subjunc(index='my_index', readfile1='rnaseq-reads.txt.gz', output_file='subjunc_results.bam')

featureCounts(files="mapping_results_SE.sam", nthreads=5)
# 提供参考注释的话
featureCounts(files="mapping_results_SE.sam",annot.ext="annotation.gtf",
isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id")

# paired-end
featureCounts(files=c('read1_mapped.bam', 'read2_mapped.bam'))
featureCounts(files='PE_mapped.bam', isPairedEnd=T,
			 requireBothEndsMapped=TRUE, # 要两个reads都比对上才计数
			 countChimericFragments=FALSE) # 不考虑chimeric fragment
```

## 其他小功能

+ repair: 把paired-end.bam文件里面的reads成对放置
+ coverageCount: 计算基因组中某个区域的read coverage
+ flattenGTF: 把GTF里面的某个feature的行提取成SAF文件.如提取所有的exons注释
+ promoterRegions:只在Rsubread里面有,提取出每个基因的启动子区域位置.
+ removeDunp: 去除SAM文件里面的重复reads
+ subread-fullscan：提取目标序列在染色体上的高度同源序列位置
