---
title: 质控软件 raw_data质控
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 生信基础
---


测序数据的质控是必须的，也是主要的。

<!--more-->

## 1. fastqc


```bash
fastqc input_path/*.fastq(.gz) -o output_path -p 10
```


## 2. cutadapt

https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage

基本的用法如下，你需要将序列改成对应的3‘端adapter。支持压缩文件作为输入和输出。如果你的输出是fasta的话，cutadapt不会进行去除adapter，而是进行转化。


```bash
	# 对于pair-end
	cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq

	# 对于single-end
	cutadapt -a AACCGGTT -o output.fastq input.fastq

	cutadapt -a AACCGGTT input.fastq > output.fastq 2><report class="txt"></report>
	# 使用管道的话，输入文件位置以-代替
	tail -n 4 input.fastq | cutadapt -a AACCGGTT - > output.fastq

	# 双端时，-u指前一个，-U指后一个。指定切除特定长度的碱基, -5 表示切除3‘端（反向切除5个碱基），5表示切除5’端（正向切除5个碱基）
	cutadapt -u 5 -o trimmed.fastq input.fastq

	# 单双端都可用. 切除碱基质量低于10的碱基（3‘端），使用的是phred quality+33.如果要使用+64版本，--quality-base=64指定
	cutadapt -q 10 -o output.fastq input.fastq
	# 切除5‘端低于15的，3’端低于10的碱基
	cutadapt -q 15,10 -o output.fastq input.fastq

	# 把read切成指定长度，注意是从3‘端切掉.如果要从两端切掉，使用-u参数
	cutadapt -l 120 -o output.fastq input.fastq

	# 给read name添加内容，下方指令把read1 改成read1 we found adapter1如果ACGT在该read1找到的话
	cutadapt -a adapter1=ACGT -y ' we found {names}' input.fastq

	# 如果要添加上length字符的话. output.fastq => ">read1 length=120"
	cutadapt -l 120 -o output.fastq --length-tag 'length=' input.fastq

	# -m length: 只输出长度大于length的reads
	# --too-short-output FILE: 把-m删除的reads输入到FILE
	# -M length：只输出短于length的reads
	# --too-long-output FILE: 把-M删除的reads出入到FILE
	# --untrimmed-output FILE
	# --discard-trimmed:去掉含有接头的reads
	# --discard-untrimmed：去掉不含接头的reads
```

除了去除adapter，cutadapt可以做的事情有很多，比如3‘端低质量碱基的切除，reads过滤等等。下表是指定adapter类型的对应指令。

接头类型 | 命令 | 接头类型 | 命令
---|---|---|---
3’adapter | -a ADAPTER | 5'adapter | -g ADAPTER
Anchored 3'adapter | -a ADAPTER$ | Anchored 5'adapter | -g ^ADAPTER
5'or3'(both possible) | -b ADAPTER | Linked adapter | -a ADAPTER1...ADAPTER2


```bash
#	1. Unconditional base removal with --cut
#	2. Quality trimming (-q)
#	3. Adapter trimming (-a, -b, -g and uppercase versions)
#	4. Read shortening (--length)
#	5. N-end trimming (--trim-n)
#	6. Length tag modification (--length-tag)
#	7. Read name suffix removal (--strip-suffix)
#	8. Addition of prefix and suffix to read name (-x/--prefix and -y/--suffix)
#	9. Double-encode the sequence (only colorspace)
#	10.Replace negative quality values with zero (zero capping, only colorspace)
```

## 3. Fastx-toolkit

用途 | 用途 | 用途
---|---|---
01.去除接头 | 02.去除低质量碱基 | 03.转换成fasta  
04.碱基质量统计 | 05.生成碱基质量箱图 | 06.生成碱基质量分布图
07.重命名read | 08.输出指定长度read(前后切除) | 09.去重reads
10.去除人工reads(artifact) | 11.生成反向互补序列 | 12.对fasta文件进行格式化
13.DNA序列与RNA互换 | 14.生成reads长度分布图 | 15.依据barcode进行reads分类成不同文件

### 3.1.去除接头


```bash
fastx_clipper -a AGATCGGAAGAGCACACG -l 25 -d 0 -Q 33 -i SRR306394_1.fastq -o 

 # [-a ADAPTER] =接头序列（默认为CCTTAAGG）
 # [-l N]       = 忽略那些碱基数目少于N的reads，默认为5
 # [-d N]       = 保留接头序列后的N个碱基默认  -d 0
 # [-c]         = 放弃那些没有接头的序列.
 # [-C]         = 只保留没有接头的序列.
 # [-k]         = 报告只有接头的序列.
 # [-n]         = 保留有N多序列，默认不保留
 # [-v]         =详细-报告序列编号
 # [-z]         =压缩输出.
 # [-D]       = 输出调试结果.
 # [-M N]   =要求最小能匹配到接头的长度N，如果和接头匹配的长度小于N不修剪
 # [-i INFILE]  = 输入文件
 # [-o OUTFILE] = 输出文件
```

### 3.2. 去除低质量碱基

```bash
 fastq_quality_filter -q 20 -p 80 -Q 33 -i input.fastq -o output.fastq

# [-q N]       = 最小的需要留下的质量值
# [-p N]       = 每个reads中最少有百分之多少的碱基需要有-q的质量值
# [-z]         =压缩输出
# [-v]         =详细-报告序列编号，如果使用了-o则报告会直接在STDOUT，如果没有则输入到STDERR
```

### 3.3. 转换成fasta

```bash
fastq_to_fasta [-h] [-r] [-n] [-v] [-z] [-i INFILE] [-o OUTFILE]
[-h]         显示帮助信息
[-r]         使用连续的数字来重命名序列
[-n]         保持序列中的未知碱基N, 默认为丢弃这样的序列
[-v]         打印出转换的reads数
[-z]         用GZIP压缩输出
[-i INFILE]  指定输入fastq
[-o OUTFILE] 指定输出fasta
```

### 3.4. 碱基质量统计

```bash
fastx_quality_stats [-h] [-i INFILE] [-o OUTFILE]
[-h]            显示帮助信息
[-i INFILE]     指定输入文件, 如果输入为fasta, 只统计碱基分布
[-o OUTFILE]    指定输出文件
```

### 3.5. 生成碱基质量箱图

```bash
fastq_quality_boxplot_graph.sh [-i INPUT.TXT] [-t TITLE] [-p] [-o OUTPUT]
[-p]            生成图片格式为Generate PostScript (.PS), 默认PNG图片
[-i INPUT.TXT]  输入文件名, 应为fastx_quality_stats的结果
[-o OUTPUT]     输出文件名
[-t TITLE]      箱线图的标题
```

### 3.6. 生成碱基质量分布图

```bash
fastx_nucleotide_distribution_graph.sh [-i INPUT.TXT] [-t TITLE] [-p] [-o OUTPUT]
[-p]            生成图片格式为Generate PostScript (.PS), 默认PNG图片
[-i INPUT.TXT]  输入文件名, 应为fastx_quality_stats的结果
[-o OUTPUT]     输出文件名
[-t TITLE]      柱状图的标题
```

### 3.7.重命名read

```bash
fastx_renamer [-n TYPE] [-h] [-z] [-v] [-i INFILE] [-o OUTFILE]
[-n TYPE]    重命名类型:
             SEQ - 使用核苷酸序列作为名称
             COUNT - 使用简单的计数重命名
[-h]         显示帮助信息
[-z]         用GZIP压缩输出
[-i INFILE]  指定输入文件, 默认为标准输入
[-o OUTFILE] 指定输出文件, 默认为标准输出
```

### 3.8. 输出指定长度read（前后切除）

```bash
fastx_trimmer [-h] [-f N] [-l N] [-z] [-v] [-i INFILE] [-o OUTFILE]
[-h]            显示帮助信息
[-f N]          从前面第几个碱基开始保留, 默认值为1
[-l N]          保留几个碱基, 默认值为全部保留
[-z]            用GZIP压缩输出
[-i INFILE]     指定输入文件, 默认为标准输入
[-o OUTFILE]    指定输出文件, 默认为标准输出
```

### 3.9.去重reads

```bash
fastx_collapser [-h] [-v] [-i INFILE] [-o OUTFILE]
[-h]         显示帮助信息
[-v]         打印出输入/输出的统计摘要
[-i INFILE]  指定输入文件, 默认为标准输入
[-o OUTFILE] 指定输出文件, 默认为标准输出
```

### 3.10.去除人工reads（artifact）

```bash
fastx_artifacts_filter [-h] [-v] [-z] [-i INFILE] [-o OUTFILE]
[-h]         显示帮助信息
[-v]         报告处理的序列数目
                 如果有'-o', 在屏幕上标准输出
                 如果无'-o', 信息在屏幕上输出, 则报告在标准错误上输出
[-z]         用GZIP压缩输出
[-i INFILE]  指定输入文件, 默认为标准输入
[-o OUTFILE] 指定输出文件, 默认为标准输出
```

### 3.11. 生成反向互补序列

```bash
 fastx_reverse_complement -i input.fastq -o output.fastq
```

### 3.12.对fasta文件进行格式化

```bash
fasta_formatter [-h] [-i INFILE] [-o OUTFILE] [-w N] [-t] [-e]
[-w N]          输出fasta每行最大的碱基数, 默认值为0, 也就是每条fasta占据一行序列, 容易脚本操作
[-t]            输出为制表符格式, 而不是fasta, 第一列为序列ID, 第二列为序列的单行显示
[-e]            输出空白序列, 默认值为丢弃, 具有序列ID但并没有任何碱基的为空序列
```

### 3.13.DNA序列与RNA互换

```bash
fasta_nucleotide_changer [-h] [-z] [-v] [-i INFILE] [-o OUTFILE] [-r] [-d]
[-r]         DNA to RNA, T to U
[-d]         RNA to DNA, U to T
[-v]         报告序列的数目
             如果有'-o', 在屏幕上标准输出
             如果无'-o', 则在标准错误输出
[-z]         用GZIP压缩输出
```

### 3.14.生成reads长度分布图

```bash
fasta_clipping_histogram.pl INPUT_FILE.FA OUTPUT_FILE.PNG
INPUT_FILE.FA   指定输入fasta文件, 可以为GZIP格式
OUTPUT_FILE.PNG 生成的柱状图
```

### 3.15. 依据barcode进行reads分类成不同文件

```bash
fastx_barcode_splitter.pl --bcfile FILE --prefix PREFIX [--suffix SUFFIX] [--bol|--eol] [--mismatches N] [--exact] [--partial N] [--quiet] [--debug]
--bcfile FILE       Barcodes文件名, 解释如下
--prefix PREFIX     文件前缀, 会添加到输出文件, 可以用来辨别输出目录
--suffix            文件后缀, 可选, 可以用来辨别文件扩展名
--bol               尝试在序列的开始匹配barcodes, 也就是5‘端, 程序从0开始
--eol               尝试在序列的末尾匹配barcodes, 也就是3‘端, 程序从最后开始
                    --bol 与 --eo l必须有其一
--mismatches N      允许的最大错配数, 默认为1
--exact             等同于'--mismatches 0', 如果与'--mismatches'同时存在, 那么'--exact'具有优先级
--partial N         允许部分barcodes匹配, 默认不允许, 解释如下
--quiet             运行结束不打印统计信息, 默认为打印
--debug             在标准错误上打印大量有用的debug信息
```

## 4.Trimmomatic

> 你需要指定adapter的path,线程数,SRR的共有部分,及循环体.还要各种path

+ Paired end

```bash
adapter=
threads=
SRR_base=
for ((id=4;id<=8;id++));
do
time trimmomatic PE ../00rawdata/${SRR_base}${id}_1.fastq.gz  ../00rawdata/${SRR_base}${id}_2.fq.gz \
    ./trim/${SRR_base}${id}_F.fastq.gz ./trim/${SRR_base}${id}_UF.fastq.gz \
    ./trim/${SRR_base}${id}_R.fastq.gz ./trim/${SRR_base}${id}_UR.fastq.gz \
    ILLUMINACLIP:${adapter}:2:30:10 \
    -threads ${threads} \
    LEADING:4 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25
done
```

+ Single end


```bash
adapter=
threads=
SRR_base=
for ((id=4;id<=8;id++));
do
time trimmomatic SE ../00rawdata/${SRR_base}${id}.fastq.gz \
    ./trim/${SRR_base}${id}_F.fastq.gz \
    ILLUMINACLIP:${adapter}:2:30:10 \
    -threads ${threads} \
    LEADING:4 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25
done
```
