---
title: 扩增子测序 QIIME1
date: 2018-10-15
toc: true
categories: Bioinformatics
tags:
  - 软件和包
  - 扩增子
---

## 1 构建mapping file,并验证是否有误
mapping file记录着样品对应的barcode、primer、treatment group等信息，以tab作为列分割符。各列名如下：

column name | Description 
---|---
SampleID | 样品名，数字、字母、点号
BarcodeSequence | barcode序列，区分样本
LinkerPrimerSequence | 5'端引物
ReversePrimerSequence | 3'端反向引物,如果有的话
Treatment | 分组信息
Description | 样品详细注释
DOB| 日期或其他信息

下机reads组成：AdapterA - BarcodeSequence - LinkerPrimerSequence - Target - ReversePrimer - AdapterB

<!--more-->

```bash
# mapping.tsv
# SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	ReversePrimer 	Description
# Example mapping file for the QIIME analysis package.
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control 	GCGCACGGGTGAGTA		Control_mouse__I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control 	GCGCACGGGTGAGTA		Control_mouse__I.D._355
```


```bash
validate_mapping_file -m mapping.tsv -o 01_map_out
```

## 2 拆分和质控(demultiplex and QC)

扩增子测序通常是混池测序，所以我们要对序列拆分到各自所属样品，并进行质控。执行这一步的脚本为split_libaries_fastq.py，质控在里面默认进行。输入的文件格式可以是 fastq.gz 也可以是 fastq。

### 2.0 只进行质控
把split_libaries_fastq.py的参数--barcode_type 设置为'not-barcoded'即可。

### 2.1 Illumina paired end fastq
使用joined_paired_ends.py合并双端数据，如果后续进行拆分序列的话，你需要提供barcodes.fastq文件。这个由extract_barcodes来实现，见3.2.8.
joined_paired_ends内部使用fast-join或SeqPrep来进行，默认为前者，-m SeqPrep可以指定后者。
后续拆分同Illumina single end fastq，见3.2.2.

```bash
join_paired_ends.py -f forward_reads.fastq -r reverse_reads.fastq -b barcodes.fastq -o join_fastq/
```

### 2.2 Illumina single end fastq

```bash
# 如果有多个的话，用逗号分隔输入, 指定--rev_comp_mapping_barcodes是使用barcode.fastq.gz与mapping里面的方向反向互补了(3‘端的单端barcode)。
split_libraries_fastq.py -i read1.fastq.gz -b barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt -q 19
```

### 2.3 Illumina fastq：multiple lane

若想同时拆分Illumina的多个lane的数据，可以按下列命令进行，但是要指定每个lane对应的barcode和mapping file.所有的拆分结果将写入seqs.fna中，但log和histogram file是按lane区分. split_libraries_log.txt记录了运行的各项参数，如果结果不满意可调整参数重新运行。


```bash
split_libraries_fastq.py -i s_6_sequence.fastq,s_7_sequence.fastq,s_8_sequence.fastq \
        -o sl_out/ -b s_6_barcode.fastq,s_7_barcode.fastq,s_8_barcode.fastq \
        -m s_6_map.txt,s_7_map.txt,s_8_map.txt
```

### 2.4 Illumina qseq

首先你需要知道你的barcode位于哪个reads上，长度多少。可以查看一个样品的不同reads文件。
然后使用process_qseq.py生成fastq文件，每个lane生成一个single fastq文件。
最后才是使用split_libraries_fastq.py进行拆分，类似于2.2.3。


```bash
# 假设文件为s_1_1_0001_qseq.txt，s_1_2_0001_qseq.txt
head -n 1 s_1_1_0001_qseq.txt
# 如果结果类似如下，第八列表示是read2, 且barcode长度为12bases。
# 如果不是，查看s_1_2_0001_qseq.txt
M10     68      1       1       28680   29475   0       2       ACTCACGGTATTA   \_J\Sa^Y[ZYK`   0
M10     68      1       1       19607   29475   0       2       AGACTGAGTACTA   PP\JJ\JQ`\RK^   1
```


```bash
# 使用process_qseq.py生成fastq
# read1
process_qseq.py -i input_dir/ -o output_dir -r 1
# read2,并且截断了12个bases,12是因为你的barcode长12 bases；只处理barcode序列
process_qseq.py -i input_dir/ -o output_dir -r 2 -b 12
```

### 2.5 iseq
使用process_iseq.py进行处理, 每个lane将产生一个single read；后续使用split_libraries_fastq.py进行拆分，类似于2.2.3。

+ 如果barcode与序列在一起(较为常见)

```bash
# s_6_1_sequences.txt,s_7_1_sequences.txt都是含barcode的；barcode长度为12
# HWI-ST753_50:6:1101:15435:9071#0/1:ACCAGACGATGCTACGGAGGGAGCTAG...
process_iseq.py -i s_6_1_sequences.txt,s_7_1_sequences.txt -o out_dir/ -b 12
```

+ 如果barcode在header里面

```bash
# ACAGCTA长为7,但是你的barcode为ACAGCT为6，可以指定--barcode-length
# HWI-6X_9267:1:1:12:410#ACAGCTA/1:TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCG
process_iseq.py -i s_6_1_sequences.txt,s_7_1_sequences.txt -o out_dir/ --barcode-length 6 --barcode_in_header
```

### 2.6 LEA-Seq数据
使用split_libraries_lea_seq.py来处理.

### 2.7 处理已经拆分的数据

+ multiple_split_libraries_fastq.py: 进行质控
+ multiple_join_paired_ends.py:	合并已经拆分的paired ends
+ multiple_extract_barcodes.py: 去除已经拆分的fastq内部的barcode

### 2.8 处理不同的barcode策略

+ single fastq starts with the barcode sequence(单个fastq单个barcode)

```bash
# 假设in_seqs.fastq的前10个碱基是barcode
# 如果你想反向互补barcode的话，设置--rev_comp_bc1
extract_barcodes.py -f in_seqs.fastq --bc1_len 10 -o parsed_barcodes/ --input_type barcode_single_end
```

+ two fastqs each starts with part of barcode: paired-end(双端fasta,单个barcode)

```bash
# 假设barcode部分分布在reads1,6 bases；部分分布在reads2,8 bases。
extract_barcodes.py --input_type barcode_paired_end -f reads1.fastq -r reads2.fastq --bc1_len 6 --bc2_len 8 -o parsed_barcodes/
# 如果你想反向互补barcode的话，--rev_comp_bc1, --rev_comp_bc2

# 有时候，reads的方向不确定，你需要根据mapping file里面的primer来确定reads的方向
extract_barcodes.py --input_type barcode_paired_end -m mapping_file.txt \
        -a -f reads1.fastq -r reads2.fastq \
        --bc1_len 6 --bc2_len 8 -o parsed_barcodes/
```

+ two index/barcode reads and two fastq reads (双端fastq, 2个barcode)

```bash
# 这里已经得到分离的barcode了，所以用不上read文件
# index1.fastq,index2.fastq是barcode文件。每个barcode长度为6.
extract_barcodes.py --input_type barcode_paired_end -f index1.fastq \
        -r index2.fastq --bc1_len 6 --bc2_len 6 -o parsed_barcodes/
```

+ single stitched read with barcodes on each end (单个fastq, 双端barcode)

```bash
extract_barcodes.py --input_type barcode_paired_stitched -f reads.fastq --bc1_len 6 --bc2_len 6 -o parsed_barcodes/
```

+ barcode在fastq文件的描述行内

```bash
# @MCIC-SOLEXA_0051_FC:1:1:4065:1039#CGATGT/1
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
# +MCIC-SOLEXA_0051_FC:1:1:4065:1039#CGATGT/1
# KPPPQWWWWWQQ________BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# 这里barcode开头是#，长度为6
extract_barcodes.py --input_type barcode_in_label --char_delineator "#" -f in_seqs.fastq --bc1_len 6 -o parsed_barcodes/
```

### 2.9 处理Joint Genome Institute fastq 文件

```bash
# @MISEQ03:64:000000000-A2H3D:1:1101:14358:1530 1:N:0:TCCACAGGAGT
# TNCAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTA
# + 
# .....

# TCCACAGGAGT为barcode，长为11
extract_reads_from_interleaved_file.py -i reads_to_extract.fastq -o extracted_reads
extract_barcodes.py -f extracted_reads/forward_reads.fastq -o barcode_reads -c barcode_in_label --char_delineator '1:N:0:' -l 11
```

## 3 去除嵌合体和OTU聚类

[pick_otus.py](http://qiime.org/scripts/pick_otus.html?highlight=pick_otus)整合十多种算法用于聚类OTU。
输出结果的第一列为OTU ID，第二列为该OTU的种子序列。以usearch为例,程序内部执行的基本过程如下：

```bash
#1. 根据序列长度排序(rank)          2. 去除重复序列(去冗余)
#3. 根据序列丰度进行排序            4. 对序列进行去噪
#5. 去除嵌合体,先是de novo, 然后是ref-based
#6. 对序列进行合并                  7. 再根据丰度进行拍醋
#8. 基于排序进行聚类                9. 为聚类的结果添加对应的序列id
#10.对序列完成分类
```


```bash
pick_otus.py -i split_library_output/seqs.fna -o picked_otus_default
# 较为重要的参数
# -i <input_fasta>
# -m <method>		使用blast, uclust_ref, usearch_ref, usearch61_ref时需要指定-r参数.默认为uclust.
# -o <out_dir>, -r <ref_seq>
# -s <similarity>	默认为0.97
# -j <ErrorRate>
# -g <min_count_in_cluster> 当cluster的序列count数低于设定值时，该cluster被忽略
# -f <chimeric_ref_seq>		当-m usearch时需要指定

# uclust-ref
pick_otus.py -i seqs.fna -r refseqs.fasta -m uclust_ref --denovo_otu_id_prefix qiime_otu_

# cdhit
# -n 100将会执行预聚类，前100个bp完全相同的序列将被预聚类;使用-t参数则不指定-n num
pick_otus.py -i seqs.fna -m cdhit -o cdhit_picked_otus/ -n 100

# blast
# 如果参考序列时预构建的blast序列,使用-b代替-r进行指定; -s指定相似度; -e指定错误发现率(blast专有)
pick_otus.py -i seqs.fna -o blast_picked_otus_90_percent/ -m blast -r refseqs.fasta -s 0.90 -e 1e-30

# prefix-suffix
# 前50个bp后25个bp完全相同的序列将会被预聚类
pick_otus.py -i seqs.fna -o prefix_suffix_picked_otus/ -m prefix_suffix -p 50 -u 25
```


```bash
# mothur
# -c 指定聚类算法, mothur专有,取值furtherest, neareast, average
pick_otus.py -i seqs.aligned.fna -o mothur_picked_otus_nn/ -m mothur -c nearest

# usearch
# --suppress_reference_chimera_detection禁止基于参考序列的嵌合体发现
pick_otus.py -i seqs.fna -m usearch --word_length 64 --db_filepath refseqs.fasta -o usearch_qf_results/ --minsize 2

# swarm
pick_otus.py -i $PWD/seqs.fna -m swarm -o $PWD/picked_otus_swarm
```

## 4 选取OTU代表性序列

输出为代表序列out_rep_seq.fasta。


```bash
pick_rep_set.py -f split_library_output/seqs.fna -m most_abundant -i seqs_outs.txt -o out_rep_seq.fasta

# -i <input_otus_file>
# -f <all_seqs_fna>
# -m <method> 	random, longest, most_abundant, first(default)
```

## 5 进行物种分类

使用assign_taxonomy.py进行物种分类；对于给定的序列集，该脚本会尝试对每条序列确定其分类，当前支持方法有BLAST、RDP分类器、RTAX、mothur和uclust。输出目录默认为methodname_assigned_taxonomy, 输出的结果文件第一列是序列ID，第二列是该序列所属物种分类，第三列是质量分数，以及其他信息，依据方法不同而不同。使用的参考序列最新版本可在qiime的resources页面下载。


```bash
# default：uclust
assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt
# SortMeRNA
assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt -m sortmerna
# Blast
# -e maximum e-value to record an assignment
assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt -e 0.01 -m blast
# RDP classifier
# 训练集默认为-t和-r指定的taxonomy和otus；如未指定，使用默认设置
# -c minimum confidence to record an assignment
assign_taxonomy.py -i repr_set_seqs.fasta -m rdp -c 0.80
# RTAX
assign_taxonomy.py -i rtax_repr_set_seqs.fasta -m rtax \
 		--read_1_seqs_fp read_1.seqs.fna \
 		--read_2_seqs_fp read_2.seqs.fna \
 		-r rtax_ref_seq_set.fna \
 		-t rtax_id_to_taxonomy.txt
# mothur
assign_taxonomy.py -i mothur_repr_set_seqs.fasta -m mothur -r mothur_ref_seq_set.fna -t mothur_id_to_taxonomy.txt
```


```bash
# -o <out_dir> 				输出结果路径
# -i <otu_rep_set.fasta>	输入的代表序列
# -t <id_to_taxonomy> 		分类时参考序列ID对应的物种分类信息。
# -r <ref_otu.fasta>		分类时的参考序列
# -p <rdp_train_set> 		指定RDP的训练数据集，会覆盖-r,-t指定的参数
# -m <method> 				选取某个方法进行分类
# -b <blast_db> 			指定使用blast时的数据库, 与-r同义
# -e <max_E_value> 			使用blast时的最大错误值
# -c <min_confidence>		RDP分类时的最小置信值,超过该值才会认为属于该分类
```

---

## 6 生成OTU表(biom文件)

make_otu_table.py将对OTU进行汇总和注释，生成biom格式的OTU表。

### 6.1 生成OTU表


```bash
make_otu_table.py -i otu_map.txt \
 		-t tax_assignments.txt \
 		-o otu_table.biom \
 		-m mapping_file.txt

# -i <input_otu_map>   	指定输入的OTU map，通常是pick_otus.py的结果文件
# -o <out_biom_dir> 	输出结果文件
# -m <mapping_file> 	指定mapping file,如果有的话
# -e <id_file>			指定排除OTU的identifiers文件
```

### 6.2 生成OTU表同时进行过滤

为了便于后续操作，一般会对OTU表进行过滤等操作.有两种方式，第一种方式是在构建的时候进行过滤，这通常是过滤嵌合体或未比对序列.另一种方式是生成otu表后，再进行过滤。其他过滤操作见3.7。


```bash
# identify_chimeric_seqs.py, align_seqs.py 生成不要的seqs对应的文件exclude.txt，然后在make_otu_table.py中以-e参数指定exclude.txt
```

## 7 对OTU表(biom)进行操作: 过滤、拆分等

### 7.1 过滤

+ Even sampling, rarefaction (稀释)

这里的抽样方式有单次抽样、多次等差抽样(稀释抽样)、多次平均抽样。single_rareefaction.py只进行一次抽样, multiple_rarefactions.py进行多次抽样(每次抽样深度不同), multiple_rarefactions_even_depth.py进行多次抽样(每次抽样深度相同)。


```bash
# 从otu-table中按每个样品100条序列随机抽取，并输出
single_rarefaction.py -i otu_table.biom -o otu_table_even100.biom -d 100
# -k <bool>		指定是否保留全0的OTU,default:fasle
# --subsample_multinomial <bool> 是否是放回抽样,default:fasle

multiple_rarefactions.py -i otu_table.biom -m 10 -x 140 -s 10 -n 2 -o rarefied_otu_tables/
# -m <min_num>		指定多次抽样的最低抽样数,seqs/sample, 抽样深度
# -x <max_num>		指定多次抽样的最大抽样数
# -s <step_num>		指定2次抽样的抽样次数的步长,差值
# -n <num_reps> 	每一步抽样的重复次数
# -k，--subsample_nultinomial 同单次抽样

multiple_rarefactions_even_depth.py -i otu_table.biom -o rarefied_otu_tables/ -d 100 -n 10
```

+ 过滤OTUs/observations

过滤时可以根据mapping file的某列值进行过滤，也可以根据OTU的count数来过滤(singletons, doubletons...),丰度过滤。 根据run_number的切分和过滤有点相似，但是切分是所有run切分。过滤可以只过滤某个值。


```bash
# 保留mapping file中某列的值等于指定值的OTUs, "Sample_Type:*,!Control_Blank"则表示删除
filter_samples_from_otu_table.py -i split_otu_tables/otu_table_1.biom -o otu_table_run1_blank_samples.biom -m map.txt -s "Sample_Type:Control_Blank"

# 排除singleton
filter_otus_from_otu_table.py -i otu_table.biom -o filtered_otu_table.biom -n 2

# -n <min_count>					如果该OTU的total counts低于指定值,排除该OTU
# -x <max_count> 					如果该OTU的total counts高于指定值,排除该OTU
# -s <min_samples>					如果属于该OTU的样品数低于指定值,排除该OTU
# -y <max_samples> 					如果属于该OTU的样品数超过指定值,排除该OTU
# -e <out_ids_to_exclude_file>		排除指定id的OTUs
```

+ 过滤样本

过滤样本可以基于丰度过滤，也可以基于metadata过滤，还可以基于ID过滤。过滤样本与过滤OTUs是相类似的.


```bash
# 保留total OTUs counts在[100,150]的样品
filter_samples_from_otu_table.py -i otu_table.biom -o filtered_otu_table.biom -n 100 -x 150

# 保留mapping file中某列的值等于指定值的样品
filter_samples_from_otu_table.py -i otu_table.biom -o filtered_otu_table.biom -m map.txt -s 'Treatment:Control'

# 保留mapping file中某列的值==不等于==指定值的样品
filter_samples_from_otu_table.py -i otu_table.biom -o filtered_otu_table.biom -m map.txt -s 'Treatment:*,!Control'

# 保留在ids.txt的样品; 如果设定了--negate_sample_id_fp,则是删除在ids.txt的样品
filter_samples_from_otu_table.py -i otu_table.biom -o filtered_otu_table.biom --sample_id_fp ids.txt
```

+ 过滤分类群: 保留属于特定分类的OTU

过滤分类群时，你可以有3种逻辑的过滤方式：or, not or, not but


```bash
# 保留属于__Bacteroidetes 或者 p__Firmicutes
filter_taxa_from_otu_table.py -i otu_table.biom -o otu_table_bac_firm_only.biom -p p__Bacteroidetes,p__Firmicutes
# 保留不属于__Bacteroidetes和 p__Firmicutes
filter_taxa_from_otu_table.py -i otu_table.biom -o otu_table_non_bac_firm.biom -n p__Bacteroidetes,p__Firmicutes
# 保留不属于c__Clostridia但是属于p__Firmicutes
filter_taxa_from_otu_table.py -i otu_table.biom -o otu_table_all_firm_but_not_clos.biom -p p__Firmicutes -n c__Clostridia

# -p <list, comma> 保留
# -n <list, comma> 排除
# --metadata_field 指定用于过滤的metadata的列名
```

### 7.2 拆分和合并

+ 基于mapping/metadata进行拆分


```bash
split_otu_table.py -i otu_table.biom -m Fasting_Map.txt -f Treatment -o per_study_otu_tables
split_otu_table.py -i otu_table.biom -m Fasting_Map.txt -f Treatment,Color -o ./per_study_otu_tables/
```

+ 基于taxonomy(物种分类)进行拆分


```bash
# 将otu_table依据Level3进行拆分
split_otu_table_by_taxonomy.py -i otu_table.biom -L 3 -o ./L3/
```

+ 合并

在进行合并时，你要保证合并的两个OTUs表的OTU id和sample id是兼容的；如果两个表有重叠内容，将会计算重叠内容的值的和。


```bash
merge_otu_tables.py -i otu_table1.biom,otu_table2.biom -o merged_otu_table.biom
```

### 7.3 统计汇总


```bash
# 按样品进行统计汇总,取总的observations的所有counts
biom summarize-table -i rich_sparse_otu_table.biom -o rich_sparse_otu_table_summary.txt
# 按样品进行汇总，但只取唯一的observations的counts
biom summarize-table -i rich_sparse_otu_table.biom --qualitative-o rich_sparse_otu_table_summary.txt
```

### 7.4 对OTU表进行排序


```bash
# 按照sample id进行排序, 字母排序
sort_otu_table.py -i otu_table.biom -o sorted_otu_table.biom
# 按照mapping/metadata file的某个列排序
sort_otu_table.py -i otu_table.biom -o dob_sorted_otu_table.biom -m Fasting_Map.txt -s Age
# 按照给定顺序的sample id进行排序
sort_otu_table.py -i otu_table.biom -o sorted_otu_table.biom -l sample_id_list.txt
```

### 7.5 biom文件的其他处理

+ 添加metadata信息到biom文件中

见[网页tutorial](http://biom-format.org/documentation/adding_metadata.html)

+ 转换biom到其他格式


```bash
# json/hdf5
biom convert -i table.txt -o table.from_txt_json.biom --table-type="OTU table" --to-json
biom convert -i table.txt -o table.from_txt_hdf5.biom --table-type="OTU table" --to-hdf5
# tsv, 并添加taxonomy物种分类信息; 指定--output-metadata-id "ConsensusLineage"可以把添加的taxonomy列名改为ConsensusLineage
biom convert -i table.biom -o table.from_biom.txt --to-tsv --header-key taxonomy
# 把tsv转换回去
biom convert -i otu_table.txt -o new_otu_table.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
```

## 8 物种分类分析 /群落组成分析(otu)

### 8.1 OTU聚类分析和分类学分析

此处是OTU聚类表(OTU table),OTU在每个样品中的分布信息就是OTU聚类分析；后续添加的分类注释信息就是分类学信息。


```bash
#1. 生成OTU聚类表, 见2.6
#2. 把表转换为tsv格式, 见2.7.5
```

### 8.2 OTU物种注释统计、OTU物种组成图

使用summeriza_taxa_through_plot.py对OTU聚类表进行统计, 以各个分类水平进行统计每个样品在该分类水平上的OTU总数,并进行可视化. 该脚本由summarize_taxa.py和plot_taxa_summary.py组成。生成的otu_table_L\*.txt文件是各个Level的水平物种丰度统计表。taxa_summary_plots文件夹则是物种丰度绘制图。


```bash
summeriza_taxa_through_plots.py -o taxa_summary -i otu_table_w_tax.biom -m Fasting_Map.txt

# -i <input_otu_table> 	指定input文件, otu_table需要包含物种信息
# -p <param_file> 		指定特定中间过程脚本的参数
# -m <mapping_file>		如果指定了-c,-m必须指定
# -c <mapping_category>	指定进行统计的分组信息(按组统计)
# -s <bool> 			指定是否排序otu table,default: false

# summeriza_taxa_through_plots.py内部依序执行的命令
summarize_taxa.py -i otus/otu_table.biom -o taxa_summary
plot_taxa_summary.py -i taxa_summary/otu_table_L2.txt,taxa_summary/otu_table_L3.txt,taxa_summary/otu_table_L4.txt,taxa_summary/otu_table_L5.txt,taxa_summary/otu_table_L6.txt -o taxa_summary/taxa_summary_plots/
```

### 8.3 OTU进化树分析(系统发育树)

构建OTU系统发育树分为三个步骤，分别使用三个脚本进行.生成的树可以用meta、Figtree等软件进行查看。系统发育树的构建有多种方法,可见[biotrainee](http://www.biotrainee.com/thread-2253-1-1.html)

1. 序列比对

可以选择多种方法进行序列比对,生成的文件夹命名为method_aligned的文件夹里,含有method_aligned.fasta, method_failures.fasta, method_log.txt


```bash
align_seqs.py -i otu_rep_set.fasta -o pynast_aligend/

# -m <align_method> 	指定比对序列的算法, pynast(default), infernal, clustalW, muscle, mafft
# -a <pairwse_align> 	指定在PyNAST中成对比对的算法, muscle, pair_hmm, clustal, blast, uclust(default), mafft
# -e <min_length>		进行比对的序列最短长度
# -p <min_percent_id> 	blast中最低的一致性占比(default:0.75)
# -d <blast_db>			当-m pynast时, 进行blast的数据库
# -t <template_align> 	指定预比对的文件;可在http://greengenes.lbl.gov/下载获得
# --muscle_max_memory

align_seqs.py -i $PWD/unaligned.fna -t core_set_aligned.fasta.imputed -o $PWD/pynast_aligned/ -e 500 -p 95.0
```

2. 过滤比对结果中的高变异区


```bash
filter_alignment.py -i pynast_aligned/otu_rep_set_aligned.fasta -o pynast_aligned/ --remove_outliers -g 0.95

# -m <lane_mask_file> 	指定在构建树时保留位点的lanemask文件
# -s <suppress_lane_mask_filter>
# -g <allow_grap_frac> 	指定某个位点所允许的最少gap数
# -t <threshold> 		指定认定为outlier的偏离标准偏差倍数, default=3
# -r <remove_outliers>	指定过滤outlier
```

3. 生成树


```bash
make_phylogeny.py -i pynast_aligned/otu_rep_set_aligned_pfilter.fasta -o rep_set.tre

# -t <tree_method> 		指定进行建树的方法, clustalw, raxml_v730, muscle, fasttree(default), clearcut
# -l <log_fp>			log的存储位置, 不指定的话不会创建log
# -r <root_method> 		指定构建树的根(root)的方法, midpoint, tree_method_default(default)
```


### 8.4 OTU Heatmap绘制（物种分类热图)

使用make_otu_heatmap.py就对OTU表绘制热图，热图的每一行对应着一个OTU，每一列对应着一个样本。颜色代表了对应样本(列)的对应OTU(行)的丰度(counts).默认情况下，OTUs(行)将会使用UPGMA层次聚类，样品(列)将会按照在OTU表中出现的顺序进行排列。用户可以按照聚类树对行/列进行排序，也可以按照mapping file排序样品，对于后者，样品将会进行分组并聚类。


```bash
# default
make_otu_heatmap.py -i otu_table.biom -o heatmap.pdf

# -t <otu_tree> 	指定用于OTUs排序的OTU tree
# -m <mapping_file> 指定用于样品排序的meta/mapping file
# -c <category>		指定用于样品排序的mapping file 列(分组)。
# -s <sample_tree> 	指定用于样品排序的sample_tree,可来源于upgma_cluster.py的输出
# -g <image_type> 	输出热图的格式, png, svg, pdf(default)
# --no_log_transform 不进行log转换, 0值被设为最小值的1/2
# --suppress_row/column_clustering 不进行upgma聚类,-t/-m设定了,不聚类的设定被忽略
# --absolute_abundance 不进行百分比化
# --color_scheme 	指定颜色主题
# --width,--height,--dpi 指定热图图片相关格式
# --obs_md_category 	observation metadata category to plot
# --obs_md_level 	The level of observation metadata to plot for hierarchical metadata
```

### 8.5 OTU network分析

make_otu_network.py可以利用OTU和metadata/mapping file生成用于cytoscape可视化的网络文件，并对网络进行统计。进行网络分析是因为他可以展示和分析样品之间OTUs的分布，这对于展示大数据中样本的相似性和差别十分便捷。可视化结果展示的是根据样品之间共有OTUs数目的聚类结果，网络的加权值是每个OTU所含有的序列数。然后根据网络节点间的关联来统计其聚类模式，这其中用到G-test(我不懂啊)。

基于network分析的样品比较与基于tree的PCoA分析相得益彰。在多数研究中，这两个方法的结果总是相一致的，但是反映这数据的不同方面。network分析提供了系统发育的可视化呈现，PCoA-unifrac聚类提供了发育的聚类水平，这有可能在network观察到。从技术上讲，OTUs和样品是network中的两种节点, OTU-node通过edges同sample-node连接，这些edges就是samples在OTUs中的序列。


```bash
make_otu_network.py -i otu_table.biom -m Fasting_Map.txt -o otu_network
# -b <color_by> 		指定依据mapping_file的特定列进行着色,可逗号分隔指定多个
# -k <bg_color> 		指定网络的背景色
```

## 9 α多样性分析 (alpha diversity analysis)

α多样性的分析主要是通过计算[各项α多样性指数](http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html)来表征α多样性的, 比如Shannon指数、ACE指数、Chao1指数等。

### 9.1 计算指数和对应的指数稀释曲线   

alpha diversity index and rarefaction curve. 在alpha_rarefaction.py中，可以生成稀释OTU table，基于稀释表计算α多样性矩阵，然后画出稀释曲线图。可计算的α多样性指数有：ace, berger_parker_d, brillouin_d, chao1, chao1_ci, dominance, doubles, enspie, equitability, esty_ci, fisher_alpha, gini_index, goods_coverage, heip_e, kempton_taylor_q, margalef, mcintosh_d, mcintosh_e, menhinick, michaelis_menten_fit, observed_otus, observed_species, osd, simpson_reciprocal, robbins, shannon, simpson, simpson_e, singles, strong, PD_whole_tree


```bash
# 指定需要计算的α多样性指数
echo 'alpha_diversity:metrics simpson,shannon,PD_whole_tree,chao1,observed_species' > alpha_params.txt
alpha_rarefaction.py -i otu_table.biom -o arare_max100/ -t rep_set.tre -m Fasting_Map.txt -e 100 -p alpha_params.txt
# -p <parameter_fp> 		指定过程脚本对应的参数
# -n <num_steps> 			指定进行抽样时的步数,不是步长, 是生成稀释OTU表的大小
# -w 						只打印过程脚本，不运行
# -a 						指定进行多线程
# -t <tree_file>			指定计算系统发育指数的tree
# --min_rare_depth <num>	指定抽样时的最小深度
# -e,--max_rare_depth <num> 指定抽样时的最大深度
# -O <jobs_to_start> 		指定跳转到第几步开始执行,重新执行时十分有用
# -f 						指定强制覆写文件
# --retain_intermediate_files 	指定保留中间过程文件

multiple_rarefactions.py -i ../04_pick_open_ref_otus/otu_table_mc2.biom \
		-m 10 -x 100 \
		-s 9 -o ./alpha/rarefaction/
alpha_diversity.py -i ./alpha/rarefaction/ 
		-o ./alpha/alpha_div/ \
		--metrics simpson,shannon,PD_whole_tree,chao1,observed-species \
		-t ../04_pick_open_ref_otus/rep_set.tre
collate_alpha.py -i ./alpha/alpha_div/ \
		-o ./alpha/alpha_div_collated/
rm -r ./alpha/rarefaction/ ./alpha/alpha_div/
make_rarefaction_plots.py -i ./alpha/alpha_div_collated/ \
		-m ../map.tsv \
		-o ./alpha/alpha_rarefaction_plots/
```

### 9.2 Rank-abundance曲线


```bash
# 比较多个样品的排序-丰度曲线图；如果想比较一个OTU表里面的所有样品, -s '*' 即可
plot_rank_abundance_graph.py -i otu_table.biom -s 'PC.354,PC.481,PC.364' -x -v -o rank_abundance.pdf
# -s <sample_name>		画图的样品名
# -a,--absolute_counts 	指定图中使用的数据是绝对丰度,
# -n,--no_legend 		指定图中不画图例
# -x,--x_linear_scale 	指定x轴要进行线性缩放
# -y,--y_linear_scale 	指定y轴要进行线性缩放
# -f <file_type> 		指定输出图的类型, pdf(default), svg, png, eps
# -v 					输出过程信息
```


### 9.3 α多样性指数差异分析

基于样品的多样性指数，可以检验组间样品的alpha diversity是否存在显著差异，检验方法可以是Wilcoxon秩和检验(2组样品)或者Kruskal-Willis秩和检验(多于两组).



```bash
# 汇总计算的所有alpha diversity index, 如果是执行了3.9.1,则不必执行本命令;
collate_alpha.py -i alpha_div/ -o alpha_div_collated/

# 比较alpha diversity index，按SampleType作为分组，得到每个diversity index下的成对组间统计
ls alpha/alpha_div_collated/*.txt | while read div;
do
name=`basename $div .txt`
compare_alpha_diversity.py -i $div -m ../map.tsv -c SampleType -d 100 -o compare_alpha_div/$name
done
# 如果想把各个diversity index合在一张图,需要用额外的R代码
```

### 9.4 物种累积曲线  (species accumulation curves)
```
# 可以使用R的vegan包
```

## 10 β 多样性分析  (beta diversity analysis)

[微生物群落差异分析方法大揭秘](http://www.omicshare.com/forum/thread-3251-1-1.html)

### 10.0 一步法进行PCoA分析

beta_diversity_through_plot.py可以进行β多样性分析、主坐标轴分析，并画出3D PCoA图. 首先该脚本会随机对输入的otu_table.biom进行抽样，以达到抽平每个样品的序列数目的目的。然后计算各个表征β多样性的距离矩阵。最后基于这些距离矩阵进行主坐标轴分析，并对结果进行可视化。


```bash
beta_diversity_through_plots.py -i otu_table.biom -m map.tsv -o beta_diversity -e 100
# -t <tree_fp> 		指定树的文件
# -p <param_fp> 	指定中间命令的参数文件
# -f,-w 			强制覆写和只打印命令不运行
# -a 				指定多线程运行
# -e <seqs_per_sample> 进行样品抽平的覆盖深度
# --supress_emperor_plots
# -O,--jobs_to_start 跳转到第几步开始运行。
```

### 10.1 计算β多样性距离矩阵

beta_diversity.py支持计算多种距离矩阵, 使用-s参数查看支持的矩阵.由于unifrac，还需要提供发育树作为输入。一般来讲，由于unifrac有用到系统发育信息，所以一般推荐使用unifrac。定量的距离矩阵(weighted)更适合用于揭示由于群落物种丰度引起的群落差异,而定性的距离矩阵更适用于因环境不同所导致的群落差异。

**2.10.2-2.10.5等分析均基于本步骤生成的距离矩阵进行后续的分析。**


```bash
# 查看支持的矩阵类型,然后选择想要计算的距离矩阵
beta_diversity.py -s
metrics=unweighted_unifrac,weighted_unifrac,bray_curtis,binary_jaccard
beta_diversity.py -i 04_pick_open_ref_otus/otu_table_mc2.biom \
   		-o 05_beta_diversity \
   		--metrics $metrics \
   		-t 04_pick_open_ref_otus/rep_set.tre

mv 05_beta_diversity/weighted_unifrac_otu_table_mc2.txt 05_beta_diversity/weighted_unifrac_distance_matrix.txt

# -i <otu_table> 	指定输入otu_table, 可以以文件夹的形式指定多个
# -r <rows> 		指定只计算指定行的距离矩阵,应该传入样品名,例如'S1,S3'之类的
# -m <metrics> 		指定要计算的beta diversity matrix类型, 可以用逗号分隔指定多个,默认: unweighted_unifrac,weighted_unifrac
# -s <show_metrics> 显示支持的距离矩阵类型
# -t <tree_file> 	指定输入的系统发育树
# -f <full_tree> 	指定使用全部的系统树;OTU table与系统树不对应的话,使用全部整个系统树会导致很多额外的结果.
```

### 10.2 PCoA分析


```bash
# -i,-o还可以是文件夹用以进行批量计算
principal_coordinates.py -i 05_beta_diversity/weighted_unifrac_dm.txt \
		-o 05_beta_diversity/weighted_unifrac_pc.txt

# 画2d图
make_2d_plots.py -i weighted_unifrac_pc.txt -m map.tsv -o 2d_plots
# -b <color_by>			以逗号分隔的metadata分类/组列名,用以指定着色分组。可以使用'col1&&col2'指定一致着色
# -k <bg_color>
# --ellipsoid_opacity 	指定透明度,只在jackknifed beta diversity中有用
# --ellipsoid_method	指定绘制的数据类型, IQR (default) 和 sdev, 只在jackknifed beta diversity中有用
# --master_pcoa			
# --scree 				绘制碎石图

# 画3d图
make_emperor.py -i 05_beta_diversity/weighted_unifrac_pc.txt \
		-o 05_beta_diversity/weighted_unifrac_emperor_pcoa_plot/ \
		-m map.tsv

make_emperor.py -i weighted_unifrac_pc.txt -m map.tsv \
		-b 'Treatment&&DOB,Treatment' \ # 依据metadata进行着色
		-a DOB \ # 指定一个轴为DOB
		-x 'DOB:20060000' \ # 当自定义轴DOB出现缺失值或值错误时,指定为20060000
		-o emperor_colored
# 见https://biocore.github.io/emperor/build/html/scripts/make_emperor.html
```

### 10.3 PCA分析



```bash
PCA分析基于OTU表，运用方差分解，将样本间的差异反映在二维坐标图上，坐标轴是两个主成分。样本组成越相似在图中则越聚集。
待续.
```

### 10.4 NMDS分析

与 PCoA 一样， 非度量多维尺度分析（ NMDS, nonmetric multidimensional scaling） 可以基于任何类型的距离/非相似性矩阵对对象（样本） 进行排序，其区别在于 NMDS 不再是特征根排序技术，也不再以排序轴承载更多的方差为目的，因此 NMDS 排序图可以任意旋转、中心化和倒置。 NMDS 分析在多维空间内构建对象的初始结构， 并用迭代程序不断地调整对象位置， 目标是尽可能的最小化应力函数（ stress function， 取值 0~1），应力函数是排序空间内对象结果与原始距离矩阵之间相异程度的度量。如果预先设定的排序轴数量比较少（如二维、三维），则在相同轴数的情况下， NMDS往往能够获得比 PCoA 更少失真的对象之间的关系。

R包vegan可以进行NMDS分析和作图。

```
nmds.py -i unweighted_unifrac_dm.txt -d 3 -o unweighted_unifrac_nmds.txt
# -i,-o 			如果-i指定的是文件夹，那么就进行批量分析
# -d <dimensions> 	指定NMDS分析的排序轴数量
```


### 10.5 ADONIS、ANOSIM分析

通过compare_categories.py，我们可以比较不同组样品的聚类结果是否有显著性差异，该脚本提供多种方法，根据方法不同输出结果不同，但大多数是tsv文件。

**adonis**, 又叫permanova, 是一种非参数检验方法。能够给出不同分组因素对样品差异的解释度（R值）与分组显著性（P值），它首先计算数据的相关性重心，然后计算这些重心点的平方偏差(R^2);最后基于原始数据的排列平方和进行F-test; 它以距离矩阵和mapping file作为输入，比如Unifrac聚类矩阵，通过指定分组列来进行比较。

**anosim**,相似性分析，也是检验多组样本之间是否具有显著性差异，与adonis相似。但是它只对类别性变量有用，如果数据是连续型变量，推荐使用Mantel。它的计算公式为R=(rb-rw)/(N(N-1)/4)。其中rb指组间的所有距离的平均排名，rw指组内所有距离的平均排名。

adnois和anosim分析都只能衡量组与组之间整体是否有显著差异, 无法找到是哪个具体的OTU导致的该差异。


```bash
# 看HOST_SUBJECT_ID分组的组间结果是否会有差异
compare_categories.py --method adonis -i unweighted_unifrac_dm.txt -m map.txt -c HOST_SUBJECT_ID -o adonis_out -n 999
# --method <method> 		指定统计方法, adonis,anosim
# -i  <distance_matrix> 	指定输入的距离矩阵
# -m <mapping_file> 		指定输入的mapping file
# -c <group/category>		指定进行差异分析的分组方式
# -n <permutations_num>		指定进行计算p-value的排名数(前多少名)
```

### 10.6 样本聚类树图 

样本聚类树可以从整体上描述和比较样本/分组见的相似性和差异性. 采用的聚类方法是UPGMA (unweighted pair group method with arithmetic mean, http://www.nmsr.org/upgma.htm), 这里以使用unweighted_unifrac聚类矩阵为例.

```
# unweighted_unifrac_dm.txt来源于3.10.1
upgma_cluster.py -i unweighted_unifrac_dm.txt -o unweighted_unifrac_cluster.tre

# 如果要探究生成的聚类树的健壮性的话，就进行jackknifed_beta_diversity分析
# otu_table.biom是unweighted_unifrac_dm.txt的计算来源table
jackknifed_beta_diversity.py -i otu_table.biom -o bdiv_jk100 -e 100 -m map.tsv -t rep_set.tre
```

### 10.7 Lefse分析

线性判别分析（ LDA）效应量方法
http://www.omicshare.com/forum/forum.php?mod=viewthread&tid=563&extra=page%3D1%26filter%3Dtypeid%26typeid%3D31


```bash
conda install lefse -y
conda install rpy2=2.8.6 -y
# 如果出现rpy2与matplotlib有conflict,删除matplotlib,先安装rpy2,后安装matplotlib
conda remove matplotlib --force
# 如果有某个包总是断网下载不下来（得到该包的版本）,先安装这个包，再去安装rpy2
# https://anaconda.org/，下载对应包的对应版本
conda install --use-local package_name

# 1. 先转化biom为TSV文件, 注意--header-key是taxonomy,非Taxonomy
biom convert -i illumina/04_pick_open_ref_otus/otu_table_mc2_w_tax.biom -o otu_consensus.txt --to-tsv --header-key taxonomy --output-metadata-id "Consensus Lineage"

# 2. 转化为lefse的初始格式,设定分组、亚组和个体信息; 这里十分重要
qiime2lefse.py --in otu.txt \
 			   --md illumina/map.tsv \
 			   --out otu_lefse.txt \
 			   -c SampleType \
 			   -u Subject
 			   
# 3. 转换为lefse接受的格式
lefse-format_input.py otu_lefse.txt otu_lefse.in -c 1 -s -1 -u 2

# 4. 走正常的分析流程
run_lefse.py otu_lefse.in otu_lefse.res
lefse-plot_res.py otu_lefse.res otu_lefse.png
lefse-plot cladogram.py otu_lefse_res otu_lefse.cladogram.png --format png
mkdir biomarkers_raw_image
lefse-plot_features.py otu_lefse.in otu_lefse.res biomarkers_raw_images/
```
