---
title: 扩增子测序 QIIME2
date: 2018-10-15
toc: true
categories: Bioinformatics
tags:
  - 软件和包
  - 扩增子
---

## 1 qiime2中数据的导入与导出

在QIIME2中，所有的输入数据都是以qiime2 artifacts的格式(.qza)存在，该格式有利于数据传递和生成路径追踪。在QIIME2中，你可以在分析的各个阶段导入数据，不论是原始的测序数据，还是经过处理产生的中间数据(eg. biom)，你都可以直接导入，接着进行后续分析。导入数据使用qiime tools import命令，你可以使用--show-importable-types和--show-importable-formats查看支持导入的数据类型，选择与你数据相应的导入类型。

你可以使用qiime tools export导出qiime的数据，用于R或者其他软件。这个命令与qiime tools extract命令的差别在于，导出命令只导出数据，与数据相关的生成追踪信息将被丢弃；但是提取命令不会，里面包含的provenance文件保存了这些追踪信息。

<!--more-->

### 1.1 导入含有测序质量的信息序列数据

#### 1.1.1 'EMP protocol' multiplexed fastq

在EMP的格式结果中，由于它是混池测序的，有不同样品不同实验的测序结果，所以有2个fastq.gz文件，一个是测序序列文件，各个样品的测序结果混杂在一起(multiplexed)，另一个是对应的barcode文件.

这里的测序数据分别有single-end和pair-end，需要指定不同的导入类型

```bash
qiime tools import \
        --type EMPSingleEndSequneces \
		--input-path EMP_single \
		--output-path EMP_single_sequences.qza

qiime tools import \
		--type EMPPairEndSequences \
		--input-path EMP_pair \
		--output-path EMP_pair_sequences.qza
```

#### 1.1.2 'Fastq manifest' formats

在fastq manifest格式中，manifest文件指定了fastq.gz/fastq文件的绝对路径，比如TCGA中的manifest文件。你可以自己生成这类文件。它通常是逗号分隔的文件，每行的第一个域是样品标识符以供qiime2使用，第二个域是对应文件的绝对路径，第三个域是read的方向。注释行用#开头。文件的首行如果不是空行和注释行的话，则必须是表头行，内容为sample-id,absolute-filepath,direction。对于single-end reads，每个sample-id只能有1行信息。而对于pair-end reads，每个sample-id则必须有2行信息，分别对应着forward或reverse。在direction这一个域(列)的值只能是forward或者reverse。在absolute-filepath中，你可以使用环境变量$HOME/$PWD来指定路径。

在这个格式中，由于质量值不确定，我们需要选择自己数据质量值对应的输出格式。并且由于我们的manifest是根据sample-id来定义的，导入后的数据格式同EMP的数据demuplexed后的数据是一致的。

```bash
# pair-end reads
sample-id,absolute-filepath,direction
sample-1,$PWD/path/sample1_R1.fastq.gz,forward
sample-1,$PWD/path/sample1_R2.fastq.gz,reverse
# single-end reads
sample-id,absolute-filepath,direction
sample-1,$PWD/path/sample1_R1.fastq.gz,forward
```


```bash
# 注意，在旧版本中--input-format对应着--source-format
qiime tools import \
		--type 'SampleData[SequencesWithQuality]' \
		--input-format SingleEndFastqManifestPhred33 \
		--input-path ss-33-manifest \
		--output-path single-end-demux.qza

qiime tools import \
		--type 'SampleData[PairedEndSequnecesWithQuality]'
		--input-format PairedEndFastqManifestPhred64
		--input-path pe-64-manifest \
		--output-path paired-end-demux.qza
```

#### 1.1.3 Casava 1.8 demultiplexed fastq

在这种格式中，每个样品只有1个fastq.gz文件，文件名就包含了样品标识符。对于single样品来说，可能是类似于‘L2S357_15_L001_R1_001.fastq.gz’的形式，含义是sample-id_barcode_lane-number_read-number_set-number。而对于paired样品来说，就额外多了一个‘L2S357_15_L001_R2_001.fastq.gz’，其中R1，R2表明是2个reads。


```bash
# single-end
qiime tools import \
		--type 'SampleData[SequencesWithQuality]' \
		--input-path casava-1.8 \
		--input-format CasavaOneEightSingleLanePerSampleDirFmt \
		--output-path demux-single-end.qza

# paired-end
qiime tools import \
		--type 'SampleData[PairedEndSequencesWithQuality]' \
		--input-path casava-1.8 \
		--input-format CasavaOneEightSingleLanePerSampleDirFmt \
		--output-path demux-paired-end.qza
```

### 1.2 导入特征表数据

#### 1.2.1 biom

具体的格式描述见[biom v1.0.0](http://biom-format.org/documentation/format_versions/biom-1.0.html)和[biom v2.1.0](http://biom-format.org/documentation/format_versions/biom-2.1.html)。

```bash
# biom v1.0.0
qiime tools import \
		--input-path feature-table-v100.biom \
		--type 'FeatureTable[Frequency]' \
		--input-format BIOMV100Format \
		--output-path feature-table-v100.qza

# biom v2.1.0
qiime tools import \
		--input-path feature-table-v210.biom \
		--type 'FeatureTable[Frequency]' \
		--input-format BIOMV210Format \
		--output-path feature-table-v210.qza
```


### 1.3 导入特征序列数据

这一类数据通常是代表序列数据，是fasta格式的DNA序列。有些没有经过比对的；所以不会含有-的gap形式。有些则是经过比对的代表序列数据，其中就会含有-，使得各个序列的长度一致。但是这2类数据可能有N这类简并碱基。有些qiime2命令不支持这类。


```bash
# unaligned representative sequences
qiime tools import \
		--type 'FeatureData[Sequence]' \
		--input-path sequences.fna \
		--output-path sequences.qza

# aligned representative sequences
qiime tools import \
		--type 'FeatureData[AlignedSequence]' \
		--input-path aligned-sequences.fna \
		--output-path aligned-sequences.qza
```

### 1.4 导入无根系统发育树


```bash
qiime tools import \
		--type 'Phylogeny[Unrooted]'
		--input-path unrooted-tree.tre \
		--output-path unrooted-tree.qza
```

### 1.5 导出特征表

```bash
qiime tools export \
		--input-path feature-table.qza \
		--output-path exported-feature-table
```

### 1.6 导出发育树

```bash
qiime tools export \
		--input-path unrooted-tree.qza \
		--output-path exproted-tree
```


## 2 qiime2中的metadata

qiime的sample meatadata存储了关于实验的元信息，诸如各类技术细节，barcod、run、Subject、time point等。 Feature metadata通常是一个特征注释，比如能域序列变化或OTU对应的分类信息。下面是一些关于metadata的操作。

### 2.1 可视化metadata

```bash
qiime metadata tabulate \
		--m-input-file sample-metadata.tsv \
		--o-visualization tabulated-sample-metadata.qzv

qiime metadata tabulate \
		--m-input-file faith_pd_vector.qza \
		--o-visualization tabulated-faith-pd-metadata.qzv
```

### 2.2 合并两个metadata
你可以联合几个metadata进行可视化

```bash
qiime metadata tabulate \
  --m-input-file rep-seqs.qza \
  --m-input-file taxonomy.qza \
  --o-visualization tabulated-feature-metadata.qzv

qiime metadata tabulate \
		--m-input-file sample-metadata.tsv \
		--m-input-file faith_pd_vector.qza \
		--o-visualization tabulated-combined.qzv

qiime emperor plot \
		--i-pcoa unweighted_unifrac_pcoa_results.qza \
		--m-metadata-file sample-metadata.tsv \
		--m-metadata-file faith_pd_vector.qza \
		--o-visualization unweighted-unifrac-emperor-with-alpha.qzv
```

## 3 qiime2中数据的过滤

### 3.1 过滤特征表

特征表有两个维度：样品和特征，分别使用filter-samples和filter-features来实现。

+ 基于总体频率的过滤


```bash
# 比如我们要过滤frequeny低于1500的样品
qiime feature-table filter-samples \
		--i-table table.qza \
		--p-min-frequency 1500 \
		--o-filtered-table sample-frequency-filtered-table.qza

# 过滤在所有样品中frequency总计少于10的features
qiime feature-table filter-features \
		--i-table table.qza \
		--p-min-frequency 10 \
		--o-filtered-table feature-frequency-filtered-table.qza
```

+ 基于偶然性的过滤

在进行计算特征表时或者特征序列时，有些样品个数或者某个特征仅仅为1，这些是需要进行过滤，因为可能是由测序引起的假阳性。

```bash
qiime feature-table filter-features \
		--i-table table.qza \
		--p-min-samples 2 \
		--o-filtered-table sample-contigency-filtered-table.qza

qiime feature-table filter-samples \
		--i-table table.qza \
		--p-min-features 10 \
		--o-filtered-table feature-contigency-filtered-table.qza
```

+ 基于id的过滤

```bash
qiime feature-table filter-samples \
		--i-table table.qza \
		--m-metadata-file samples-to-keep.tsv \
		--o-filtered-table id-filtered-table.qza
```

+ 基于metadata的过滤

```bash
# 单个过滤
qiime feature-table filter-samples \
		--i-table table.qza \
		--m-metadata-file sample-metadata.tsv \
		--p-where "Subject='subject-1'" \
		--o-filtered-table subject-1-filtered-table.qza

# 多个过滤
qiime feature-table filter-samples \
		--i-table table.qza \
		--m-metadata-file sample-metadata.tsv \
		--p-where "BodySite IN ('left palm', 'right palm')" \
		--o-filtered-table subject-1-filtered-table.qza

# 联合过滤，支持SQL选择语句，AND，OR，AND NOT
qiime feature-table filter-samples \
		--i-table table.qza \
		--m-metadata-file sample-metadata.tsv \
		--p-where "Subject='subject-1' AND BodySite='gut'"
		--o-filtered-table subject-1-filtered-table.qza
```

+ 基于分类的过滤


```bash
qiime taxa filter-table \
		--i-table table.qza \
		--i-taxonomy taxonomy.qza \
		--p-exclude mitochondria \
		--o-filtered-table table-no-mitochondria.qza

qiime taxa filter-table \
		--i-table table.qza \
		--i-taxonomy taxonomy.qza \
		--p-include mitochondria, chloroplast \
		--o-filtered-table table-no-mitochondria.qza

# 保存有phylum-level注释的特征，从中保存的结果中去掉mitochondria,chloroplast
qiime taxa filter-table \
		--i-table table.qza \
		--i-taxonomy taxonomy.qza \
		--p-include p__ \
		--p-exclude mitochondria,chloroplast \
		--o-filtered-table table-with-phyla-no-mitochondria-no-chloroplast.qza

# 通过指定匹配模式为extact进行过滤
qiime taxa filter-table \
		--i-table table.qza \
		--i-taxonomy taxonomy.qza \
		--p-mode exact \
		--p-exclude "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__mitochondria" \
		--o-filtered-table table-no-mitochondria-exact.qza
```

### 3.2 过滤序列

在q2-feature-table中也有个filter-seqs，你可以用它进行多条件过滤。q2-quality-control里面的exclude-seqs也能进行序列过滤。

```bash
qiime taxa filter-seqs \
		--i-sequences sequences.qza \
		--i-taxonomy taxonomy.qza \
		--p-include p__ \
		--p-exclude mitochondria,chloroplast \
		--o-filtered-sequences sequences-with-phyla-no-mitochondria-no-chloroplast.qza
```

### 3.3 过滤距离矩阵


```bash
qiime diversity filter-distance-matrix \
		--i-distance-matrix distance-matrix.qza \
		--m-metadata-file samples-to-keep.tsv \
		--o-filtered-distance-matrix identifier-filtered-distance-matrix.qza

qiime diversity filter-distance-matrix \
		--i-distance-matrix distance-matrix.qza \
		--m-metadata-file sample-metadata.tsv \
		--p-where "Subject='subject-2'" \
		--o-filtered-distance-matrix subject-2-filtered-distance-matrix.qza
```

## 4 qiime2中特征分类器的训练

我们使用greengenes的参考序列对Native Bayes classifier进行训练，然后对目标代表序列进行分类。我们已经预先训练好了[多个分类器](https://docs.qiime2.org/2018.8/data-resources/)，你可以直接使用。**要想训练分类器的话，有2个必须的东西：参考序列和物种分类**。真实情况下，一般使用97%及以上的相似度的参考OTU数据集进行训练。参考序列和物种分类应该是要对应起来的额。比如你使用的是Greengenes 99% OTU sequences，那么你应该使用99% OTU taxonomy。

### 4.1 下载和导入数据


```bash
# 这里作为演示，用的是85 otu。真实情况不能使用。
wget https://data.qiime2.org/2018.8/tutorials/training-feature-classifiers/85_otus.fasta
wget https://data.qiime2.org/2018.8/tutorials/training-feature-classifiers/85_otu_taxonomy.txt
wget https://data.qiime2.org/2018.8/tutorials/training-feature-classifiers/rep-seqs.qza

qiime tools import \
		--type 'FeatureData[Sequence]' \
		--input-path 85_otus.fasta \
		--output-path 85_otus.qza

qiime tools import \
		--type 'FeatureData[Taxonomy]' \
		--input-path 85_otu_taxonomy.txt \
		--output-path 85_otu_taxonomy.qza
```

### 4.2 提取参考reads
在处理16s数据时，使用目标序列的特定区域训练朴素贝叶斯分类器会提高分类的准确度。在对ITS序列时，则不用这样，因为根据经验，在处理ITS的UNITE参考序列的预先训练时，进行提取或截断的效果并不是很好。所以对于ITS分析，使用全长更好。
在这里，我们先使用515F/806R primer把能匹配的reads提取出来，然后截断成120bp的长度。注意，只有当待分类的代表序列被截断成一致长度时，参考序列才需要截断成对应的长度(这里是120bp)。用于提取的primer应该是具有生物学意义的序列，不能是adapter、linker、barcode之类的，一般来讲，若是你的primer长于30 nt，肯定有非生物学序列在里面。

```bash
qiime feature-classifier extract-reads \
		--i-sequences 85_otus.qza \
		--p-f-primer GTGCCAGCMGCCGCGGTAA \
		--p-r-primer GGACTACHVGGGTWTCTAAT \
		--p-trunc-len 120 \
		--o-reads ref-seqs.qza
```

### 4.3 训练分类器

```bash
qiime feature-classifier fit-classifier-naive-bayes \
		--i-reference-reads ref-seqs.qza \
		--i-reference-taxonomy ref-taxonomy.qza \
		--o-classifier classifier.qza
```

### 4.4 测试分类器

```bash
qiime feature-classifier classify-sklearn \
		--i-classifier classifier.qza \
		--i-reads rep-seqs.qza \
		--o-classification taxonomy.qza

qiime metadata tabulate \
		--m-input-file taxonomy.qza \
		--o-visualization taxonomy.qzv
```

## 5 使用q2-quality-control进行质控和过滤

下面演示了q2-quality-control如何使用模拟群落/已知群落组成的已知样品来进行数据质控和过滤。

### 5.1 下载数据

```bash
wget https://data.qiime2.org/2018.8/tutorials/quality-control/query-seqs.qza
wget https://data.qiime2.org/2018.8/tutorials/quality-control/reference-seqs.qza
wget https://data.qiime2.org/2018.8/tutorials/quality-control/query-table.qza
wget https://data.qiime2.org/2018.8/tutorials/quality-control/qc-mock-3-expected.qza
wget https://data.qiime2.org/2018.8/tutorials/quality-control/qc-mock-3-observed.qza
```

### 5.2 通过比对过滤序列
exclude-seqs将查询序列和参考序列进行比对，根据多个比对标准区分能比对上和不能比对上的序列。你可以用本步骤来去除已知的污染序列、排除宿主序列、或这其他非目标序列。
exclude-seqs目前支持blast、vsearch、blastn-short作为比对方法，如果查询序列短于30nt，应当使用blastn-short。

```bash
qiime quality-control exclude-seqs \
		--i-query-sequences query-seqs.qza \
		--i-reference-sequences reference-seqs.qza \
		--p-method blast \
		--p-perc-identity 0.97 \
		--p-perc-query-aligned 0.97 \
		--o-sequence-hits hist.qza \
		--o-sequence-misses misses.qza

qiime feature-table filter-features \
		--i-table query-table.qza \
		--m-metadata-file hits.qza \
		--o-filtered-table no-hits-filtered-table.qza \
		--p-exclude-ids
```

### 5.3 用已知组成的样品进行质控
模拟群落的群落组成和丰度都是已知的，可以用于生信的基准测试(benchmarking bioinformatics methods)，比如测定某个方法或流程的处理准确度。
quality-control的evaluate-composition可以验证分类器分类结果的准确度。现有的模拟群落可在[mockrobiota查询](https://github.com/caporaso-lab/mockrobiota/blob/master/inventory.tsv).
在结果中，预期和观测的物种丰度的TAR(taxon accuracy rate)、TDR(taxon detection rate)和linear regression scores都会计算，并根据每个物种水平进行作图。

```bash
qiime quality-control evaluate-composition \
		--i-expected-features qc-mock-3-expected.qza \
		--i-observed-features qc-mock-3-observed.qza \
		--o-visualiation qc-mock-3-cmparison.qzv
```

### 5.4 序列质量的验证
evaluate-seqs将查询序列和参考序列进行比对验验证比对的质量。利用预测序列和观测序列可以验证desnoising或OTU picking结果的准确度。

```bash
qiime quality-control evaluate-seqs \
  --i-query-sequences query-seqs.qza \
  --i-reference-sequences reference-seqs.qza \
  --o-visualization eval-seqs-test.qzv
```


## 6 使用q2-sample-classifier预测metadata

正如多数统计方法一样，sample-classifier的预测功能也需要一定量的样本才能产生有意义的结果，一般来说，需要至少50个样品；样本太少产生的结果将会导致不准确的模型或者错误。用于预测的类别型元数据列的每个值都至少需要10个sample；用于回归的连续型元数据列不能有太多离群值或者分布不均匀。


```bash
wget https://data.qiime2.org/2018.8/tutorials/moving-pictures/sample_metadata.tsv
wget https://data.qiime2.org/2018.8/tutorials/sample-classifier/moving-pictures-table.qza
```

### 6.1 预测类别型样品数据

这里我们将根据样品的微生物组成来预测这个样品时来源于哪个身体部位。按照classify-samples进行,内部步骤如下：
1. 把样品随机分成训练集和测试集；通过--p-test-size来指定测试集的大小。
2. 选择一个估值器(estimator)，设置参数进行训练。
3. K-fold cross-validation进行模型优化，通过--p-cv设置K-fold值。
4. 使用测试集进行预测
5. 计算预测模型准确度

#### 6.1.1 建立预测模型

```bash
qiime sample-classfifier classify-samples \
		--i-table moving-pictures-table.qza \
		--m-metadata-file moving-pictures-sample-metadata.tsv \
		--m-metadata-column BodySite \
		--p-optimize-feature-selection \
		--p-parameter-tuning \
		--p-estimator RandomForestClassifier \
		--p-n-estimators 20 \
		--output-dir moving-pictures-classifier
# output:
# sample_estimator.qza, feature_importance.qza, predictions.qza, accuracy_result.qzv, model_summary.qzv

qiime metadata tabulate \
		--m-input-file moving-pictures-classifier/predictions.qza \
		--o-visualization moving-pictures-classifier/predictions.qzv

qiime metadata tabulate \
		--m-input-file moving-pictures-classifier/feature_importance.qza \
		--o-visualization moving-pictures-classifier/feature_importance.qzv

# 利用预测重要的feature过滤原始table
qiime feature-table filter-features \
		--i-table moving-pictures-table.qza \
		--m-metadata-file moving-pictures-classifier/feature_importance.qza \
		--o-filtered-table moving-pictures-classifier/important-feature-table.qza
```

#### 6.1.2 进行预测


```bash
# 预测新的metadata
qiime sample-classifier predict-classification \
  		--i-table new-table.qza \
  		--i-sample-estimator moving-pictures-classifier/sample_estimator.qza \
  		--o-predictions moving-pictures-classifier/new_predictions.qza


# 用新的样品重新测试模型的准确性
qiime sample-classifier confusion-matrix \
		--i-predictions moving-pictures-classifier/new_predictions.qza \
		--m-truth-file moving-pictures-sample-metadata.tsv \
		--m-truth-column BodySite \
		--o-visualization moving-pictures-classifier/new_confusion_matrix.qzv
```


### 6.2 预测连续型样品数据

```bash
wget https://data.qiime2.org/2018.8/tutorials/longitudinal/sample_metadata.tsv
wget https://data.qiime2.org/2018.8/tutorials/longitudinal/ecam_table_maturity.qza
```

#### 6.2.1 训练模型

```bash
qiime sample-classifier regress-samples \
		--i-table ecam-table.qza \
		--m-metadata-file ecam-metadata.tsv \
		--m-metadata-column month \
		--p-estimator RandomForestRegressor \
		--p-n-estimators 20 \
		--output-dir ecam-regressor

# ouput
# sample_estimator.qza, feature_importance.qza, accuracy_results.qzv, model_summary.qzv
```

### 6.3 使用nested cross-validation法进行预测

#### 6.3.1 类别型变量：训练模型

```bash
qiime sample-classifier classify-samples-ncv \
		--i-table moving-pictures-table.qza \
		--m-metadata-file moving-pictures-sample-metadata.tsv \
		--m-metadata-column BodySite \
		--p-estimator RandomForestClassifier \
		--p-n-estimators 20 \
		--o-predictions BodySite-predictions-ncv.qza \
		--o-feature-importance BodySite-importance-ncv.qza

# 查看模型准确性
qiime sample-classifier confusion-matrix \
		--i-predictions BodySite-predictions-ncv.qza \
		--m-truth-file moving-pictures-sample-metadata.tsv \
		--m-truth-column BodySite \
		--o-visualization ncv_confusion_matrix.qzv
```

#### 6.3.2 连续型变量：训练模型

```bash
qiime sample-classifier regress-samples-ncv \
		--i-table ecam-table.qzv \
		--m-metadata ecam-metadata.tsv \
		--m-metadata-column month \
		--p-n-estimators 20 \
		--p-estimator RandomForestRegressor \
		--o-predictions ecam-predictions-ncv.qza \
		--o-feature-importance ecam-importance-ncv.qza

# 查看结果准确性
qiime sample-classifier scatterplot \
		--i-predictions ecam-predictions-ncv.qza \
		--m-truth-file ecam-metadata.tsv \
		--m-truth-column month \
		--o-visualization ecam-scatter.qzv
```


### 6.4.最佳实践：你不应该做q2-sample-classifier做的事
+ 数据泄漏
+ 过拟合

## 7 使用q2-longitudinal进行纵向和成对差异比较

q2-longitudinal可以用来进行纵向研究设计和成对样品的统计和可视化，以探究样品是如何在可观测的状态之间变化的。可观测的状态通常是时间或者环境梯度，也可以是不同时间点的成对分析比如成对的距离和差异分析；比如在治疗前后的临床研究。可观测的状态也可以是方法上的，在同个时间点上采用不同方法的结果的比较。比如q2-longitudinal可以比较不同样品收集方法、存储方法、DNA提取方法对单个样品的微生物群落特征组成影响。

在q2-longitudinal中的操作多数以某个度量值作为输入，这通常是metadata里面的某个列、alpha diversity vectors、PCoA results或者是feature ID。合法的度量值名字可以通过qiime metadata tabulate命令查看得到，feature name也可以通过qiime feature-data summarize命令查看得到。

### 7.1 成对差异比较

成对差异比较探究的是成对样品之间某个特定的度量值是否具有显著性变化；目前支持的比较类型有：
+ 特征丰度的比较(comparison of feature abundance): 比如在特征表中微生物序列变体或者变化率
+ 元数据值的比较(comparison of metadata values): 比如将alpha/beta diversity values作为度量值，结合某个元数据值进行比较

下面的命令是探究以delivery mode作为分组，不同的两个时间点上(state-column)两组的shannon diversity index(alpha diversity)是否发现显著变化。

```bash
qiime longitudinal pairwise-differences \
		--m-metadata-file ecam-sample-metadata.tsv \
		--m-metadata-file shannon.qza \
		--p-metric shannon \
		--p-group-column delivery \
		--p-state-column month \
		--p-state-1 0 \
		--p-state-2 12 \
		--p-individual-id-column studyid \
		--p-replicate-handling random \
		--o-visualization pairwise-differences.qzv
```

### 7.2 成对距离比较
成对距离比较探究的也是成对样品之间某个特定的度量值是否具有显著性变化，但是，相对于使用metadata的某个列或者qiime artifact作为输入，成对距离比较以距离矩阵作为输入，来探究成对样品在pre和post下的距离变化，并显示其差异是否显著。
下面的命令是探究以delivery作为分组，正常分娩和剖腹产在12个月后对微生物组成的稳定性的影响。

```bash
qiime longitudinal pairwise-distances \
		--i-distance-matrix unweighted_unifrac_distance_matrix.qza \
		--m-metadata-file ecam-sample-metadata.tsv \
		--p-group-column delivery \
		--p-state-column month \
		--p-state-1 0 \
		--p-state-2 12 \
		--p-individual-id-column studyid \
		--p-replicate-handling random \
		--o-visualization pairwise-distances.qzv
```

### 7.3 线性混合效应模型

线性混合效应模型(linear mixed effects models, LME)探究单个因变量和一至多个自变量之间的关系；因变量的观测值通常在跨相关样本中获得，比如在重复测量的抽样实验中。该模型以至少一个数字型state-column(如时间)或者一至多个逗号分隔的group-columns作为自变量，然后绘出因变量(metric,某个度量)的回归曲线。此外，individual-id-column参数应该是个可以指明样品重复情况的metadata column。因变量可以是metadata的某列(column)，也可以是feature table的feature ID。逗号分隔的多个随机效应值可以作为该模型的输入；模型默认每个个体都有随机的截距；如果用户为每个个体提供随机的斜率(slope)的话，可以使用state-column值来作为random-effects的值。

注意，判断某个因素是为固定效应还是随机效应比较复杂。一般来说，如果一个因子的不同水平(metadata column values)可以代表所有的离散值的话，该因子应该是固定效应；例如delivery mode、sex、diet在下面的分析中就是固定效应。那相反地，如果一个因子的值是代表群体内的随机样品，那么这个因子就是随机效应；比如body-weight、daily-kcal-from-breastmilk等，这些值只能代表从人群中的某些人，并不能代表所有的情况。

这里我们使用LME来探究shannon diversity index是否会随时间、delivery mode、diet、sex发生变化。

```bash
qiime longitudinal linear-mixed-effects \
		--m-metadata-file ecam-sample-metadata.tsv \
		--m-metadata-file shannon.qza \
		--p-metric shannon \
		--p-group-columns delivery,diet,sex \
		--p-state-column month \
		--p-individual-id-column studyid \
		--o-visualization linear-mixed-effects.qzv
```

产生的可视化文件包含多个图。首先是输入参数显示在最上方。然后是model summary显示LME模型的训练信息。主要的结果是第三个，所谓model results table。它显示了每个固定效应和效应之间的互作对因变量(shannon diversity)的影响，它还显示了其他信息如参数估计、估计标准差、z-score、p-values等。从该表可以看出shannon diversity受month、diet和其他几个交互因子影响明显。最后是依据各个group column制作的散点图和线性回归曲线。如果--p-lowess指定了，还会显示加权平均值。第一类散点图只是进行了一下快速统计，可以使用volatility进行交互式地统计和作图。第二类散点图显示每个样品的观测残差值和度量预测值之间的关系：度量预测值应该在0值附近，与残差值没有相关性。如果观察到U形或者非随机分布图，说明你的预测变量(group_columns, random_effects)不能很好地表征数据。


### 7.4 波动分析(volatility anaylsis)
波动分析(volatility analysis)可产生交互式的折线图以探究因变量在单组或多组之间的连续型自变量的依赖关系。多个metadata file(alpha,beta diversity)和FeatureTable[RelativeFrequency]可以作为输入，也可以选择不同的因变量作为图的y轴。
这里我们探究一下shannon diversity和其他metadata是如何随时间变化的。

```bash
qiime longitudinal volatility \
		--m-metadata-file ecam-sample-metadata.tsv \
		--m-metadata-file shannon.qza \
		--p-default-metric shannon \
		--p-default-group-column delivery \
		--p-state-column month \
		--p-individual-id-column studyid \
		--o-visualization volatility.qzv
```

### 7.5 用一阶差分来追踪变化率
观察时间序列数据的另一个方法就是计算变化率随时间的变化，要实现这个目的就需要计算一阶差分，它表征了在连续时间点之间的变化大小。delta Yt = Y(t) -Y(t-1),这个变换在first-differences/first-distances中用到。对于生成的结果文件，你可以进行波动分析。

```bash
# alpha diversity or metadata column (metric)
qiime longitudinal first-differences \
  --m-metadata-file ecam-sample-metadata.tsv \
  --m-metadata-file shannon.qza \
  --p-state-column month \
  --p-metric shannon \
  --p-individual-id-column studyid \
  --p-replicate-handling random \
  --o-first-differences shannon-first-differences.qza
# beta diversity (matrix)
qiime longitudinal first-distances \
  --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ecam-sample-metadata.tsv \
  --p-state-column month \
  --p-individual-id-column studyid \
  --p-replicate-handling random \
  --o-first-distances first-distances.qza
# 经过计算，beta diversity也可以进行LME分析
qiime longitudinal linear-mixed-effects \
  --m-metadata-file first-distances.qza \
  --m-metadata-file ecam-sample-metadata.tsv \
  --p-metric Distance \
  --p-state-column month \
  --p-individual-id-column studyid \
  --p-group-columns delivery,diet \
  --o-visualization first-distances-LME.qzv
```

你也可以从固定时间点开始追踪变化率的大小变化。


```bash
qiime longitudinal first-distances \
  --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ecam-sample-metadata.tsv \
  --p-state-column month \
  --p-individual-id-column studyid \
  --p-replicate-handling random \
  --p-baseline 5 \
  --o-first-distances first-distances-baseline-0.qza
```



### 7.6 微生物相互依赖性非参数检验

在微生物群落中，微生物种群并不会独立存在，而是以生态互作网存在。而这些互作网是否会在同组个体之家产生临时效应将指明不同的研究方向。微生物相互依赖性非参数检验(non-parametric microbial interdependence test, NMIT)将评估同组的同个群落内的特征相互依赖性(比如种类、序列变体、OTU)如何依据时间不同而不同。MNIT评估纵向样本相似性。对于每个个体，NMIT计算成对特征的相关性；然后基于NMIT矩阵计算成对个体的距离。该方法只对纵向数据有用，例如时间序列数据。为使检验结果鲁棒，我们建议每个个体最少要测5个时间点的数据；但它不要求所有的样品数据都在同个时间点测得。

NMIT的检验结果是一个距离矩阵，你可以像处理beta diversity的结果数据一样进行处理。


```bash
wget https://data.qiime2.org/2018.8/tutorials/longitudinal/ecam_table_taxa.qza

qiime longitudinal nmit \
  --i-table ecam-table-taxa.qza \
  --m-metadata-file ecam-sample-metadata.tsv \
  --p-individual-id-column studyid \
  --p-corr-method pearson \
  --o-distance-matrix nmit-dm.qza

qiime diversity beta-group-significance \
  --i-distance-matrix nmit-dm.qza \
  --m-metadata-file ecam-sample-metadata.tsv \
  --m-metadata-column delivery \
  --o-visualization nmit.qzv

qiime diversity pcoa \
  --i-distance-matrix nmit-dm.qza \
  --o-pcoa nmit-pc.qza

qiime emperor plot \
  --i-pcoa nmit-pc.qza \
  --m-metadata-file ecam-sample-metadata.tsv \
  --o-visualization nmit-emperor.qzv
```


### 7.7 特征波动分析

特征波动分析是一个监回归方法，可以确定某个数字型metadata column、state-column的值，并绘出不同state的相对频率图。state-column通常是个时间度量值，不过它也可以时任何数字型column；但是它不能和individual-id-column相同。


```bash
wget -O ecam-table.qza https://data.qiime2.org/2018.8/tutorials/longitudinal/ecam_table_maturity.qza

qiime longitudinal feature-volatility \
  --i-table ecam-table.qza \
  --m-metadata-file ecam-sample-metadata.tsv \
  --p-state-column month \
  --p-individual-id-column studyid \
  --p-n-estimators 10 \
  --p-random-state 123 \
  --output-dir ecam-feat-volatility
```

### 7.8 成熟指数预测

这个分析在对各个时间点样品数目分布均匀时的分析结果较好; 对于某些时间点有缺失值的情况下无法进行分析；此外，本分析还需要大量数据才能正常分析，特别是对照组在各个时间点需要足够的生物学重复。这个分析基于特征数据(feature data)构建回归模型，计算微生物成熟指数，一次来预测给定连续型metadata column的值，例如给定微生物组成预测个体的年龄。这个分析域标准监督回归分析不同在于它量化了多个组之间随时间变化的相对速率。

下面我们比较在不同年龄下，剖腹产和自然分娩之间微生物群落的成熟度。

```bash
qiime longitudinal maturity-index \
  --i-table ecam-table.qza \
  --m-metadata-file ecam-sample-metadata.tsv \
  --p-state-column month \
  --p-group-by delivery \
  --p-individual-id-column studyid \
  --p-control Vaginal \
  --p-test-size 0.4 \
  --output-dir maturity
```
## 8 qiime2中嵌合体的去除

qiime2中的嵌合体检测基于FeatureTable[Frequency]和FeatureData[Sequences]，它整合了Uchime de novo和vsearch这两个软件。


```bash
wget https://data.qiime2.org/2018.8/tutorials/chimera/atacama-table.qza
wget https://data.qiime2.org/2018.8/tutorials/chimera/atacama-rep-seqs.qza
```

### 8.1 de novo chimera checking

```bash
qiime vsearch uchime-denovo \
		--i-table atacama-table.qza \
		--i-sequences atacama-rep-seqs.qza \
		--output-dir uchime-dn-out
# out
# nonchimeras.qza, chimeras.qza, stats.qza

# 可视化统计结果
qiime metadata tabulate \
		--m-input-file uchime-dn-out/stats.qza \
		--o-visualization uchime-dn-out/stats.qzv
```

### 8.2 对输入的tables和序列进行过滤

+ 去除chimeras和borderline chimeras


```bash
qiime feature-table filter-features \
  --i-table atacama-table.qza \
  --m-metadata-file uchime-dn-out/nonchimeras.qza \
  --o-filtered-table uchime-dn-out/table-nonchimeric-wo-borderline.qza

qiime feature-table filter-seqs \
  --i-data atacama-rep-seqs.qza \
  --m-metadata-file uchime-dn-out/nonchimeras.qza \
  --o-filtered-data uchime-dn-out/rep-seqs-nonchimeric-wo-borderline.qza

qiime feature-table summarize \
  --i-table uchime-dn-out/table-nonchimeric-wo-borderline.qza \
  --o-visualization uchime-dn-out/table-nonchimeric-wo-borderline.qzv
```

+ 去除嵌合体，但是保留 borderline chimeras


```bash
qiime feature-table filter-features \
  --i-table atacama-table.qza \
  --m-metadata-file uchime-dn-out/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-table uchime-dn-out/table-nonchimeric-w-borderline.qza

qiime feature-table filter-seqs \
  --i-data atacama-rep-seqs.qza \
  --m-metadata-file uchime-dn-out/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-data uchime-dn-out/rep-seqs-nonchimeric-w-borderline.qza
  
qiime feature-table summarize \
  --i-table uchime-dn-out/table-nonchimeric-w-borderline.qza \
  --o-visualization uchime-dn-out/table-nonchimeric-w-borderline.qzv
```

## 9 qiime2中的reads合并与去噪

这里描述的方法并不是说要代替DADA2中的合并序列、去噪过程。相反地，这里的方法的关注点在于分析 paired-end reads的其他方法。如果你对使用DADA2进行合并和去噪感兴趣的话，[atacama soil](https://docs.qiime2.org/2018.8/tutorials/atacama-soils/)写明了如何使用qiime dada2 denoise-paired实现；如果使用dada2的话，请不要预先进行合并序列。


```bash
wget -O demux.qza https://data.qiime2.org/2018.8/tutorials/read-joining/atacama-seqs.qza
```

### 9.1 joining reads

```bash
qiime vsearch join-pairs \
		--i-demultiplexed-seqs demux.qza \
		--o-joined-sequences demux-joined.qza

qiime demux summarize \
		--i-data demux-joined.qza \
		--o-visualization demux-joined.qzv
```

### 9.2 sequence quality control

```bash
qiime quality-filter q-score-joined \
		--i-demux demux-joined.qza \
		--o-filtered-sequences demux-joined-filtered.qza \
		--o-filter-stats demux-joined-filter-stats.qza
```

### 9.3 deblur

在第二步之后，你可以选择使用deblur进一步质控，也可以进行去重复序列然后进行聚类。

```bash
qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-joined-filtered.qza \
  --p-trim-length 250 \
  --p-sample-stats \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-stats deblur-stats.qza

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv
```

### 9.4 导入预合并的reads

```bash
wget https://data.qiime2.org/2018.8/tutorials/read-joining/fj-joined.zip
unzip fj-joined.zip

qiime tools import \
  --input-path fj-joined/manifest \
  --output-path fj-joined-demux.qza \
  --type SampleData[JoinedSequencesWithQuality] \
  --input-format SingleEndFastqManifestPhred33

qiime demux summarize \
  --i-data fj-joined-demux.qza \
  --o-visualization fj-joined-demux.qzv
```

## 10 qiime2中的OTU聚类

在qiime2中，de novo, closed-reference, open-reference clustering均已支持。在拆分序列、质控、去重复后，可以使用vsearch进行序列/特征聚类成OTUs。


```bash
wget https://data.qiime2.org/2018.8/tutorials/otu-clustering/seqs.fna
wget https://data.qiime2.org/2018.8/tutorials/otu-clustering/85_otus.qza

qiime tools import \
		--input-path seqs.fna \
		--output-path seqs.qza \
		--type 'SampleData[Sequences]'
```

### 10.1 对序列进行去重复


```bash
qiime vsearch dereplicate-sequences \
		--i-sequences seqs.qza \
		--o-dereplicated-table table.qza \
		--o-dereplicated-sequences rep-seqs.qza
```

### 10.2 clustering

+ de novo clustering


```bash
qiime vsearch cluster-features-de-novo \
		--i-table table.qza \
		--i-sequences rep-seqs.qza \
		--p-perc-identity 0.99 \
		--o-clustered-table table-dn-99.qza \
		--o-clustered-sequences rep-seqs-dn-99.qza
```

+ closed-reference clustering


```bash
qiime vsearch cluster-features-closed-reference \
		--i-table tablq.qza \
		--i-sequences rep-seqs.qza \
		--i-reference-sequences 85_otus.qza \
		--p-perc-identity 0.85 \
		--o-clustered-table table-cr-85.qza \
		--o-clustered-sequences  rep-seqs-cr-85.qza \
		--o-unmatched-sequences unmatched-cr-85.qza
```

+ open-reference clustering


```bash
qiime vsearch cluster-features-open-reference \
		--i-table tablq.qza \
		--i-sequences rep-seqs.qza \
		--i-reference-sequences 85_otus.qza \
		--p-perc-identity 0.85 \
		--o-clustered-table table-or-85.qza \
		--o-clustered-sequences  rep-seqs-or-85.qza \
		--o-new-reference-sequences new-ref-seqs-or-85.qza
```
