---
title: 生信基础 常用gene_ID的转换
date: 2018-10-12
toc: true
categories: Bioinformatics
tags: 生信基础
---

## ID 类型

ID 示例 | ID 来源
---|---
ENSG00000116717	| Ensemble ID
GA45A_HUMAN	| UniProtKB/Swiss-Prot, entry name
A5PJB2_BOVIN	| UniProtKB/TrEMBL, entry name
A2BC19, P12345, | A0A022YWF9	UniProt, accession number
GLA, GLB, UGT1A1	| HGNC Gene Symbol
U12345, AF123456	| GenBank, NCBI, accession number
NT_123456, NM_123456, NP_123456	| RefSeq, NCBI, accession number
10598, 717	| Entrez ID, NCBI
uc001ett, uc031tla.1	| UCSC ID

<!--more-->

## 1. Ensembl stable ID
Ensembl stable ID 的结构是根据不同物种设置的前缀, 加上数据所指的类型, 如基因蛋白质, 再加上一系列的数字. 有的时候可以有不同的版本, 则在 Ensembl ID 后面加上小数点和版本号.

Ensembl_gene_identifier是Ensembl ID里的一种, Enseml ID包括exon, protein family, gene, gene tree, protein, regulatory feature 和 transcript.

Ensembl ID的由5部分构成: ENS(species)(object type)(identifier).(version)  
第一部分ENS代表这是一个Ensembl ID  
第二部分代表物种, 如MUS代表小鼠(如果物种是人则此处为空)  
第三部分代表ID的类型, 如G代表基因, T代表转录本, P代表蛋白, E代表外显子, S代表  
第四部分是一个特殊的数字标志  
第五部分代表版本号  

如:ENSMUSG00000017167.6
我们知道这是一个Ensembl ID (ENS), 物种为小鼠(MUS), 代表一个基因(G), 并且这是第6个版本(.6).

### 物种前缀

物种前缀	| 学名
---|---
ENSCEL	| Caenorhabditis elegans (Caenorhabditis elegans)
ENSCAF	| Canis lupus familiaris (Dog)
ENSDAR	| Danio rerio (Zebrafish)
FB	| Drosophila melanogaster (Fruitfly)
ENS	| Homo sapiens (Human)
ENSMUS	| Mus musculus (Mouse)
ENSRNO	| Rattus norvegicus (Rat)
ENSXET	| Xenopus tropicalis (Xenopus)

### 类型前缀
类型前缀	| 类型
---|---
E	| exon
FM	| Ensembl protein family
G	| gene
GT	| gene tree
P	| protein
R	| regulatory feature
T	| transcript

## 2. UniProt

UniProt中录入的数据都被分配了一个唯一的entry name.

- UniProtKB/Swiss-Prot entry name 是最多有 11 位包含大写字母的字符串, 一般有着 "X_Y" 的形式, 其中 "X" 是最多五个便于记忆的蛋白质编号, "_" 是下划线, "Y" 是最多五个便于记忆的物种编号.

- UniProtKB/TrEMBL entry name 是最多 16 位包含大写字母的字符串, 一般有着 "X_Y" 的形式, 其中 "X" 是 6 到 10 个字符组成的 accession number, "_" 是下划线, "Y" 是最多五个便于记忆的物种编号.

+ UniProtKB 的 Accession Number 相当于数据库的主键, 由 6 到 10 个大写字母或者数字组成. 其构成规律为: [OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}

### 蛋白质编号示例
Code(X)	| Recommended protein name	| Gene name
---|---|---
B2MG	| Beta-2-microglobulin	| B2M
HBA	| Hemoglobin subunit alpha	| HBA1
INS	| Insulin	| INS
CAD17	| Cadherin-17	| CDH17

### 物种编号
编号 | 物种
---|---
BOVIN	| Bovine
CHICK	| Chicken
ECOLI	| Escherichia coli
HORSE	| Horse
HUMAN	| Homo sapiens
MAIZE	| Maize (Zea mays)
MOUSE	| Mouse

## 3. Gene Symbol

Gene Symbol 是用来表示基因的编码, 由大写字母构成, 或由大写字母和数字构成, 首字母均应该是字母.如: GLA "galactosidase, alpha"; GLB "galactosidase, beta"; UGT1A1 "UDP glycosyltransferase 1 family, polypeptide A1" 再到 UGT1A13 代表了 13 个不同的 gene symbol.

## 4. NCBI
+ GenBank的通用accession number通常是由一个大写字母加上 5 个数字的组合, 或者两个大写字母加上6个数字的组合.
+ RefSeq有一套特殊的Accesion Number.形式是: [A-Z]{2}[_][0-9]{6:},两个大写字母,一个下划线,6个或更多的数字.

Accession | 前缀	| 类型	| 说明
---|---|---|---
AC_	| Genomic	| Complete genomic molecule, usually alternate assembly
NC_	| Genomic	| Complete genomic molecule, usually reference assembly
NG_	| Genomic	| Incomplete genomic region
NT_	| Genomic	| Contig or scaffold, clone-based or WGS
NW_	| Genomic	| Contig or scaffold, primarily WGS
NS_	| Genomic	| Environmental sequence
NZ_	| Genomic	| Unfinished WGS
NM_	| mRNA	| none
NR_	| RNA	| none
XM_	| mRNA	| Predicted model
XR_	| RNA	| Predicted model
AP_	| Protein	| Annotated on AC_ alternate assembly
NP_	| Protein	| Associated with an NM_ or NC_ accession
YP_	| Protein	| none
XP_	| Protein	| Predicted model, associated with an XM_ accession
ZP_	| Protein	| Predicted model, annotated on NZ_ genomic records

## 5. Entrez ID

Entrez 是 NCBI 使用的能够对众多数据库进行联合搜索的搜索引擎, 其对不同的 Gene 进行了编号, 每个 gene 的编号就是 entrez gene id. 由于 entrez id 相对稳定, 所以也被众多其他数据库, 如 KEGG 等采用. Entrez Gene ID 就是一系列数字, 也比较容易辨识. R 或网站都有众多的工具可以帮助从不同的 ID 转换为 entrez id 或者反向转换.

## 6.UCSC ID

UCSC ID由小写字母和数字构成, 起始均为 uc, 然后是三位数字, 接着又是三位小写字母, 最后有小数点和数字构成版本号.
如: uc010qfk.3, uc010qfk.3.

## 7.使用R包org.Hs.eg.db转换ID(org.xx.eg.db)

```r
## 类似数据库查询一样
library(org.Hs.eg.db)
# 显示可转换的ID类型
keytypes(org.Hs.eg.db)
ids = yourid_vector
# 你要转换的选项,由keytypes里面选
cols <- c("SYMBOL", "GENENAME", 'ENTREZID')
getIDs = select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="ENSEMBL")
```
## 8.使用biomaRt包转换

```r
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
listMarts()

ensembl<-useMart("ENSEMBL_MART_ENSEMBL")
all_datasets <- listDatasets(ensembl)

library(DT)
datatable(all_datasets,
        options = list(searching = FALSE,pageLength = 5,lengthMenu = c(5, 10, 15, 20)))
        
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
filters = listFilters(ensembl)
datatable(filters, options = list(searching = FALSE,pageLength = 5,lengthMenu = c(5, 10, 15, 20)))
attributes = listAttributes(ensembl)
datatable(attributes, options = list(searching = FALSE,pageLength = 5,lengthMenu = c(5, 10, 15, 20)))
entrez.mart<-getBM(attributes=c('ensembl_gene_id', 'entrezgene'), filters = 'ensembl_gene_id', values = rownames(fpkmtpens), mart = ensembl)
```

## 9. 其他数据库整合包

http://www.bio-info-trainee.com/1399.html
