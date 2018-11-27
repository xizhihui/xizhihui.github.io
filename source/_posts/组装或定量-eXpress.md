---
title: 组装或定量 eXpress
date: 2018-10-12
toc: true
categories: Bioinformatics
tags:
  - 软件和包
  - 定量
---

<!--more-->

> https://pachterlab.github.io/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz

## 使用方法


```bash
express [options]* <target_seqs.fasta> <aligned_reads.(sam/bam)>

# <target_seqs.fasta>:	参考序列的fasta文件
# <aligned_reads.(sam/bam)>:	比对结果文件
# -o path: 	指定输出路径
# -B n:	指定进行batch EM rounds计算的次数,以时间为代价提高定量准确性
# -O n:	指定进行online EM rounds计算的次数,以时间为代价提高定量准确性
# -m n:	指定fragment的平均长度
# -s n：指定fragment长度的标准差
# -H str:	以逗号分隔的位置坐标,指定haplotypes的位置。利于等位基因的特异表达定量
# --output-align-samp:	
# --output-align-prob:
# --fr-stranded: 指定链特异性比对的方向,read1比对到forward链,read2比对到reverse链
# --rf-stranded：指定链特异性比对的方向,read1比对到reverse链,read2比对到forward链
# --f-stranded:	 single-end, read比对到forward链
# --r-stranded:	 single-end, read比对到reverse链
# 
```

## 输入文件

eXpress需要一个multi-fasta格式的输入文件,计算该文件序列的转录丰度。
