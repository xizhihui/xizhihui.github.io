---
title: 生信基础 reads 重复和比对重复
date: 2018-10-12
toc: true
description: 呦呦呦
categories: Bioinformatics
tags: 生信基础
---

在说duplication之前,我们有必要说一下PCR bias和unique reads. PCR bias是指在PCR扩增的过程中,PCR引物会偏好性地同某条DNA链结合,导致的结果就是这条DNA链被扩增的数目要更多;而不是所有的DNA链被平行扩增.那么什么又是unique reads呢?在得到测序结果后,对于任意两条reads,只要其reads的起点,中间序列,终点这三点有一点不同,这两条reads就是互为unique reads. 所以我们来看, 建库后PCR扩增的测序结果里面,肯定有很多条reads不是unique reads,这些reads就是所谓的duplication, 这其中肯定也包含了PCR bias引起的额外重复.

<!--more-->

> duplication rate = 1 - unique reads / total reads.

## 1. duplication的来源

+ PCR duplicate
+ cluster duplicate
+ optical duplicate: 单个cluster的光信号重影导致2个信号产生

## 2. 影响重复率(duplication rate)的因素

+ 扩增模版的多样性
+ 模版序列的多样性: PCR bias
+ adapter与fragment的连接效率
+ fragment的一致性: 片段长度不同,扩增效率也是不同的
+ cluster PCR: cluster的大小影响测序信号的强弱
+ 其他影响测序过程的方面

### 3. 去重
如果是原始数据的重复在质控阶段就可以了,类似的包Trimmomatic和fastp都能去重.如果是比对后的重复,使用samtools可以去重,也可以使用picard去重.(去重的原理是把bam文件里面的比对结果改成1024,表示duplicate)
