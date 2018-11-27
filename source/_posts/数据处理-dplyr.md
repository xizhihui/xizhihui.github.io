---
title: 数据处理 dplyr
date: 2018-11-02
toc: true
categories: Bioinformatics
tags:
  - 软件和包
  - 数据处理
---

以前就看过一些 dplyr 的使用，记得当初说是和 tplyr，reshape2 并称 R 数据处理三剑客，想想和 web 开发的 HTML，JavaScript，CSS 的三剑客有些类似。这两天又在《R 数据科学》看到 dplyr 的使用，还挺详细的。现在把书里的相关用法记录在此，刚好也把书上的习题在这里做个回答，以加强使用。

<!--more-->

## 准备工作

dplyr 的这些函数的第一个参数是一个数据框，随后的参数是对数框的操作，使用列名（不带引号）引用列，输出一个新的数据框，不会更改源数据。这里使用的数据是 nycflights13 包里的 flights，不是 tidyverse 里的，需要我们通过 *install.packages* 安装。具体的列名代表了什么意思，可以通过 *?flights* 查看。


```r
suppressPackageStartupMessages(library(nycflights13))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
flights[1:4,]
```

## 使用 filter 筛选行

我们使用 filter 对数据框的行进行筛选，对于数据框列的操作有以下几种。

*   比较运算：>, >=, <=, <, !=, ==（比较浮点数使用 near 函数）
*   逻辑运算：&, |, !, xor(x,y), isTRUE(x), isFALSE(x), %in%
*   缺失值：is.na(x)

问题是我们要找到“在夏季 7、8、9 月出发的航班”和“到达时间延误 2 小时，但出发时间没有延误的航班”。


```r
# 问题 1，有多种解法, 看下方
filter(flights, month == 7 | month == 8 | month == 9)
filter(flights, month %in% c(7,8,9))
filter(flights, between(momtn, 7, 9))
# 问题 2,
filter(flights, arr_delay == 120 & dep_delay == 0)
```

## 使用 arrange 排列行

arrange 的作用是改变行的顺序，如果有多个参数，那么就按输入参数的顺序对数据框依次进行排序，注意排序的效果是嵌套的。如果要使用降序，可以使用 *des()* 函数。一般来说，NA 值总是排在最后的。

这里的问题是“如何将缺失值排序在前面”和“在所有航班中，哪个航班的飞行时间最长，哪个最短”。


```r
# 问题 1，比如把到达时间的 NA 排在前面
arrange(flights, desc(is.na(flights$arr_time)))
# 问题 2,
arrange(flights, air_time)
arrange(flights, desc(air_time))
```

## 使用 select 选择列

select 选择列的时候，可以直接指定待选列的列名，也可以使用切片语法（闭区间），和排除语法“-”。此外，还有以下一些辅助函数。对于 everything，通常在对待选列重排到列首有奇效。在指定一个列名多次时，只会选取一次该列，一是因为列名不允许重复，二是我们从 everything 提取目标列到前面也能猜测到。

*   starts_with("abc"), ends_with("abc")
*   contains("abc")：大小写不敏感的
*   matches("(.)\\1")
*   num_range("x", 1:3): 获得 x1，x2，x3，两个参数不能同时为多个
*   one_of(c("abc", "def"))
*   everything(): 全选

在 dplyr 中有个与 select 相近的函数 rename，他可以对列进行重命名。

这里的问题是“使用尽可能多的方式来选择 dep_time, dep_delay, arr_time, arr_delay ”。


```r
select(flights, starts_with("dep"), starts_with("arr"))
select(flights, matches("^(dep|arr)"))
select(flights, num_range("dep", c("_time", "_delay")), num_range("arr", c("_time", "_delay")))
select(flights, num_range(c("dep", "arr"), "_delay"), num_range(c("dep", "arr"), "_time"))
select(flights, one_of(c("dep_time", "dep_delay", "arr_time", "arr_delay")))
```

## 使用 mutate 添加新列

mutate 可以对数据框添加新列，但总是放在最后面。新列在添加时就可以使用。如果你只是想保留新列，则可以使用 transmute 代替。与 select 类似，有多种辅助函数来帮助创建新列，这些函数的一个共同点就是输入和输出都是一个向量数据结构，一维的，长度与数据框的行数相同。常见的辅助函数有下面这些。

*   运算符：+, -, \*, /, ^, %/%, %%
*   运算函数: log, commin, cummax, [sum, min, max, lead, lag], 对于摘要性质的，需要与其他函数联用
*   逻辑比较：<, <=, >, >=, !=, ==
*   排秩：  min_rank, row_number, dense_rank, percent_rank, cume_dist, ntile

这里的问题是“比较 air_time 和 arr_time - dep_time，你期望看到什么？实际结果是什么？”


```r
transmute(flights, air_time = air_time, arr_dep = arr_time - dep_time)
# 因为 arr_tim 和 dep_time 是以“ HHMM ”的格式录入的，所以应该如下
HMM2min = function(HMM) {
    H = HMM %/% 100
    min = HMM %% 100 + H * 60
    min
}
transmute(flights, air_time = air_time, arr_dep = HMM2min(arr_time) - HMM2min(dep_time))
```


## 使用 group_by 分组

group_by 可以根据使用的列名个数对数据框进行逐级分组，得到分组的统计总量。每进行一次摘要统计，就会使用掉一个分组。这里值得注意的是，进行摘要统计的函数必须是分组层级一致的，比如分组求均值和总体求均值结果一致，但是分组求中位数却不是这样。ungroup 函数可以取消分组，回到未分组的数据。单独使用 group_by 函数的情况较少，通常是与 summrise 联用。


```r
group_by(flights, year, month, day)
```

## 使用 summarise 摘要数据

summarise 函数可以对数据框进行摘要统计，把数据框整合成一行，它通常是同groupby联合使用，在 dplyr 中最常用的操作就是进行分组摘要。%>% 的引入让分组摘要的流程更加明确、清晰、自然。在进行分组摘要的时候，由于 NA 具有扩撒性，na.rm 要记得指定为真，n() 与 sum(!is_na()) 总是会不一样。常用的摘要函数如下面所示。

*   位置度量：mean, median
*   分散程度度量：sd, IQR, mad
*   秩度量：min, max, quantile(x, 0.25)
*   定位度量：first, nth(x, 2), last
*   数量度量：n, sum(!is.na(x)), distinct, count

这里的问题是“每天取消航班的比例与平均延误时间有什么关系”。


```r
# 定义取消航班为 is.na(dep_time) & is.na(dep_delay)
filter(flights, is.na(dep_time) & is.na(dep_delay)) %>%
    mutate(flight = flight, mean_time = HMM2min(sched_arr_time) - HMM2min(sched_dep_time)) %>%
    ggplot(mapping = aes(x=flight, y=abs(mean_time))) +
        geom_point() +
        geom_smooth(se=F)
```

## 参考

*   《R 数据科学》
*   vignette("dplyr")
