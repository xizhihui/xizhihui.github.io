---
title: 可视化 ggplot2入门
date: 2018-11-02
toc: true
description: 听说你要画图？试试 ggplot2 吧，可以搭积木的。
categories: Bioinformatics
tags:
  - 软件和包
  - 可视化
---




久闻 ggplot2 大名，它的出图也在各个生信分析包中随处可见。今天恰好买来的新书《R 数据科学》里面有讲解，遂随着它一起学习一下，也将以前的相关学习一齐记录于此。

<!--more-->

这里我们使用的是 ggplot2 自带的 mpg 和 diamonds 数据集，还有 nlme 的Oxboys 数据集。

```r
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
```

```r
mpg[1:5, ]
diamonds[1:5, ]
data(Oxboys, package = "nlme")
Oxboys[1:5, ]
```

## ggplot2 初探
### 用 ggplot2 画个图

```r
# step 1 生成 ggplot 画板，指定数据对象
myplot = ggplot(data=mpg)

# step 2 添加几何对象 (geom_point)，添加图形映射 (aes)
ggplot(data=mpg) + 
    geom_point(mapping=aes(x=displ, y=hwy))

# step 3 添加其他图形映射, 如 color, size, alpha, shape
# 颜色属性 color 以 mpg$class 进行映射
# 大小属性 size 以 mpg$class 进行映射
ggplot(data=mpg) + 
    geom_point(mapping = aes(x=displ, y=hwy, color=class, size=class))

# step 4 尝试改变几何对象的图形属性, 注意与图形映射属性的差别
ggplot(data=mpg) + 
    geom_point(mapping = aes(x=displ, y=hwy), color="blue")

# step 5 尝试使用连续变量来映射图形属性, 
# 这时候分类变量就会变成连续的图形属性，比如 color 是颜色过渡而非颜色分类; 
# 对于本就是连续变量的图形属性来讲，就没啥影响，比如 size
# 这里的 dispal < 5 生成了一个分类变量
ggplot(data = mpg) +
    geom_point(mapping = aes(x=displ, y=hwy, size=cty, color=displ < 5))

# step 6 分面的引入
# 前面通过添加图形属性映射来添加新的变量；
# 除此之外，也可以通过分面实现，但是要注意，二者的表现形式不同
# 分面会对变量进行 level，根据 level 设置对应的分面；即便使用的是连续变量
# 单变量分面, 我们得到的是网格形式排列（行列）的分面，ncol 和 nrow 指定行列数
ggplot(data=mpg) +
    geom_point(mapping=aes(x=displ, y=hwy)) +
    facet_wrap(~ class)
# 双变量分面
ggplot(data=mpg) +
    geom_point(mapping = aes(x=displ, y=hwy)) +
    facet_grid(drv ~ cyl)
ggplot(data=mpg) +
    geom_point(mapping = aes(x=displ, y=hwy)) +
    facet_grid(cyl ~ drv)
# 单变量分面，不做网格形式，仅做行排列或列排列
ggplot(data=mpg) +
    geom_point(mapping = aes(x=displ, y=hwy)) +
    facet_grid(class ~ .)
ggplot(data=mpg) +
    geom_point(mapping = aes(x=displ, y=hwy)) +
    facet_grid(. ~ class)

# step 7 使用其他几何对象来表现数据
ggplot(data=mpg) +
    geom_point(mapping=aes(x=displ, y=hwy, color=class)) +
    geom_smooth(mapping=aes(x=displ, y=hwy, linetype=class, color=class))

# step 8 不同的几何对象的公有和私有图形属性映射
ggplot(data=mpg, mapping=aes(x=displ, y=hwy)) +
    geom_point(mapping=aes(color=class)) +
    geom_smooth()
# step 9 不同几何对象的私有数据对象, geom_smooth 只对 subcompact 生成拟合曲线
ggplot(data=mpg, mapping=aes(x=displ, y=hwy)) +
    geom_point(mapping=aes(color=class)) +
    geom_smooth(data=filter(mpg, class=="subcompact"), se=F)
```

### 尝试一个数据的不同画法


```r
# step 10 尝试一个数据的不同画法
ggplot(data=mpg, mapping=aes(x=displ, y=hwy)) +
    geom_point() +
    geom_smooth(se=F)
```

```r
ggplot(data=mpg, mapping=aes(x=displ, y=hwy)) +
    geom_point() +
    geom_smooth(se=F, mapping=aes(group=drv))
```

```r
ggplot(data=mpg, mapping=aes(x=displ, y=hwy, color=drv)) +
    geom_point() +
    geom_smooth(se=F)
```

```r
ggplot(data=mpg, mapping=aes(x=displ, y=hwy)) +
    geom_point(mapping=aes(color=drv)) +
    geom_smooth(se=F)
```

```r
ggplot(data=mpg, mapping=aes(x=displ, y=hwy, color=drv)) +
    geom_point() +
    geom_smooth(se=F, mapping=aes(linetype=drv))
```

```r
ggplot(data=mpg, mapping=aes(x=displ, y=hwy, color=drv)) +
    geom_point()
```

### 引入其他新元素如统计变换、坐标变换


```r
# step 10 引入统计变换
ggplot(data=diamonds) +
    stat_summary(mapping=aes(x=cut, y=depth), fun.ymin=min, fun.ymax=max, fun.y=median)
ggplot(data=diamonds, mapping=aes(x=cut, y=depth)) + geom_boxplot()

ggplot(data=diamonds) +
    geom_bar(mapping=aes(x=cut, y=..prop.., group=1))
ggplot(data=diamonds) +
    geom_bar(mapping=aes(x=cut, y=..prop..))

ggplot(data=mpg, mapping=aes(x=cty, y=hwy)) + geom_point(position = "jitter")

## step 11 引入坐标变换
mybar = ggplot(data=diamonds, aes(x=cut, fill=clarity)) +
    geom_bar()
mybar + coord_flip()
mybar + coord_polar()

ggplot(data=mpg, mapping=aes(x=cty, y=hwy)) +
    geom_point() + geom_abline()
ggplot(data=mpg, mapping=aes(x=cty, y=hwy)) +
    geom_point() + geom_abline() + coord_fixed()
```

## ggplot2 进阶

> 在for循环和函数里, ggplot的图形需要显示调用print(plot)才能画出来.

### ggplot2 图形分层语法总结


```r
ggplot(data=<DATA>) +
    <geom_function>(
        mapping = aes(<MAPPINGS>),
        stat = <STAT>,
        position = <POSITION>
    ) +
    <scale_function> +
    <coord_function> +
    <facet_function>
```

+ \<DATA\>, 数据：必须是一个数据框
+ \<MAPPINGS\>, 一组图形属性映射，设定数据集中的变量如何映射到该图层的图形属性
+ \<geom_function\>, 几何对象，指定在图层中进行绘图的几何对象类型
+ \<STAT\>, 统计变换，指定时候对元数据进行统计变换
+ \<POSITION\>, 位置调整，指定图形元素的位置，避免图形重合
+ \<coord_function\>, 坐标变换: coordinate system, 描述数据与图形所在平面的映射关系
+ \<facet_function\>, 分面设定: 条件作图或网格作图，分解数据成为子集并对后者作图、联合展示
+ \<scale_function\>, 标度调整指定数据-图形属性、数据-位置、数据-坐标等的映射关系( f(x) )

### 数据和图形属性映射
#### 数据与映射

ggplot函数对数据的类型要求为数据框(data.frame),它不会直接修改原数据，而是创建一个副本。你可以使用 “%+%” 来改变数据集。你要保证aes的变量(这里是carat,price,cut)都是来自于data(这里是diamonds)，以保证ggplot对象是自含型的。在定义映射时，我们要注意区分设定和映射。例如下方darkblue，在定义为映射时，“darkblue”作为只有一个元素的变量，然后对该变量进行映射，所以你会发现画出的图形上会出现分组信息（单个组darkblue）而图形颜色不是darkblue。而设定则是直接把图形颜色设定为darkblue。这在前面 gglot2 初探的 step4 有涉及。


```r
p = ggplot(data=diamonds, aes(carat, price, colour=cut))
p %+% mtcars		# 进行映射的数据集被改为mtcars

p + geom_point(aes(colour='darkblue'))	# 映射
p + geom_point(colour='darkblue') # 设定
```
#### 分组、群组

有时候，我们想对数据进行分组展示。就需要对映射指定分组变量。分组变量的指定与否会产生很明显的差别（群组）。

```r
ggplot(Oxboys, aes(age, height, group=Subject)) + geom_line()
ggplot(Oxboys, aes(age, height, group=1)) + geom_line()
```

群组几何对象还要考虑的是如何将个体的图形属性映射到整体的图形属性。线条和路径遵循差一原则：观测点比线段数目多一，第一条线段将使用第一条观测的图形属性，第二条使用第一条的，最后一条观测的图形属性不会被用到。

### 几何对象

geom，几何对象执行着图层的实际渲染，控制着生成的图像类型。每个几何对象都有一个默认的统计变换，且每个统计变换都有一个默认的几何对象。

<table><thead><tr><th>几何对象</th><th>描述</th><th>默认统计变换</th><th>图形属性</th></tr></thead><thead></thead><tbody><tr><td>abline</td><td>线,由斜率和截距决定</td><td>abline</td><td>colour,linetype,size</td></tr><tr><td>area</td><td>面积图</td><td>identity</td><td>colour,fill,linetype,size,x,y</td></tr><tr><td>bar</td><td>条形图, 以x轴为底的矩形</td><td>bin</td><td>colour,fill,linetype,size,weight,x</td></tr><tr><td>bin2d</td><td>2维热图</td><td>bin2d</td><td>colour,fill,linetype,size,weight,xmax,xmin,ymax,ymin</td></tr><tr><td>blank</td><td>空白,什么也不画</td><td>identity</td><td></td></tr><tr><td>boxplot</td><td>箱线图</td><td>boxplot</td><td>colour,fill,lower,middle,size,upper,weight,x,ymax,ymin</td></tr><tr><td>contour</td><td>等高线图</td><td>contour</td><td>colour,linetype,size,weight,x,y</td></tr><tr><td>crossbar</td><td>带有水平中心线的盒子图</td><td>identity</td><td>colour,fill,linetype,size,x,y,ymax,ymin</td></tr><tr><td>density</td><td>光滑密度曲线图</td><td>density</td><td>colour,fill,linetype,size,weight,x,y</td></tr><tr><td>density2d</td><td>二维密度等高线图</td><td>density2d</td><td>colour,linetype,size,weight,x,y</td></tr><tr><td>dotplot</td><td>点直方图,用点表示观测值的个数</td><td>bindot</td><td>colour, fill, x, y</td></tr><tr><td>errorbar</td><td>误差棒</td><td>identity</td><td>colour,linetype,size,width,x,ymax,ymin</td></tr><tr><td>errorbarh</td><td>水平的误差行</td><td>identity</td><td>colour,linetype,size,width,x,ymax,ymin</td></tr><tr><td>freqpoly</td><td>频率多边形图</td><td>bin</td><td>colour,linetype,size</td></tr><tr><td>hex</td><td>用六边形表示的二维热图</td><td>binhex</td><td>colour,fill,size,x,y</td></tr><tr><td>histogram</td><td>直方图</td><td>bin</td><td>colour,fill,linetype,size,weight,x</td></tr><tr><td>hline</td><td>水平线</td><td>hline</td><td>colour,linetype,size</td></tr><tr><td>jitter</td><td>给点添加扰动，减轻图形重叠问题</td><td>identity</td><td>colour,fill,shape,size,x,y</td></tr><tr><td>line</td><td>按照x坐标的大小顺序依次连接各个观测值</td><td>identity</td><td>colour,linetype,size,x,y</td></tr><tr><td>linerange</td><td>一条代表一个区间的竖直线</td><td>identity</td><td>colour,linetype,size,x,ymax,ymin</td></tr><tr><td>map</td><td>基准地图里的多边形</td><td>identity</td><td>colour,fill,linetype,size,x,y,map_id</td></tr><tr><td>path</td><td>按数据的原始顺序连接各个观测值</td><td>identity</td><td>colour,linetype,size,x,y</td></tr><tr><td>point</td><td>点,用来绘制散点图</td><td>identity</td><td>colour,shape,fill,size,x,y</td></tr><tr><td>pointrange</td><td>用一条中间带点的竖直线代表一个区间</td><td>identity</td><td>colour,fill,linetype,shape,size,x,y,ymax,ymin</td></tr><tr><td>polygon</td><td>多边形,相当于一个有填充的路径</td><td>identity</td><td>colour,fill,linetype,size,x,y</td></tr><tr><td>quantile</td><td>添加分位数回归线</td><td>quantile</td><td>colour,linetype,size,weight,x,y</td></tr><tr><td>raster</td><td>高效的矩形瓦片图</td><td>identity</td><td>colour,fill,linetype,size,x,y</td></tr><tr><td>rect</td><td>2维的矩形图</td><td>identity</td><td>colour,fill,linetype,size,xmax,xmin,ymax,ymin</td></tr><tr><td>ribbon</td><td>色带图,连续的x值所对应的y的范围</td><td>identity</td><td>colour,fill,linetype,size,x,ymax,ymin</td></tr><tr><td>rug</td><td>边际地毯图</td><td>identity</td><td>colour,linetype,size</td></tr><tr><td>segment</td><td>添加线段或箭头</td><td>identity</td><td>colour,linetype,size,x,xend,y,yend</td></tr><tr><td>smooth</td><td>添加光滑的条件均值线</td><td>smooth</td><td>alpha,colour,fill,linetype,size,weight,x,y</td></tr><tr><td>step</td><td>以阶梯形式连接各个观测值</td><td>identity</td><td>colour,linetype,size,x,y</td></tr><tr><td>text</td><td>文本注释</td><td>identity</td><td>angle,colour, hjust,label,size, vjust,x,y</td></tr><tr><td>tile</td><td>瓦片图</td><td>identity</td><td>colour,fill,linetype,size,x,y</td></tr><tr><td>violin</td><td>小提琴图</td><td>yidentity</td><td>weight,colour,fill,size,linetype,x,y</td></tr><tr><td>vline</td><td>竖直线</td><td>vline</td><td>colour,linetype,size</td></tr></tbody></table>

### 统计变换

stat，对数据进行统计变换，通常是对数据信息进行汇总。一个统计变换必须是一个位置尺度不变量，f(x+a)=f(x)+a && f(b\*x)=b\*f(x).统计变量通常会向原数据集中插入新的变量,你可以使用..xxx..的方式直接调用这些变量.


```r
ggplot(diamonds, aes(carat)) + geom_histogram(aes(y=..density..), binwidth=0.1)
```

![stats in ggplot](./stat_in_ggplot.png)


### 位置调整

位置调整参数 | 描述
---|---
dodge | 避免重叠，并排放置
fill | 堆叠图形元素并将高度标准化为1
identity | 不做任何调整
jitter | 添加扰动避免重合
stack | 把图形元素堆叠起来

### 标度
标度（scale）控制着数据到图形属性的映射，每一种标度都是从数据空间的某个区域（标度定义域）到图形属性空间的某个区域（标度值域）的一个函数.标度的执行过程分为3步，变换、训练和映射。变换是对数据进行汇总（统计变换），训练基于变换后的数据得到标度的定义域，最后映射到图形上去。使用时通过scale_aesType_scaleType来实现改变标度.
```
p = qplot(sleep_total, sleep_cycle, data=msleep, colour=vore)
p +	scale_colour_hue('what does it eat?',
					breaks=c('herbi', 'carni', 'omni', NA),
					labels=c('plants', 'meants', 'both', 'dont know'))
p + scale_colour_brewer(palette='Set1')
```

#### 标度分类

标度可以粗略分为4类：位置标度、颜色标度、手动离散型标度和同一型标度。

+ 位置标度：连续型、离散型、日期-时间型变量的映射到绘图区域，构造坐标轴  
+ 颜色标度：连续型、离散型变量映射到颜色  
+ 手动标度：离散型变量映射到我们选择的符号大小、线条类型、形状或颜色，创建对应图例  
+ 同一型标度：直接用变量值进行绘制图形属性，不执行映射。如颜色变量值  


```r
# name: 使用xlab(name),ylab(name),labs()这三个辅助函数进行设置
p + xlab('City mpg')
p + labs(x='City mpg', y='Highway', colour='Displacement')

# limits: 固定标度的定义域, 可以移除不想在图形上展示的数据/保持不同范围数据的范围一致
# breaks&labels: 控制坐标轴显示的断点和在该断点上的标签，二者相匹配
p + scale_x_continuous(breaks=c(5.5, 6.5))		# 坐标轴在5.5、6.5处进行断点
p + scale_x_continuous(limits=c(5.5, 6.5))		# 整个图x轴的范围是5.5-6.5

# formatter：未指定标签时，在断点处调用formatter来生成格式化的标签。可用的标签刷为：
# 对于连续型：comma, percent, dollar, scientific
# 对于离散型标度：abbreviate
```

#### 位置标度


```r
p + scale_x_log10() + scale_y_log10()
p + scale_x_date(
	limits=as.Date(c('2004-01-01', '2005-01-01')),
	labels=date_format('%Y-%m-%d')
)
```

#### 颜色标度

```r
# 连续型
p + scale_fill_gradient(limits=c(0,0.04))
p + scale_fill_gradient(limits=c(0,0.04), low='white', high='black')	# 双色标度
p + scale_fill_gradient2(limits=c(-0.04,0.04), midpoint=0.02)			# 三色标度, gradientn是多色标度

# 离散型
p + scale_fill_brewer(palette='Set2')
p + scale_colour_brewer(palette='Pastel1')
```

#### 手动离散标度

```r
# 离散型标度
p + scale_linetype()
p + scale_size_discrete()
p + scale_shape()

# 如果定制这些离散型标度,则手动创建
p + scale_colour_manual(values=c('red','orange', 'yellow', 'green', 'blue'))
p + scale_shape_manual(values=c(1,2,6,0,23))
# 实际应用
ggplot(huron, aes(year)) + 
	geom_line(aes(y=level-5, colour='below')) +
	geom_line(aes(y=level+5, colour='above')) +
	scale_colour_manual('Direction', values=c('below'='blue', 'above'='red'))
```

#### 同一型标度
当你的数据能被R中的绘图函数理解（数据空间和图形属性空间相同时），可以使用同一型标度。

## 常见图表代码


```r
attach(msleep)
```

```r
mdata = msleep
```

```r
festival.data = read.table('DownloadFestival.dat', sep='\t', header=T)
```

```
## Warning in file(file, "rt"): cannot open file 'DownloadFestival.dat': No
## such file or directory
```

```r
head(mdata)
```

```r
head(festival.data)
```

### 常见几何对象及参数

+ geom_bar(x,color,size,fill,linetype,alpha): 创建代表不同统计性质的条形图图层
+ geom_point(x,y,shape,color,size,fill,alpha): 创建数据点图层
+ geom_line(x,y,color,size,linetype,alpha): 创建直线图层
+ geom_smooth(x,y,color,size,linetype,alpha): 创建平滑曲线图层(类似于回归曲线/点连线,但是更为平滑的曲线)
+ geom_histogram(x,color,size,fill,linetype,alpha): 创建柱状图图层
+ geom_boxplot(x,color,siz,fill,alpha): 创建box-whisker图图层
+ geom_text(x,y,color,size,angle,hjust,vjust,alpha): 创建文本图层
+ geom_errorbar(x,ymin,ymax,color,linetype,width,alpha): 创建误差线图层
+ geom_hline(yintercept)/vline(xintercept,color,size,linetype,alpha): 创建水平线或垂直线图层

### 点图

```r
{
    scatterplot = ggplot(data = mdata, aes(x=bodywt, y=sleep_total)) + geom_point()
    scatterplot
    # 添加颜色
    scatterplot = ggplot(data = mdata, aes(x=bodywt, y=sleep_total, col=vore)) + geom_point()
    scatterplot
    # 数据变换
    scatterplot = ggplot(data = mdata, aes(x=log(bodywt), y=sleep_total, col=vore)) + geom_point()
    scatterplot
    # 改变各种外观
    scatterplot = scatterplot + geom_point(size=5) + 
                                xlab('Log Body Weight') + 
                                ylab('Total Hours Sleep') +
                                ggtitle('Some Sleep Data')
    # 修改x轴刻度线，title位置
    scatterplot = scatterplot + scale_color_brewer(palette = 'Set2')
    scatterplot = scatterplot + theme(plot.title=element_text(vjust=+2))+scale_x_continuous(breaks=-5:10)
    # 调整背景颜色
    scatterplot = scatterplot + theme_set(theme_bw(base_size=18))
    scatterplot
}
```

### 带状点图

```r
{
    my.stripchart = ggplot(data=mdata, aes(vore, sleep_total)) + geom_point()
    # 改变大小，位置(左右偏移一些)
    my.stripchart = ggplot(data=mdata, aes(vore, sleep_total)) + geom_point(size=5, position='jitter')
    # 通过引入jitter图层来改变位置,默认用于point，所以不用在geom_point中指定
    my.stripchart = ggplot(data=mdata, aes(vore, sleep_total, col=vore)) +
                    geom_jitter(position=position_jitter(width = 0.2), size=5)
    # 添加误差线
    my.stripchart = my.stripchart + 
                    ylab('Total Hours Sleep') + xlab('Trophic Level') + 
                    stat_summary(fun.data = mdata, geom = "errorbar", width = .5) +
                    stat_summary(geom = "errorbar", fun.y = mean, aes(ymin = ..y.., ymax = ..y..))
    # 添加标题，去掉NA，调整y轴
    my.stripchart = my.stripchart +
                    scale_x_discrete(limits=c('carni', 'herbi', 'insecti', 'omni')) +
                    ggtitle('Some Sleep Data') +
                    theme(plot.title=element_text(vjust = +2)) +
                    scale_y_continuous(breaks=seq(0,20,2))
    my.stripchart
}
```

### 柱状图

```r
{
    my.hist = ggplot(data=festival.data, aes(x=day1)) + geom_histogram()
    # 上点颜色咯
    my.hist = ggplot(data=festival.data, aes(x=day1)) + 
                geom_histogram(binwidth=0.3, color='black', fill='yellow') +
                labs(x='Score', y='Counts') + ggtitle('Hygiene at Day 1')
    my.hist
    
    # 构造一下数据
    festival.data.stack = melt(festival.data, id=c('ticknumb', 'gender'))
    colnames(festival.data.stack)[3:4]<-c('day','score')
    head(festival.data.stack)

    # 运用facet进行分组,有两种方式 facet_grid(), facet_wrap()
    my.hist.day3 = ggplot(data=festival.data.stack,aes(score)) + 
                geom_histogram(binwidth=0.3, color='black', fill='yellow') +
                labs(x='Score', y='Counts') + ggtitle('Hygiene at Day 1') +
                facet_grid(gender~day)
    my.hist.day3
}
```

### 画个stripchart;

```r
{
    # color用于分亚组，gender用于分大组；所以color分组的话与x是一致的分组
    # facet_gripd用于分大组
    my.stripchart.day3 = ggplot(data=na.omit(festival.data.stack), aes(day,score, color=day)) +
                geom_point(position='jitter') + facet_grid(~gender)
    my.stripchart.day3 = my.stripchart.day3 + 
                stat_summary(fun.data = na.omit(festival.data.stack), geom = "errorbar", width = .5) +
                stat_summary(geom = "errorbar", fun.y = mean, aes(ymin = ..y.., ymax = ..y..))
# 带箱图的stripchart
    my.stripchart.day3 = ggplot(data=na.omit(festival.data.stack), aes(day,score, color=day)) +
                geom_point(position='jitter') + facet_grid(~gender) +
                scale_colour_manual(values=c("darkorange", "darkorchid4",'blue')) +
                geom_boxplot(alpha=0, colour="black")
    my.stripchart.day3
}
```

### 箱图

```r
{
    my.boxplot = ggplot(data = festival.data.stack, aes(gender,score, fill=gender)) + geom_boxplot() + facet_grid(~day)
    # scale_fill_manual要配合fill=gender来使用
    my.boxplot = my.boxplot + scale_fill_manual(values=c('orange', 'purple')) + stat_boxplot(geom='errorbar')
    # 改变x-轴标签顺序
    my.boxplot = my.boxplot + scale_x_discrete(limits=c('Male', 'Female'))
    my.boxplot    
}
```

### 小提琴图

```r
{
    my.violinplot = ggplot(festival.data.stack, aes(gender, score, fill=gender)) + 
                    geom_violin(trim=FALSE) +
                    facet_grid(~day) + 
                    stat_summary(fun.data = festival.data.stack, geom = "errorbar", width = .5) +
                    stat_summary(geom = "errorbar", fun.y = mean, aes(ymin = ..y.., ymax = ..y..)) +
                    scale_fill_manual(values=c('orange', 'yellow'))
    my.violinplot
}
```


### 条形图


```r
score.sem = data.frame(gender=rep(c('Female', 'Male'), each=3),
                       day=rep(c('day1', 'day2', 'day3'), 2),
                       mean=c(1.88, 1.08, 1.10, 1.60, 0.77, 0.83),
                       sem=c(0.032, 0.061, 0.099, 0.036, 0.058, 0.073))
head(score.sem)
```

```
##   gender  day mean   sem
## 1 Female day1 1.88 0.032
## 2 Female day2 1.08 0.061
## 3 Female day3 1.10 0.099
## 4   Male day1 1.60 0.036
## 5   Male day2 0.77 0.058
## 6   Male day3 0.83 0.073
```

```r
{
    # geom_errorbar要和stat='identity'一起使用,不然无法知道该进行何类统计
    my.barplot = ggplot(data = score.sem, aes(day, mean, fill=gender)) +
                 geom_bar(stat='identity', position='dodge', color='black', size=1) +
                 geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position=position_dodge(width=0.8), size=1, width=.5) +
                 ylab('Mean Scores') + ggtitle('Levels of hygiene over 3 days of concert') +
                # 去掉x轴标签
                 theme(axis.title.x=element_blank()) + 
                 theme(plot.title=element_text(vjust=+2))
    my.barplot
}
```

### 堆积条形图/百分比堆积图

```r
Changing = read.csv('Changing.csv')
```

```
## Warning in file(file, "rt"): cannot open file 'Changing.csv': No such file
## or directory
```

```r
head(Changing)
```

```r
{
    my.stackbarchart = ggplot(Changing, aes(Type.of.Behaviour, Sample.Size, fill=Stage.of.Change)) + 
                geom_bar(stat='identity') +
                coord_flip() # 水平和垂直堆积变换
    
    # 由于各个堆积图内部的排序都不一致，把stage.of.change改为因子就可以
    Changing$Stage.of.Change = factor(Changing$Stage.of.Change, 
                                      levels=c("Precontemplation", "Contemplation", "Preparation", "Action",
                                                                         "Maintenance"))
    my.stackbarchart = ggplot(Changing, aes(Type.of.Behaviour, Sample.Size, fill=factor(Stage.of.Change))) + 
        geom_bar(stat='identity') +
        coord_flip() + 
        scale_fill_brewer(palette=3) + 
        labs(title='Stages of Each of the 12 Problems Behaviours',
             y='Type of Behaviour',
             x='Sample Size',
             fill='Stage of Change')
    my.stackbarchart
    
    # 画成百分比堆积图
    contingency.table = xtabs(Sample.Size~Type.of.Behaviour+Stage.of.Change, Changing)
    contingency.tab100 = prop.table(contingency.table, 1)
    contingency.tab100 = contingency.tab100*100
    contingency.percent = as.data.frame(contingency.tab100)
    contingency.percent$Stage.of.Change = factor(contingency.percent$Stage.of.Change, 
                                                 levels=c("Precontemplation", "Contemplation", "Preparation", "Action",
                                                          "Maintenance"))
    my.stackbarchart.percent = ggplot(contingency.percent, aes(Type.of.Behaviour, Freq, fill=factor(Stage.of.Change))) + 
        geom_bar(stat='identity') +
        coord_flip() + 
        scale_fill_brewer(palette='RdYlGn') + 
        labs(title='Freq of Each of the 12 Problems Behaviours',
             y='Freq',
             x='Type of Behaviour',
             fill='Stage of Change')
    
    my.stackbarchart.percent
}
```

### 线图

```r
{
    my.linegraph = ggplot(data = score.sem, aes(day, mean, group=gender, color=gender)) + 
                    geom_line(size=1.5) + geom_point() + 
                    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem))
    
    my.linegraph = ggplot(data = score.sem, aes(day, mean, group=gender, color=gender)) + 
        geom_line(size=1.5) + geom_point(size=5) + 
        geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2, size=1.5) +
        ylab('Mean Scores') + ggtitle('Levels of hygiene over 3 days concert') + 
        theme(axis.title.x = element_blank()) + 
        scale_color_manual(values=c('lightblue', 'steelblue')) + 
        theme(legend.justification = c(1,1), legend.position = c(0.95,0.95)) + 
        theme(legend.text = element_text(size=16, face='bold')) + 
        theme(legend.title=element_blank(), plot.title=element_text(vjust=+2))
    my.linegraph
}
```

## 参考

+ 《R 数据科学》
+ [ggplot reference](https://ggplot2.tidyverse.org/reference/)
+ 《ggplot2：数据分析与图形艺术》
+ Tutorials from Anne Segonds-Pichon in Babraham Bioinformatics (links not found)
