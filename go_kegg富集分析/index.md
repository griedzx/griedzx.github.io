# Go_kegg富集分析

# GO_KEGG富集分析

## 富集分析

富集分析（ ***Enrichment Analysis***），是一种识别基因集合与已知*生物过程*、*细胞组分*和*通路*之间关联的统计方法。这些工具通过使用数据库中的注释信息来找到对应的基因集合。

### 富集分析的效果、目的

1. 把差异基因或者物质根据其功能进行归类，使具有相似功能的基因或者物质就被放在一起
2. 实现功能和表型相关联
3. 解读一组基因背后所代表的生物学知识，揭示其在细胞内或细胞外扮演了什么样的角色。

### 富集分析的原理

富集分析通常是分析一组基因在某个功能节点上是否相比于随机水平过于出现(over-presentation)

### 常用方法

1. 目前最常用的方法是基于GO和KEGG的富集分析。 首先通过多种方法得到大量的感兴趣的基因，例如差异表达基因集、共表达基因模块、蛋白质复合物基因簇等，然后寻找这些感兴趣基因集显著富集的GO节点或者KEGG通路，这有助于进一步深入细致的实验研究。
2. 依据富集分析过程中基因选择、注释数据库的不同，常用的富集分析可以分为以下四种类型：GO term功能富集、KEGG pathway通路富集、MSigDB基因集富集和单基因富集等等

## GO & KEGG是什么

对于每个基因而言，其基本的功能基于他们的蛋白结构域以及研究的文献已经可以大致的知道一个基因具有什么样子的功能了。**GO和KEGG就是基于不同的分类思想而储存的基因相关功能的数据库。**

### GO数据库

全称是Gene Ontology(基因本体)，他们把基因的功能分成了三个部分分别是： **细胞组分（cellular component, CC）、分子功能（molecular function, MF）** 、 **生物过程（biological process, BP）** 。利用GO数据库，我们就可以得到我们的目标基因在CC, MF和BP三个层面上，主要和什么有关。

### KEGG数据库

除了对基因本身功能的注释，我们也知道基因会参与人体的各个通路，基于人体通路而形成的数据库就是通路相关的数据库。

京都基因与基因组百科全书（Kyoto encyclopedia of genes and genomes, KEGG）是系统分析基因功能、基因组信息的数据库，整合了基因组学、生物化学及系统功能组学的信息，有助于研究者把基因及表达信息作为一个整体进行研究。目前KEGG共包含了19个子数据库，富集分析常用在**KEGG Pathway**通路中。

## GO_KEGG_GSEA分析

### 下载加载需要的包

#### 下载

```R
install.packages("包名称")
```

或者

```R
biocManager::install("包名称")
```

1. `BiocManager`：是一个 R 包，其中包含了一些用于管理 Bioconductor 包的函数，如 `install()`、`update()` 等。使用 `::` 表示法，你可以指定调用特定包中的函数。通过 BiocManager，你可以轻松地访问和管理 Bioconductor 存储库中的生物信息学和生物统计学相关的软件包。它简化了安装和维护 Bioconductor 包的过程。
2. `Bioconductor`：是一个专门用于生物信息学和生物统计学研究的 R 软件包存储库。它是一个巨大的代码资源库，包含了许多用于基因表达数据分析、基因组学、蛋白质组学等生物学领域的软件包。Bioconductor 提供了丰富的分析工具和算法，能够帮助生物学家处理和解释生物学数据。BiocManager 用于访问和安装 Bioconductor 中的软件包。
3. `::` 不仅可以用于指定调用的包，还可以在不加载整个包的情况下调用该包中的函数。这样可以在特定代码行中临时使用特定包中的函数，而不需要在整个会话中加载该包。`::` 只在当前代码行中生效

#### 加载

```R
library(clusterProfiler)#GO&KEGG
library(enrichplot)#GO&KEGG
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(forcats)#绘图中对因子处理
```

### 读取数据，将基因ID从GENE_SYMBOL转换为ENTREZ_ID

```R
#载入数据，只需要基因ID(GO,KEGG,GSEA需要)和log2FOldChange(GSEA(基因集富集分析)需要)
data <- read.csv("DEG.csv",header = T)
#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种
KEGG_database <- 'hsa' #KEGG分析指定物种
```

GO分析指定物种，详见 [GO物种缩写索引表](http://bioconductor.org/packages/release/BiocViews.html#___OrgDb)

KEGG分析指定物种，详见 [keggle物种缩写索引表](http://www.genome.jp/kegg/catalog/org_list.html)

```R
#安装注释数据库
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

```

有问题参见官网 [Bioconductor - org.Hs.eg.db](https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

#### OrgDb

org数据库是一个用于存储特定物种的基因和基因注释信息的数据库。每个org数据库专门用于一个特定的生物物种，例如人类、小鼠、大鼠等。org数据库包含了与该物种相关的基因、基因产物以及它们的注释信息，如基因的功能、组织表达、基因本体注释等。在生物学研究中，研究人员经常需要将基因与GO术语相关联，从而了解基因的功能和参与的生物过程。为此，他们需要使用org数据库中的基因注释信息，将基因与GO术语进行映射，以便进行GO富集分析、基因集对比等功能。

1. 在 `http://bioconductor.org/packages/release/BiocViews.html#___OrgDb`可以找到对应物种

```R
#gene ID转换
gene <- bitr(data$gene,fromType = 'SYMBOL', toType = 'ENTREZID',OrgDb = GO_database)
```

2. 如果自己研究的物种不在 `http://bioconductor.org/packages/release/BiocViews.html#___OrgDb` 之列，很大可能就需要自己构建OrgDb，然后用clusterProfiler分析。非模式生物要想找到自己的注释包，又分成两类：

* 一类是在 [AnnotationHub](https://links.jianshu.com/go?to=https%3A%2F%2Fbioconductor.org%2Fpackages%2Frelease%2Fbioc%2Fhtml%2FAnnotationHub.html) 中存在的，例如玉米
* 另一类是在AnnotationHub也不存在相应物种，就需要用 [AnnotationForge](https://links.jianshu.com/go?to=https%3A%2F%2Fbioconductor.org%2Fpackages%2Frelease%2Fbioc%2Fhtml%2FAnnotationForge.html) 来自己构建

若出现

```R
> gene <- bitr(data$gene,fromType = 'SYMBOL', toType = 'ENTREZID',OrgDb = GO_database)
'select()' returned 1:many mapping between keys and columns
Warning message:
In bitr(data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = GO_database) :
  0.7% of input gene IDs are fail to map...
```

一对多映射：某些基因 ID（SYMBOL）可能对应多个基因 ID（ENTREZID）。这种情况通常是因为一个基因可能有多个变体或不同的数据库中存在多个记录。

### GO分析

```R
#GO分析
GO <- enrichGO(gene$ENTREZID,
               OrgDb = GO_database,
               keyType = "ENTREZID",#设定读取的gene ID类型
               ont = 'ALL', #ont=all即包括GO数据库中的三个部分CC、MF、BP
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = T)
```

`readable` 参数用于控制 GO 富集分析结果中的基因 ID 是否转换成易读的描述性名称,即ENTREZ_ID -> SYMBOL_ID,默认为 `False`

### KEGG分析 （参数和GO差不多）

```R
KEGG <- enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
```

GO & KEGG 都可以保存为csv文件查看

`write.csv(GO,"./result.csv")`

`GO_result <- GO@result` 所得变量存储格式和输出文件格式相同，可以查看选择感兴趣的定制画图

## GO & KEGG结果可视化

### 富集柱状图 + 点状图

```R
barplot(GO, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale = "free")  +
  scale_y_discrete(labels=function(x) str_wrap(x, width = 100))

dotplot(GO, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale = "free")
```

* `barplot()` 函数用于创建柱状图
* 参数 `split = "ONTOLOGY"` 表示按照 "ONTOLOGY" 列的不同取值（GO_db的三个层面）将数据拆分为不同的柱状图
* `facet_grid()` 函数用于将拆分后的柱状图按照 "ONTOLOGY" 列的不同取值排列在一个网格中。每一行代表一个 "ONTOLOGY" 类别
* 子图根据数据的分布自动调整高度，即参数 `scale = "free"`，以便更好地展示数据的差异
* `barplot(GO, showCategory=20)` 参数 showCategory 可以控制展示的数量
* scale_y_discrete(labels=function(x) str_wrap(x, width = 100)) #调整y轴标签长度，使其放在一行

### 棒棒糖图

```R
ggplot(GO, showCategory = 20, 
       aes(GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  ylab(NULL) 
```

`aes(GeneRatio, fct_reorder(Description, GeneRatio))`: 这里设置了图形的"美学映射"（aesthetic mappings）

* GeneRatio 是x轴的变量，fct_reorder(Description, GeneRatio) 是y轴的变量
* fct_reorder() 函数可以根据 GeneRatio 的值对 Description 进行重新排序，这样可以根据 GeneRatio 的大小对y轴的标签进行排列。

### 富集基因与所在功能集/通路集的关联网络图

```R
enrichplot::cnetplot(GO, circular = F, colorEdge = T)
```

* `circluar`为指定是否环化，基因过多建议设置成 `F`

### 富集到的功能集/通路集之间的关联网络图

```
GO2 <- pairwise_termsim(GO)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")#通路间关联网络图
```

热图展现关联关系

```
enrichplot::heatplot(GO,showCategory = 50)#基因-通路关联热图
enrichplot::heatplot(KEGG,showCategory = 50)
```

### GO富集功能网络图

```R
GO_BP<-enrichGO( gene$ENTREZID,#GO富集分析BP模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_BP)#GO-BP功能网络图
GO_CC<-enrichGO( gene$ENTREZID,#GO富集分析CC模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "CC",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_CC)#GO-CC功能网络图
GO_MF<-enrichGO( gene$ENTREZID,#GO富集分析MF模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "MF",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_MF)#GO-MF功能网络图
```

### GO富集弦图

先要将所有GO富集到的基因集所对应的类型写入本地文件从而得到BP/CC/MF各自的起始位置如我的数据里是1，2103，2410.

```R
genedata<-data.frame(ID=info$gene_symbol,logFC=info$log2FoldChange)
write.table(GO$ONTOLOGY, file = "/Users/ZYP/Downloads/KEGG_GO/GO_ONTOLOGYs.txt", #将所有GO富集到的基因集所对应的类型写入本地文件从而得到BP/CC/MF各自的起始位置如我的数据里是1，2103，2410
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

```

```R
GOplotIn_BP<-GO[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC<-GO[2103:2112,c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MF<-GO[2410:2419,c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
GOplotIn_BP$Category = "BP"#分类信息
GOplotIn_CC$Category = "CC"
GOplotIn_MF$Category = "MF"
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) #GOplot导入数据格式整理
circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata) 
circ_MF<-GOplot::circle_dat(GOplotIn_MF,genedata) 
chord_BP<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
chord_CC<-chord_dat(data = circ_CC,genes = genedata) 
chord_MF<-chord_dat(data = circ_MF,genes = genedata) 
GOChord(data = chord_BP,#弦图
        title = 'GO-Biological Process',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 10) #GO Term字体大小
GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10) 
GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10)
```

## 参考文章

1. [GO分析和KEGG分析都是啥？ - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/135410211)
2. [最全的GO, KEGG, GSEA分析教程(R),你要的高端可视化都在这啦！[包含富集圈图] - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/377356510)

