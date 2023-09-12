# clusterProfiler包GSEA富集

GESA富集分析可以使用 java 版本的 `GSEA`软件进行富集。但是这个软件有一个问题就是自带的物种数据库有限（主要是人、小鼠的 PMT 文件、且不支持 KEGG 库），如果想要分析一些其他物种，需要上传自己准备的 PMT 文件。

这时候有一个可选的方案就是用 `clusterProfiler`包进行 GSEA 富集分析，该软件提供了两种方式进行富集：

1. `gseGO`、`gseKEGG`等函数分别从 `OrgDb`、KEGG 官网读取或者下载 PMT 基因集
2. `GSEA`函数，一个通用的 GSEA 富集框架，支持从本地读取自己已经准备好的 PMT 基因集

## gseGO

```R
gseGO(
  geneList,
  ont = "BP",
  OrgDb,
  keyType = "ENTREZID",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
  ...
)
```

` geneList	order ranked geneList`

根据 `logFC`(log2folderchange)列进行 `降序`排列（上调基因在顶部，下调在底部），制作 `gene_list`集(logFC对应元素名为对应基因id):

```R
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
geneList <- sort(GSEA_input, decreasing = TRUE)
```

gseGO:GO 本体富集，可以进行GO本体的GSEA富集

```text
gsea_go <- gseGO(
  gene_list,    # 根据logFC排序的基因集
  ont = "BP",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
  OrgDb = maize, #orgdb数据库  
  keyType = "ENTREZID",    # 基因id类型
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",    # p值校正方法
)
```

查看分析结果数据和对应列名 `head(gsea_go,2)`

* ID ：GO term 或 KEGG 通路
* Description ：GO term 描述信息
* setSize ：富集到该 term 的基因个数
* enrichmentScore ：富集分数，也就是 ES
* NES ：标准化以后的 ES，全称 normalized enrichment score
* pvalue：富集的 P 值
* p.adjust ：校正后的 P 值
* qvalues ：FDR （false discovery rate）错误发现率
* rank ：当 ES 最大时，对应基因所在排序好的基因列表中所处的位置
* leading_edge：tags 表示核心基因占该通路基因集的百分比；list 表示核心基因占所有基因的百分比；signal，将前 2 项统计值结合在一起计算出的富集信号强度
* core_enrichment：核心或者 leading 基因列表。

### 上、下调的 GO term 分开展示：

```R
dotplot(
  gsea_go,
  showCategory=10,
  split=".sign") + facet_grid(.~.sign)
```

