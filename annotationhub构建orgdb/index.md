# Annotationhub构建orgdb

# AnnotationHub构建OrgDb

要进行GO或者KEGG富集分析，就需要知道每个基因对应什么样的GO/KEGG分类，OrgDb就是存储不同数据库基因ID之间对应关系，以及基因与GO等注释的对应关系的 R 软件包

如果自己研究的物种不在 `http://bioconductor.org/packages/release/BiocViews.html#___OrgDb` 之列，很大可能就需要自己构建OrgDb，然后用clusterProfiler分析

## 利用AnnotationHub得到org.db

其中一种情况是在（[AnnotationHub](https://links.jianshu.com/go?to=https%3A%2F%2Fbioconductor.org%2Fpackages%2Frelease%2Fbioc%2Fhtml%2FAnnotationHub.html)）中存在对应的注释包

```R
require(AnnotationHub)
hub <- AnnotationHub()#下载失败多试几次
query(hub,"zea mays")

AnnotationHub with 8 records
# snapshotDate(): 2023-04-24
# $dataprovider: NCBI,DBCLS, ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/, WikiPathways
# $species: Zea mays, Zea mays_var._japonica
# $rdataclass: SQLiteFile, OrgDb, Tibble
# additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
#   rdatapath, sourceurl, sourcetype 
# retrieve records with, e.g., 'object[["AH91642"]]' 

             title                    
  AH91642  | MeSHDb for Zea mays (Corn, v001)   
  AH91817  | wikipathways_Zea_mays_metabolites.rda
  AH97909  | MeSHDb for Zea mays (Corn, v002)   
  AH100374 | MeSHDb for Zea mays (Corn, v003)   
  AH107139 | MeSHDb for Zea mays (Corn, v004)   
  AH111528 | MeSHDb for Zea mays (Corn, v005)   
  AH111691 | org.Zea_mays.eg.sqlite   
  AH111692 | org.Zea_mays_var._japonica.eg.sqlite 
```

阅读命令输出信息，可以看到

* 数据库来源不同，选择需要的（org）以及确定物种（zea_mays）
* `retrieve records with, e.g., 'object[["AH111691"]]'` 检索下载db用指定格式

