# maize_gene_id转换

## 问题描述

做GO富集分析时，师姐给的数据是Genome assembly B73 RefGen_v4版本的数据

* 对应基因id没有需要使用的symbol或者entrez_id形式
* 并且大多数在线GO分析网站和id转换网站不支持v4版本基因id与symbol或者entrez_id相互转换

Genome assembly B73 RefGen_v4

Genome assembly B73 RefGen_v4是玉米的第四版参考基因组序列，它是由美国农业部农业研究局（USDA-ARS）和美国能源部联合基因组研究所（DOE-JGI）的科学家合作完成的。它是目前最完整、最准确的玉米基因组序列，它包含了21条染色体的端粒到端粒的连续序列，共有2.06亿个碱基对，覆盖了99.8%的可编码区域，注释了39,656个基因。

## 解决方法

### 网页在线转换

[MaizeGDB gene Search Page](https://chinese.maizegdb.org/gene_center/gene#translate)

MaizeGDB数据库有一个[Translate Gene Model IDs](https://chinese.maizegdb.org/gene_center/gene#translate)工具，可以识别各种类型的玉米的gene_id，并转换为比较常用的几种id(如)

### 爬虫技术解决

[R语言爬取NCBI大豆基因Locus tag数据 - 简书 (jianshu.com)](https://www.jianshu.com/p/4cf02b79f574)

[生信笔记01：Locus tag转换为Entrez Gene ID - 简书 (jianshu.com)](https://www.jianshu.com/p/6513f0a3ceb7)

通过

