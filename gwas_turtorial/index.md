# Gwas_turtorial

## GWAS主要步骤

1. 选择一个感兴趣的性状或疾病，以及一个适当的样本群体，如病例和对照组，或者具有连续性状值的个体。
2. 对样本进行全基因组分型，即测定数百万个SNP位点的基因型。
3. 对每个SNP位点进行关联分析，即计算其与性状或疾病的相关性强度和显著性水平，通常使用回归模型或卡方检验等统计方法。
4. 根据设定的显著性阈值，筛选出具有显著关联信号的SNP位点，并绘制曼哈顿图（manhattan plot）来展示全基因组关联结果。
5. 对显著关联信号进行进一步的验证和解释，如进行基因型插补（genotype imputation），精细定位（fine-mapping），功能注释（functional annotation），共定位分析（colocalization analysis）等，以确定最可能的因果变异和相关基因。

## GWAS分析基本流程及分析思路

[GWAS分析基本流程及分析思路 (qq.com)](https://mp.weixin.qq.com/s?__biz=Mzg2MDA2MDQzMQ==&mid=2247483663&idx=1&sn=acfacdf0a0ee6df2c003875c0db06476&chksm=ce2d6f6ff95ae67902343fd81768fa0e3cf48949d8ea01e1c3d50f4e09763e779906fda8de54&scene=21#wechat_redirect)

## GWAS的数据QC

qc,(quality  control)质量控制，比起表观学研究，GWAS研究很少有引起偏差的来源，一般来说，一个人的基因型终其一生几乎不会改变的，因此很少存在同时影响表型又影响基因型的变异。但即便这样，我们在做GWAS时也要去除一些可能引起偏差的因素。

这种因素主要有：群体结构、个体间存在血缘关系、技术性操作。

### **怎么进行质量控制（PLINK）**

质量控制包括两个方向，一个是样本的质量控制，一个是SNP的质量控制

**必须先进行SNP过滤，才能进行个体过滤 ???**

#### PLINK

对于PLINK来说，它既可以处理文本格式的文件，也可以处理二进制格式的文件。但是大文本的文件处理起来十分消耗计算资源，所以我们一般 **推使用二进制格式的输入文件** 。

文本格式的PLINK数据包括两份文件 .ped文件 和 .map文件

* ped文件 包含个体信息（例如个体标识符ID，性别等等）以及他们的基因型信息
* map文件 包含遗传标记的信息（染色体号，snp号等等）

**二进制格式的PLINK数据则包括三份文件** ：.bed文件，.fam文件 和 .bim文件

* **bed文件** 含有  *每个个体的识别符（ID）* 和*每个个体相对应的基因型*
* **fam文件** 含有  *个体信息* （例如性别之类的）
* **bim文件** 含有  *遗传标记的信息* （染色体号，snp号等等）


PLIINK基础使用命令

`plink --bfile MY_DATA --assoc --out gwas_results`

‐‐file {your_file}输入文本格式文件

**‐‐bfile {your_file} 输入二进制格式文件**

**--assoc 关联分析,这一步会对每个SNP和研究者感兴趣的性状进行卡方检验**

**--out {outfile} 输出文件**


#### 样本质量控制

样本的质量控制包括：缺失率、杂合性、基因型性别和记录的性别是否一致。

检测缺失率(个体缺失)，一般将样本缺失率大于5%的个体去除

`plink --bfile file --mind 0.05 --make-bed --out file_mind`

检测杂合性

`plink --bfile file --het --make-bed --out file_het`

* 阈值设置：偏差可能表明样品污染，近亲繁殖。我们**建议去除偏离样本的杂合率平均值±3 SD的个体**

检测性别不一致的个体 --

受试者信息填写的性别和遗传性别不一致

`plink --bfile file --check-sex --make-bed --out file_checksex`

将上述筛选出来的不符合的样本去除

`plink --bfile file --remove removesample.txt --make-bed --out file_qcsample`

## 处理群体分层

## 参考

1. [一篇手把手教你做GWAS的Guideline文献解读 - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/148905500)
2. [全基因组关联分析学习资料（GWAS ](https://zhuanlan.zhihu.com/p/90414014)[tutorial） - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/90414014)
3. [easyGWAS - Running GWAS easily over the web (mpg.de)](https://easygwas.biochem.mpg.de/)

