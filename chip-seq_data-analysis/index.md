# Chip-seq_data-analysis

# Chip-seq数据分析

## 背景知识

* ChIP-seq:染色质免疫共沉淀技术，也称作结合位点分析法。用于转录因子结合位点或蛋白特异性修饰位点的研究
* Chip-seq数据分析的使用工具和对应的流程
  1.bowtie aligber   原始数据与比对参考基因组(hg19)比对
  2.picard  tools      pcr数据去重 (pcr duplicating)
  3.MACS2              call peak
  4.deeptools          可视化
* 4.homer annotatePeaks.pl   peak注释到对应基因上
* 染色质开放区：(OCRs) 基因组DNA不与组蛋白缠绕而形成的区域，可与蛋白分子(TF)直接结合调控下游基因表达。

## Chip-seq实验顺序

![1687001002814](image/Chip-seq_data-analysis/1687001002814.png)

## Chip-seq数据处理过程

1. Chip-seq原始数据
2. （人类）参考基因组:用于建立index
3. bowtie2比对：在建库和测序（我认为就是上述实际实验中）后，read之间无位置关系，需要read与参考基因组比较，在参考基因组上定位
4. PCR去重复
5. call peak:peak意味着有序列，在chip-seq中指对应区域蛋白结合在dna上
   5.1  某些特殊的chip-seq实验可以用rose筛选超级增强子
6. 下游分析（注释）

![1687000979344](image/Chip-seq_data-analysis/1687000979344.png)

## 软件介绍以及安装

### conda

conda 强大的包管理器和环境管理器

因为在家使用个人云服务器所以从头安装一遍:

`wget -c http://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh`

`bash Miniconda3-latest-Linux-x86_64.sh`

#### **激活刚安装完成的软件**

一般安装软件完成后需要重启，在Linux叫激活，有两种方式，**第一种**是重新登录服务器，**第二种**是输入以下命令：

`source ~/.bashrc `

使用conda安装相应的软件

#### 一次性安装所需要的所有软件

`conda install -y -c bioconda  bowtie2 samtools picard macs2 deeptools`

最好还是逐个安装,不同的软件可能适配的python环境不同，可能需要conda新的虚拟环境

安装代码出现问题可以上[anaconda.org](https://anaconda.org/bioconda "安装出错")查询，会有推荐代码

#### 安装中出现的问题

![1687000960404](image/Chip-seq_data-analysis/1687000960404.png)

存在两个冲突：

Python 版本不符合要求：deeptools 要求的 Python 版本是 2.7.x 或 3.x（其中包括了 2.7.x、3.5.x、3.6.x 等），而你的环境中安装的是 python=3.11。这会导致 deeptools 无法与你的 Python 版本兼容。

glibc 版本不符合要求：deeptools 还依赖于 libgcc-ng，而 libgcc-ng 又依赖于 __glibc[version='>=2.17']。然而，你的系统中安装的 glibc 版本是 2.35，与要求的版本不兼容。

当使用 conda 创建python=2.7 环境，这两个问题就被解决了。这是因为 Python 2.7 符合 deeptools 对 Python 版本的要求，并且 conda 会自动处理软件包之间的依赖关系，安装适配于 Python 2.7 的 deeptools 版本。

![1687000940528](image/Chip-seq_data-analysis/1687000940528.png)

这个似乎只显示了glibc 版本的不适配，直接更改可能会对系统的稳定性产生影响，因此需要谨慎操作。
这里似乎可以使用***[docker](https://www.docker.com/)***进行折腾，**挖个坑放假来填**

#### macs2&deeptools

针对python3.11不适配，和macs一样适配与python2.7，故安装在环境macs2中

#### 软件功能

![1687028485214](image/Chip-seq_data-analysis/1687028485214.png)![1687028495430](image/Chip-seq_data-analysis/1687028495430.png)

