# ATAC-Seq数据处理

## ATAC-seq数据peak calling

### fastq数据trim_galore软件质控

对双端测序的左右端结果文件进行处理

```shell
$ pwd
/home/yuanhx/dzx/ATAC_seq/clean
```

```shell
#shell脚本
bin_trim_galore="trim_galore"
ls ../raw/*gz | while read -r fq1 && read -r fq2;
do 
$bin_trim_galore -q 10 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o ./ $fq1 $fq2 
done
```

### bowtie2比对

用bowtie2进行比对和统计比对率,需要提前下载参考基因组然后使用命令构建索引，或者直接就下载索引文件：

这里使用下载好的参考基因组构建索引：

```shell
$ pwd
/home/ljx/yuanh/work_data/230926_Wheat_Ref_genome/index
```

```shell
bowtie2-build --threads 25 ../Wheat_Ref_genome.fasta ./Wheat
```

双端测序数据的比对：

```shell
cd /home/yuanhx/dzx/ATAC_seq/align

bin_bowtie2='/home/ljx/yuanh/bin/bowtie2'
bin_samtools='/home/ljx/yuanh/bin/samtools'
index="/home/ljx/yuanh/work_data/230926_Wheat_Ref_genome/index/Wheat"
ls ../clean/*gz |while read -r fq1 && read -r fq2;
do
sample=$(basename $fq1 | cut -d '_' -f1)
$bin_bowtie2 -p 10 -X 1000 -x $index -1 $fq1 -2 $fq2 |$bin_samtools sort -O bam -@ 20 -o ->${sample}.bam
done
```

```shell
#slurm提交脚本
#!/bin/bash
#SBATCH -J dzx
#SBATCH -p GPU-3090-1
#SBATCH -N 1
#SBATCH -o /home/yuanhx/dzx/ATAC_seq/align/out.txt
#SBATCH -e /home/yuanhx/dzx/ATAC_seq/align/err.txt

# 记录开始时间
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Job started at $start_time"

# 运行任务
bash run.sh

# 记录结束时间
end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Job finished at $end_time"

# 计算运行时间
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
runtime=$((end_seconds - start_seconds))
echo "Job took $runtime seconds to complete"

```

对bam文件过滤

```shell
mkdir  /home/yuanhx/dzx/ATAC_seq/mark_duplicate; cd $_
picard='/home/ljx/yuanh/bin/picard'
samtools='/home/ljx/yuanh/bin/samtools'
#去除PCR重复
$picard MarkDuplicates -I ../align/${sample}.bam \
			-O ${sample}.rmdup.bam -M ${sample}.rmdup.metric --REMOVE_DUPLICATES true

#去除低质量reads(-q 30)以及未必对到同一条染色体(-f 2)的数据
$samtools view -h -f 2 -q 30 ${sample}.rmdup.bam | grep -v chrM | $samtools sort -O bam -@ 20 -o - >${sample}.last.bam
bedtools bamtobed -i ${sample}.last.bam >${sample}.last.bed
```

## MACS2进行call peak

```shell
cd /home/yuanhx/dzx/ATAC_seq/mark_duplicate
ls *.bed | while read id;
do
macs2 callpeak -t $id -g 14300719022 --nomodel --shift -100 --extsize 200 -n ${id%%.*} --outdir ../peaks
done
```

后续脚本整合

`sbatch --dependency=afterok:21520 run.slurm`

run.slurm:

```shell
#!/bin/bash
#SBATCH -J dzx
#SBATCH -p GPU-3090-1
#SBATCH -N 1
#SBATCH -o /home/yuanhx/dzx/ATAC_seq/out.txt
#SBATCH -e /home/yuanhx/dzx/ATAC_seq/err.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=2719323380@qq.com

bash run.sh
```

```shell
#run.sh
mkdir  /home/yuanhx/dzx/ATAC_seq/mark_duplicate; cd $_
picard='/home/ljx/yuanh/bin/picard'
samtools='/home/ljx/yuanh/bin/samtools'
#去除PCR重复
$picard MarkDuplicates -I ../align/${sample}.bam \
			-O ${sample}.rmdup.bam -M ${sample}.rmdup.metric --REMOVE_DUPLICATES true

#去除低质量reads(-q 30)以及未必对到同一条染色体(-f 2)的数据
$samtools view -h -f 2 -q 30 ${sample}.rmdup.bam | grep -v chrM | $samtools sort -O bam -@ 20 -o - >${sample}.last.bam
bedtools bamtobed -i ${sample}.last.bam >${sample}.last.bed

cd /home/yuanhx/dzx/ATAC_seq/mark_duplicate
ls *.bed | while read id;
do
macs2 callpeak -t $id -g 14300719022 --nomodel --shift -100 --extsize 200 -n ${id%%.*} --outdir ../peaks
done
```

