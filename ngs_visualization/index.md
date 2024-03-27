# NGS可视化

数据准备

创建合适的文件结构，使用软连接连接对应bam文件

```shell
#all_bam.txt
/home/data/yuanh/work_data/230612_encher_data/00.shoot_1/result/00.DNase/bam/shoot_DNase-seq.bam
/home/data/yuanh/work_data/230612_encher_data/00.shoot_1/result/01.ATAC_seq/bam/shoot_ATAC-seq.bam
/home/data/yuanh/work_data/230612_encher_data/00.shoot_1/result/02.Chip_seq/bam/H3K9ac_Chip-seq.bam
/home/data/yuanh/work_data/230612_encher_data/00.shoot_1/result/02.Chip_seq/bam/shoot_1_H3K27ac.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/00.shoot_1/result/02.Chip_seq/bam/shoot_1_H3K4me3.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/00.shoot_1/result/03.methylation/bam/output-prefix.bsmap.mkdup.bam
/home/data/yuanh/work_data/230612_encher_data/01.ear_1/result/01.ATAC_seq/bam/Zmays_ear_small.bam
/home/data/yuanh/work_data/230612_encher_data/01.ear_1/result/02.Chip_seq/bam/ear_1_H3K27ac.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/01.ear_1/result/02.Chip_seq/bam/ear_1_H3K4me3.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/02_ear_2/result/01.ATAC_seq/bam/ear_ATAC-seq.bam
/home/data/yuanh/work_data/230612_encher_data/02_ear_2/result/02.Chip_seq/bam/ear_2_H3K27me3.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/02_ear_2/result/02.Chip_seq/bam/ear_2_H3K4me3.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/02_ear_2/result/02.Chip_seq/bam/ear_2_H3K9ac.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/02_ear_2/result/03.methylation/bam/output-prefix.bsmap.mkdup.bam
/home/data/yuanh/work_data/230612_encher_data/03.tassel/result/01.ATAC_seq/bam/tassel_ATAC-seq.bam
/home/data/yuanh/work_data/230612_encher_data/03.tassel/result/02.Chip_seq/bam/tassel_H3K27me3.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/03.tassel/result/02.Chip_seq/bam/tassel_H3K4me3.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/03.tassel/result/02.Chip_seq/bam/tassel_H3K9ac.rmdup.bam
/home/data/yuanh/work_data/230612_encher_data/03.tassel/result/03.methylation/bam/output-prefix.bsmap.mkdup.bam
/home/data/yuanh/work_data/230612_encher_data/04.shoot_2/result/00.DNase/bam/DNase.bam
/home/data/yuanh/work_data/230612_encher_data/04.shoot_2/result/01.ATAC_seq/bam/Zmays_young_leaf_2c.bam
/home/data/yuanh/work_data/230612_encher_data/04.shoot_2/result/02.Chip_seq/bam/H3K9ac_Chip-seq.bam
/home/data/yuanh/work_data/230612_encher_data/04.shoot_2/result/02.Chip_seq/bam/SRR9016266_H3K27ac_ChIP-seq_rep1_picard_.sort.bam.bam
/home/data/yuanh/work_data/230612_encher_data/04.shoot_2/result/02.Chip_seq/bam/SRR9016267_H3K4me1_ChIP-seq_Rep1_picard_.sort.bam.bam
/home/data/yuanh/work_data/230612_encher_data/04.shoot_2/result/02.Chip_seq/bam/SRR9016268_RNA11_ChIP-seq_Rep1_picard_.sort.bam.bam
/home/data/yuanh/work_data/230612_encher_data/04.shoot_2/result/02.Chip_seq/bam/SRR9016269_H3k27me3_ChIP-seq_Rep1_picard_.sort.bam.bam
/home/data/yuanh/work_data/230612_encher_data/04.shoot_2/result/02.Chip_seq/bam/SRR9016270_H3K4me3_ChIP-seq_Rep1_picard_.sort.bam.bam
/home/data/yuanh/work_data/230612_encher_data/04.shoot_2/result/03.methylation/bam/shoot_2_methylation.sorted.bam
```

```shell
#! /bin/bash

while IFS= read -r line; do
  # tissue、omic、file_path
  file_path="$line"
  omic=$(echo "$file_path" | awk -F'/' '{gsub(/[0-9]+\.|_/, "", $(NF-2)); print $(NF-2)}')
  tissue=$(echo "$file_path" | awk -F'/' '{gsub(/[0-9]+\.|_/, "", $(NF-4)); print $(NF-4)}')

  # mkdir
  target_dir="./$tissue/$omic/bam"
  echo "$target_dir"
  mkdir -p "$target_dir"

  # ln -s
  if [ "$omic" == "Chip_seq" ]; then
    # Chip_seq
    link_name=$(basename "$file_path" | grep -oP 'H[123]K[1-9]*(me[123]|ac)|RNA11')
  else
    #omic
    link_name="$omic.bam"
  fi
  ln -s "$file_path" "$target_dir/$link_name"
done < all_bam.txt
```

递归查询有多少软链接建立，和txt中文件路径相符

```shell
(base) yuanhx@T640server 23:08:13 ~/dzx/ultimate/1.input_circumstance/data
$ find ./ -type l |wc -l
28
```

**deeptools**去处理链特异性的数据

**我们要分析测序数据在基因上的分布特征，显然不能将来自于负链基因数据统计到同一区域的正链基因上。** 一个简单的思路就是：

* 首先将BAM文件根据链拆分成正链和负链的BAM文件；
* 再用正链的BAM文件去和正链基因的BED文件去 `computeMatrix`，负链同理；
* 最后再将正负链 `computeMatrix`的结果进行合并。

bam文件区分正负链

```shell
# 选择比对到正链的reads
samtools view -F 16 -b -o positive.bam example.bam

# 选择比对到负链的reads
samtools view -f 16 -b -o negative.bam example.bam 
```

"-f"参数用于选择具有特定标志的读取，"-F"参数用于排除具有特定标志的读取

SAM/BAM文件中，每个读取都有一个FLAG字段，这个字段是一个位标志，用于表示读取的各种属性，这些属性包括读取是否配对，是否比对到参考序列，是否是配对中的第一读取，等等

对于链特异性，我们关注的是第5位，如果这一位是1，那么读取比对到参考序列的负链，如果是0，那么读取比对到参考序列的正链。这一位的值是16，所以使用16作为参数

```shell
#divide the bam files into positive and negative
bam_files=$(find ./ -type l)
for bam in $bam_files; do
    #remove the .bam which is the last 4 characters
    file_name=${bam::-4}
    samtools view -@ 9 -f 16 -b -o ${file_name}_negative.bam $bam
    samtools view -@ 9 -F 16 -b -o ${file_name}_positive.bam $bam
    samtools index ${file_name}_negative.bam ${file_name}_negative.bam.bai
    samtools index ${file_name}_positive.bam ${file_name}_positive.bam.bai
    bamCoverage -b ${file_name}_negative.bam -o ${file_name}_negative.bw --binSize 10 --normalizeUsing RPKM --scaleFactor 1 --numberOfProcessors 10
    bamCoverage -b ${file_name}_positive.bam -o ${file_name}_positive.bw --binSize 10 --normalizeUsing RPKM --scaleFactor 1 --numberOfProcessors 10
done
```

bed基因文件区分正负链

```shell
#divide gene.bed into positive and negative
bed="./gene.bed"
awk '{if($5=="+") print $0}' $bed > gene_positive.bed
awk '{if($5=="-") print $0}' $bed > gene_negative.bed
```

分正负链基因计算reads富集分数

```shell
#! /bin/bash
# Path: static/NGS_visualization/script

#sbatch --dependency=afterok:31482 ./run.slurm

#using computeMatrix to calculate the matrix
bam_files=$(find ./ -type l| grep "bam$")
for bam in $bam_files; do
    file_path=${dirname(dirname($bam))}
    mkdir -p $file_path/matrix
    omic=${basename(bam::-4)}
    file_name=${bam::-4}
    computeMatrix scale-regions -p 12 \
      -S ${file_name}_positive.bw \
      -R gene_positive.bed \
      -b 1000 -a 1000 \
      -o $file_path/matrix/${omic}_positive.gz \
      --binSize 10 --missingDataAsZero \
      --skipZeros \
      --outFileNameMatrix $file_path/matrix/${omic}_positive.tab
    computeMatrix scale-regions -p 12 \
      -S ${file_name}_negative.bw \
      -R gene_negative.bed \
      -b 1000 -a 1000 \
      -o $file_path/matrix/${omic}_negative.gz
      --binSize 10 --missingDataAsZero \
      --skipZeros \
      --outFileNameMatrix $file_path/matrix/${omic}_negative.tab
done
```

合并正负链tab文件

```shell
#合并正负链tab文件
bam_files=$(find ./ -type l| grep "bam$")
for bam in $bam_files; do
    file_path=$(dirname "$(dirname "$bam")")
    omic=$(basename "${bam::-4}")
    #根据两个文件表头确定新表头
    #first line
    group_boundaries_positive=$(zcat $file_path/matrix/${omic}_positive.gz | head -n 1 | grep -oP '"group_boundaries":\[\K[^\]]*')
    group_boundaries_negative=$(zcat $file_path/matrix/${omic}_negative.gz | head -n 1 | grep -oP '"group_boundaries":\[\K[^\]]*')
    # 计算新的 group_boundaries 的值
    num1=$(echo $group_boundaries_positive | awk -F ',' '{print $2}')
    num2=$(echo $group_boundaries_negative | awk -F ',' '{print $2}')
    new_group_boundaries="0,$((num1+num2))"
    # 从 positive 文件中提取 JSON 数据，并将 group_boundaries 和 sample_labels 的值进行修改
    new_json=$(zcat $file_path/matrix/${omic}_positive.gz | head -n 1 | sed "s/\"group_boundaries\":\[[^]]*\]/\"group_boundaries\":[$new_group_boundaries]/g" | sed "s/\"sample_labels\":\[[^]]*\]/\"sample_labels\":[\"$omic\"]/g")
    # 输出新的 JSON 数据
    echo $new_json > $file_path/matrix/${omic}.tab

    #每个gz文件去除第一行行合并
    zcat $file_path/matrix/${omic}_positive.gz | tail -n +2 >> $file_path/matrix/${omic}.tab
    zcat $file_path/matrix/${omic}_negative.gz | tail -n +2 >> $file_path/matrix/${omic}.tab
    #gzip the tab file
    gzip -c $file_path/matrix/${omic}.tab > $file_path/matrix/${omic}.tab.gz
done
```

画图

```shell
#plotHeatmap
#! /bin/bash

bam_files=$(find ./ -type l| grep "bam$")
for bam in $bam_files; do
    file_path=$(dirname "$(dirname "$bam")")
    omic=$(basename "${bam::-4}")
    mkdir -p $file_path/plotheatmap
    plotHeatmap -m $file_path/matrix/${omic}.tab.gz \
      -out $file_path/plotheatmap/${omic}.pdf \
      --dpi 300 \
      --colorMap YlGnBu  \
      --missingDataColor "#FFF6EB"\
      --heatmapHeight 21 \
      --startLabel "TSS" \
      --endLabel "TTS" \
      --regionsLabel "gene" \
      --legendLocation "none" \
      --plotTitle $omic 
done
```

