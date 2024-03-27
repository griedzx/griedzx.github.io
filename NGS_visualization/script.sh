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

#divide gene.bed into positive and negative
bed="./gene.bed"
awk '{if($5=="+") print $0}' $bed > gene_positive.bed
awk '{if($5=="-") print $0}' $bed > gene_negative.bed

#! /bin/bash
# Path: static/NGS_visualization/script

#sbatch --dependency=afterok:31482 ./run.slurm

#using computeMatrix to calculate the matrix
bam_files=$(find ./ -type l| grep "bam$")
bed_path="./bed"
for bam in $bam_files; do
    file_path=$(dirname "$(dirname "$bam")")
    mkdir -p $file_path/matrix
    omic=$(basename "${bam::-4}")
    file_name=${bam::-4}
    computeMatrix scale-regions -p 12 \
      -S ${file_name}_positive.bw \
      -R $bed_path/gene_positive.bed \
      -b 1000 -a 1000 \
      -o $file_path/matrix/${omic}_positive.gz \
      --binSize 10 --missingDataAsZero \
      --skipZeros \
      --outFileNameMatrix $file_path/matrix/${omic}_positive.tab
    computeMatrix scale-regions -p 12 \
      -S ${file_name}_negative.bw \
      -R $bed_path/gene_negative.bed \
      -b 1000 -a 1000 \
      -o $file_path/matrix/${omic}_negative.gz \
      --binSize 10 --missingDataAsZero \
      --skipZeros \
      --outFileNameMatrix $file_path/matrix/${omic}_negative.tab
done

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

# 找到所有匹配的文件，下载到本地
tab_files=$(find ./ -type f -name "*.tab" | grep -Pv "/.*(_negative|_positive)\.tab$")
tar -czf myfiles.tar.gz --transform 's:^./::' $tab_files

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