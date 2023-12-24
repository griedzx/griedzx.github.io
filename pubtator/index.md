# PubTator

## API调用

```shell
#! /bin/bash
id_list='/home/yh/dzx/work/BioNLP/EntityAnnotation/pmid.txt'

#读取id时跳过pmid.txt第一行列名
while read line
do
    if [ $line != 'pmid' ]
    then
    curl  https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=$line >> /home/yh/dzx/work/BioNLP/EntityAnnotation/abstract_pubtator.txt
    echo >> /home/yh/dzx/work/BioNLP/EntityAnnotation/abstract_pubtator.txt
    fi
    sleep 2
done < $id_list
```

```shell
cat abstract_pubtator.txt|grep -v '|t|' |grep -v '|a|'|awk -F"[[:space:]]+" '{print $5}'|head
```

