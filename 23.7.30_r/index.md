# 23.7.30_R

## 运用代码

```R
    keywords <- c("calcium", "calmodulin")
    matched_indices <- sapply(keywords, function(keyword){
    str_detect(GO@result$Description, fixed(keyword, ignore_case = T))
    })
    matched_indices <- apply(matched_indices, 1, any)
    matched_terms <- GO@result[matched_indices,]
```

