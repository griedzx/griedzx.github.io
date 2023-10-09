# R语言笔记2

```
ifelse(test, yes, no)


all <- inter_data %>% 
  filter(gene1 %in% tf_data | gene2 %in% tf_data) %>% 
  distinct()

all$param_1 <- ifelse(all$gene1 %in% tf_data, 1, 0)
all$param_2 <- ifelse(all$gene2 %in% tf_data, 2, 0)
all$sum <- all$param_1 + all$param_2

inter_data <- all %>% 
  filter(sum != 0)

write.csv(all, "e:/1/TF_inter.csv", row.names = FALSE)
write.csv(inter_data, "e:/1/TF_inter.1.0.csv", row.names = FALSE)
```

```
library(dplyr)
inter_data <- read.csv("E:/1/ear_1.csv", header = FALSE, sep = ",") %>%
    select(1,5) %>%
    mutate(param_1 = 0, param_2 = 0, sum = 0) %>%
    rename(gene1 = V1, gene2 = V5)
TF <- read.csv("e:/1/TF_summary.csv", header = TRUE, sep = ",")
tf_data <- TF$TF_zong

all_1 <- inter_data %>% 
    filter(inter_data$gene1 %in% tf_data) 
all_2 <- inter_data %>% 
    filter(inter_data$gene2 %in% tf_data)

all <- rbind(all_1, all_2) %>% 
    distinct()
write.csv(all, "e:/1/TF_inter.csv", row.names = FALSE)

for(i in 1:nrow(inter_data)){
    if(inter_data$gene1[i] %in% tf_data){
        inter_data$param_1[i] <- 1
    }
    if(inter_data$gene2[i] %in% tf_data){
        inter_data$param_2[i] <- 2
    }
    inter_data$sum[i] <- inter_data$param_1[i] + inter_data$param_2[i]
}
inter_data <- inter_data %>% 
    filter(inter_data$sum != 0)
write.csv(inter_data, "e:/1/TF_inter.1.0.csv", row.names = FALSE)


#两列同时存在的
same <- inter_data %>% 
    filter(inter_data$gene1 %in% tf_data & inter_data$gene2 %in% tf_data)
write.csv(same, "e:/1/TF_same.csv", row.names = FALSE)
```

