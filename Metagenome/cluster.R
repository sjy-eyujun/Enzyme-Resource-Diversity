rm(list = ls())
library(arrow)
library(dplyr)
library(tidyr)

dat <- read.table(file = '/mnt/SJY/E_working/life_history/mmseq2/354_noanno_clu.tsv',header = F)


superfamily_list <- dat %>%
  group_by(V1) %>%
  summarize(superfamily_list = list(V2),num_merged = n())

superfamily_list$superfamily_list <- sapply(superfamily_list$superfamily_list,function(x) paste(x, collapse = ', '))

dat1 <- merge(dat,superfamily_list,by = 'V1',all = T,sort = F)

dat1 <- dat1 %>% distinct(V1, .keep_all = T)



load(file = "gmgc_cluster.Rdata")

index3 <- which(dat1$num_merged>=3)

dat3 <- dat1[index3,]

cluster3 <- data.frame(dat3$V1)

index10 <- which(dat1$num_merged>=10)

dat10 <- dat1[index10,]

cluster10 <- data.frame(dat10$V1)


index100 <- which(dat1$num_merged>=100)

dat100 <- dat1[index100,]

cluster100 <- data.frame(dat100$V1)

write.table(cluster3,file = '/mnt/SJY/E_working/life_history/mmseq2/354_noanno_cluster3.txt',row.names = F,col.names = F,quote = FALSE)
write.table(cluster10,file = '/mnt/SJY/E_working/life_history/mmseq2/354_noanno_cluster10.txt',row.names = F,col.names = F,quote = FALSE)
write.table(cluster100,file = '/mnt/SJY/E_working/life_history/mmseq2/354_noanno_cluster100.txt',row.names = F,col.names = F,quote = FALSE)


dat_uni <- unique(dat$V1)
write.table(dat_uni,'/mnt/SJY/E_working/life_history/mmseq2/354_noanno_uni.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
