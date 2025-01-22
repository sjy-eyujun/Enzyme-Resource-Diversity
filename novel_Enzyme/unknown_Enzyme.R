rm(list=ls())
setwd("/mnt/user_group3/SUNJY/E_working/life_history/Uniprot/diamond")
library(dplyr)
diamond <- read.table("/mnt/user_group3/SUNJY/E_working/life_history/Uniprot/diamond/diamond_2.7kw_27.txt", 
                      sep = "\t",
                      quote = "\"", fill = TRUE, comment.char = "", na.strings = c("NA", ""))
colnames(diamond) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

diamond_uni <- diamond %>%
  group_by(qseqid) %>%
  slice(which.max(pident)) %>%
  ungroup()


#clean <- read.csv("0.2.csv")
#uniprot_ec <- read.table("/mnt/user_group3/SUNJY/E_working/life_history/Uniprot/uniprot.tsv", sep="\t", header=TRUE, fill=TRUE, quote="")

#merge1 <- merge(diamond_uni,clean,by.x = "qseqid",by.y="ref_clu")
#merge2 <- merge(merge1,uniprot_ec,by.x="sseqid",by.y = "Entry")

ident_80 <- diamond_uni[diamond_uni$pident >= 80, ]



clean_ref_allclu <- read.csv("/mnt/user_group3/SUNJY/E_working/life_history/Uniprot/diamond/2.7kw_ref_clu.csv")
ident_80_name <- ident_80[,c(1)]
merge1 <- merge(ident_80,clean_ref_allclu,by.x="qseqid",by.y="clu_member")

write.csv(merge1,"ident_2.7kw_80.csv")

merge1 <- read.csv("ident_2.7kw_80.csv")
ref_member_number <- merge1 %>%
  group_by(ref_clu) %>%
  summarise(member_count = n())


write.csv(ref_member_number,"ident_2.7kw_80_ref_number.csv")




ref_15 <- read.csv("~/E_working/life_history/Uniprot/diamond/ident80_ref_allclu.txt", header=FALSE)
ref_15 <- ref_15[,c(2,3)]
colnames(ref_15) <- c("clu_member","ref_clu")

uni_15 <- as.data.frame(unique(ref_15$ref_clu))
colnames(uni_15) <- c("ref_clu")
colnames(uni_ref) <- c("ref_clu")

in2.7kw <- as.data.frame(uni_ref[!(uni_ref$ref_clu %in% uni_15$ref_clu), ])
colnames(in2.7kw) <- c("ref_clu")
in2.7kw <- merge(in2.7kw,ref_member_number,by.x="ref_clu",by.y="ref_clu")



diamond_15w <- read.csv("/mnt/user_group3/SUNJY/E_working/life_history/Uniprot/diamond/diamond_uni.csv")
in2.7kw <- merge(in2.7kw,diamond_15w,by.x="ref_clu",by.y="qseqid")

clu_2.7 <- read.csv("2.7kw_ref_clu.csv")

clu_2.7_member_number <- clu_2.7 %>%
  group_by(ref_clu) %>%
  summarise(member_count = n())
colnames(clu_2.7_member_number)[2] <- "clu_all_member"

clu_all_number <- merge(ref_member_number,clu_2.7_member_number,by.x="ref_clu",by.y="ref_clu")
total_sum <- sum(clu_all_number$clu_all_member, na.rm = TRUE)



clean <- read.csv("0.2.csv")
clean <- clean[,c(1,2,3)]
colnames(clean)[3] <- "clean_ec"



uniprot <- read.table("/mnt/user_group3/SUNJY/E_working/life_history/Uniprot/uniprot.tsv", 
                      sep = "\t", 
                      header = TRUE, 
                      fill = TRUE, 
                      quote = "",
                      comment.char = "")
uniprot <- uniprot[,c(1,8)]
colnames(uniprot)[2] <- "uniprot_ec"


clu <- merge1[,c(2,3,4,15)]

clu_clean <- merge(clu,clean,by.x="ref_clu",by.y="ref_clu")
clu_clean_uniprot <- merge(clu_clean,uniprot,by.x="sseqid",by.y="Entry")


compare_ec <- function(clean_ec, uniprot_ec) {
  uniprot_list <- strsplit(uniprot_ec, split = ",")[[1]]
  return(clean_ec %in% uniprot_list)
}
clu_clean_uniprot$if_same <- apply(clu_clean_uniprot, 1, function(row) {
  compare_ec(row["clean_ec"], row["uniprot_ec"])
})


if_same <- clu_clean_uniprot %>%
  group_by(ref_clu) %>%
  summarize(if_same = any(if_same))


num_true <- sum(if_same$if_same)

num_true 


length <- read.csv("2.7kw_length.csv",header = F)
colnames(length)[1] <- "clu_name"
colnames(length)[2] <- "clu_length"
clu_clean_uni_length <- merge(clu_clean_uniprot,length,by.x="qseqid",by.y="clu_name")
clu_clean_uni_length <- merge(clu_clean_uni_length,merge1[,c(2,5)],by.x = "qseqid",by.y = "qseqid")

clu_clean_uni_length$coverage <- clu_clean_uni_length$length / clu_clean_uni_length$clu_length


know_clu <- clu_clean_uni_length[clu_clean_uni_length$coverage >= 0.8, ]
colnames(know_clu)[8] <- "same"

if_same2 <- know_clu %>%
  group_by(ref_clu) %>%
  summarize(same = any(same))


num_true2 <- sum(if_same2$same)

num_true2 


same_member <- if_same2[if_same2$same == TRUE, ]
same_member <- merge(same_member,clu_2.7_member_number,by.x="ref_clu",by.y="ref_clu")

total_sum <- sum(same_member$clu_all_member, na.rm = TRUE)
total_sum


write.csv(same_member,"/mnt/user_group3/SUNJY/E_working/life_history/Uniprot/know_clu/know_clu_member_number.csv")



