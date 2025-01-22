rm(list = ls())
setwd("/mnt/user_group3/SUNJY/E_working/life_history/accumulation_curves/")
library(vegan)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)



#读取文件
data_all_t <- read.csv("/mnt/user_group3/SUNJY/E_working/life_history/accumulation_curves/new/t_15w_table_new.csv",header = T,row.names = 1)
colnames(data_all_t) <-data_all_t[1, ]
data_all_t <- data_all_t[-1, ]
know_clu <- read.csv("/mnt/user_group3/SUNJY/E_working/life_history/Uniprot/know_clu/know_clu_member_number.csv")[,c(2,3,4)]
data_all_t[data_all_t != 0] <- 1


# 设置循环次数和样本数量
set.seed(1015)
known_clusters <- know_clu$ref_clu
num_iterations <- 100  #重复100次随机抽样
sample_sizes <- seq(5, 354, by = 5)  # 样本数量从5到354，步长为5

# 初始化数据框
random_novel <- data.frame(sample_number = integer(0))


# 初始化存储所有抽样结果的列表
random_novel_list <- list()
mean_ratios <- numeric(length(sample_sizes))  # 存储每个样本数量的均值
all_unknown_ratios <- list()  # 用于存储每个样本数量的未知簇占比
all_unknown_counts <- list()  # 用于存储每次抽样的未知簇数量
all_total_counts <- list()  # 用于存储每次抽样的整体簇数量

# 循环进行每个样本数量的随机抽样
for (sample_number in sample_sizes) {
  unknown_ratios <- numeric(num_iterations)
  unknown_counts <- numeric(num_iterations)
  total_counts <- numeric(num_iterations)
  
  for (i in 1:num_iterations) {
    # 随机抽样 sample_number 个样本
    random_samples <- sample(colnames(data_all_t), sample_number)
    data_all_t_subset <- as.data.frame(lapply(data_all_t[, random_samples], as.numeric), row.names = rownames(data_all_t))
    
    # 合并为二进制矩阵（1 表示该簇至少出现一次）
    combined_data <- rowSums(data_all_t_subset)
    combined_data[combined_data > 0] <- 1
    
    # 计算未知簇占比
    total_unknown <- sum(!(rownames(data_all_t_subset) %in% known_clusters) & combined_data == 1)
    total_clusters <- sum(combined_data == 1)
    unknown_ratios[i] <- total_unknown / total_clusters
    
    # 记录每次抽样的未知簇数量和整体簇数量
    unknown_counts[i] <- total_unknown
    total_counts[i] <- total_clusters
  }
  
  # 存储每次抽样的未知簇占比及数量
  random_novel_list[[as.character(sample_number)]] <- data.frame(
    sample_number = rep(sample_number, num_iterations), 
    unknown_ratios = unknown_ratios,
    unknown_counts = unknown_counts,
    total_counts = total_counts
  )
  # 每次循环结束后合并数据
  random_novel1 <- do.call(rbind, random_novel_list)
  
  # 计算当前样本数量的均值
  mean_ratios[which(sample_sizes == sample_number)] <- mean(unknown_ratios)
  
  # 存储所有未知簇占比
  all_unknown_ratios[[as.character(sample_number)]] <- unknown_ratios
  all_unknown_counts[[as.character(sample_number)]] <- unknown_counts
  all_total_counts[[as.character(sample_number)]] <- total_counts
  write.csv(random_novel1,"/mnt/user_group3/SUNJY/E_working/life_history/accumulation_curves/new/random_novel_15w_clu_method2_1.csv", row.names = FALSE)
}

# 将所有抽样结果合并为一个数据框
random_novel <- do.call(rbind, random_novel_list)
random_novel$iteration <- as.numeric(gsub(".*\\.", "", rownames(random_novel)))


write.csv(random_novel,"/mnt/user_group3/SUNJY/E_working/life_history/accumulation_curves/new/random_novel_15w_clu_method2.csv", row.names = FALSE)

