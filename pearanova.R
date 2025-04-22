library(vegan)
library(ggpubr)
library(reshape2)
library(ggsci)
library ("phyloseq")
library("ggplot2")
library(permute)
library(lattice)
library(vegan) 
library(tidyverse)  
library(ggpubr)
library(dplyr)
options(stringsAsFactors=F)

df=read.delim("/Users/dangyan/Desktop/metagenomics/filtered_data_norm_beta.txt", sep="\t", header=TRUE)
# 统计每列非零值的数量
non_zero_counts <- colSums(df != 0)

# 去掉非零值小于6的列
df_filtered <- df[, non_zero_counts >= 5]

# 计算每列（除了第一列）的总和
column_sums <- colSums(df[,-1])

# 保留第一列，并去掉总和小于10的列
df_filtered <- df[, c(TRUE, column_sums >= 10)]

df<-df_filtered


sd <- read.table("/Users/dangyan/Desktop/covariable.txt", sep="\t",header=TRUE)


data<-inner_join(df, sd, by="sampleID")

site = data.frame(sample = data$sampleID, group = data$group, age =data$age,
                 gender = data$gender)
colnames(site)

row.names(data)<-data[,1]
data=data[,-1]
data<-subset(data,select=-c(group,gender,age, MES))
site$group <- as.factor(site$group)
site$gender <- as.factor(site$gender)
site$mes <- as.factor(site$mes)
#dataT<-as.data.frame(t(data))
dist <- vegdist(data, method = "bray")
dist_matrix <- as.matrix(dist)
# 检查对角线元素是否为 0
diag(dist_matrix)

adonis_result_dis = adonis2(dist~group, site, permutations = 9999)
adonis_result <- adonis2(dist ~ group * age * gender, data = site, permutations = 9999)

adonis_result_dis

# 计算组内离散度
betadisper_result <- betadisper(dist, group)

# 进行置换检验
permutest_result <- permutest(betadisper_result, permutations = 9999)

# 打印结果
print(permutest_result)


grouping_vector <- sd$group
anosim_result <- anosim(dist, grouping_vector, permutations = 9999)
print(anosim_result)

library(vegan)

#计算香农指数
shannon <- diversity(df, index = "shannon")

# 计算辛普森指数
simpson <- diversity(df, index = "simpson")
chao1_values <- estimateR(df)["S.chao1", ]

# 打印 Chao1 指数
print(chao1_values)

# 计算观测物种数（OTU数）
observed <- specnumber(df)

# 将所有的Alpha多样性指数合并为一个数据框
alpha_diversity <- data.frame(SampleID = rownames(df),
                              Shannon = shannon,
                              Simpson = simpson,
                              Observed = observed)

# 打印结果
print(alpha_diversity)

# 假设 alpha_diversity 是你要保存的数据框

# 使用 write.table 将数据写入文件，制表符作为分隔符
write.table(alpha_diversity, file = "alpha_diversity_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
getwd()
