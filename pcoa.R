library(Rtsne)
unifrac_dist <- read.table("/taxonomy_selected_aitchison.tsv")
# 使用 PCoA 进行降维分析
ordination <- cmdscale(unifrac_dist, k = 2, eig = TRUE)

# 创建一个数据框来存储 PCoA 的结果
pcoa_df <- data.frame(PC1 = ordination$points[, 1],
                      PC2 = ordination$points[, 2],
                      SampleID = rownames(ordination$points))

# 如果你有分组信息，添加到数据框中
# 假设 sample_data 是你的样本元数据，包含分组信息
sd <- read.table("/covariable.txt", sep="\t",header=TRUE)
pcoa_df <- merge(pcoa_df, sd, by.x = "SampleID", by.y = "sampleID")

library(ggplot2)

# 使用 ggplot2 绘制 PCoA 图
pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  labs(title = "PCoA of UniFrac distances",
       x = paste("PC1 - ", round(ordination$eig[1] / sum(ordination$eig) * 100, 1), "%", sep = ""),
       y = paste("PC2 - ", round(ordination$eig[2] / sum(ordination$eig) * 100, 1), "%", sep = "")) +
  theme_minimal()
# 如果没有安装 ggrepel 包，请先安装
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)
# 显示图形
print(pcoa_plot)
pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = SampleID), size = 3, box.padding = 0.5) +  # 添加样本号标签，并避免标签重叠
  labs(title = "PCoA of UniFrac distances",
       x = paste("PC1 - ", round(ordination$eig[1] / sum(ordination$eig) * 100, 1), "%", sep = ""),
       y = paste("PC2 - ", round(ordination$eig[2] / sum(ordination$eig) * 100, 1), "%", sep = "")) +
  theme_minimal()

# 显示图形
print(pcoa_plot)

# 使用 PCoA 进行降维分析
ordination <- cmdscale(unifrac_dist, k = 2, eig = TRUE)

# 创建一个数据框来存储 PCoA 的结果
pcoa_df <- data.frame(PC1 = ordination$points[, 1],
                      PC2 = ordination$points[, 2],
                      SampleID = rownames(ordination$points))

# 读取样本元数据，并合并分组信息
sd <- read.table("/covariable.txt", sep="\t", header = TRUE)
pcoa_df <- merge(pcoa_df, sd, by.x = "SampleID", by.y = "sampleID")
# 自定义颜色向量，假设分组 0 是蓝色，分组 1 是红色
group_colors <- c("0" = "#add9ee", "1" = "#ffd47f")

# 使用 ggplot2 绘制 PCoA 图，并添加星状放射线
pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = as.factor(group))) +
  geom_point(size = 3) +
  geom_segment(aes(xend = PC1.centroid, yend = PC2.centroid), linetype = "dashed", alpha = 0.5) +  # 添加星状放射线
  geom_point(aes(x = PC1.centroid, y = PC2.centroid), shape = 4, size = 5, color = "black") +  # 绘制质心
  scale_color_manual(values = group_colors, labels = c("Group 0", "Group 1")) +  # 设置分组颜色
  labs(title = "PCoA of UniFrac distances",
       x = paste("PC1 - ", round(ordination$eig[1] / sum(ordination$eig) * 100, 1), "%", sep = ""),
       y = paste("PC2 - ", round(ordination$eig[2] / sum(ordination$eig) * 100, 1), "%", sep = ""),
       color = "Group") +  # 设置图例标题
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # 添加黑色边框
  )

# 显示图形
print(pcoa_plot)

# 导出图形为 PDF 格式
ggsave("pcoa_plot_aitchison.pdf", plot = pcoa_plot, width = 8, height = 6)

# 计算每个组的质心
centroids <- aggregate(cbind(PC1, PC2) ~ group, data = pcoa_df, FUN = mean)

# 为每个样本点找到相应组的质心
pcoa_df <- merge(pcoa_df, centroids, by = "group", suffixes = c("", ".centroid"))

# 使用 ggplot2 绘制 PCoA 图，并添加星状放射线
library(ggplot2)

pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_segment(aes(xend = PC1.centroid, yend = PC2.centroid), linetype = "dashed", alpha = 0.5) +  # 添加星状放射线
  geom_point(aes(x = PC1.centroid, y = PC2.centroid), shape = 4, size = 5, color = "black") +  # 绘制质心
  labs(title = "PCoA of UniFrac distances",
       x = paste("PC1 - ", round(ordination$eig[1] / sum(ordination$eig) * 100, 1), "%", sep = ""),
       y = paste("PC2 - ", round(ordination$eig[2] / sum(ordination$eig) * 100, 1), "%", sep = "")) +
  theme_minimal()

# 显示图形
print(pcoa_plot)


library(vegan)

# 使用 NMDS 进行降维分析
nmds <- metaMDS(unifrac_dist, k = 2)

# 创建一个数据框来存储 NMDS 的结果
nmds_df <- data.frame(NMDS1 = nmds$points[, 1],
                      NMDS2 = nmds$points[, 2],
                      SampleID = rownames(nmds$points))

# 如果你有分组信息，添加到数据框中
nmds_df <- merge(nmds_df, sd, by.x = "SampleID", by.y = "sampleID")

# 使用 ggplot2 绘制 NMDS 图
nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = group)) +
  geom_point(size = 3) +
  labs(title = "NMDS of UniFrac distances") +
  theme_minimal()

# 计算每个组的质心
centroids <- nmds_df %>%
  group_by(group) %>%
  summarise(
    NMDS1_centroid = mean(NMDS1),
    NMDS2_centroid = mean(NMDS2)
  )

# 合并质心信息到原始数据框
nmds_df <- merge(nmds_df, centroids, by = "group")

# 绘制星状放射图
nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = group)) +
  # 绘制从质心到每个样本点的连线
  geom_segment(aes(x = NMDS1_centroid, y = NMDS2_centroid, xend = NMDS1, yend = NMDS2), 
               color = "grey80", linetype = "dashed") +
  # 绘制样本点
  geom_point(size = 3) +
  # 绘制质心
  geom_point(aes(x = NMDS1_centroid, y = NMDS2_centroid), size = 4, shape = 4, color = "black") +
  # 添加标题和其他图形美化
  labs(title = "NMDS of UniFrac distances with Star Plot from Centroids") +
  theme_minimal()

# 显示图形
print(nmds_plot)

# 显示图形
print(nmds_plot)

