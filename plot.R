library(ggplot2)

mydata<- read.csv("/maaslin2-pwy-baseline-coef.txt", header = T)


# 根据 Group 列对数据进行排序
mydata <- mydata[order(mydata$coef), ]

p <- ggplot(mydata, aes(x = coef, y = reorder(feature, coef, fill = coef)) +
  geom_bar(stat = 'identity') +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Linear model coefficient", y = "Pathway associated with Ineffective therapy (p<0.05)") +
  theme_minimal()

library(ggplot2)
  
  
# Plotting
p <- ggplot(mydata, aes(x = coef, y = reorder(feature, coef), fill = coef)) +
    geom_bar(stat = 'identity') +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Linear model coefficient", y = "Pathway associated with ineffective therapy (p<0.05)") +
    theme_minimal()
  
# Display the plot
print(p)

  
# 获取不同 Group 的唯一值，并按照需要的顺序排列
group_order <- unique(mydata$coef)

# 设置 x 轴顺序为按照 Group 排序
p <- ggplot(mydata, aes(factor(feature, levels = unique(mydata$coef)), fill = Group, label = Gene_number)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c('#FB8072', '#80B1D3', '#FDB462', '#8DD3C7', '#FFFFB3', '#BEBADA','#B3DE69', '#FCCDE5')) +
  coord_flip() +
  labs(title = 'KEGG annotation', y = 'Number of annotated genes', x = 'Function class', fill = 'Group') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  geom_text(position = position_dodge(width = 0.9), hjust = -0.3, size = 2)

# 显示图形
print(p)
