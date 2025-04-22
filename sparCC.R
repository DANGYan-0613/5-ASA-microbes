#观测值的相关矩阵
cor_sparcc <- read.delim('X.txt', row.names = 1, sep = '\t', check.names = FALSE)
#伪 p 值矩阵
pvals <- read.delim('X.txt', row.names = 1, sep = '\t', check.names = FALSE)

#保留 |相关性|≥0.8 且 p<0.01的值
cor_sparcc[abs(cor_sparcc) < 0.22] <- 0

pvals[pvals>0.08] <- -1
pvals[pvals<=0.08 & pvals>=0] <- 1
pvals[pvals==-1] <- 0

#筛选后的邻接矩阵
adj <- as.matrix(cor_sparcc) * as.matrix(pvals)
diag(adj) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
write.table(data.frame(adj, check.names = FALSE), 'neetwork.adj_baseline_effective.txt', col.names = NA, sep = '\t', quote = FALSE)

##网络格式转换
library(igraph)

#输入数据，邻接矩阵
neetwork_adj <- read.delim('X.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(neetwork_adj)[1:6]    #邻接矩阵类型的网络文件

#邻接矩阵 -> igraph 的邻接列表，获得含权的无向网络
g <- graph_from_adjacency_matrix(as.matrix(neetwork_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
g    #igraph 的邻接列表

#自相关也可以通过该式去除
g <- simplify(g)

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
E(g)$color <- ifelse(E(g)$correlation > 0, "red", "blue")

# 生成子图
subgraph <- induced_subgraph(g, subgraph_nodes)

# 为特定节点及其邻居添加标签
V(subgraph)$label <- V(subgraph)$name

# 查看子图
plot(subgraph, vertex.label.cex = 0.8, vertex.label.dist = 0.5, vertex.label.color = "black", vertex.size = 5)
