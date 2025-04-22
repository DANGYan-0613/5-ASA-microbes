#观测值的相关矩阵
cor_sparcc <- read.delim('/Users/dangyan/Desktop/metagenomics/cor_sparcc.out_effective_baseline.txt', row.names = 1, sep = '\t', check.names = FALSE)
#伪 p 值矩阵
pvals <- read.delim('/Users/dangyan/Desktop/metagenomics/pvals.two_sided_effective_baseline.txt', row.names = 1, sep = '\t', check.names = FALSE)

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
neetwork_adj <- read.delim('/Users/dangyan/Desktop/metagenomics/neetwork.adj_baseline_effective.txt', row.names = 1, sep = '\t', check.names = FALSE)
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


# 找到 F.prau 及其相关节点
target_node <- "s__Faecalibacterium_prausnitzii"
neighbors <- neighbors(g, target_node)
nodes_to_label <- c(target_node, names(neighbors))

# 为特定节点及其邻居添加标签
V(g)$label <- ifelse(V(g)$name %in% nodes_to_label, V(g)$name, NA)

# 查看网络图
plot(g, vertex.label.cex = 0.8, vertex.label.dist = 0.5, vertex.label.color = "black",vertex.size = 0.5)
g
plot(g)


# 找到 F.prau 及其相关节点
target_node <- "s__Faecalibacterium_prausnitzii"
neighbors <- neighbors(g, target_node)
subgraph_nodes <- c(target_node, names(neighbors))

# 生成子图
subgraph <- induced_subgraph(g, subgraph_nodes)

# 为特定节点及其邻居添加标签
V(subgraph)$label <- V(subgraph)$name

# 查看子图
plot(subgraph, vertex.label.cex = 0.8, vertex.label.dist = 0.5, vertex.label.color = "black", vertex.size = 5)


# 读取相关性矩阵
cor_sparcc <- read.delim('/Users/dangyan/Desktop/metagenomics/cor_sparcc.out_effective_baseline.txt', row.names = 1, sep = '\t', check.names = FALSE)
# 读取P矩阵
pvals <- read.delim('/Users/dangyan/Desktop/metagenomics/pvals.two_sided_effective_baseline.txt', row.names = 1, sep = '\t', check.names = FALSE)
pvals[pvals>=0.01] <- -1   #防止下一步数据重合，目的是有效保留>=0.01的然后下下步设置为零
pvals[pvals<0.01 & pvals>=0] <- 1
pvals[pvals==-1] <- 0
#筛选后的邻接矩阵
adj <- as.matrix(cor_sparcc) * as.matrix(pvals)
#将相关矩阵中对角线中的值（代表了自相关）转为 0
diag(adj) <- 0    
# 将矩阵中非0的值变为1
adj[adj != 0] <- 1

# R语言下加载igraph包
library(igraph)
# convert adjacency matrix into a graph
sparcc.graph=graph.adjacency(adj,mode="undirected")
# display the graph
layout=layout.spring
plot(sparcc.graph, layout=layout)

#再转为其它类型的网络文件，例如
#再由 igraph 的邻接列表转换回邻接矩阵
adj_matrix <- as.matrix(get.adjacency(g, attr = 'sparcc'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(g, 'network_baseline_effective.graphml', format = 'graphml')
getwd()

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(g, 'network.gml', format = 'gml')

#边列表，也可以直接导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
edge <- data.frame(as_edgelist(g))

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  sparcc = E(g)$sparcc
)
head(edge_list)

write.table(edge_list, 'network.edge_list_effective_baseline.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
getwd()

#节点属性列表，对应边列表，记录节点属性，例如
node_list <- data.frame(
  nodes_id = V(g)$name,    #节点名称
  degree = degree(g)    #节点度
)
head(node_list)

write.table(node_list, 'network.node_list_effective_baseline.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
