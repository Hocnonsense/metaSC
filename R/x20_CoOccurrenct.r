###
#* @Date: 2021-01-20 13:59:22
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-02-03 11:55:19
#* @FilePath: /Analyze/2021/01/x20_CoOccurrenct.r
#* @Description:
###
rm(list = ls())

library(vegan)
library(corrplot)
library(Hmisc)
library(igraph)
library(psych)

file="Analyze/level-6-small.csv"
data0=read.csv(file,header = T,row.names = 1)
data1=data0[,-which(names(data0)%in%c("name"))]
data2=as.data.frame(
    apply(data1, 2,
          function(x){return(x/rowSums(data1))}
    )
)
data3=as.data.frame(t(data2))
data=data3
rm(file,data0,data1,data2,data3)

# 分析 http://www.biotrainee.com/thread-542-1-1.html
occor = corr.test(t(data),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occor.r = occor$r  # 相关系数
occor.p = occor$p  # 显著性水平
occor.r[occor.p>=0.01|abs(occor.r)<=0.6] = 0
write.csv(occor.r,"/home/hwrn/Work/2021_01-AOA16S/Analyze/level-6-small-occor.csv")

# my_igraph
my_igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs = V(my_igraph)[degree(my_igraph) == 0]
my_igraph = delete.vertices(my_igraph, bad.vs)
my_igraph.weight = E(my_igraph)$weight
my_igraph.color=ifelse(
    my_igraph.weight>0,
    "red",ifelse(my_igraph.weight<0,
                 "blue","gray"))

comps = membership(cluster_fast_greedy(my_igraph,weights =NULL))
# 简单画图
colbar = rainbow(max(comps))
V(my_igraph)$color = colbar[comps]

set.seed(599)
p=plot(my_igraph,
       main="Co-occurrence network",
       vertex.frame.color=NA,
       edge.width=1, vertex.size=5, edge.lty=1,
       #vertex.size=8, edge.lty=2,
       edge.curved=TRUE,
       margin=c(0,0,0,0),
       vertex.label=NA,#V(my_igraph)$name,
       #vertex.label.size=0.01,
       edge.color=my_igraph.color,
       layout=c(
           layout_with_graphopt,
           layout_nicely,
           layout_randomly,
           layout_with_dh,
           layout_with_fr,
           layout_with_gem,
           layout_with_kk,
           layout_with_lgl,
           layout_with_mds,
           layout_with_sugiyama
       )[[2]]
)

# 选择 GraphML 格式导出
write_graph(my_igraph,
            "Analyze/level-6-tp.graphml","graphml")

quit()


# 统计网络属性
summerize_graph <- function(i,an_igraph){
    cc=transitivity(my_igraph)
    apl=average.path.length(my_igraph)
    ad=mean(igraph::degree(my_igraph))
    gd=edge_density(my_igraph,loops = FALSE)
    diameter=diameter(my_igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
    comps = membership(cluster_fast_greedy(my_igraph,weights =NULL))
    md = modularity(my_igraph,comps)
    return(cbind(cc,apl,ad,gd,diameter,md))
}

# 可选
c_score=oecosimu(data,nestedchecker,method='swap',nsimul = 10000)

summerize_graph(my_igraph)

#计算多少个边和多少个点
e_num=length(E(my_igraph)) #430
v_num=length(V(my_igraph)) #96
#ER随机网络，生成文件，节点和节点的频数，随机网络的属性
df_out=matrix(NA,0,7)
for (i in 1:10000){
    g=erdos.renyi.game(v_num,e_num,'gnm',weight=T,mode="undirected")
    g=simplify(g)
    write.graph(g,"random.txt",format = 'edgelist')

    df_out=rbind(df_out,cbind(i,summerize_graph(g)))
}
write.csv(as.data.frame(df_out), "Analyze/ER_random_network_attributes.csv")

d <- degree(my_igraph, mode="in")
fit <- fit_power_law(d, impelementation = "R.mle")
##　这里可以作出拟合的图出来
curve(x^fit$alpha)
