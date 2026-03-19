setwd("D:/老师的生信处理/蔡老师文件夹/ljy/动物")

library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library ("org.Mm.eg.db") #人是Hs；小鼠 Mm；大鼠 Rn；
library(patchwork)
library(WGCNA)
library(GSEABase)
library(GSVA)

exp <- read.table(file = "matrix.txt",header = T, quote="\"")
#exp <-read.csv("GSE275110_raw_counts.csv",row.names = 1)
exp <- as.matrix(exp)
rownames(exp)<-exp[,1]
exp=exp[,2:ncol(exp)]
mode(exp)
dimnames=list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames) ##变成数字形式
exp <- avereps(exp) ##重复的基因名 平均化
exp[exp==0]<-NA
exp=na.omit(exp) 
exp1<-log2(exp+1)
range(exp1)
boxplot(exp1,col=c("skyblue","skyblue","skyblue","skyblue","skyblue",
                  "lightgreen","lightgreen","lightgreen","lightgreen","lightgreen",
                  "orange","orange","orange","orange","orange"),
        las = 2)
exp=normalizeBetweenArrays(exp)
###PCA
dat=as.data.frame(t(exp1))
#BiocManager::install("factoextra")#按照包
library(FactoMineR)#画主成分分析图需要加载这两个包
library(factoextra)

# pca的统一操作走起
group=c(rep("A",5),rep("B",5),rep("C",5))
dat.pca <- PCA(dat, graph = FALSE)
??fviz_pca_ind()
pca_plot <- fviz_pca_ind(dat.pca,
                         
                         geom.ind = c("point"), # show points only (nbut not "text")
                         
                         col.ind = group, # color by groups
                         
                         palette = c("skyblue","lightgreen","orange"),
                         
                         addEllipses = T, # Concentration ellipses
                         
                         legend.title = "Groups",
)+theme(panel.grid=element_blank())+theme(panel.border = element_rect(fill=NA,color="black", 
                                                                      size=1, linetype="solid"));pca_plot

pca_plot

##热图
library(pheatmap)

# 提取 voom 转换后的表达量（log2CPM）
exp_aim <- read.table(file = "Heatmap2.txt",header = T, quote="\"")
exp_aim <- as.matrix(exp_aim)
rownames(exp_aim)<-exp_aim[,1]
exp_aim=exp_aim[,2:ncol(exp_aim)]
mode(exp_aim)
dimnames=list(rownames(exp_aim),colnames(exp_aim))
exp_aim<- matrix(as.numeric(as.matrix(exp_aim)),nrow=nrow(exp_aim),dimnames=dimnames) ##变成数字形式
exp_aim <- avereps(exp_aim) ##重复的基因名 平均化
exp_aim[exp_aim==0]<-NA
exp_aim=na.omit(exp_aim) 
exp_aim<-log2(exp_aim+1)
# z-score 标准化（按基因 row）
heat_expr <- t(scale(t(exp_aim)))  # 每行 z-score

# 可选：设置分组颜色
ann_col <- data.frame(Group = group)
rownames(ann_col) <- colnames(heat_expr)

# 画图
pheatmap(heat_expr,
         annotation_col = ann_col,
         color = colorRampPalette(c("#016DDB", "white", "#E43401"))(100),
         cluster_rows = TRUE,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "none",
         fontsize_row = 8,
         fontsize_col = 10,
         main = "")

#co=cor(exp)；pheatmap(co)这两句代码可做相关性热图

#rm(list = ls())

##差异分析
rownames(design)
colnames(exp)
expAB<-exp[,c(1:10)]
expBC<-exp[,c(6:15)]
group <- factor(c(rep("B", 5), rep("C", 5)), levels = c("B", "C"))

design=model.matrix(~group)

v <- voom(expBC, design, plot = TRUE)

fit=lmFit(v,design)

fit=eBayes(fit)

deg=topTable(fit,coef=2,number = Inf)

#为deg数据框添加几列

#1.加probe_id列，把行名变成一列

library(dplyr)

deg <- mutate(deg,SYMBOL=rownames(deg))

#或者 deg$probe_id=rownames(deg)

head(deg)


#3.加change列,标记上下调基因
#######
logFC_t=1
P.Value_t = 0.05

k1 = (deg$adj.P.Val< P.Value_t)&(deg$logFC < -logFC_t)#可以用table(k1)看看具体有几个
table(k1)
k2 = (deg$adj.P.Val < P.Value_t)&(deg$logFC > logFC_t)
table(k2)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))#或者用case_when()
table(change)
deg <- mutate(deg,change)
write.csv(deg, file = "DEGsAB_change.csv", row.names = FALSE)

mycolor<- c( "#016DDB","gray80","#E43401")
ggplot(
  deg, 
  aes(x = logFC, 
      y = -log10(adj.P.Val), 
      colour=change)) +
  geom_point(alpha=0.65, size=3) +
  scale_color_manual(values=mycolor)+
  
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=1.0) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=1.0) +
  
  labs(x="log2(Fold Change)",
       y="-log10(FDR)")+
  theme_bw()+
  
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )+theme(panel.grid=element_blank())
dev.off()

#######GSEA富集分析
fix(deg)
gene.df <- bitr(deg$SYMBOL,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Mm.eg.db) 
Enrich <- gene.df %>%
  inner_join(deg, by = "SYMBOL")
head(Enrich)
Enrichment <- Enrich %>%
  as.data.frame()%>%
  dplyr::select(ENTREZID,logFC)
Enrichment <- Enrichment[!duplicated(Enrichment[,1]),]
geneList <- Enrichment$logFC
names(geneList) <- Enrichment$ENTREZID
geneList <-sort(geneList,decreasing = T)
###GSEA_GO
GSEA_GO<- gseGO(geneList = geneList,
                        ont = "BP",   # 可以选择 "BP"（生物过程）、"CC"（细胞组分）或 "MF"（分子功能）
                        OrgDb = org.Mm.eg.db,  # 大鼠的基因注释数据库
                        keyType = "ENTREZID", # 基因ID类型
                        nPerm = 10000,        # 自定义置换次数，越大结果越稳健
                        pvalueCutoff = 1,   # P值显著性阈值
                        pAdjustMethod = "BH",  # p值多重比较校正方法
                        verbose = TRUE)
write.csv(GSEA_GO,file="GSEA_GOBP_AB.csv",quote=F,row.names = F) 

GoCategory=c("positive regulation of cytosolic calcium ion concentration",
"calcium ion transport",
"regulation of calcium ion transport",
"response to calcium ion",
"calcium ion homeostasis",
"intracellular calcium ion homeostasis",
"negative regulation of cytosolic calcium ion concentration",
"positive regulation of calcium ion import",
"regulation of calcium ion import",
"calcium ion transmembrane import into cytosol"
)
ridgeplot(GSEA_GO,
          showCategory = GoCategory,##如果写数字
          fill = "qvalue", #填充色 "pvalue", "p.adjust", "qvalue" 
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "NES",
          decreasing = F
)+theme(axis.text.y = element_text(size=8))+theme(panel.grid=element_blank())

###GSEA_KEGG
options(timeout = 3000)
GSEA1 <- gseKEGG(geneList     = geneList,
                 organism     = 'mmu',
                 nPerm        = 10000,
                 minGSSize    = 1,
                 maxGSSize    = 10000,
                 pvalueCutoff = 1,
                 pAdjustMethod = "none" ) # hsa人 mmu鼠
write.csv(GSEA1,file="GSEA_BC.csv",quote=F,row.names = F) #
###山峦图展示
library(enrichplot)
library(ggplot2)

Category=c("TNF signaling pathway",
           "NF-kappa B signaling pathway",
           "Neutrophil extracellular trap formation",
           "IL-17 signaling pathway",
           "Toll-like receptor signaling pathway",
           "JAK-STAT signaling pathway",
           "Cytokine-cytokine receptor interaction",
           "MAPK signaling pathway",
           "Chemical carcinogenesis - reactive oxygen species",
           "NOD-like receptor signaling pathway"
)

##选择展示部分
#Category=c("TNF signaling pathway - Mus musculus (house mouse)",
#           "Bacterial invasion of epithelial cells",
#           "PPAR signaling pathway",
#           "ECM-receptor interaction",
#           "Focal adhesion")

ridgeplot(GSEA1,
          showCategory = Category,##如果写数字
          fill = "qvalue", #填充色 "pvalue", "p.adjust", "qvalue" 
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "NES",
          decreasing = T
)+theme(axis.text.y = element_text(size=10))+theme(panel.grid=element_blank())

## Picking joint bandwidth of 0.212
??ridgeplot()
p1=gseaplot2(GSEA_GO, geneSetID = "GO:0034599", title = "cellular response to oxidative stress", color = "black");p1
p1=gseaplot2(GSEA1, geneSetID = "mmu00190", title = "OXPHS Signaling Pathway", color = "black");p1

###GO and KEGG
deg2 = deg %>% 
  filter( abs( logFC ) > 1 & adj.P.Val < 0.05 )# &：和；|：或。
  #%>% # 按 logFC 绝对值筛选
# filter( avg_log2FC > 0 ) # 按 logFC 正负筛选
gene.df <- bitr(deg2$SYMBOL,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Mm.eg.db)

gene <- gene.df$ENTREZID
go<- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05,ont="all",readable =T)

GOcategoryBP<-c("negative regulation of neurogenesis",
              "neural precursor cell proliferation",
              "synapse assembly",
              "axon guidance",
              "learning",
              "cognition",
              "memory",
              "memory",
              "regulation of synapse structure or activity",
              "regulation of neurogenesis"
)

GOcategoryCC<-c("presynaptic membrane",
"postsynaptic membrane",
"receptor complex",
"presynaptic active zone",
"GABA-ergic synapse",
"presynaptic active zone membrane",
"cytosolic ribosome",
"postsynaptic specialization",
"transmembrane transporter complex",
"monoatomic ion channel complex")

GOcategoryMF<-c("channel activity",
"passive transmembrane transporter activity",
"monoatomic ion channel activity",
"monoatomic ion-gated channel activity",
"gated channel activity",
"monoatomic cation channel activity",
"metal ion transmembrane transporter activity",
"ligand-gated monoatomic ion channel activity",
"ligand-gated channel activity",
"voltage-gated monoatomic ion channel activity")

write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F) #
write.csv(go,file="GO.csv",quote=F,row.names = F) #
##可视化
##条形图
pdf(file="GO-barplot.pdf",width = 10,height = 15)
barplot(go, drop = TRUE, showCategory =10,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')+theme(panel.grid=element_blank())
dev.off()

##气泡图
pdf(file="GO-bubble.pdf",width = 10,height = 15)
dotplot(go,showCategory = 5,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')+theme(panel.grid=element_blank())
dev.off()

#kegg分析
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism = "mmu", pvalueCutoff =1, qvalueCutoff =1, pAdjustMethod = "fdr")   
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                         
write.csv(kk,file="KEGG_AB.csv",quote=F,row.names = F) #

##可视化
##条形图
KEGGcategory<-c("TNF signaling pathway",
"Cytokine-cytokine receptor interaction",
"NF-kappa B signaling pathway",
"JAK-STAT signaling pathway",
"ECM-receptor interaction"
)
keggcategory<-c("Chemical carcinogenesis - reactive oxygen species","Cellular senescence",
                "MAPK signaling pathway","Oxidative phosphorylation","Cell cycle","Apoptosis",
                "p53 signaling pathway","Motor proteins","FoxO signaling pathway","Cell adhesion molecules")

pdf(file="KEGG-barplot.pdf",width = 10,height = 13)
barplot(kk, drop = TRUE, showCategory = 10,label_format=1000)+theme(panel.grid=element_blank())
dev.off()

##气泡图
pdf(file="KEGG-bubble.pdf",width = 10,height = 13)
dotplot(kk, showCategory = KEGGcategory,label_format=1000)+theme(panel.grid=element_blank())
dev.off()

####提取KEGG/GSEA通路的DEGs
kure<-GSEA_GO@result
View(kure)
GS<-GSEA_GO@geneSets
View(GS)
a<-GS$"GO:0001666"
PI3K_symbol <- bitr(a,
                    fromType = "ENTREZID",
                    toType = c("SYMBOL"),
                    OrgDb = org.Mm.eg.db)

#View(PI3K_symbol)
PI3K_symbol<-PI3K_symbol$SYMBOL
PI3K_symbol<-as.data.frame(PI3K_symbol)
PI3K_symbolexp<-deg[deg$SYMBOL %in% PI3K_symbol$PI3K_symbol, ]###提取通路的差异分析
write.csv(PI3K_symbolexp, file = "HIF_DEGs.csv", row.names = FALSE) ##导出数据
PI3K_symbolexpress<-rt[rt$SYMBOL %in% PI3K_symbol$PI3K_symbol, ]###表达矩阵
rownames(PI3K_symbolexp)<-NULL
rownames(PI3K_symbolexpress)<-NULL
rownames(PI3K_symbolexpress)=PI3K_symbolexpress[,1]
exp=PI3K_symbolexpress[,2:ncol(PI3K_symbolexpress)]
dimnames=list(rownames(exp),colnames(exp))
PI3K_symbolexpress=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
PI3K_symbolexpress=avereps(PI3K_symbolexpress)
rownames(PI3K_symbolexp)=PI3K_symbolexp[,1]

####PPI可视化
if (!requireNamespace("STRINGdb", quietly = TRUE)) {
  BiocManager::install("STRINGdb")
}

# 加载包
library(STRINGdb)
string_db <- STRINGdb$new(version = "11.5", species = 10116, score_threshold = 400, input_directory = "")
mapped_genes <- string_db$map(data.frame(gene = a), "gene", removeUnmappedRows = TRUE)
head(mapped_genes)
ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)
head(ppi_network)
library(igraph)

# 构建 igraph 对象
ppi_graph <- graph_from_data_frame(ppi_network, directed = FALSE)

# 绘制网络
plot(ppi_graph, vertex.label.cex = 0.8, vertex.size = 10)
library(ggraph)
library(ggplot2)

# 绘制网络图
ggraph(ppi_graph, layout = "fr") +
  geom_edge_link() +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()

ggraph(ppi_graph, layout = "fr") + # layout 可以替换为更合适的布局
  geom_edge_link(aes(edge_alpha = 0.5, edge_width = 0.8), color = "gray") + # 边的透明度和颜色
  geom_node_point(aes(color = degree(ppi_graph)), size = 6) +  # 节点大小和颜色
  geom_node_text(aes(label = name), size = 3, repel = TRUE) +  # 标签大小和位置
  scale_color_viridis_c() +  # 节点颜色梯度
  theme_void() +  # 去掉背景
  theme(legend.position = "right")  # 调整图例位置

##铁死亡
library(GSVA)
library(GSEABase)
library(ggplot2)

# 2. 读取表达矩阵（行为基因名，列为样本）
expAB <- exp1[,c(6:15)]
expr <- exp1  # 替换为你的路径
expr <- as.matrix(expr)

# 3. 定义铁/铜死亡基因集
FE_positive <- c("Fdx1", "Lias", "Dld", "Pdha1", "Pdhb", "Dlat", "Dbt", "Gcsh")
FE_negative <- c("Mtf1", "Cdkn2a", "Glrx5", "Nfe2l1", "Nfe2l2", "Atox1", "Cox17")
  #c("ACSL4", "NCOA4", "SAT1", "TFRC", "ALOX15", "LPCAT3", "NOX1", "STEAP3", "SLC11A2")

FE_positive <- c("Acsl4", "Ncoa4", "Sat1", "Tfrc", "Alox15", "Lpcat3", "Nox1", "Steap3", "Slc11a2")
FE_negative <- c("Gpx4", "Slc7a11", "Fth1", "Ftl", "Gclm", "Gss", "Slc3a2", "Prnp", "Pcbp1")
  #c("GPX4", "SLC7A11", "FTH1", "FTL", "GCLM", "GSS", "SLC3A2", "PRNP", "PCBP1")
#Cupro_positive <- c("Fdx1", "Lias", "Dld", "Pdha1", "Pdhb", "Dlat", "Dbt", "Gcsh")
#Cupro_negative <- c("Mtf1", "Cdkn2a", "Glrx5", "Nfe2l1", "Nfe2l2", "Atox1", "Cox17")
# 4.定义双硫死亡基因集
FE_positive <- c("Slc7a11", "Nckap1", "Rpn1", "Cyb5r1", "Rac1", "Actb")
FE_negative <- c("G6pd", "Prdx1", "Txnrd1", "Glrx", "Sod1", "Cat", "Mtf1")

# 5.定义经典凋亡基因集
FE_positive <- c("Bax", "Bak1", "Casp3", "Casp7", "Casp8", "Fas", "Fasl", "Trp53", "Bbc3", "Bad", "Bid")
FE_negative <- c("Bcl2", "Bcl2l1", "Mcl1", "Bcl2a1")

# 6. 定义程序性坏死Necroptosis基因集
FE_positive <- c("Ripk1", "Ripk3", "Mlkl", "Tnf", "Tnfrsf1a", "Zbp1", "Fadd", "Tradd")
FE_negative <- c("Cflar", "Birc2", "Birc3")

# 7. 定义钠死亡Sodioptosis基因集
FE_positive <- c("Scnn1a", "Scnn1b", "Scnn1g","Scn1a", "Scn2a", "Scn3a", "Scn8a", "Slc9a1", "Slc9a2", "Slc9a3", 
                 "Slc12a1", "Slc12a2", "Slc12a3", "Slc5a1", "Slc5a2", "Slc5a3", "Aqp1", "Aqp3", "Aqp5", "Trpv4", "Trpm7",                   
                 "Slc8a1", "Slc8a2", "Slc8a3", "Atp1a1", "Atp1a2", "Atp1b1")
FE_negative <-c("Hmox1", "Nfe2l2","Fxyd2", "Fxyd3", "Fxyd7", "Kcnj2", "Kcnj10", "Kcna1",
                "Kcnq1", "Kcnq2", "Kcnq3","Tmem206", "Clcn3", "Best1", "Sgk1","Hsp90aa1", "Hspa1a", "Hspa8")
# 定义Parthanatos死亡基因集
FE_positive  <- c(
  "Parp1", "Parp2","Aifm1","Mif", "H2afx","Xrcc1", "Xrcc5", "Xrcc6", "Neil1", "Neil2", "ParpBp","Nampt", "Nmrk1", 
  "Nmrk2","Sarm1","Nos2", "Cybb", "Nox1", "Nox4", "Endog","Dffb"
)

FE_negative  <- c(
  "Parp3", "Parp9","Nmnat1", "Nmnat2", "Nmnat3", "Qprt", "Idh1", "Idh2", "Mdh1", "Ppargc1a", "Ppargc1b", "Nfe2l2", "Hmox1", 
  "Gpx4", "Prdx1", "Sod2","Bcl2", "Bcl2l1","CypA")


gs1 <- GeneSet(FE_positive, setName = "FE_pos")
gs2 <- GeneSet(FE_negative, setName = "FE_neg")

gsc <- GeneSetCollection(list(gs1, gs2))
# 4. ssGSEA计算
parm<-gsvaParam(expr, gsc)
ssgsea_scores <- gsva(parm)

# 5. 计算FPI
FPI <- ssgsea_scores["FE_pos", ] - ssgsea_scores["FE_neg", ]

# 6. 可视化：箱线图（可根据你的分组替换 group 向量）
group <- c(rep("A", 5), rep("B", 5), rep("C", 5))  # 替换为你的分组
FPI_df <- data.frame(Sample = names(FPI), FPI = FPI, Group = group)
FPI_df$Group <- factor(FPI_df$Group, levels = c("A", "B", "C"))
library(ggpubr)
##单代码
ggplot(FPI_df, aes(x = Group, y = FPI, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  stat_compare_means(method = "t.test", label.y = max(FPI_df$FPI) * 1.1) +
  scale_fill_manual(values = c("A"="skyblue","B" = "lightgreen", "C" = "orange")) +
  labs(title = "Parthanatos Potential Index (SPI)", y = "SPI score") +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 加边框
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
library(ggpubr)

ggplot(FPI_df, aes(x = Group, y = FPI, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  
  # 添加全局检验：ANOVA
  #stat_compare_means(method = "anova", label.y = max(FPI_df$FPI) * 1.15, label = "p.format") +
  
  # 添加两两比较（pairwise t-test）
  stat_compare_means(
    comparisons = list(c("A", "B"), c("A", "C"), c("B", "C")),
    method = "t.test",
    label = "p.signif",
    tip.length = 0.01
  ) +
  
  scale_fill_manual(values = c("A"="skyblue", "B"="lightgreen", "C"="orange")) +
  labs(title = "Parthanatos Potential Index (PPI)", y = "PPI score") +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
# 7. 可选保存
write.table(FPI_df, "Parthanatos_scores.txt", sep = "\t", row.names = FALSE)

#多种死亡热图
library(pheatmap)
library(dplyr)

# 一、设置三种细胞死亡通路的基因（小鼠格式）
ferroptosis_genes <- c("Acsl4", "Ncoa4", "Sat1", "Tfrc", "Alox15", "Lpcat3", "Nox1", "Steap3", "Slc11a2",
                       "Gpx4", "Slc7a11", "Fth1", "Ftl", "Gclm", "Gss", "Slc3a2", "Prnp", "Pcbp1")

cuproptosis_genes <- c("Fdx1", "Lias", "Dld", "Pdha1", "Pdhb", "Dlat", "Dbt", "Gcsh",
                       "Mtf1", "Cdkn2a", "Glrx5", "Nfe2l1", "Nfe2l2", "Atox1", "Cox17")

disulfidptosis_genes <- c("Slc7a11", "Nckap1", "Rpn1", "Cyb5r1", "Rac1", "Actb",
                          "G6pd", "Prdx1", "Txnrd1", "Glrx", "Sod1", "Cat", "Mtf1")

# 二、样本分组信息（每组5个样本）
group <- rep(c("A", "B", "C"), each = 5)
names(group) <- colnames(expr)  # 确保你的expr列名为 A1–A5, B1–B5, C1–C5

# 三、定义分组颜色
group_colors <- list(Group = c("A" = "skyblue", "B" = "lightgreen", "C" = "orange"))

# 四、热图绘图函数（加入颜色设置）
plot_heatmap <- function(gene_list, expr, group_vector, title_text) {
  # 过滤存在于表达矩阵的基因
  gene_list_filtered <- gene_list[gene_list %in% rownames(expr)]
  expr_sub <- expr[gene_list_filtered, ]
  
  # 样本顺序按分组排序
  expr_sub <- expr_sub[, order(group_vector)]
  annotation_col <- data.frame(Group = group_vector[order(group_vector)])
  rownames(annotation_col) <- colnames(expr_sub)
  
  # 画热图
  pheatmap(expr_sub,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_col = annotation_col,
           annotation_colors = group_colors,
           fontsize_row = 10,
           fontsize_col = 10,
           color = colorRampPalette(c("#4874CB", "white", "#D10000"))(50),
           main = title_text)
}
#"#3EC1D3", "white", "#D10000"
#
# 五、绘制三类死亡通路热图
plot_heatmap(ferroptosis_genes, expr, group, "Ferroptosis Gene Expression")
plot_heatmap(cuproptosis_genes, expr, group, "Cuproptosis Gene Expression")
plot_heatmap(disulfidptosis_genes, expr, group, "Disulfidptosis Gene Expression")
