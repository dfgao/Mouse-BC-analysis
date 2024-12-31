
########################################################################## bulk RNA-seq section ############################################################

# 01.data import -----

# gene info 
mus.gene.info <- read.table('mouse.110.gtf',header = F,sep = '\t')
colnames(mus.gene.info) <- c('geneid','symbol','biotype')
mus.gene.info[mus.gene.info$symbol == 'ensembl',]$symbol <- mus.gene.info[mus.gene.info$symbol == 'ensembl',]$geneid
mus.gene.info[mus.gene.info$symbol == 'ensembl_havana',]$symbol <- mus.gene.info[mus.gene.info$symbol == 'ensembl_havana',]$geneid
# counts
rna.counts <- read.table('../04.rnaseq/02.quant/all_counts.txt',header = T,sep = '\t',row.names = NULL)
colnames(rna.counts)

samp.name <- c(paste0('Young+PBS',c('10',seq(1,9))), 
               paste0('Aged+PBS',c('10',seq(1,9))),
               paste0('Aged+NR',c('10',seq(1,9))),
               paste0('Aged+Cd38i',c('10',seq(1,9))))

colnames(rna.counts) <- c('geneid',samp.name)
rna.counts <- rna.counts %>% 
  left_join(y = mus.gene.info, by = 'geneid') %>% 
  na.omit() %>% 
  dplyr::filter(biotype == 'protein_coding') %>% 
  dplyr::group_by(`symbol`) %>% slice_sample(n = 1) %>% 
  column_to_rownames(var = 'symbol') %>% 
  dplyr::select(!c('geneid','biotype')) %>% 
  round(0)

# TPM 
rna.tpm <- read.table('../04.rnaseq/02.quant/all_quant_TPM.txt',header = T,sep = '\t',row.names = NULL)
colnames(rna.tpm) <- c('geneid',samp.name)
rna.tpm <- rna.tpm %>% 
  left_join(y = mus.gene.info, by = 'geneid') %>% 
  na.omit() %>% 
  dplyr::filter(biotype == 'protein_coding') %>% 
  dplyr::group_by(`symbol`) %>% slice_sample(n = 1) %>% 
  column_to_rownames(var = 'symbol') %>% 
  dplyr::select(!c('geneid','biotype')) %>% round(2)

# pick genes 
pick.genes <- unique(c(rownames(rna.tpm[matrixStats::rowMedians(as.matrix(rna.tpm),cols = 1:10) > 0.1, ]),
                       rownames(rna.tpm[matrixStats::rowMedians(as.matrix(rna.tpm),cols = 11:20) > 0.1, ]),
                       rownames(rna.tpm[matrixStats::rowMedians(as.matrix(rna.tpm),cols = 21:30) > 0.1, ]),
                       rownames(rna.tpm[matrixStats::rowMedians(as.matrix(rna.tpm),cols = 31:40) > 0.1, ])
))

# clean dataset 
rna.counts.clean <- rna.counts[pick.genes,]
rna.tpm.clean <- rna.tpm[pick.genes,]

# 02.PCA plot -----

# annotation info
ann <- data.frame(Treatment = factor(c(rep('Young+PBS',10),rep('Aged+PBS',10),rep('Aged+NR',10),rep('Aged+Cd38i',10)), levels = c('Young+PBS','Aged+PBS','Aged+NR','Aged+Cd38i')))
rownames(ann) <- colnames(rna.counts)
group.cols <- pal_npg(palette = c("nrc"), alpha = 1)(4)
group_list_large <- ann$Treatment

# gene exp section 
d_sec <- function(tpm){
  gene.exp.sec <- data.frame(TPM_100 = apply(tpm,2,function(x) {table(x>=100)["TRUE"]}),
                             TPM_100_50 = apply(tpm,2,function(x) {table(x>=50 & x< 100)["TRUE"]}),
                             TPM_50_5 = apply(tpm,2,function(x) {table(x>=5 & x< 50)["TRUE"]}),
                             TPM_5_0.5 = apply(tpm,2,function(x) {table(x>=0.5 & x< 5)["TRUE"]}),
                             TPM_0.5 = apply(tpm,2,function(x) {table(x< 0.5)["TRUE"]})
  )
  gene.exp.sec <- gene.exp.sec %>% rownames_to_column(var = 'Sample')
  
  test <- melt(gene.exp.sec,id.vars = 'Sample',variable.name = 'Range',value.name = 'genenum')
  test$genenum <- as.integer(test$genenum)
  test$Sample <- factor(test$Sample, levels = gene.exp.sec$Sample)
  ggplot(test) +
    theme_bw() +
    geom_bar(stat = "identity",aes(x=Sample,y=genenum,fill=Range),alpha=.8,position = 'stack') +
    scale_fill_brewer(palette = "RdBu",direction = 1,) +
    labs(y='Gene number',x=NULL) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 15)
    )
}
d_sec(rna.tpm.clean)

# gene complex 
d_comp <- function(tpm,group_list){
  all.tpm.sort <- data.frame(apply(tpm, 2, function(x){sort(x,decreasing = T)}))
  gene.comp <- all.tpm.sort
  
  for (i in colnames(all.tpm.sort)) {
    gene.comp[,i] <- all.tpm.sort[,i]/colSums(all.tpm.sort)[i]
  }
  
  gene.comp <- data.frame(apply(gene.comp, 2, cumsum))
  
  plot.data <- data.frame(t(gene.comp),Group=group_list) %>% rownames_to_column(var = 'Sample')
  colnames(plot.data) <- c('Sample',(seq(1:nrow(tpm))),'Group')
  plot.data.log <- melt(plot.data,id.vars = c('Group','Sample'),variable.name = 'Genenum',value.name = 'Percent')
  plot.data.log$Genenum <- as.numeric(plot.data.log$Genenum)
  plot.data.log$Sample <- factor(plot.data.log$Sample, levels = rownames(ann))
  plot.data.log$Percent <- plot.data.log$Percent*100
  
  cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=1)(6))
  cols <- as.matrix(apply(data.frame(cols),1,function(x){rep(x,12)}))
  cols <- as.vector(cols)
  ggplot(data = plot.data.log) +
    theme_bw() +
    geom_line(aes(x=Genenum, y=Percent,color=Sample),size=.7 ) + 
    scale_x_log10() +
    scale_color_manual(values = group.cols[as.numeric(group_list)]) +
    labs(x = 'Genes number',y = 'Gene accumulation ratio (%)') +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.key.height=unit(.85,"line"))
}
d_comp(rna.counts.clean, group_list = rownames(ann))

# tpm cor 
tpm_cor <- cor(rna.tpm.clean,method = 'spearman',use = 'pairwise.complete.ob')
ha_left.col <- list(Treatment = c('#00afb9','#eae2b7','#f77f00','#d62828'))
names(ha_left.col$Treatment) <- factor(unique(ann$Treatment))

ComplexHeatmap::pheatmap(tpm_cor,
                         name = "Spearman's coeff",
                         border_color = NA, 
                         clustering_method = 'ward.D',
                         show_colnames = F,
                         # color = scales::alpha(colorRampPalette(colors = c('#00509d','gray80','#e63946'),alpha=T,bias=1)(256),alpha = 1),
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1)(256),alpha = 1),
                         angle_col = '45',
                         annotation_col = ann,
                         annotation_colors = ha_left.col)

# TPM PCA
library(factoextra)
library(FactoMineR)
d_p <- function(tpm,group_list,name){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 40),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               axes = c(1, 2),
               mean.point=F,
               addEllipses = T,
               ellipse.type = "convex", #convex norm confidence t euclid
               # repel = T,
               habillage = group_list,
               label = "none",
               geom.ind = c("point"),
               fill.ind = tpm.t$group_list,
               palette = c('#00afb9','#eae2b7','#f77f00','#d62828'),
               # palette = c("#7C6DAF",'#C7AED4','#0B996F','#83c5be','#D6570D','#e29578'),
               legend.title = "Groups",
               pointsize = 3,
               pointshape = 21,
               col.ind = 'dark',
               title = name
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}

d_p(rna.tpm.clean,group_list_large,'tpm')
d_p(rna.counts.clean,group_list_large,'data_raw')

# rlog PCA 
library(DESeq2)
library(EDASeq)
colData = data.frame(row.names = colnames(rna.counts.clean),group_list = group_list_large)
dds <- DESeqDataSetFromMatrix(countData = rna.counts.clean,colData = colData,design = ~group_list)

set.tmp <- newSeqExpressionSet(as.matrix(rna.counts.clean),
                               phenoData = colData)
EDASeq::plotRLE(set.tmp, outline=FALSE, ylim=c(-4, 4))
DESeq2::plotPCA(set.tmp,col=group.cols, cex=1.2)

rld_rlog <- rlog(dds, blind=TRUE,fitType = "glmGamPoi")
rld_vst <- vst(dds, blind=TRUE)
DESeq2::plotPCA(rld_rlog, intgroup = "group_list",ntop = 500)
vst_mat <- assay(rld_vst) %>% as.data.frame()
rlog_mat <- assay(rld_rlog) %>% as.data.frame()

rlog_cor <- cor(rlog_mat, method = 'spearman')
vst_cor <- cor(vst_mat, method = 'spearman')

ComplexHeatmap::pheatmap(vst_cor,
                         name = "Spearman's coeff",
                         border_color = NA, 
                         clustering_method = 'ward.D2',
                         show_colnames = F,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1.2)(256),alpha = 1),
                         angle_col = '45',
                         annotation_col = ann,
                         annotation_colors = ha_left.col
)

p1 <- d_p(rlog_mat,group_list = group_list_large,'RLOG')
p2 <- d_p(vst_mat,group_list = group_list_large,'VST')
p3 <- d_p(rna.tpm.clean,group_list = group_list_large,'TPM')
p1 + p2 + p3

# 03.Euclidean Distance ----
d_p_pca <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
}
tpm.pca <- d_p_pca(vst_mat,group_list_large)

dist.pca <- tpm.pca$ind$coord[,1:2] %>% as.data.frame()
dist.pca.euc <- dist(dist.pca,method = 'euclidean') %>% as.matrix()
dist.young.pbs_age.pbs <- dist.pca.euc[1:20,1:20] %>% as.data.frame()
dist.young.pbs_age.nr <- dist.pca.euc[c(1:10,21:30),c(1:10,21:30)] %>% as.data.frame()
dist.young.pbs_age.cd38i <- dist.pca.euc[c(1:10,31:40),c(1:10,31:40)] %>% as.data.frame()

dist.young.pbs_age.pbs.long <- melt(dist.young.pbs_age.pbs)[1:c(dim(dist.young.pbs_age.pbs)[1] * dim(dist.young.pbs_age.pbs)[2] / 2),]
dist.young.pbs_age.nr.long <- melt(dist.young.pbs_age.nr)[1:c(dim(dist.young.pbs_age.nr)[1] * dim(dist.young.pbs_age.nr)[2] / 2),]
dist.young.pbs_age.cd38i.long <- melt(dist.young.pbs_age.cd38i)[1:c(dim(dist.young.pbs_age.cd38i)[1] * dim(dist.young.pbs_age.cd38i)[2] / 2),]

dist.group <- data.frame(dist = c(dist.young.pbs_age.pbs.long$value, dist.young.pbs_age.nr.long$value, dist.young.pbs_age.cd38i.long$value),
                         group = c(rep('Young.pbs_Age.pbs',200),rep('Young.pbs_Age.NR',200), rep('Young.pbs_Age.cd38i',200)))

p1 <- ggpubr::ggboxplot(dist.group, x="group", y="dist", width = 0.5, 
                        color = "black",#轮廓颜色 
                        fill="group",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细 
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) + 
  geom_jitter(color="black", size=0.4, alpha=1,width = 0.1, height = 0.5)+
  ylab('Euclidean distance')  + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),legend.position = 'none')

my_comparisons <- list(c("Young.pbs_Age.pbs", "Young.pbs_Age.NR"),c('Young.pbs_Age.pbs','Young.pbs_Age.cd38i'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

# 04.key genes ------
Elovl3 <- rna.tpm.clean['Cd38',] %>% t() %>% as.data.frame()
Elovl3$group <- rep(c(rep('MC+PBS',6), rep('MO+FGF21',6), rep('MO+PBS',6)),2)
Elovl3$Elovl3 <- log2(Elovl3$Elovl3+1)

p1 <- ggpubr::ggboxplot(Elovl3, x="group", y="Elovl3", width = 0.5, 
                        color = "black",#轮廓颜色 
                        fill="group",#填充
                        # palette = 
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=0.7, #箱型图边线的粗细 
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.6) + 
  geom_jitter(color="black", size=0.4, alpha=1,width = 0.1, height = 0.5) +
  ylab('Log2(TPM+1)')  + 
  theme_bw() + scale_fill_manual(values = c('#0B996F',"#7C6DAF",'#D6570D'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("MO+FGF21","MO+PBS"),c("MC+PBS", "MO+PBS"))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

# key genes -

key.genes <- c('Ucp1','Prdm16','Cidea','Cidea','Cidea','Cidea')

colData.clean <- data.frame(group_list = colData)
rownames(colData.clean) <- rownames(colData)
dds <- DESeqDataSetFromMatrix(countData = rna.counts.clean,colData = colData.clean,design = ~group_list)

rld <- rlog(dds, blind=TRUE,fitType = "glmGamPoi")
rld_mat <- assay(rld) %>% as.data.frame()
# rld <- vst(dds, blind=TRUE)
ha_left.col <- list(Group = c('#D6570D','#7C6DAF','#0B996F'),
                    Gender = c('#f47068','#0e606b'))
names(ha_left.col$Group) <- factor(unique(ann$Group))
names(ha_left.col$Gender) <- factor(unique(ann$Gender))

p1 <- ComplexHeatmap::pheatmap(hl.tpm.clean[key.genes,19:36],
                               name = "Z-score",
                               border_color = NA, 
                               clustering_method = 'ward.D',cluster_cols = F,
                               show_colnames = F,
                               scale = 'row',
                               color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1.2)(256),alpha = 1),
                               angle_col = '45',
                               annotation_col = ann[19:36,],
                               annotation_colors = ha_left.col
)
p2 <- ComplexHeatmap::pheatmap(hl.tpm.clean[key.genes,c(1:6,8,9,11:18)],
                               name = "Z-score",
                               border_color = NA, 
                               clustering_method = 'ward.D',cluster_cols = F,
                               show_colnames = F,
                               scale = 'row',
                               color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1.2)(256),alpha = 1),
                               angle_col = '45',
                               annotation_col = ann[c(1:6,8,9,11:18),],
                               annotation_colors = ha_left.col
)
p2 <- ComplexHeatmap::pheatmap(hl.tpm.clean[key.genes,c(1:18)],
                               name = "Z-score",
                               border_color = NA, 
                               clustering_method = 'ward.D',cluster_cols = F,
                               show_colnames = F,
                               scale = 'row',
                               color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1.2)(256),alpha = 1),
                               angle_col = '45',
                               annotation_col = ann[c(1:18),],
                               annotation_colors = ha_left.col)
p1+p2

p3 <- ComplexHeatmap::pheatmap(rld_mat[key.genes,c(1:6,8,9,11:18)],
                               name = "Z-score",
                               scale = 'row',
                               border_color = NA, 
                               clustering_method = 'complete',cluster_cols = F,
                               show_colnames = F,
                               color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1.2)(256),alpha = 1),
                               angle_col = '45',
                               annotation_col = ann[c(1:6,8,9,11:18),],
                               annotation_colors = ha_left.col
)

p4 <- ComplexHeatmap::pheatmap(rld_mat[key.genes,19:36],
                               name = "Z-score",
                               scale = 'row',
                               border_color = NA, 
                               clustering_method = 'complete',cluster_cols = F,
                               show_colnames = F,
                               color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1.2)(256),alpha = 1),
                               angle_col = '45',
                               annotation_col = ann[19:36,],
                               annotation_colors = ha_left.col
)
p4+p3

Cidea <- hl.tpm.clean['Cidea',]
# Cidea[1,] <- log2(Cidea[1,])
Cidea <- reshape2::melt(Cidea)
# Cidea$variable <- rownames(Cidea)


ggplot(Cidea[1:18,],aes(x=variable,y=value, fill=variable)) +
  geom_bar(stat = "identity",width = .6) +
  theme_bw() +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10),
        axis.text.y = element_text(vjust = 1, hjust=1,size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  ggtitle('Cidea-female') +
  xlab('Samples') + ylab('TPM') +
  # scale_fill_manual(values = rep(c("#7C6DAF",'#C7AED4','#0B996F','#83c5be','#D6570D','#e29578'),each = 6))
  scale_fill_manual(values = c(rep('#0B996F',6), rep('#7C6DAF',6),rep('#D6570D',6)))

ggplot(Cidea[19:36,],aes(x=variable,y=value, fill=variable)) +
  geom_bar(stat = "identity",width = .6) +
  theme_bw() +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10),
        axis.text.y = element_text(vjust = 1, hjust=1,size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  ggtitle('Cidea-male') +
  xlab('Samples') + ylab('TPM') +
  # scale_fill_manual(values = rep(c("#7C6DAF",'#C7AED4','#0B996F','#83c5be','#D6570D','#e29578'),each = 6))
  scale_fill_manual(values = c(rep('#83c5be',6), rep('#C7AED4',6),rep('#e29578',6)))


# 05.DEGs ----
# BiocManager::install("IHW") # v1.31
# BiocManager::install("ashr")
# BiocManager::install('apeglm')
options(scipen = 10)

library(IHW)
library(ashr)
library(tidyverse)
#fat: obe vs nor

## data perpare
comp_fat <- rna.counts.clean[,c(13:18,31:36,1:6,19:24)]
colData_fat <- data.frame(group_list = colData[c(13:18,31:36,1:6,19:24),])
colData_fat$group_list_large <- c(rep('Obesity',12),rep('Normal',12))
rownames(colData_fat) <- colnames(comp_fat)

dds <- DESeqDataSetFromMatrix(countData = comp_fat,colData = colData_fat,design = ~group_list_large)

## pearson and PCA
rld <- rlog(dds, blind=TRUE,fitType = "glmGamPoi")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat, method = 'pearson')

ha_left.col <- list(Group = c('#D6570D','#0B996F'),
                    Gender = c('#f47068','#0e606b'))
names(ha_left.col$Group) <- factor(unique(ann[c(13:18,31:36,1:6,19:24),]$Group))
names(ha_left.col$Gender) <- factor(unique(ann[c(13:18,31:36,1:6,19:24),]$Gender))

ComplexHeatmap::pheatmap(rld_cor,
                         name = "Pearson's coeff",
                         border_color = NA, 
                         clustering_method = 'ward.D',
                         show_colnames = F,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1)(256),alpha = 1),
                         angle_col = '45',
                         annotation_col = ann[c(13:18,31:36,1:6,19:24),],
                         annotation_colors = ha_left.col
)

d_p_2g <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 40),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               mean.point=F,
               addEllipses = T,
               ellipse.type = "convex", #convex norm confidence t euclid
               # repel = T,
               habillage = group_list,
               label = "none",
               geom.ind = c("point"),
               fill.ind = tpm.t$group_list,
               palette = c("#0B996F",'#D6570D'),
               legend.title = "Groups",
               pointsize = 3,
               pointshape = 21,
               col.ind = 'dark',
               title = 'PCA norm ellipse'
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}
d_p_2g(rld_mat,group_list_large[c(13:18,31:36,1:6,19:24)])

## DEGs
dds <- DESeq(dds,minReplicatesForReplace = 7)
plotDispEsts(dds)
# fat_res_ihw <- results(dds, contrast = c('group_list_large','Obesity','Normal'), filterFun = ihw) # adjust P value
# fat_res_ihw <- results(dds, contrast = c('group_list_large','Obesity','Normal')) # adjust P value
# 
# summary(fat_res_ihw)
# fat_res_ihw <- fat_res_ihw %>%
#   data.frame() %>%
#   rownames_to_column(var="geneid") %>%
#   as_tibble() %>%
#   arrange('padj')

fat_res_ashr <- lfcShrink(dds, contrast = c('group_list_large','Normal','Obesity'), type = 'ashr') # adjust LFC
summary(fat_res_ashr)
fat_res_ashr <- fat_res_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')
# fat_res_ihw_ashr <- fat_res_ihw %>% 
#   left_join(y = fat_res_ashr[,c(1,3)], by = 'geneid')

# rescue fgf vs nor

comp_resc <- rna.counts.clean[,c(7:12,25:30,1:6,19:24)]
colData_resc <- data.frame(group_list = colData[c(7:12,25:30,1:6,19:24),])
colData_resc$group_list_large <- c(rep('FGF',12),rep('Normal',12))
rownames(colData_resc) <- colnames(comp_resc)

dds <- DESeqDataSetFromMatrix(countData = comp_resc,colData = colData_resc,design = ~group_list_large)

## pearson and PCA
rld <- rlog(dds, blind=TRUE,fitType = "glmGamPoi")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat, method = 'pearson')

ha_left.col <- list(Group = c('#0B996F','#7C6DAF'),
                    Gender = c('#f47068','#0e606b'))
names(ha_left.col$Group) <- factor(unique(ann[c(7:12,25:30,1:6,19:24),]$Group))
names(ha_left.col$Gender) <- factor(unique(ann[c(7:12,25:30,1:6,19:24),]$Gender))

ComplexHeatmap::pheatmap(rld_cor,
                         name = "Pearson's coeff",
                         border_color = NA, 
                         clustering_method = 'ward.D',
                         show_colnames = F,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1)(256),alpha = 1),
                         angle_col = '45',
                         annotation_col = ann[c(7:12,25:30,1:6,19:24),],
                         annotation_colors = ha_left.col
)

d_p_2g <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 40),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               mean.point=F,
               addEllipses = T,
               ellipse.type = "convex", #convex norm confidence t euclid
               # repel = T,
               habillage = group_list,
               label = "none",
               geom.ind = c("point"),
               fill.ind = tpm.t$group_list,
               palette = c("#7C6DAF",'#0B996F'),
               legend.title = "Groups",
               pointsize = 3,
               pointshape = 21,
               col.ind = 'dark',
               title = 'PCA norm ellipse'
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}
d_p_2g(rld_mat,group_list_large[c(7:12,25:30,1:6,19:24)])

## DEGs
dds <- DESeq(dds,minReplicatesForReplace = 7)
plotDispEsts(dds)
# resc_res_ihw <- results(dds, contrast = c('group_list_large','FGF','Normal')) # adjust P value
# summary(resc_res_ihw)
# resc_res_ihw <- resc_res_ihw %>%
#   data.frame() %>%
#   rownames_to_column(var="geneid") %>%
#   as_tibble() %>%
#   arrange('padj')

resc_res_ashr <- lfcShrink(dds, contrast = c('group_list_large','FGF','Normal'), type = 'ashr') # adjust LFC
summary(resc_res_ashr)
resc_res_ashr <- resc_res_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')
# resc_res_ihw_ashr <- resc_res_ihw %>% 
#   left_join(y = resc_res_ashr[,c(1,3)], by = 'geneid')

# effect obe vs fgf 

comp_efct <- rna.counts.clean[,c(13:18,31:36,7:12,25:30)]
colData_efct <- data.frame(group_list = colData[c(13:18,31:36,7:12,25:30),])
colData_efct$group_list_large <- c(rep('Obesity',12),rep('FGF',12))
rownames(colData_efct) <- colnames(comp_efct)

dds <- DESeqDataSetFromMatrix(countData = comp_efct,colData = colData_efct,design = ~group_list_large)

## pearson and PCA
rld <- rlog(dds, blind=TRUE,fitType = "glmGamPoi")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat, method = 'pearson')

ha_left.col <- list(Group = c('#D6570D','#7C6DAF'),
                    Gender = c('#f47068','#0e606b'))
names(ha_left.col$Group) <- factor(unique(ann[c(13:18,31:36,7:12,25:30),]$Group))
names(ha_left.col$Gender) <- factor(unique(ann[c(13:18,31:36,7:12,25:30),]$Gender))

ComplexHeatmap::pheatmap(rld_cor,
                         name = "Pearson's coeff",
                         border_color = NA, 
                         clustering_method = 'average',
                         show_colnames = F,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1)(256),alpha = 1),
                         angle_col = '45',
                         annotation_col = ann[c(13:18,31:36,7:12,25:30),],
                         annotation_colors = ha_left.col
)

d_p_2g <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 40),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               mean.point=F,
               addEllipses = T,
               ellipse.type = "norm", #convex norm confidence t euclid
               # repel = T,
               habillage = group_list,
               label = "none",
               geom.ind = c("point"),
               fill.ind = tpm.t$group_list,
               palette = c("#7C6DAF",'#D6570D'),
               legend.title = "Groups",
               pointsize = 3,
               pointshape = 21,
               col.ind = 'dark',
               title = 'PCA norm ellipse'
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}
d_p_2g(rld_mat,group_list_large[c(13:18,31:36,7:12,25:30)])

## DEGs
dds <- DESeq(dds,minReplicatesForReplace = 7)
plotDispEsts(dds)
# efct_res_ihw <- results(dds, contrast = c('group_list_large','Obesity','FGF')) # adjust P value
# summary(efct_res_ihw)
# efct_res_ihw <- efct_res_ihw %>%
#   data.frame() %>%
#   rownames_to_column(var="geneid") %>%
#   as_tibble() %>%
#   arrange('padj')

efct_res_ashr <- lfcShrink(dds, contrast = c('group_list_large','FGF','Obesity'), type = 'ashr') # adjust LFC
summary(efct_res_ashr)
efct_res_ashr <- efct_res_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')
# efct_res_ihw_ashr <- efct_res_ihw %>% 
#   left_join(y = efct_res_ashr[,c(1,3)], by = 'geneid')

# 06.GSEA -----
library(clusterProfiler)
library(enrichplot)
library(GSEABase)
library(org.Mm.eg.db)

gene_for_gsea <- fat_res_ashr$log2FoldChange 
names(gene_for_gsea) <- fat_res_ashr$geneid

geneList=sort(gene_for_gsea,decreasing = T)

mus.gmt <- read.gmt('E:/GSEA/v2023.mouse/msigdb_v2023.2.Mm_files_to_download_locally/msigdb_v2023.2.Mm_GMTs/m5.go.bp.v2023.2.Mm.symbols.gmt')
mus.gmt.kegg <- read.gmt('E:/GSEA/v2023.mouse/msigdb_v2023.2.Mm_files_to_download_locally/msigdb_v2023.2.Mm_GMTs/m2.cp.v2023.2.Mm.symbols.gmt')

length(unique(mus.gmt$term))
egmt <- GSEA(geneList, TERM2GENE=mus.gmt, 
             minGSSize = 1,
             pvalueCutoff = 0.1,
             verbose=FALSE)

egmt.order <- egmt[order(egmt$enrichmentScore, decreasing = T),]

gseaplot2(egmt,geneSetID = 'GOBP_REGULATION_OF_FATTY_ACID_BETA_OXIDATION', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)
gseaplot2(egmt,geneSetID = 'GOBP_ADAPTIVE_THERMOGENESIS', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)
gseaplot2(egmt,geneSetID = 'GOBP_ADAPTIVE_IMMUNE_RESPONSE', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)
egmt.order$core_enrichment[1]

egmt.wp <- GSEA(geneList, TERM2GENE=mus.gmt.kegg, 
                minGSSize = 1,
                pvalueCutoff = 0.1,
                verbose=FALSE)

egmt.wp.order <- egmt.wp[order(egmt.wp$enrichmentScore, decreasing = T),]

gseaplot2(egmt.wp,geneSetID = 'WP_FATTY_ACID_BETA_OXIDATION', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)
gseaplot2(egmt.wp,geneSetID = 'REACTOME_FATTY_ACID_METABOLISM', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)


# efct
gene_for_gsea <- efct_res_ashr$log2FoldChange 
names(gene_for_gsea) <- fat_res_ashr$geneid

geneList=sort(gene_for_gsea,decreasing = T)

egmt.efct <- GSEA(geneList, TERM2GENE=mus.gmt, 
                  minGSSize = 1,
                  pvalueCutoff = 0.1,
                  verbose=FALSE)

egmt.efct.order <- egmt[order(egmt.efct$enrichmentScore, decreasing = T),]

gseaplot2(egmt.efct,geneSetID = 'GOBP_REGULATION_OF_FATTY_ACID_BETA_OXIDATION', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)
gseaplot2(egmt.efct,geneSetID = 'GOBP_ADAPTIVE_THERMOGENESIS', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)
gseaplot2(egmt.efct,geneSetID = 'GOBP_ADAPTIVE_IMMUNE_RESPONSE', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)

egmt.efct.wp <- GSEA(geneList, TERM2GENE=mus.gmt.kegg, 
                     minGSSize = 1,
                     pvalueCutoff = 0.1,
                     verbose=FALSE)

egmt.efct.wp.order <- egmt.efct.wp[order(egmt.efct.wp$enrichmentScore, decreasing = T),]

gseaplot2(egmt.efct.wp,geneSetID = 'WP_FATTY_ACID_BETA_OXIDATION', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)
gseaplot2(egmt.efct.wp,geneSetID = 'REACTOME_FATTY_ACID_METABOLISM', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)

library(clusterProfiler)
library(org.Mm.eg.db)


gene.ent <- clusterProfiler::bitr(fat_res_ashr[fat_res_ashr$change == 'mND+PBS up',]$geneid, fromType = 'SYMBOL', 
                                  toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
mND <- clusterProfiler::enrichGO(gene = gene.ent$ENTREZID, 
                                 keyType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', ont = 'BP', 
                                 pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                 readable = T)
mND.df <- mND@result[mND@result$pvalue < 0.05,]


gene.ent <- clusterProfiler::bitr(fat_res_ashr[fat_res_ashr$change == 'mHFD+PBS up',]$geneid, fromType = 'SYMBOL', 
                                  toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
mHFD <- clusterProfiler::enrichGO(gene = gene.ent$ENTREZID, 
                                  keyType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', ont = 'BP', 
                                  pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                  readable = T)
mHFD.df <- mHFD@result[mHFD@result$pvalue < 0.05,]

gene.ent <- clusterProfiler::bitr(efct_res_ashr[efct_res_ashr$change == 'mHFD+FGF21 up',]$geneid, fromType = 'SYMBOL', 
                                  toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
mFGF <- clusterProfiler::enrichGO(gene = gene.ent$ENTREZID, 
                                  keyType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', ont = 'BP', 
                                  pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                  readable = T)
mFGF.df <- mFGF@result[mFGF@result$pvalue < 0.05,]


gene.ent <- clusterProfiler::bitr(efct_res_ashr[efct_res_ashr$change == 'mHFD+PBS up',]$geneid, fromType = 'SYMBOL', 
                                  toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
mHFD.2 <- clusterProfiler::enrichGO(gene = gene.ent$ENTREZID, 
                                    keyType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', ont = 'BP', 
                                    pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                    readable = T)
mHFD2.df <- mHFD.2@result[mHFD.2@result$pvalue < 0.05,]

bulk.go <- list(MC_MO.MC.up.GO = mND.df,
                MC_MO.MO.up.GO = mHFD.df,
                FGF_MO.FGF.up.GO = mFGF.df,
                FGF_MO.MO.up.GO = mHFD2.df)
rio::export(bulk.go, file = '../../../all.fig&table.gdf_24.5.11/add.5.21/bulk.GO.xlsx')


# 07.decon RNAseq ----
library(DeconRNASeq)
cell.type.mean <- data.frame(AverageExpression(Clean_sct.inte.rm.lpt, layer = 'counts', assays = 'SCT', group.by = 'cell.type.percise.new'))
colnames(cell.type.mean) <- levels(Clean_sct.inte.rm.lpt$cell.type.percise.new)

cell.type.mean.cpm <- data.frame(edgeR::cpm(cell.type.mean))
rna.counts.clean.cpm <- data.frame(edgeR::cpm(rna.counts.clean))
rna.counts.clean.cpm <- rna.counts.clean.cpm[,-c(8,9,
                                                 17,20,
                                                 22,26,
                                                 32,40
)]

cell.fraction <- DeconRNASeq(rna.counts.clean.cpm, cell.type.mean.cpm,checksig=FALSE,
                             known.prop = F, use.scale = TRUE, fig = TRUE)
# cell.fraction <- DeconRNASeq(rna.counts.clean.cpm, cell.type.mean.cpm[unique(rm.lpt.ct.deg$gene),],checksig=FALSE,
#                              known.prop = F, use.scale = TRUE, fig = TRUE)
cell.fraction <- DeconRNASeq(rna.counts.clean.cpm, cell.type.mean.cpm[,c(1:3,8:15)],checksig=FALSE,
                             known.prop = F, use.scale = TRUE, fig = TRUE)

cell.fraction.plot <- data.frame(cell.fraction$out.all)
cell.fraction.plot$Samples <- colnames(rna.counts.clean.cpm)
cell.fraction.plot <- data.table::melt(cell.fraction.plot,variable.name = 'cell.type',id.vars = 'Samples',value.name = 'Fraction')
# cell.fraction.plot$cell.type <- factor(cell.fraction.plot$cell.type, levels = c('GCs','S-TGCs','SpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells'))

ggplot(cell.fraction.plot, aes(fill=cell.type, y=Fraction, x=Samples)) + 
  geom_bar( stat="identity") + 
  coord_flip() + 
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_bw() +
  scale_fill_manual( values = c(use.cols,npg.cols))

# cibersortx web
randoms <- sample(1:ncol(Clean_sct.inte.rm.lpt),10000)
sc.mat <- GetAssayData(Clean_sct.inte.rm.lpt,assay = 'SCT',layer = 'counts')
colnames(sc.mat) <- Clean_sct.inte.rm.lpt@meta.data$cell.type.percise.new
sc.mat <- sc.mat[,randoms] %>% as.data.frame()
write.table(sc.mat, file = 'sc.mat.txt',sep = '\t',quote = F,row.names = T,col.names = T)
bulk.mat <- rna.counts.clean
write.table(bulk.mat, file = 'bulk.mat.txt',sep = '\t',quote = F,row.names = T,col.names = T)
