
# 01.umap plot ----
Clean_sct.inte.rm.lpt$cell.type.percise.new <- as.character(Clean_sct.inte.rm.lpt$cell.type.percise.raw)

Clean_sct.inte.rm.lpt$cell.type.percise.new <- factor(Clean_sct.inte.rm.lpt$cell.type.percise.new, levels = c('Top2a + cancer cells', 'Padi4 + cancer cells', 'Cdh3 + Epithelial cells', 'Erbb4 + Epithelial cells', 'Prlr + Epithelial cells','Endothelial cells', 'Pericytes', 'MSCs','Macrophages','DCs', 'Monocytes', 'Neutrophils','T cells'))

Clean_sct.inte.rm.lpt$condition <- as.character(Clean_sct.inte.rm.lpt$condition)
Clean_sct.inte.rm.lpt$condition[which(str_detect(Clean_sct.inte.rm.lpt$condition, "WT"))] <- "Untreated"
Clean_sct.inte.rm.lpt$condition <- factor(Clean_sct.inte.rm.lpt$condition, levels = c('Untreated','Targeted'))

DimPlot(Clean_sct.inte.rm.lpt, group.by = 'cell.type.percise.new', cols = c(use.cols,npg.cols),split.by = 'condition', label = F,repel = T) + ggtitle("") + NoAxes()

p1 <- DimPlot(Clean_sct.inte.rm.lpt, group.by = 'cell.type.percise.new', cols = c( use.cols,npg.cols), label = F) + ggtitle("") + NoAxes()
p2 <- DimPlot(Clean_sct.inte.rm.lpt, group.by = 'cell.type.percise.new', cols = c(use.cols,npg.cols), label = T,repel = T) + ggtitle("") + NoAxes() + NoLegend()
p3 <- DimPlot(Clean_sct.inte.rm.lpt, group.by = 'cell.type.percise.new', cols = c(use.cols,npg.cols), label = T,repel = T) + ggtitle("") + NoAxes()
p1 + p2 +p3

# 02.miloR:composition of neighbourhood of cells ----
library(miloR)
library(scales)
library(SingleCellExperiment)
library(ggbeeswarm)

milo.all.seu <- subset(Clean_sct.inte.rm.lpt, cells = c(sample( rownames(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$condition == 'Untreated',]), 14761),
                                                        rownames(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$condition == 'Targeted',] )))
milo.all.seu <- as.SingleCellExperiment(milo.all.seu,assay = 'SCT') %>%
  Milo() %>%
  miloR::buildGraph(k = 30, d = 50) %>% 
  makeNhoods(
    prop = 0.2,                 
    k = 30,       
    d=50,                   
    refined = T)

milo.all.seu <- countCells(milo.all.seu, 
                           meta.data = data.frame(colData(milo.all.seu)),                     
                           sample="orig.ident")
milo.traj_design <- data.frame(colData(milo.all.seu))[,c("orig.ident", "condition")]
milo.traj_design$orig.ident <- as.factor(milo.traj_design$orig.ident)
milo.traj_design <- distinct(milo.traj_design)
rownames(milo.traj_design) <- milo.traj_design$orig.ident

milo.all.seu <- calcNhoodDistance(milo.all.seu, d=50)
milo.da_results <- testNhoods(milo.all.seu,                 
                              design = ~ condition,            
                              design.df = milo.traj_design)
milo.all.seu <- buildNhoodGraph(milo.all.seu)
ggplot(milo.da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(milo.da_results, aes(logFC, -log10(SpatialFDR))) +  
  geom_point() + 
  geom_hline(yintercept = 1)
scater::plotReducedDim(milo.all.seu, dimred = "UMAP", colour_by="cell.type.percise.new", text_size = 3, point_size=0.1)

milo.da_results$logFC <- -milo.da_results$logFC
miloR::plotNhoodGraphDA(milo.all.seu, milo.da_results, alpha=0.1,) +  
  scale_fill_gradient2(low="#4cc9f0",#修改颜色            
                       mid="#F2F2F2",
                       high="#FF5E5B",                  
                       name="log2FC",                   
                       limits=c(-5,5),                  
                       oob=squish)

milo.da_results <- annotateNhoods(milo.all.seu, milo.da_results, coldata_col = "cell.type.percise.new")
milo.da_results$cell.type.percise.new <- factor(milo.da_results$cell.type.percise.new, levels = rev(c('Cancer cells', 'CAFs', 'Cdh3 + Epithelial cells', 'Erbb4 + Epithelial cells', 'Prlr + Epithelial cells','Endothelial cells', 'Pericytes', 'MSCs','Macrophages','DCs', 'Monocytes', 'Neutrophils','T cells')))

plotbee <- function (da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL) 
{
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", 
           group.by, ")?")
    }
    if (is.numeric(da.res[, group.by])) {
    }
    da.res <- mutate(da.res, group_by = da.res[, group.by])
  }
  else {
    da.res <- mutate(da.res, group_by = "g1")
  }
  if (!is.factor(da.res[, "group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, 
                                               levels = unique(group_by)))
  }
  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods, ]
  }
  beeswarm_pos <- ggplot_build(da.res %>% mutate(is_signif = ifelse(SpatialFDR < 
                                                                      alpha, 1, 0)) %>% arrange(group_by) %>% ggplot(aes(group_by, 
                                                                                                                         logFC)) + geom_quasirandom())
  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y
  n_groups <- unique(da.res$group_by) %>% length()
  da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 
                                       1, 0)) %>% mutate(logFC_color = ifelse(is_signif == 
                                                                                1, logFC, NA)) %>% arrange(group_by) %>% mutate(Nhood = factor(Nhood, 
                                                                                                                                               levels = unique(Nhood))) %>% mutate(pos_x = pos_x, pos_y = pos_y) %>% 
    ggplot(aes(pos_x, pos_y, color = logFC_color)) + scale_color_gradient2() + 
    guides(color = "none") + xlab(group.by) + ylab("Log Fold Change") + 
    scale_x_continuous(breaks = seq(1, n_groups), labels = setNames(levels(da.res$group_by), 
                                                                    seq(1, n_groups))) + geom_point(size = .5) + coord_flip() + 
    theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0))
}

plotbee (milo.da_results, alpha = .1,
         group.by = "cell.type.percise.new") +  
  scale_color_gradient2(low="#4cc9f0",
                        mid="#F2F2F2",                  
                        high="#FF5E5B",                  
                        limits=c(-5,5),                  
                        oob=squish) + 
  labs(x="", y="Log2 Fold Change") +  
  theme_bw(base_size=10)+  
  theme(axis.text = element_text(colour = 'black')) + 
  scale_size_continuous(range = .1)


# 03.LTSR composition of cell types ----
library(magrittr)
library(lme4)
library(numDeriv)
## data perpare
Clean_sct.inte.rm.lpt$orig.ident <- factor(Clean_sct.inte.rm.lpt$orig.ident, levels = c('A3_T','A5_T','C2_T','C3_T'))
Clean_sct.inte.rm.lpt$condition <- factor(Clean_sct.inte.rm.lpt$condition, levels = c('Untreated','Targeted'))
cell.number <- FetchData(Clean_sct.inte.rm.lpt, 
                         vars = c("orig.ident", "cell.type.percise.new")) %>%
  dplyr::count(orig.ident, cell.type.percise.new) %>% 
  tidyr::spread(orig.ident, n) 
cell.number[is.na(cell.number)] <- 0

sample_ids <- colnames(cell.number)[-1]
cell_types <- cell.number$cell.type.percise.new
n_cells_per_sample <- colSums(cell.number[,-1])
n_var_cats <- 2 

sample_cats <- tibble(
  Sample_ID = sample_ids,
  Treatment = c(rep('Targeted',2),rep('Untreated',2)),
  Rep = c(rep(c('one','two'),2)),
  # Tumor.size.change = c(108.3/85.3, 152.8/99.7, 875.7/88.4, 891.4/113.2),
  cell.num = n_cells_per_sample
)

sample_num1_values<- rep(1,4)
obs_tbl <- data.frame(
  Sample_ID = rep(sample_ids, c(n_cells_per_sample)),
  Treatment = rep(sample_cats$Treatment, c(n_cells_per_sample)),
  Rep = rep(sample_cats$Rep, c(n_cells_per_sample)),
  # Tumor.size.change = rep(sample_cats$Tumor.size.change, n_cells_per_sample),
  Var_Num1 = rep(sample_num1_values, c(n_cells_per_sample))
)

obs_tbl$Cell_type <- c(rep(cell.number$cell.type.percise.new,c(cell.number$A3_T)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$A5_T)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$C2_T)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$C3_T))
)

## RUN LTSR 
source('./LTSR.raw.code.R')
results <- CellTypeCompositionAnalysis(obs_tbl, "Sample_ID", "Cell_type", c("Treatment",'Rep'), "Var_Num1")
ranef_tbl <- results$ranef
sdse_tbl <- results$sdse

vars1 <- list(Treatment = c('Untreated','Targeted'))
ranef_plot <- plot_ranef(ranef_tbl, vars = vars1, celltypes = cell_types, celltype_order = rev(cell_types),
                         maxFC = 1.5, LTSR2p = FALSE) + xlab('Condition')
sdse_plot <- plot_sdse(sdse_tbl, "Sample_ID", ci = 0.95, xlim = c(0, 1))

ranef_plot + sdse_plot

# 04.trans noise ----
library(MASS)
library(ggpubr)
celltypes <- unique(Clean_sct.inte.rm.lpt@meta.data$cell.type.percise.new)
celltypes <- celltypes[which(!is.na(celltypes))]
Idents(Clean_sct.inte.rm.lpt) <- Clean_sct.inte.rm.lpt$cell.type.percise.new
set.seed(1234)

getEuclideanDistance <- function(celltype, obj, assay, slot, ident1, ident2, group.by, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp <- subset(obj, cells = WhichCells(obj, idents = celltype))
  
  counts <- GetAssayData(object = tmp, slot = slot, assay = assay)
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) > 0
  expr <- counts[keep_genes, ]
  
  ifelse(min(table(tmp@meta.data[[group.by]])) > 300,
         expr <- expr[,c(rownames(tmp@meta.data[tmp@meta.data[[group.by]] == ident1,])[sample(1:nrow(tmp@meta.data[tmp@meta.data[[group.by]] == ident1,]),300)],
                         rownames(tmp@meta.data[tmp@meta.data[[group.by]] == ident2,])[sample(1:nrow(tmp@meta.data[tmp@meta.data[[group.by]] == ident2,]),300)])
         ],
         expr <- expr)
  tmp <- subset(tmp,cells = colnames(expr))
  
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  
  nsample <- min(table(tmp@meta.data[[group.by]])[c(ident1,ident2)])
  
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  print(nsample)
  ident2_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data[[group.by]] == ident2)], nsample)
  ident1_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data[[group.by]] == ident1)], nsample)
  ds_expr_r <- ds_expr[, c(ident1_r, ident2_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
  
  calcEuclDist <- function(matr, ident1, ident2){
    tmp <- data.matrix(sqrt(matr[genes, ident1]))
    mean <- rowMeans(sqrt(matr[genes, ident1]))
    d_ident1 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident1) <- ident1
    
    tmp <- data.matrix(sqrt(matr[genes, ident2]))
    mean <- rowMeans(sqrt(matr[genes, ident2]))
    d_ident2 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident2) <- ident2
    
    list(ident1 = d_ident1, ident2 = d_ident2)
  }
  ds <- calcEuclDist(matr = ds_expr_r, ident2 = ident2_r, ident1 = ident1_r)
  ds
}

res <- lapply(celltypes, function(x) getEuclideanDistance(x, 
                                                          obj = Clean_sct.inte.rm.lpt,
                                                          assay = 'SCT',
                                                          slot = 'counts',
                                                          group.by = 'condition',
                                                          ident1 = 'Untreated',
                                                          ident2 = 'Targeted',
                                                          lowcv = F))
names(res) <- celltypes
res.df <- data.frame(TN.value = unlist(do.call(c, res))) %>% 
  rownames_to_column(var = 'info') %>% 
  separate(col = info,
           into = c('Celltype','Condition','Sample_cells'),
           sep = '\\.',
           remove = T,
           extra = "merge")
res.df$Condition <- factor(ifelse(res.df$Condition == 'ident1','Untreated','Targeted'), levels = c('Untreated','Targeted'))
res.df$Celltype <- factor(res.df$Celltype, levels = c('Cancer cells', 'CAFs', 'Cdh3 + Epithelial cells', 'Erbb4 + Epithelial cells', 'Prlr + Epithelial cells','Endothelial cells', 'Pericytes', 'MSCs','Macrophages','DCs', 'Monocytes', 'Neutrophils','T cells'))

my_comparisons <- list(c("Untreated","Targeted"))
ggplot(res.df, aes(x=Condition, y=TN.value, fill=Condition)) + 
  geom_boxplot(alpha = .7) + 
  # geom_violin(alpha = .5) + 
  theme_bw() +
  labs(x = '', y = 'Transcriptional \nheterogeneity') +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#4cc9f0',"#FF5E5B")) +
  facet_wrap(~Celltype, scale="free",nrow = 2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

# TN.cts <- data.frame(div = 0, pvalue = 0)
# for (ct in celltypes) {
#   tmp <- res.df[res.df$Celltype == ct,]
#   div <- sum(tmp[tmp$Condition == 'Targeted',]$TN.value) / sum(tmp[tmp$Condition == 'Untreated',]$TN.value)
#   pvalue <- wilcox.test(tmp[tmp$Condition == 'Targeted',]$TN.value, tmp[tmp$Condition == 'Untreated',]$TN.value)
#   pvalue <- pvalue[['p.value']]
#   tmp.df <- data.frame(div = div, pvalue = pvalue)
#   TN.cts <- rbind(TN.cts,tmp.df)
# }
# TN.cts <- TN.cts[-1,]
# rownames(TN.cts) <- celltypes

# 05.dotplot of key genes -----
library(plot1cell)

complex_dotplot_single <- function (seu_obj, feature, celltypes = NULL, groups, splitby = NULL, 
                                    color.palette = NULL, font.size = 12, strip.color = NULL, 
                                    do.scale = T, scale.by = "radius") 
{
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(c("grey80", "lemonchiffon1", 
                                        "indianred1", "darkred"))(255)
  }
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (is.null(celltypes)) {
    celltypes <- levels(seu_obj)
  }
  if (length(groups) == 1) {
    groups_level <- levels(seu_obj@meta.data[, groups])
    if (is.null(groups_level)) {
      seu_obj@meta.data[, groups] <- factor(seu_obj@meta.data[, 
                                                              groups], levels = names(table(seu_obj@meta.data[, 
                                                                                                              groups])))
      groups_level <- levels(seu_obj@meta.data[, groups])
    }
    if (!is.null(splitby)) {
      if (is.null(levels(seu_obj@meta.data[, splitby]))) {
        seu_obj@meta.data[, splitby] <- factor(seu_obj@meta.data[, 
                                                                 splitby], levels = names(table(seu_obj@meta.data[, 
                                                                                                                  splitby])))
      }
      splitby_level <- levels(seu_obj@meta.data[, splitby])
      count_df <- extract_gene_count(seu_obj, features = feature, 
                                     cell.types = celltypes, meta.groups = c(groups, 
                                                                             splitby))
      count_df$new_group <- paste(count_df[, groups], 
                                  count_df[, "celltype"], count_df[, splitby], 
                                  sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by = "new_group")
      data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[1]]
                                                           }))
      data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                             split = "___"), FUN = function(x) {
                                                               x[[2]]
                                                             }))
      data_plot$splitby <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                            split = "___"), FUN = function(x) {
                                                              x[[3]]
                                                            }))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$splitby <- factor(data_plot$splitby, levels = splitby_level)
      data_plot$celltype <- factor(data_plot$celltype, 
                                   levels = rev(celltypes))
    }
    else {
      count_df <- extract_gene_count(seu_obj, features = feature, 
                                     cell.types = celltypes, meta.groups = groups)
      count_df$new_group <- paste(count_df[, groups], 
                                  count_df[, "celltype"], sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by = "new_group")
      data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[1]]
                                                           }))
      data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                             split = "___"), FUN = function(x) {
                                                               x[[2]]
                                                             }))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$celltype <- factor(data_plot$celltype, 
                                   levels = rev(celltypes))
    }
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 
                               2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    p <- ggplot(data_plot, aes(y = celltype, x = groups)) + 
      geom_tile(fill = "white", color = "white") + geom_point(aes(colour = avg.exp, 
                                                                  size = pct.exp)) + scale_color_gradientn(colours = color.palette) + 
      theme(panel.background = element_rect(fill = "white", 
                                            colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                          hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                            2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
            legend.text = element_text(size = (font.size - 
                                                 2)), legend.title = element_text(size = (font.size)), 
            strip.text = element_text(size = font.size), 
            legend.position = "right") + ylab("") + xlab("") + 
      ggtitle(feature)
    if (do.scale) {
      p = p + scale_size(range = c(0, 10))
    }
    else {
      if (max(data_plot$pct.exp) >= 20) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                        20))
      }
    }
    if (!is.null(splitby)) {
      p <- p + facet_wrap(~splitby, scales = "free_x")
      g <- change_strip_background(p, type = "top", strip.color = strip.color)
      print(grid.draw(g))
    }
    else {
      p
    }
  }
  else {
    gene_count <- extract_gene_count(seu_obj = seu_obj, 
                                     features = feature, cell.types = celltypes, meta.groups = c(groups, 
                                                                                                 splitby))
    allgroups <- c(groups, splitby)
    for (i in 1:length(allgroups)) {
      if (is.null(levels(seu_obj@meta.data[, allgroups[i]]))) {
        seu_obj@meta.data[, allgroups[i]] <- factor(seu_obj@meta.data[, 
                                                                      allgroups[i]], levels = names(table(seu_obj@meta.data[, 
                                                                                                                            allgroups[i]])))
      }
      group_level <- levels(seu_obj@meta.data[, allgroups[i]])
      gene_count[, allgroups[i]] <- factor(gene_count[, 
                                                      allgroups[i]], levels = group_level)
    }
    gene_count$celltype <- factor(gene_count$celltype, levels = celltypes)
    all_levels <- list()
    for (i in 1:length(groups)) {
      if (is.null(levels(seu_obj@meta.data[, groups[i]]))) {
        seu_obj@meta.data[, groups[i]] <- factor(seu_obj@meta.data[, 
                                                                   groups[i]], levels = names(table(seu_obj@meta.data[, 
                                                                                                                      groups[i]])))
      }
      group_level <- levels(seu_obj@meta.data[, groups[i]])
      all_levels[[i]] <- group_level
    }
    all_levels <- as.character(unlist(all_levels))
    data_plot <- list()
    for (i in 1:length(groups)) {
      count_df <- gene_count
      count_df$new_group <- paste(gene_count[, groups[i]], 
                                  gene_count[, "celltype"], sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      df1 <- merge(exp_df, pct_df, by = "new_group")
      df1$groupID <- groups[i]
      data_plot[[i]] <- df1
    }
    data_plot <- do.call("rbind", data_plot)
    data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                         split = "___"), FUN = function(x) {
                                                           x[[1]]
                                                         }))
    data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[2]]
                                                           }))
    data_plot$groups <- factor(data_plot$groups, levels = all_levels)
    data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    data_plot$groupID <- factor(data_plot$groupID, levels = groups)
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 
                               2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    if (is.null(splitby)) {
      p <- ggplot(data_plot, aes(y = celltype, x = groups)) + 
        geom_tile(fill = "white", color = "white") + 
        geom_point(aes(colour = avg.exp, size = pct.exp)) + 
        scale_color_gradientn(colours = color.palette) + 
        theme(panel.background = element_rect(fill = "white", 
                                              colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                            hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                              2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
              legend.text = element_text(size = (font.size - 
                                                   2)), legend.title = element_text(size = (font.size)), 
              strip.text = element_text(size = font.size), 
              legend.position = "right") + ylab("") + xlab("") + 
        ggtitle(feature) + facet_wrap(~groupID, scales = "free_x")
      if (do.scale) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        if (max(data_plot$pct.exp) >= 20) {
          p = p + scale_size(range = c(0, 10))
        }
        else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                          20))
        }
      }
      g <- change_strip_background(p, type = "top", strip.color = strip.color)
      print(grid::grid.draw(g))
    }
    else {
      df2 <- reshape2::melt(gene_count[, c(groups, splitby)], 
                            measure.vars = groups)
      df2 <- df2[!duplicated(df2$value), ]
      colnames(df2)[colnames(df2) == "value"] <- "groups"
      data_plot2 <- list()
      for (i in 1:length(groups)) {
        df3 <- data_plot[data_plot$groupID == groups[i], 
        ]
        df4 <- df2[df2$variable == groups[i], c("groups", 
                                                splitby[i])]
        colnames(df4)[2] <- "split"
        df5 <- merge(df3, df4, by = "groups")
        data_plot2[[i]] <- df5
      }
      data_plot2 <- do.call("rbind", data_plot2)
      fill_x1 <- grDevices::rainbow(length(groups), alpha = 0.5)
      fill_x2 <- list()
      for (i in 1:length(splitby)) {
        n_col <- unique(gene_count[, splitby[i]])
        fill_x2[[i]] <- (scales::hue_pal(l = 90))(length(n_col))
      }
      fill_x2 <- as.character(unlist(fill_x2))
      fill_x <- c(fill_x1, fill_x2)
      p <- ggplot(data_plot2, aes(y = celltype, x = groups)) + 
        geom_tile(fill = "white", color = "white") + 
        geom_point(aes(colour = avg.exp, size = pct.exp)) + 
        scale_color_gradientn(colours = color.palette) + 
        theme(panel.background = element_rect(fill = "white", 
                                              colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                            hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                              2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
              legend.text = element_text(size = (font.size - 
                                                   2)), legend.title = element_text(size = (font.size)), 
              strip.text = element_text(size = font.size), 
              legend.position = "right") + ylab("") + xlab("") + 
        ggtitle(feature) + facet_nested(~groupID + split, 
                                        scales = "free_x", strip = strip_nested(background_x = elem_list_rect(fill = fill_x)))
      if (do.scale) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        if (max(data_plot$pct.exp) >= 20) {
          p = p + scale_size(range = c(0, 10))
        }
        else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                          20))
        }
      }
      p
    }
  }
}
com.dot.new <- com.dot.new <- function (seu_obj, features, celltypes = NULL, groups, color.palette = NULL, 
                                        strip.color = NULL) 
{
  pb <- progress_bar$new(format = "  Ploting [:bar] :percent eta: :eta", 
                         clear = FALSE, total = length(features), width = 100)
  plot_list <- list()
  for (i in 1:length(features)) {
    pp <- invisible(complex_dotplot_single(seu_obj = seu_obj, 
                                           feature = features[i], groups = groups, celltypes = celltypes))
    pp <- pp$data
    pp$gene <- features[i]
    plot_list[[i]] <- pp
    pb$tick()
    Sys.sleep(1/length(features))
  }
  all_data <- do.call("rbind", plot_list)
  all_data$gene <- factor(all_data$gene, levels = rev(features))
  all_data$celltype <- factor(all_data$celltype, levels = levels(seu_obj))
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(c("grey80", "lemonchiffon1",
                                        "indianred1", "darkred"))(255)
  }
  p <- invisible(ggplot(all_data, aes(x = groups, y = gene)) + 
                   geom_tile(fill = "white", color = "white") + 
                   geom_point(aes(colour = avg.exp, size = pct.exp), alpha = 0.9) + 
                   scale_color_gradientn(colours = color.palette) + 
                   scale_size(range = c(0, 5)) +
                   theme(
                     panel.background = element_rect(fill = "white", colour = "black"), 
                     axis.text.x = element_text(angle = 45,hjust = 1),
                     axis.text.y = element_text(face = 'italic'),
                     plot.title = element_text(size = 10, hjust = 0.5,face = "bold"), 
                     axis.text = element_text(size = 12), 
                     axis.title = element_text(size = 8), 
                     legend.text = element_text(size = 8), 
                     legend.title = element_text(size = 12),
                     legend.position = "right", 
                     strip.text = element_text(size = 8, colour = "black",face = "bold")) + 
                   ylab("") + xlab("") + ggtitle("") + 
                   facet_wrap(~celltype, ncol = length(levels(seu_obj))))
  g <- change_strip_background(p, type = "top", strip.color = strip.color)
  print(grid.draw(g))
}
extract_gene_count <- function (seu_obj, features, cell.types = NULL, data.type = "data", 
                                meta.groups = NULL) 
{
  if (is.null(cell.types)) {
    cell.types = levels(seu_obj)
  }
  seu_obj@meta.data$celltype <- as.character(seu_obj@active.ident)
  if (is.null(meta.groups)) {
    meta.groups = colnames(seu_obj@meta.data)
  }
  if (!is.null(cell.types)) {
    new_seu <- subset(seu_obj, idents = cell.types)
  }
  feature_count <- Seurat::FetchData(new_seu, slot = data.type, 
                                     vars = c(features, meta.groups, "celltype"))
  umap_data <- data.frame(new_seu[["umap"]]@cell.embeddings)
  feature_count$UMAP1 <- umap_data$UMAP_1
  feature_count$UMAP2 <- umap_data$UMAP_2
  feature_count
}
change_strip_background <- function (ggplt_obj, type = "top", strip.color = NULL) 
{
  g <- ggplot_gtable(ggplot_build(ggplt_obj))
  if (type == "top") {
    strip_both <- which(grepl("strip-t", g$layout$name))
    fills <- strip.color
    if (is.null(fills)) {
      fills <- (scales::hue_pal(l = 90))(length(strip_both))
    }
  }
  else if (type == "right") {
    strip_both <- which(grepl("strip-r", g$layout$name))
    fills <- strip.color
    if (is.null(fills)) {
      fills <- (scales::hue_pal(l = 90))(length(strip_both))
    }
  }
  else {
    strip_t <- which(grepl("strip-t", g$layout$name))
    strip_r <- which(grepl("strip-r", g$layout$name))
    strip_both <- c(strip_t, strip_r)
    fills <- strip.color
    if (is.null(fills)) {
      fills <- c((scales::hue_pal(l = 90))(length(strip_t)), 
                 (scales::hue_pal(l = 90))(length(strip_r)))
    }
  }
  k <- 1
  for (i in strip_both) {
    j <- which(grepl("rect", g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  g
}

Idents(Clean_sct.inte.rm.lpt) <- Clean_sct.inte.rm.lpt$cell.type.percise.new
Clean_sct.inte.rm.lpt <- PrepSCTFindMarkers(Clean_sct.inte.rm.lpt)
rm.lpt.ct.deg <- FindAllMarkers(Clean_sct.inte.rm.lpt, only.pos = T,logfc.threshold = 0.25, min.pct = 0.2) %>% dplyr::filter(p_val_adj < 0.05)
rm.lpt.ct.deg.top10 <- rm.lpt.ct.deg %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
Clean_sct.inte.rm.lpt$cell.type.percise.new <- factor(Clean_sct.inte.rm.lpt$cell.type.percise.new, levels = c('Top2a + cancer cells', 'Padi4 + cancer cells', 'Cdh3 + Epithelial cells', 'Erbb4 + Epithelial cells', 'Prlr + Epithelial cells','Endothelial cells', 'Pericytes', 'MSCs','Macrophages','DCs', 'Monocytes', 'Neutrophils','T cells'))
Idents(Clean_sct.inte.rm.lpt) <- Clean_sct.inte.rm.lpt$cell.type.percise.new
com.dot.new(Clean_sct.inte.rm.lpt, feature = c(
                                               'Top2a','Brca1','Kif15',
                                               'Rerg','Padi4','Cdon',
                                               'Trp63','Cdh3','Brinp1',
                                               'Erbb4','Apod','Epcam',
                                               'Prlr','Ptn',
                                               'Pecam1','Kdr','Tek',
                                               'Rgs5','Notch3','Abcc9',
                                               'Twist2','Col3a1','Slit3',
                                               'Mrc1','Ms4a7','C1qc',
                                               'Flt3','Btla','Itgae',
                                               'Ccr2','Plac8','F10',
                                               'S100a9','Cxcl2','Slfn4',
                                               'Cd3e','Cd3g','Cd8b1')
            ,groups = "condition",strip.color = c(use.cols,npg.cols)[1:15])

# 06.umap of key genes ----
p1 <- FeaturePlot(Clean_sct.inte.rm.lpt, features = 'Top2a',split.by = 'condition',order = T ,combine = F, cols = c('gray90','red3'))
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

# 07.Augur:prioritize the cell types most responsive to biological perturbations -----
library(Augur)
augur <- calculate_auc(Clean_sct.inte.rm.lpt, cell_type_col = 'cell.type.percise.new', label_col = 'condition', n_threads = 50)
head(augur$AUC,5)

plot_umap(augur,
          Clean_sct.inte.rm.lpt,
          mode = "default",
          reduction = "umap",
          palette = "inferno", #  "viridis", "plasma", "magma", "inferno"
          # augur_mode = "default",
          cell_type_col = "cell.type.percise.new")

plot_loli <- function (augur) {
  aucs = augur$AUC
  size_sm = 10
  size_lg = 10
  range = range(aucs$auc)
  expand = abs(diff(range)) * 0.1
  p = aucs %>% ggplot(aes(x = reorder(cell_type, auc), y = auc)) + 
    geom_hline(aes(yintercept = 0.5), linetype = "dotted", 
               size = 0.8) +
    geom_point(size = 2) + 
    geom_text(aes(label = format(auc, digits = 3), 
                  y = ifelse(auc < 0.5, 0.5, auc)), 
              size = 4, nudge_y = expand, hjust = 0.5) + 
    geom_segment(aes(xend = cell_type, yend = 0.5)) + 
    scale_y_continuous("AUC", limits = c(min(range[1] - expand, 0.5), range[2] + expand * 1.5)) + 
    coord_flip() + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = size_sm), 
          axis.text.y = element_text(size = size_sm + 2),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_blank(), panel.grid = element_blank(), 
          strip.text = element_text(size = size_lg), strip.background = element_blank(),
          axis.line.y = element_blank(), axis.line.x = element_blank(), 
          legend.position = "top", legend.text = element_text(size = size_sm), 
          legend.title = element_text(size = size_sm), 
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          legend.background = element_blank(), 
          plot.title = element_text(size = size_lg, hjust = 0.5))
  p
}

plot_loli(augur)+
  geom_segment(aes(xend=cell_type,yend=0.5),size=1)+
  geom_point(size=3,aes(color=cell_type))+
  scale_color_manual(values = c('Top2a + cancer cells' = use.cols[1],
                                'Padi4 + cancer cells' = use.cols[2],
                                'Cdh3 + Epithelial cells' = use.cols[3], 
                                'Erbb4 + Epithelial cells' = use.cols[4],
                                'Prlr + Epithelial cells' = use.cols[5],
                                'Endothelial cells' = use.cols[6],
                                'Pericytes' = use.cols[7],
                                'MSCs' = use.cols[8],
                                'Macrophages' = use.cols[9],
                                'DCs' = use.cols[10],
                                'Monocytes' = use.cols[11],
                                'Neutrophils' = use.cols[12],
                                'T cells' = use.cols[13])) + 
  theme(legend.position = "")

## scDist
library(scDist)
sim <- list(Y = Clean_sct.inte.rm.lpt@assays$SCT@data ,
            meta.data = Clean_sct.inte.rm.lpt@meta.data %>% as.data.frame())
scdist <- scDist(normalized_counts = sim$Y, 
                   meta.data = sim$meta.data, 
                   d = 20, 
                   fixed.effects = "condition", 
                   random.effects = 'orig.ident', 
                   clusters="cell.type.percise.new" 
)
DistPlot(scdist) + theme_bw()

# 08.EMT score -----
orth <- rio::import('hsa2mmu.110.txt')
orth <- orth[orth$`Mouse homology type` == 'ortholog_one2one',] %>% na.omit()
EMT.genenset <- rio::import('EMT.genes.txt')
EMT.genes <- EMT.genenset %>% left_join(y = orth[,1:2], by = c('GeneSymbol' = 'Gene name')) %>% na.omit()
EMT.genes.use <- list(c(EMT.genes$`Mouse gene name`))
Clean_sct.inte.rm.lpt <- AddModuleScore(Clean_sct.inte.rm.lpt, features = EMT.genes.use,name = 'EMT')

p1 <- FeaturePlot(Clean_sct.inte.rm.lpt, features = 'EMT1', split.by = 'condition', order = T, min.cutoff = 'q50', cols = c('grey90', "red3"), combine = F)
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

EMT.score <- Clean_sct.inte.rm.lpt@meta.data[,c('cell.type.percise.new','EMT1','condition')] %>% na.omit()
EMT.score$condition <- factor(EMT.score$condition, levels = c('Untreated','Targeted'))

ggplot(EMT.score, aes(x=condition, y=EMT1, fill=condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'EMT score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#4cc9f0',"#FF5E5B")) +
  facet_wrap(~cell.type.percise.new, scale="free",nrow = 2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")
median(EMT.score[EMT.score$condition == 'Targeted',]$EMT1)

ggdensity(EMT.score, 
          x = "EMT1",
          add = "median", rug = F,
          color = "condition", fill = "condition",
          palette = c('#4cc9f0',"#FF5E5B")) +
  ylab("Density") + 
  xlab('EMT score') + 
  ggtitle('EMT') + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) 

FeaturePlot(Clean_sct.inte.rm.lpt, features = c('Twist1','Snai1','Zeb1','Vim','Sox2'),split.by = 'condition',order = T,min.cutoff = 0.4)

# 09.infercnv -----
library(infercnv)
library(AnnoProbe)
dat <- GetAssayData(Clean_sct.inte.rm.lpt,layer = 'counts',assay = 'SCT')
# groupinfo <- data.frame(v1 = colnames(dat),
#                         v2 = paste0(Clean_sct.inte.rm.lpt@meta.data$condition,"_",Clean_sct.inte.rm.lpt@meta.data$cell.type.percise.new))
groupinfo <- data.frame(v1 = colnames(dat),
                        v2 = Clean_sct.inte.rm.lpt$cell.type.percise.new)

geneInfor <- annoGene(rownames(dat),"SYMBOL","mouse") 
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrM", "chrX", "chrY"), ]
geneInfor$chr_num <- as.numeric(sub("chr", "", geneInfor$chr))
colnames(geneInfor)
geneInfor <- geneInfor[with(geneInfor,order(chr_num,start)),c(1,4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]

dat <- dat[rownames(dat) %in% geneInfor[,1],]
dat <- dat[match(geneInfor[,1],rownames(dat)),]

expFile <- 'expFile.txt'
colnames(dat) <- gsub("-", "_", colnames(dat))
write.table(dat, file = expFile, sep = '\t', quote = F)

groupFiles <- 'groupFiles.txt'
groupinfo$v1 <- gsub("-", "_", groupinfo$v1)
write.table(groupinfo,file = groupFiles, sep = '\t',
            quote = F, col.names = F, row.names = F)
geneFile <- 'geneFile.txt'
write.table(geneInfor, file = geneFile, sep = '\t',
            quote = F, col.names = F, row.names = F)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = expFile,
                                     annotations_file = groupFiles,
                                     delim = "\t",
                                     gene_order_file = geneFile,
                                     ref_group_names = c("Endothelial cells","MSCs")
)

infercnv_obj2 <- infercnv::run(infercnv_obj,
                               cutoff =  0.1, 
                               out_dir = "infercnv_output.group", 
                               cluster_by_groups =  T,
                               hclust_method = "ward.D2",
                               # analysis_mode = "samples", 
                               denoise = TRUE, 
                               HMM = F,  
                               plot_steps = F, 
                               leiden_resolution = "auto",
                               num_threads = 80, 
                               output_format = "pdf",
                               write_expr_matrix=T
)

# cnv vlnplot
data <- read.table("infercnv_output.group/infercnv.observations.txt", header=T)
expr <- data %>% as.matrix()
expr.scale <- scale(t(expr))
tmp1 <- sweep(expr.scale, 2, apply(expr.scale, 2, min),'-')
tmp2 <- apply(expr.scale, 2, max) - apply(expr.scale,2,min)
expr_1 <- t(2*sweep(tmp1, 2, tmp2, "/")-1)
cnv_score <- as.data.frame(colSums(expr_1 * expr_1))
colnames(cnv_score)="cnv_score"
cnv_score <- rownames_to_column(cnv_score, var='cell')
colnames(groupinfo) <- c('cell','cell.type')
test <- cnv_score %>% left_join(y = groupinfo, by = 'cell') %>% na.omit()
ggplot2::ggplot(test,aes(x=cell.type,y=cnv_score))+
  geom_violin(aes(fill=cell.type),cex=0.5)+
  scale_fill_manual(values = use.cols)+
  # geom_boxplot(width=0.1,cex=1.2)+
  # theme_classic(base_size = 20)+
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(color = 'black',angle = 45, vjust = 1, hjust = 1),
        legend.position = 'none') + 
  xlab("") + 
  ylab("CNV score")

# cnv in Padi4+ & Top2a+ cells
choose.cnv <- test[test$cell.type == 'Padi4+ cancer cells' | test$cell.type == 'Top2a+ cancer cells',]
choose.cnv <- choose.cnv %>% left_join(y = Clean_sct.inte.rm.lpt@meta.data[,c('cells','condition')], by = c('cell' = 'cells'))
median(choose.cnv[choose.cnv$cell.type == 'Top2a+ cancer cells' & choose.cnv$condition == 'WT',]$cnv_score)
median(choose.cnv[choose.cnv$cell.type == 'Top2a+ cancer cells' & choose.cnv$condition == 'Targeted',]$cnv_score)
median(choose.cnv[choose.cnv$cell.type == 'Padi4+ cancer cells' & choose.cnv$condition == 'WT',]$cnv_score)
median(choose.cnv[choose.cnv$cell.type == 'Padi4+ cancer cells' & choose.cnv$condition == 'Targeted',]$cnv_score)

my_comparisons <- list(c("Untreated","Targeted"))
ggplot(choose.cnv, aes(x=condition, y=cnv_score, fill=condition)) + 
  geom_boxplot(alpha = .7) + 
  # geom_violin(alpha = .5) + 
  theme_bw() +
  labs(x = '', y = 'CNV score') +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#4cc9f0',"#FF5E5B")) +
  facet_wrap(~cell.type, scale="free",nrow = 1) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

plot.cnv.seu <- subset(Clean_sct.inte.rm.lpt, cells = colnames(data))
plot.cnv.seu$cnv.score <- scale(colMeans(data))
FeaturePlot(plot.cnv.seu, features = 'cnv.score', cols = c('#4cc9f0',"#FF5E5B"))

# 10.cytotrace -----
library(CytoTRACE2)

cytotrace2_result <- cytotrace2(Clean_sct.inte.rm.lpt,is_seurat = T,ncores = 80)
annotation <- data.frame(phenotype = Clean_sct.inte.rm.lpt@meta.data$cell.type.percise.new) %>% 
  set_rownames(., colnames(Clean_sct.inte.rm.lpt))

plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = T)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
(p1+p2+p3+p4) + plot_layout(ncol = 2)

my_comparisons <- list(c("WT","Targeted"))
ggboxplot(cytotrace2_result@meta.data, x="condition", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",
                fill="condition",
                palette = "npg",
                xlab = F, 
                bxp.errorbar=T,
                bxp.errorbar.width=0.5, 
                size=.1, 
                outlier.shape=NA,
                legend = "right",
                alpha = 0.8) + 
  ylab('Potency score')  + 
  scale_fill_manual(values = c('#4cc9f0',"#FF5E5B")) +
  theme_bw() + 
  # ylim(0,0.5) + 
  facet_wrap(~ cell.type.percise.new,nrow = 2) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),legend.position = 'none') +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

median(cytotrace2_result@meta.data[cytotrace2_result@meta.data$cell.type.percise.new == Top2a+ cancer cells' & cytotrace2_result@meta.data$condition == 'WT',]$CytoTRACE2_Score)
median(cytotrace2_result@meta.data[cytotrace2_result@meta.data$cell.type.percise.new == 'Top2a+ cancer cells' & cytotrace2_result@meta.data$condition == 'Targeted',]$CytoTRACE2_Score)
median(cytotrace2_result@meta.data[cytotrace2_result@meta.data$cell.type.percise.new == 'Padi4+ cancer cells' & cytotrace2_result@meta.data$condition == 'WT',]$CytoTRACE2_Score)
median(cytotrace2_result@meta.data[cytotrace2_result@meta.data$cell.type.percise.new == 'Padi4+ cancer cells' & cytotrace2_result@meta.data$condition == 'Targeted',]$CytoTRACE2_Score)
# 11.T cells -----

## 01.UMAP
T.inte$cell.type.new <- factor(T.inte$cell.type.new, levels = c('Cd8 + T', 'Cd4 + T', 'Memory T','Treg','NKT'))
p1 <- DimPlot(T.inte, group.by = 'cell.type.new', cols = c( npg.cols[1:8]), label = F) + ggtitle("") + NoAxes()
p2 <- DimPlot(T.inte, group.by = 'cell.type.new', cols = c(npg.cols), label = T,repel = T) + ggtitle("") + NoAxes() + NoLegend()
p3 <- DimPlot(T.inte, group.by = 'cell.type.new', cols = c(npg.cols), label = T,repel = T) + ggtitle("") + NoAxes()
p1 + p2 +p3

DimPlot(T.inte, group.by = 'cell.type.new', cols = c( npg.cols[1:8]), label = F,split.by = 'condition') + ggtitle("") + NoAxes() 

## 02.dotplot
T.ct.markers <- FindAllMarkers(T.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .2) %>% dplyr::filter(p_val_adj < 0.05)
T.ct.markers.top10 <- T.ct.markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
T.inte$condition <- factor(T.inte$condition, levels = c('Untreated', 'Targeted'))
Idents(T.inte) <- T.inte$cell.type.percise.new
com.dot.new(T.inte, feature = c('Cd8a','Filip1','Srgap3',
                                'Cd4', 'Emb','Slamf6',
                                'Sidt1','Bcl2','Ccl5',
                                'Foxp3','Nckap5','Atrnl1',
                                'Ncr1','Aoah','Gas7')
  ,groups = "condition",strip.color = c(npg.cols))

## 03.milo fail because cell number insufficient
# milo.T.seu <- subset(T.inte, cells = c(sample( rownames(T.inte@meta.data[T.inte$condition == 'Untreated',]), 332), rownames(T.inte@meta.data[T.inte$condition == 'Targeted',] )))
# milo.T.seu <- as.SingleCellExperiment(milo.T.seu,assay = 'SCT') %>%
#   Milo() %>%
#   miloR::buildGraph(k = 20, d = 30) %>% 
#   makeNhoods(
#     prop = 0.5,                  
#     k = 20,         
#     d=30,                    
#     refined = T)
# 
# milo.T.seu <- countCells(milo.T.seu, 
#                            meta.data = data.frame(colData(milo.T.seu)),                     
#                            sample="orig.ident")
# milo.traj_design <- data.frame(colData(milo.T.seu))[,c("orig.ident", "condition")]
# milo.traj_design$orig.ident <- as.factor(milo.traj_design$orig.ident)
# milo.traj_design <- distinct(milo.traj_design)
# rownames(milo.traj_design) <- milo.traj_design$orig.ident
# 
# milo.T.seu <- calcNhoodDistance(milo.T.seu, d=30)
# milo.da_results <- testNhoods(milo.T.seu,                 
#                               design = ~ condition,            
#                               design.df = milo.traj_design)
# milo.T.seu <- buildNhoodGraph(milo.T.seu)
# milo.da_results$logFC <- -milo.da_results$logFC
# miloR::plotNhoodGraphDA(milo.T.seu, milo.da_results, alpha=0.1,) +
#   scale_fill_gradient2(low="#4cc9f0",
#                        mid="#F2F2F2",
#                        high="#FF5E5B",
#                        name="log2FC",
#                        limits=c(-1.5,1.5),
#                        oob=squish)
# # 
# milo.da_results <- annotateNhoods(milo.T.seu, milo.da_results, coldata_col = "cell.type.new")
# milo.da_results$cell.type.new <- factor(milo.da_results$cell.type.new, levels = rev(c('Cd8 + T', 'Cd4 + T', 'Memory T','Treg','NKT')))
# plotDAbeeswarm (milo.da_results, alpha = 1,
#          group.by = "cell.type.new") +
#   scale_color_gradient2(low="#4cc9f0",
#                         mid="#F2F2F2",
#                         high="#FF5E5B",
#                         limits=c(-1.5,1.5),
#                         oob=squish) +
#   labs(x="", y="Log2 Fold Change") +
#   theme_bw(base_size=10)+
#   theme(axis.text = element_text(colour = 'black')) +
#   scale_size_continuous(range = .1)

## 04.LTSR
cell.number <- FetchData(T.inte, 
                         vars = c("orig.ident", "cell.type.new")) %>%
  dplyr::count(orig.ident, cell.type.new) %>% 
  tidyr::spread(orig.ident, n) 
cell.number[is.na(cell.number)] <- 0

sample_ids <- colnames(cell.number)[-1]
cell_types <- cell.number$cell.type.new
n_cells_per_sample <- colSums(cell.number[,-1])
n_var_cats <- 2 

sample_cats <- tibble(
  Sample_ID = sample_ids,
  Treatment = c(rep('Targeted',2),rep('Untreated',2)),
  Rep = c(rep(c('one','two'),2)),
  cell.num = n_cells_per_sample
)

sample_num1_values<- rep(1,4)
obs_tbl <- data.frame(
  Sample_ID = rep(sample_ids, c(n_cells_per_sample)),
  Treatment = rep(sample_cats$Treatment, c(n_cells_per_sample)),
  Rep = rep(sample_cats$Rep, c(n_cells_per_sample)),
  Var_Num1 = rep(sample_num1_values, c(n_cells_per_sample))
)

obs_tbl$Cell_type <- c(rep(cell.number$cell.type.new,c(cell.number$A3_T)),
                       rep(cell.number$cell.type.new,c(cell.number$A5_T)),
                       rep(cell.number$cell.type.new,c(cell.number$C2_T)),
                       rep(cell.number$cell.type.new,c(cell.number$C3_T))
)

results <- CellTypeCompositionAnalysis(obs_tbl, "Sample_ID", "Cell_type", c("Treatment",'Rep'), "Var_Num1")
ranef_tbl <- results$ranef
sdse_tbl <- results$sdse

vars1 <- list(Treatment = c('Untreated','Targeted'))
ranef_plot <- plot_ranef(ranef_tbl, 
                         vars = vars1, 
                         celltypes = cell_types,
                         celltype_order = rev(cell_types),
                         maxFC = 0.9, LTSR2p = FALSE) + xlab('Condition')
sdse_plot <- plot_sdse(sdse_tbl, "Sample_ID", ci = 0.95, xlim = c(0, 1))
ranef_plot + sdse_plot

## 05.cytotrace
my_comparisons <- list(c("WT","Targeted"))
cytotrace2.T <- cytotrace2(T.inte,is_seurat = T,ncores = 10)
annotation <- data.frame(phenotype = T.inte@meta.data$cell.type.new) %>% 
  set_rownames(., colnames(T.inte))

plots <- plotData(cytotrace2_result = cytotrace2.T, 
                  annotation = annotation, 
                  is_seurat = T)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
(p1+p2+p3+p4) + patchwork::plot_layout(ncol = 2)

ggboxplot(cytotrace2.T@meta.data, x="condition", y="CytoTRACE2_Score", width = 0.6, 
          color = "black",
          fill="condition",
          palette = "npg",
          xlab = F, 
          bxp.errorbar=T,
          bxp.errorbar.width=0.5,
          size=.1, 
          outlier.shape=NA, 
          legend = "right",
          alpha = 0.8) + 
  ylab('Potency score')  + 
  scale_fill_manual(values = c('#4cc9f0',"#FF5E5B")) +
  theme_bw() + 
  # ylim(0,0.5) + 
  facet_wrap(~ cell.type.new,nrow = 1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),legend.position = 'none') +
  stat_compare_means(comparisons = my_comparisons,              
                     method = "wilcox.test")
median(cytotrace2.T@meta.data[cytotrace2.T@meta.data$cell.type.new == 'Cd8 + T' & cytotrace2.T@meta.data$condition == 'WT',]$CytoTRACE2_Score)
median(cytotrace2.T@meta.data[cytotrace2.T@meta.data$cell.type.new == 'Cd8 + T' & cytotrace2.T@meta.data$condition == 'Targeted',]$CytoTRACE2_Score)
median(cytotrace2.T@meta.data[cytotrace2.T@meta.data$cell.type.new == 'Cd4 + T' & cytotrace2.T@meta.data$condition == 'WT',]$CytoTRACE2_Score)
median(cytotrace2.T@meta.data[cytotrace2.T@meta.data$cell.type.new == 'Cd4 + T' & cytotrace2.T@meta.data$condition == 'Targeted',]$CytoTRACE2_Score)

## 06.augur
T.augur <- calculate_auc(T.inte, cell_type_col = 'cell.type.new', label_col = 'condition', n_threads = 50)
head(T.augur$AUC,5)

plot_umap(T.augur,
          T.inte,
          mode = "default",
          reduction = "umap",
          palette = "inferno", #  "viridis", "plasma", "magma", "inferno"
          # augur_mode = "default", 
          cell_type_col = "cell.type.new")

plot_loli(T.augur)+
  geom_segment(aes(xend=cell_type,yend=0.5),size=1)+
  geom_point(size=3,aes(color=cell_type))+
  scale_color_manual(values = c('Cd8 + T' = npg.cols[1],
                                'Cd4 + T' = npg.cols[2],
                                'Memory T' = npg.cols[3], 
                                'Treg' = npg.cols[4],
                                'NKT' = npg.cols[5])) + 
  theme(legend.position = "")

### 06.1. scDist
library(scDist)
sim <- list(Y = T.inte@assays$SCT@data ,
            meta.data = T.inte@meta.data %>% as.data.frame())
T.scdist <- scDist(normalized_counts = sim$Y,
                          meta.data = sim$meta.data,
                          d = 20,
                          fixed.effects = "condition",
                          random.effects = 'orig.ident',
                          clusters="cell.type.new"
)
DistPlot(T.scdist) + theme_bw()

## 07.T deg & pathways
library(clusterProfiler)
comlist <- t(combn(unique(T.inte$condition), 2))

for (ct in unique(T.inte$cell.type.new)) {
  print(ct)
  Idents(T.inte) <- 'cell.type.new'
  sub.seu.1 <- subset(T.inte, idents = ct)
  for (gp in nrow(comlist)) {
    treat <- comlist[gp,1]
    ctrl <- comlist[gp,2]
    DEGs <- FindMarkers(sub.seu.1,
                        ident.1 = treat,
                        ident.2 = ctrl, 
                        logfc.threshold = 0.25,
                        recorrect_umi = F,
                        min.pct = .2,
                        group.by = 'condition') %>% dplyr::filter(p_val < 0.01) %>% rownames_to_column(var = 'gene')
    DEGs.name <- paste('T.cells.',ct,treat,ctrl,sep = '_')
    assign(DEGs.name,DEGs)
  }
}
rio::export(`T.cells._Cd8 + T_Targeted_WT`, file = './new.analysis/add.11.8/T.cells._Cd8 + T_Targeted_WT.xlsx')

hub<-AnnotationHub::AnnotationHub()
AnnotationHub::query(hub, "Musculus.v110")
ensdb110<-hub[["AH113713"]]

plot.list <- list()
for (ct in levels(T.inte$cell.type.new)) {
  tmp <- base::get(paste('T.cells.',ct,'Targeted_WT',sep = "_"))
  gene2entrzid <- bitr(tmp[tmp$avg_log2FC > 0,]$gene, fromType = 'GENENAME', toType = "ENTREZID", OrgDb = ensdb110)
  erich.go.BP <- enrichGO(gene=gene2entrzid$ENTREZID,
                          # OrgDb = ,
                          'org.Mm.eg.db',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP")
  er.plot <- dotplot(erich.go.BP,showCategory = 15) + ggtitle(paste0(ct,'_Targeted up'))
  erich.name <- paste0(ct,'_plot')
  # assign(erich.name, er.plot)
  plot.list[[erich.name]] <- er.plot
}

patchwork::wrap_plots(c(plot.list),nrow = 2)

plot.list.down <- list()
for (ct in levels(T.inte$cell.type.new)) {
  tmp <- base::get(paste('T.cells.',ct,'Targeted_Untreated',sep = "_"))
  gene2entrzid <- bitr(tmp[tmp$avg_log2FC < 0,]$gene, fromType = 'GENENAME', toType = "ENTREZID", OrgDb = ensdb110)
  erich.go.BP <- enrichGO(gene=gene2entrzid$ENTREZID,
                          # OrgDb = ,
                          'org.Mm.eg.db',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP")
  er.plot <- dotplot(erich.go.BP,showCategory = 15) + ggtitle(paste0(ct,'_Targeted down'))
  erich.name <- paste0(ct,'_plot')
  # assign(erich.name, er.plot)
  plot.list.down[[erich.name]] <- er.plot
}

patchwork::wrap_plots(c(plot.list.down),nrow = 2)

### T genes analysis & moudle score
exhaustion.gene <- list(c('Havcr2','Pdcd1','Tigit','Lag3','Ctla4','Rbpj','Vcam1','Gzmb','Tox'))
effect.mem.gene <- list(c('Prf1','Ifng','Ccl4','Gzmk','Gzma','Cd44','Dusp2','Klrd1','Ctsw','Bcl2'))
ifn.gene <- list(c('Stat1','Ifit1','Isg15','Ccr1'))
exhau_ifn.gene <- list(c('Hspa1a','Bag3','Dnajb1'))

Idents(T.inte) <- Clean_sct.inte.rm.lpt$cell.type.percise.new
com.dot.new(T.inte, feature = c(exhaustion.gene[[1]]),groups = "condition",strip.color = c(use.cols,npg.cols)[1:15])
com.dot.new(T.inte, feature = c(effect.mem.gene),groups = "condition",strip.color = c(use.cols,npg.cols)[1:15])
com.dot.new(T.inte, feature = c(ifn.gene),groups = "condition",strip.color = c(use.cols,npg.cols)[1:15])
com.dot.new(T.inte, feature = c('Gzmk','Gzma','Gzmf','Gzmc','Gzmd','Gzme'),groups = "condition",strip.color = c(use.cols,npg.cols)[1:15])

T.inte <- AddModuleScore(T.inte,features = exhaustion.gene,name = 'exhaustion.score')
T.inte <- AddModuleScore(T.inte,features = effect.mem.gene,name = 'effect.mem.score')
T.inte <- AddModuleScore(T.inte,features = ifn.gene,name = 'ifn.score')
T.inte <- AddModuleScore(T.inte,features = exhau_ifn.gene,name = 'exhau_ifn.score')

com.dot.new(T.inte, feature = c('exhaustion.score1','effect.mem.score1','ifn.score1','exhau_ifn.score1')
            ,groups = "condition",strip.color = c(npg.cols)[1:15])

FeaturePlot(T.inte, 
            features = c('Gzmk','Gzma','Gzmh'),split.by = 'condition',
            order = T)
VlnPlot(T.inte, 
        features = c('H2-K1'),split.by = 'condition')

eff.mem.score <- T.inte@meta.data[,c('cell.type.new','effect.mem.score1','condition')] %>% na.omit()
eff.mem.score$condition <- factor(eff.mem.score$condition, levels = c('WT','Targeted'))

ggplot(eff.mem.score, aes(x=condition, y=effect.mem.score1, fill=condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'Effect memory score') +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#4cc9f0',"#FF5E5B")) +
  facet_wrap(~cell.type.new, scale="free",nrow = 1) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

median(eff.mem.score[eff.mem.score$cell.type.new == 'Cd8 + T' & eff.mem.score$condition == 'WT',]$effect.mem.score1)
median(eff.mem.score[eff.mem.score$cell.type.new == 'Cd8 + T' & eff.mem.score$condition == 'Targeted',]$effect.mem.score1)
median(eff.mem.score[eff.mem.score$cell.type.new == 'Cd4 + T' & eff.mem.score$condition == 'WT',]$effect.mem.score1)
median(eff.mem.score[eff.mem.score$cell.type.new == 'Cd4 + T' & eff.mem.score$condition == 'Targeted',]$effect.mem.score1)

### T. splice mRNA
t.sp <- sp_unsp[sp_unsp$cell.type.percise.new == 'T cells',]
t.meta <- T.inte@meta.data[,c(17,26)]
t.sp <- t.sp %>% left_join(y = t.meta, by = c( 'my_sample' = 'cells' ))
my_comparisons <- list(c('Untreated','Targeted'))
ggplot(t.sp, aes(x=condition, y=unsp_ratio, fill=condition)) + 
  geom_boxplot(alpha = .7) + 
  # geom_violin(alpha = .5) + 
  theme_bw() +
  labs(x = '', y = 'Unspliced RNA ratio') +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) + 
  scale_fill_manual(values = c('#4cc9f0',"#FF5E5B")) +
  facet_wrap(~cell.type.new, scale="free",nrow = 1) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")


# 12.Padi4 + cancer cells ----

## integrate
CAF.sub <- subset(Clean_qc.merge.filtered, cells = rownames(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$cell.type.percise.new == 'Padi4 + cancer cells',]))

CAF.inte <- CAF.sub %>%
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)

DimPlot(CAF.inte, reduction = "umap", label = T,group.by = 'condition')

## non integrate
CAF.merge <- CAF.sub %>%
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  FindNeighbors( reduction = "pca", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "pca", 
           dims = 1:30)
CAF.merge[['RNA']] <- JoinLayers(CAF.merge[['RNA']])

DimPlot(CAF.merge, reduction = "umap", label = T,group.by = 'condition') + ggtitle('CAFs')
FeaturePlot(CAF.merge, features = 'Mki67')

# 13.Top2a + Cancer cells ----

## integrate
Cancer.sub <- subset(Clean_qc.merge.filtered, cells = rownames(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$cell.type.percise.new == 'Top2a + cancer cells',]))

Cancer.inte <- Cancer.sub %>%
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)

DimPlot(Cancer.inte, reduction = "umap", label = T,group.by = 'condition')

## non integrate
Cancer.merge <- Cancer.sub %>%
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  FindNeighbors( reduction = "pca", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "pca", 
           dims = 1:30)
Cancer.merge[['RNA']] <- JoinLayers(Cancer.merge[['RNA']])

DimPlot(Cancer.merge, reduction = "umap", label = T,group.by = 'condition') + ggtitle('Cancers')
FeaturePlot(Cancer.merge, features = c('Mki67','Top2a'))

# 14.velocyte -----
library(velocyto.R)
library(tidyverse)
loom_dir <- '/data/02.project/00.other.lab/10.LHH.sc.tumor/01.velocity/03.loom.result/'
merge_rds <- "/data/02.project/00.other.lab/10.LHH.sc.tumor/01.velocity/CAF.merge.rds"
my_condition_all <- c("Targeted", "WT")
my_pdf <- paste0('../lhh.gland/new.analysis/add.11.8/CAF.velocity.pdf')

rawdata <- read_rds(merge_rds)
Idents(rawdata) <- rawdata$condition
DefaultAssay(rawdata) <- 'RNA'
#rawdata <- rawdata[, sample(1:ncol(rawdata), 2000)]
# 
all_emat <- list()
all_nmat <- list()
my_seq <- 1
for (ii in list.files(loom_dir, full.names = T)) {
  ii <- list.files(loom_dir, full.names = T)
  ldat <- read.loom.matrices(ii[1])
  emat <- ldat$spliced
  nmat <- ldat$unspliced
  
  old_name <- colnames(emat)
  # old_name <- 'sorted_A3_T_S1_L001_Aligned_ATNK4:AACGAGTATCGCTCCGATAGx'
  my_sample <- str_split(old_name, ':', simplify = T)[, 1] %>% 
    str_remove('sorted_') %>% str_remove('........_Aligned_.....') %>% unique()
  my_cell <- str_split(old_name, ':', simplify = T)[, 2] %>% str_remove('x')
  
  new_name <- paste0(my_sample, "_",my_cell)
  
  colnames(emat) <- new_name
  colnames(nmat) <- new_name
  emat <- emat[rownames(emat)%in%rownames(rawdata), colnames(emat) %in% colnames(rawdata)]
  nmat <- nmat[rownames(nmat)%in%rownames(rawdata), colnames(nmat) %in% colnames(rawdata)]
  
  all_emat[[my_seq]] <- emat
  all_nmat[[my_seq]] <- nmat
  my_seq <- my_seq+1
  print(my_sample)
  print(dim(emat))
}

all_emat <- as.data.frame(all_emat) %>% as.matrix()
all_nmat <- as.data.frame(all_nmat) %>% as.matrix()

rawdata <- rawdata[, colnames(all_nmat)]

emb <- rawdata@reductions$umap@cell.embeddings
cell.dist <- as.dist(1-cor(t(rawdata@reductions$pca@cell.embeddings[,1:30]
)))

rvel.cd <- gene.relative.velocity.estimates(all_emat,
                                            all_nmat,
                                            deltaT=1,
                                            kCells=20,
                                            cell.dist=cell.dist,
                                            fit.quantile=0.02)
print("ok")
pdf(my_pdf, width = 8, height = 8)
show.velocity.on.embedding.cor(emb,
                               rvel.cd,
                               n=300,
                               scale='sqrt',
                               #cell.colors=my_color,
                               cex=0.8,arrow.scale=5,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,
                               grid.n=40,arrow.lwd=1,
                               do.par=F,
                               cell.border.alpha = 0.1)
dev.off()
print('OK')

all_emat <- as.data.frame(all_emat)
all_nmat <- as.data.frame(all_nmat)
sp <- all_emat %>% colSums() %>% as.data.frame()
unsp <- all_nmat %>% colSums() %>% as.data.frame()

colnames(sp)[1] <- 'sp'
colnames(unsp)[1] <- 'unsp'

cell_type <- Clean_sct.inte.rm.lpt@meta.data %>% rownames_to_column('my_sample')
cell_type <- cell_type[,c('cell.type.percise.new', 'my_sample','condition')]
sp_1 <- sp %>%
  rownames_to_column('my_sample') %>%
  left_join(cell_type, by = 'my_sample') %>%
  filter()

sp_unsp <- cbind(sp, unsp) %>%  rownames_to_column('my_sample') %>%
  left_join(cell_type, by = 'my_sample') %>% na.omit()
sp_unsp$ unsp_ratio <- sp_unsp$unsp/(sp_unsp$sp + sp_unsp$unsp)
sp_unsp$sp_ratio <- sp_unsp$sp/(sp_unsp$sp + sp_unsp$unsp)
saveRDS(sp_unsp, file = 'scvelo_sp_unsp_ratio.rds')

median(sp_unsp[sp_unsp$cell.type.percise.new == 'Top2a + cancer cells' & sp_unsp$condition == 'Untreated',]$unsp_ratio)
median(sp_unsp[sp_unsp$cell.type.percise.new == 'Top2a + cancer cells' & sp_unsp$condition == 'Targeted',]$unsp_ratio)
median(sp_unsp[sp_unsp$cell.type.percise.new == 'Padi4 + cancer cells' & sp_unsp$condition == 'Untreated',]$unsp_ratio)
median(sp_unsp[sp_unsp$cell.type.percise.new == 'Padi4 + cancer cells' & sp_unsp$condition == 'Targeted',]$unsp_ratio)


# 15.SCENIC run by pySCENIC, DISCARD -----

# 16.metabolism for T cells discard DISCARD-----
dir.create('/data/02.project/00.other.lab/07.lihaohuan/03.analysis/compass')

T.counts <- GetAssayData(T.inte,assay = 'RNA',slot = 'counts') %>% as.data.frame()
write.table(T.counts, file = '/data/02.project/00.other.lab/07.lihaohuan/03.analysis/compass/T.counts.tsv',row.names=T,sep = "\t",quote=F)

# then we run compass on linux server
# bash: compass --data T.counts.tsv --num-processes 80 --microcluster-size 8 --species mus_musculus

T.reaction <- rio::import('/data/02.project/00.other.lab/07.lihaohuan/03.analysis/compass/reactions.tsv')
T.reaction <- T.reaction %>% separate(col = V1, into = c('rxn_code_nodirection','direction'), sep = '_', remove = T) %>% dplyr::filter(direction == 'pos') %>% dplyr::select(!c('direction'))
T.micro <- rio::import('/data/02.project/00.other.lab/07.lihaohuan/03.analysis/compass/micropools.tsv')
T.micro$microcluster <- paste0('cluster_',T.micro$microcluster)
T.meta <- T.inte@meta.data
rxn_md <- rio::import('/data/02.project/00.other.lab/07.lihaohuan/03.analysis/compass/rxn_md.csv')

T.cell.reaction <- T.reaction %>% 
  left_join(y = rxn_md[,c(1,3)], by = c('rxn_code_nodirection')) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = 'microcluster') %>% 
  left_join(y = T.micro, by = c('microcluster'))
  

# 17.pseudotime ----
library(monocle3)
library(ggpmisc)
## meta
all.cell.meta <- Clean_sct.inte.rm.lpt@meta.data
all.cell.meta$umap1 <- Embeddings(Clean_sct.inte.rm.lpt,reduction = 'umap')[1]
all.cell.meta$umap2 <- Embeddings(Clean_sct.inte.rm.lpt,reduction = 'umap')[2]

CAF.meta <- CAF.merge@meta.data
CAF.meta$umap1 <- Embeddings(CAF.merge,reduction = 'umap')[1]
CAF.meta$umap2 <- Embeddings(CAF.merge,reduction = 'umap')[2]

Cancer.meta <- Cancer.merge@meta.data
Cancer.meta$umap1 <- Embeddings(Cancer.merge,reduction = 'umap')[1]
Cancer.meta$umap2 <- Embeddings(Cancer.merge,reduction = 'umap')[2]

## run mono
use.cell <- 'Cancer'

mono <- GetAssayData(get(paste0(use.cell,'.merge')),assay = 'RNA', slot = 'counts')
mono.cell_meta <- get(paste0(use.cell,'.meta'))

mono.gene_annotation <- data.frame(gene_short_name = rownames(mono))
rownames(mono.gene_annotation) <- rownames(mono)

cds <- new_cell_data_set(mono,
                         cell_metadata = mono.cell_meta,
                         gene_metadata = mono.gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
plot_pc_variance_explained(cds)

## mono umap
cds <- reduce_dimension(cds, preprocess_method = "PCA",cores = 80)
p1 <- plot_cells(cds,
                 reduction_method="UMAP", 
                 color_cells_by="condition",
                 label_cell_groups = F,
                 show_trajectory_graph = F,
                 cell_stroke = .2,
                 group_label_size = 5,
                 cell_size = .1) +
  ggtitle(use.cell) +
  scale_color_manual(values = c("#709AE1FF","#FED439FF")) + 
  NoAxes()

# seurat umap
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(bac.1w, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed
# plot_cells(cds, reduction_method="UMAP", color_cells_by="condition") + ggtitle('bov.umap')
# cds <- cluster_cells(cds,reduction_method = 'UMAP',resolution = .001)

cds <- cluster_cells(cds,reduction_method = 'UMAP',resolution = 0.001)
## get trajectory
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, 
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F,
           color_cells_by="condition",
           group_label_size = 0,
           cell_size = 1) + 
  scale_color_npg() +
  ggtitle("Pseudotime") +
  NoAxes() + 
  theme(panel.background=element_rect(fill='transparent', color='black'),
        title = element_text(size = 10),
        legend.text = element_text(size = 8), 
        legend.key.height=unit(.8,"line"))

plot_cells(cds,
           genes=c("Top2a",'Mki67'),
           show_trajectory_graph=F,
           label_cell_groups=F,
           label_leaves=F,
           cell_size = .5) 

cds <- order_cells(cds)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", 
                 label_cell_groups = F, 
                 label_leaves = FALSE,  
                 label_branch_points = F,
                 group_label_size = 1,
                 cell_size = .1,
                 label_roots = F,
                 cell_stroke = .2,
                 trajectory_graph_segment_size = 1
) + NoAxes() 

p1 + p2

## get pseudotime genes
library(ClusterGVis)
mono.modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 50)
mono.genes <- c(row.names(subset(mono.modulated_genes, q_value < 0.001 & morans_I > 0.2)))
mono.mat <- pre_pseudotime_matrix(cds_obj = cds,assays = 'counts',
                             gene_list = mono.genes)

mat.plot <- mono.mat[,seq(1,10000,100)]
ck <- clusterData(exp = mat.plot,
                  cluster.method = "kmeans",
                  cluster.num = 5)

pdf(paste0('../analysis.data/',use.cell,'pseudotime.heatmap.pdf'),height = 10,width = 8,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           markGenes = c('Ucp1'))
dev.off()

## bin plot 
pseudotime.df <- data.frame(pseudotime = pseudotime(cds)) %>% 
  rownames_to_column(var = 'cells') %>% 
  left_join(y = mono.cell_meta[,c(17,18)], by = 'cells')

pseudotime.df$bin <- cut(pseudotime.df$pseudotime, breaks = seq(0, ceiling(max(pseudotime.df$pseudotime)), by = 1), include.lowest = TRUE, right = FALSE)
result <- pseudotime.df %>%
  group_by(bin, condition) %>%
  summarise(count = n()) %>%
  spread(key = condition, value = count, fill = 0) %>%
  mutate(treat_ratio = Targeted / (WT + Targeted)) 

result$bin_numeric <- as.numeric(result$bin)

model <- lm(treat_ratio ~ bin_numeric, data = result)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4] 

ggplot(result, aes(x = bin_numeric, y = treat_ratio)) +
  geom_point(size = 2.5,color = '#FF5E5B') + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +  
  annotate("text", x = 5, y = 0.8, label = paste("R² = ", round(r_squared, 3), "\np = ", format(p_value, digits = 3)),
           color = "black", size = 5) +  
  labs(x = "Pseudotime", y = "Targted nuclei ratio") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18))

## bin plot for velocyto

### padi4 + cancer
pseudotime.velo.df.caf <- pseudotime.df %>% left_join(y = sp_unsp[,c(1,6)], by = c('cells' = 'my_sample'))
result <- pseudotime.velo.df.caf %>%
  group_by(bin, condition) %>%
  summarise(ratio = median(unsp_ratio)) %>%
  spread(key = condition, value = ratio, fill = 0) 
result <- result[result$Targeted !=0 & result$WT !=0,]

result$bin_numeric <- as.numeric(result$bin)
result$ratio <- result$Targeted/result$WT

model <- lm(ratio ~ bin_numeric, data = result)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4] 

ggplot(result, aes(x = bin_numeric, y = ratio)) +
  geom_point(size = 2.5,color = '#FF5E5B') + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +  
  annotate("text", x = 5, y = 1.2, label = paste("R² = ", round(r_squared, 3), "\np = ", format(p_value, digits = 3)),
           color = "black", size = 5) +  
  labs(x = "Pseudotime", y = "Targeted/untreated unspliced RNA ratio") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 15))

### top2a + cancer
pseudotime.velo.df.cancer <- pseudotime.df %>% left_join(y = sp_unsp[,c(1,6)], by = c('cells' = 'my_sample'))
result <- pseudotime.velo.df.cancer %>%
  group_by(bin, condition) %>%
  summarise(ratio = median(unsp_ratio)) %>%
  spread(key = condition, value = ratio, fill = 0) 
result <- result[result$Targeted !=0 & result$WT !=0,]

result$bin_numeric <- as.numeric(result$bin)
result$ratio <- result$Targeted/result$WT

model <- lm(ratio ~ bin_numeric, data = result)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4] 

ggplot(result, aes(x = bin_numeric, y = ratio)) +
  geom_point(size = 2.5,color = '#FF5E5B') + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +  
  annotate("text", x = 5, y = 1.4, label = paste("R² = ", round(r_squared, 3), "\np = ", format(p_value, digits = 3)),
           color = "black", size = 5) +  
  labs(x = "Pseudotime", y = "Targeted/untreated unspliced RNA ratio") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 15))

# 18.Pi & Ti NOT FIT as samples is 2 vs 2 DISCARD----
cell.number <- FetchData(Clean_sct.inte.rm.lpt, 
                         vars = c("condition", "cell.type.percise.new")) %>%
  dplyr::count(condition, cell.type.percise.new) %>% 
  tidyr::spread(condition, n) 

cell.ratio <- FetchData(Clean_sct.inte.rm.lpt, 
                       vars = c("orig.ident", "cell.type.percise.new")) %>%
  dplyr::count(orig.ident, cell.type.percise.new) %>% 
  group_by(orig.ident) %>% 
  summarise(Prop = n / sum(n)) 

cell.ratio$cell.type <- factor(rep(levels(Clean_sct.inte.rm.lpt$cell.type.percise.new), 4),
                               levels = levels(Clean_sct.inte.rm.lpt$cell.type.percise.new))

tumor.size <- data.frame(targeted = c(99.7,85.3,108.3,152.8),
                         untreated = c(88.4,113.2,875.7,891.4),
                         condition = c(rep(c('pre','post'),each = 2)))

Pi.macro.data <- data.frame(macro = c(0.0231076558, 0.0146311221),tumor.size = c(152.8-99.7,108.3-85.3)) 
fit <- lm(tumor.size ~ macro, data = Pi.macro.data)
slope <- coef(fit)["macro"]
R_squared <- summary(fit)$r.squared
Pi <- (-slope / abs(slope)) * R_squared

# 19.DEGs in cell types ----
library(ClusterGVis)

ht.data <- prepareDataFromscRNA(Clean_sct.inte.rm.lpt,
                                diffData = rm.lpt.ct.deg,
                                showAverage = T,
                                assays = 'SCT',slot = 'data',
                                group.by = 'cell.type.percise.new',keep.uniqGene = F,
                                scale.data = T)

enrich.go <- enrichCluster(object = ht.data,
                             OrgDb = org.Mm.eg.db,
                             type = "BP",
                             organism = "mmu",
                             pvalueCutoff = 0.05,
                             topn = 5,
                             seed = 1234)

pdf('../lhh.gland/new.analysis/add.11.8/sc.go.pdf',height = 12,width = 14,onefile = F)

visCluster(object = ht.data,
        ht.col.list = list(col_range = c(-4, 0, 4)),
        plot.type = "both",
        column_names_rot = 45,
        show_row_dend = F,
        # markGenes = markGenes,
        markGenes.side = "left",
        annoTerm.data = enrich.go,
        line.side = "left",
        cluster.order = c(1:13),
        go.col = rep(use.cols,each = 5),
        sample.col = use.cols,
        ctAnno.col = use.cols,
        go.size = 8,
        add.bar = T)
dev.off()
