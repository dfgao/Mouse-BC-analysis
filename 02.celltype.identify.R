# intergrate data ------
plan(multisession, workers=30)
Clean_sct.inte <- Clean_qc.merge.filtered %>%
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = .5) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)
Clean_sct.inte <- FindClusters(Clean_sct.inte, resolution = c(.2,.4,.6,1))
DimPlot(Clean_sct.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.2',split.by = 'condition')
VlnPlot(Clean_sct.inte, group.by = 'SCT_snn_res.1', pt.size = 0,features = c('nFeature_RNA','nCount_RNA'))

Clean_sct.inte <- RunUMAP(Clean_sct.inte, reduction = "integrated.dr", 
                             dims = 1:30)
DimPlot(Clean_sct.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.2',split.by = 'condition')

# cell types identify----
Idents(Clean_sct.inte) <- Clean_sct.inte$SCT_snn_res.0.2
DefaultAssay(Clean_sct.inte) <- 'SCT'
Clean_sct.inte <- PrepSCTFindMarkers(Clean_sct.inte)
deg.res.0.2 <- FindAllMarkers(Clean_sct.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .2) %>% dplyr::filter(p_val_adj < 0.05)
deg.res.0.2.top10 <- deg.res.0.2 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

Clean_sct.inte <- FindSubCluster(Clean_sct.inte, cluster = '10',graph.name = 'SCT_snn',resolution = .1)
DimPlot(Clean_sct.inte, reduction = "umap", label = T,group.by = 'sub.cluster',split.by = 'condition')
Idents(Clean_sct.inte) <- Clean_sct.inte$sub.cluster
deg.res.0.2 <- FindAllMarkers(Clean_sct.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .2) %>% dplyr::filter(p_val_adj < 0.05)
deg.res.0.2.top10 <- deg.res.0.2 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

## major cell type

# epi 9
FeaturePlot(Clean_sct.inte, 
            features = c('Epcam','Krt19','Car3'),
            order = T)
FeaturePlot(Clean_sct.inte.rm.lpt,split.by = 'condition', 
            features = c('Mki67'),
            order = T)
FeaturePlot(Clean_sct.inte, 
            features = c('Adipoq','Adipor1','Krt14','Cdh1','Muc1'),
            order = T)
# endo 8 LEC-res1-25
FeaturePlot(Clean_sct.inte, 
            features = c('Pecam1','Kdr','Tek','Vwf','Cldn5','Reln'),
            order = T)
# pericyte 13 
FeaturePlot(Clean_sct.inte, 
            features = c('Rgs5','Notch3','Abcc9','Myh3','Mybpc1'),
            order = T)
# MSC
FeaturePlot(Clean_sct.inte, 
            features = c('Pdgfrb','Col3a1','Fap','Acta2'),
            order = T)
# immune-myeloid 2 4 7 11 12 14
FeaturePlot(Clean_sct.inte, 
            features = c('Ptprc','Itgam'),
            order = T)
# Blood NA
FeaturePlot(Clean_sct.inte, 
            features = c('Hbb-bt','Hba-a1'),
            order = T)
# test epi cancer 0 1 3 5 
FeaturePlot(Clean_sct.inte, 
            features = c('Krt8','Krt18'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)
# test breast cancer ESR+
FeaturePlot(Clean_sct.inte, 
            features = c('Esr1','Pgr','Erbb2','Prlr'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)
# test CSCs NA
FeaturePlot(Clean_sct.inte,
            features = c('Cd44','Cd24a','Aldh1a1'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)
# test CCs 1 5 
FeaturePlot(Clean_sct.inte,
            features = c('Mki67','Top2a','Ccnb1'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)
# test EMT all VIM
FeaturePlot(Clean_sct.inte,
            features = c('Vim','Snai1','Twist1','Zeb1'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)
# test MEM
FeaturePlot(Clean_sct.inte,
            features = c('Vegfa','Csf1'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)
# test transfer
FeaturePlot(Clean_sct.inte,
            features = c('Mmp9','S100a4'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)
# test special gene
FeaturePlot(Clean_sct.inte,
            features = c('Myc','Platr22'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)

Clean_sct.inte <- RenameIdents(Clean_sct.inte,
                               '0' = 'Padi4 + cancer cells',
                               '3' = 'Padi4 + cancer cells',
                               '1' = 'Top2a + cancer cells',
                               '5' = 'Top2a + cancer cells',
                               '6' = 'MSCs',
                               '8' = 'Endothelial cells',
                               '13' = 'Pericytes',
                               '9' = 'Erbb4 + Epithelial cells',
                               '12' = 'Prlr + Epithelial cells',
                               '10_0' = 'Cdh3 + Epithelial cells',
                               '10_1' = 'Adipocytes',
                               '2' = 'ICs',
                               '4' = 'ICs',
                               '7' = 'ICs',
                               '11' = 'ICs',
                               '14' = 'ICs'
)
Clean_sct.inte$cell.type.major <- Idents(Clean_sct.inte)
DimPlot(Clean_sct.inte, group.by = 'cell.type.major')

## ICs subcluster ----

### integrate 
ICs.sub <- subset(Clean_qc.merge.filtered, cells = rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.major == 'ICs',]))

ICs.inte <- ICs.sub %>%
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = .5) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)

ICs.inte <- FindClusters(ICs.inte, resolution = c(.2,.5,1))
DimPlot(ICs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.1')

### HVGs
Idents(ICs.inte) <- ICs.inte$SCT_snn_res.1
DefaultAssay(ICs.inte) <- 'SCT'
ICs.inte <- PrepSCTFindMarkers(ICs.inte)
ICs.deg.res.1 <- FindAllMarkers(ICs.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .1) %>% dplyr::filter(p_val_adj < 0.05)
ICs.deg.res.1.top10 <- ICs.deg.res.1 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

## major markers
FeaturePlot(ICs.inte, 
            features = c('Nkg7','Cd14','Cd8a','Cd3e'),
            order = T)

#NKT NA
FeaturePlot(ICs.inte, 
            features = c('Gzmb','Klrb1c','Tbx21','Ncr1','Foxp3','Il2ra','Cd4','Cd8a','Cd3e'),
            order = T)
#NK NA
FeaturePlot(ICs.inte, 
            features = c('Klrb1c','Itga2','Klrk1'),
            order = T)
# macro 0 2 3 4 6 10 11 14
FeaturePlot(ICs.inte, 
            features = c('Cd68','Pf4','Adgre1','C1qc','Apoe',"Ctsd","Lpl","Csf1r","Atp6v0d2",'Spp1'))

# mono 9
FeaturePlot(ICs.inte, 
            features = c('Ly6c2','Gbp2','Plac8',"Fscn1","Traf1","Epsti1","F13a1","Pou2f2","Ifi27l2a",'Ifitm3'), order = T)

# B NA
FeaturePlot(ICs.inte, 
            features = c('Cd79a','Cd79b','Ms4a1','Prdm1'),order = T)

# T 7 8 15
FeaturePlot(ICs.inte,
            features = c('Ccl5','Cd3g','Cd8b1','Cd3e'),order = T)

# neut 1 5 16
FeaturePlot(ICs.inte,
            features = c('S100a9','S100a8','Asprv1','G0s2','Lcn2','Rsad2'),order = T)

# mast NA
FeaturePlot(ICs.inte,
            features = c('Ifitm1','Ccr3','Cd200r3','Cpa3'),order = T)

# mast na
FeaturePlot(ICs.inte,
            features = c('Kit','Fcer1g','Fcer2','Kit'),order = F)

# GMP sub of 5
FeaturePlot(ICs.inte,
            features = c('Stfa1','Ngp','Mpo'),order = T)

# Immature dendritic cell 13 18 19
FeaturePlot(ICs.inte,
            features = c('Cd14','Cd83'),order = T,cols = c('gray90','red3'), split.by = 'condition')

FeaturePlot(ICs.inte,
            features = c('Esr1','Esr2','Gper1'),order = T)

DimPlot(Clean_sct.inte, cells.highlight = rownames(ICs.inte@meta.data[ICs.inte$SCT_snn_res.1 == '18',]))

ICs.inte <- RenameIdents(ICs.inte,
                         '1' = 'Neutrophils',
                         '16' = 'Neutrophils',
                         '5' = 'Neutrophils',
                         '9' = 'Monocytes',
                         '0' = 'Macrophages',
                         '2' = 'Macrophages',
                         '3' = 'Macrophages',
                         '4' = 'Macrophages',
                         '6' = 'Macrophages',
                         '11' = 'Macrophages',
                         '14' = 'Macrophages',
                         '10' = 'Macrophages',
                         '17' = 'Macrophages',
                         '12' = 'Macrophages',
                         '13' = 'DCs',
                         '18' = 'DCs',
                         '19' = 'DCs',
                         '15' = 'T cells',
                         '7' = 'T cells',
                         '8' = 'T cells')
ICs.inte$cell.type <- Idents(ICs.inte)
DimPlot(ICs.inte, group.by = 'cell.type',label = T,repel = T)

## T subcluster ----
T.sub <- subset(Clean_qc.merge.filtered, cells = rownames(ICs.inte@meta.data[ICs.inte$cell.type == 'T cells',]))

## iontegrate
T.inte <- T.sub %>%
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = .5) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)

T.inte <- FindClusters(T.inte, resolution = c(.2,.5,1))
DimPlot(T.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.1')

## HVGs
Idents(T.inte) <- T.inte$SCT_snn_res.1
DefaultAssay(T.inte) <- 'SCT'
T.inte <- PrepSCTFindMarkers(T.inte)
T.deg.res.1 <- FindAllMarkers(T.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .2) %>% dplyr::filter(p_val_adj < 0.05)
T.deg.res.1.top10 <- T.deg.res.1 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
VlnPlot(T.inte, group.by = 'SCT_snn_res.1', pt.size = 0,features = c('nFeature_RNA','nCount_RNA'))

# T 
FeaturePlot(T.inte, 
            features = c('Gzmb','Klrb1c','Tbx21','Ncr1','Foxp3','Il2ra','Cd4','Cd8a','Cd3e'),
            order = T)
DimPlot(Clean_sct.inte, cells.highlight = rownames(T.inte@meta.data[T.inte$SCT_snn_res.0.2 == '3',]))
FeaturePlot(T.inte, 
            features = c('Foxp3','Cd4','Cd8a','Tbx21','Bcl2','Il7r'),
            order = T)
T.num <- FetchData(T.inte, 
                         vars = c("condition", "SCT_snn_res.0.2")) %>%
  dplyr::count(condition, SCT_snn_res.0.2) %>% 
  tidyr::spread(condition, n) 
T.inte <- RenameIdents(T.inte,
                       '0' = 'Cd8 + T',
                       '1' = 'Cd8 + T',
                       '3' = 'Cd4 + T',
                       '2' = 'Memory T',
                       '4' = 'Treg',
                       '5' = 'NKT',
                       '6' = 'NKT')
T.inte$cell.type.new <- Idents(T.inte)
DimPlot(T.inte, group.by = 'cell.type.new')

## return cell types ----
Clean_sct.inte$cell.type.minor <- as.character( Clean_sct.inte$cell.type.major) 
Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.minor == 'ICs',]$cell.type.minor <- as.character(ICs.inte$cell.type)

DimPlot(Clean_sct.inte, group.by = 'cell.type.minor',label = T,cols = c(npg.cols, use.cols), repel = T) + NoAxes() + ggtitle("") 
Idents(Clean_sct.inte) <- Clean_sct.inte$cell.type.minor
table(Clean_sct.inte$cell.type.minor)

# remove adipocytes -----
Clean_sct.inte.rm.lpt <- subset(Clean_sct.inte, cells = rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.minor == 'Adipocytes',]), invert = T)
Clean_sct.inte.rm.lpt <- RunUMAP(Clean_sct.inte.rm.lpt,dims = 1:30)
DimPlot(Clean_sct.inte.rm.lpt,group.by = 'cell.type.minor',label = T, cols = c(npg.cols, use.cols))

Clean_sct.inte.rm.lpt <- subset(Clean_qc.merge.filtered, cells = rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.minor != 'Adipocytes',]))
Clean_sct.inte.rm.lpt <- Clean_sct.inte.rm.lpt %>%
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = .5) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)
Clean_sct.inte.rm.lpt <- FindClusters(Clean_sct.inte.rm.lpt, resolution = c(.2,.3,.5,1))
DimPlot(Clean_sct.inte.rm.lpt, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.3') + scale_color_igv()

Clean_sct.inte.rm.lpt$cell.type.percise.raw <- Clean_sct.inte@meta.data[match(colnames(Clean_sct.inte.rm.lpt), colnames(Clean_sct.inte)),]$cell.type.minor
DimPlot(Clean_sct.inte.rm.lpt, reduction = "umap", label = T,group.by = 'cell.type.percise.raw',split.by = 'condition')
