library(Seurat)

s1 <- Read10X_h5("C://Bioinf/Penny_project/filtered_feature_bc_matrix_s1.h5", use.names = TRUE, unique.features = TRUE)
s1 <- CreateSeuratObject(counts = s1, project = "s1", min.cells = 3, min.features = 200)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat)
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20)
  Seurat <- FindClusters(Seurat, resolution = 0.5)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20 )
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

s1 [["percent.rb"]] <- PercentageFeatureSet(s1 , pattern = "^Rps|^Rpl|^Mrps|^Mrpl", assay = 'RNA')

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

s1 <- CellCycleScoring(s1, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^mt-")
VlnPlot(s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)

s1 <- subset(s1, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20 & percent.rb < 30)

s1 <- ScaleData(s1, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))

s1 <- ProcessSeu(s1)

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters) 
  nExp_poi <- round(0.1*length(colnames(tmp)))  ## Assuming 10% doublet formation rate 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp) 
}

library(DoubletFinder)

s1 <- RDoublet(s1)

s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][12] == 'Singlet')])
s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][13] == 'Singlet')])

s1 <- ProcessSeu(s1)

s1_final <- s1
s1_final$condition <- 's1'
s1_final$group <- 'Naive'
s1_final$celltype <- 'Microglia'

s1 <- Read10X_h5("C://Bioinf/Penny_project/filtered_feature_bc_matrix_s2.h5", use.names = TRUE, unique.features = TRUE)
s1 <- CreateSeuratObject(counts = s1, project = "s1", min.cells = 3, min.features = 200)
s1 [["percent.rb"]] <- PercentageFeatureSet(s1 , pattern = "^Rps|^Rpl|^Mrps|^Mrpl", assay = 'RNA')
s1 <- CellCycleScoring(s1, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^mt-")
VlnPlot(s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
s1 <- subset(s1, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20 & percent.rb < 30)
s1 <- ScaleData(s1, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))
s1 <- ProcessSeu(s1)

s1 <- RDoublet(s1)

s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][12] == 'Singlet')])
s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][13] == 'Singlet')])

s1 <- ProcessSeu(s1)

s2_final <- s1
s2_final$condition <- 's2'
s2_final$group <- 'Naive'
s2_final$celltype <- 'RGC'

s1 <- Read10X_h5("C://Bioinf/Penny_project/filtered_feature_bc_matrix_s3.h5", use.names = TRUE, unique.features = TRUE)
s1 <- CreateSeuratObject(counts = s1, project = "s1", min.cells = 3, min.features = 200)
s1 [["percent.rb"]] <- PercentageFeatureSet(s1 , pattern = "^Rps|^Rpl|^Mrps|^Mrpl", assay = 'RNA')
s1 <- CellCycleScoring(s1, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^mt-")
VlnPlot(s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
s1 <- subset(s1, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20 & percent.rb < 30)
s1 <- ScaleData(s1, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))
s1 <- ProcessSeu(s1)

s1 <- RDoublet(s1)

s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][12] == 'Singlet')])
s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][13] == 'Singlet')])

s1 <- ProcessSeu(s1)

s3_final <- s1
s3_final$condition <- 's3'
s3_final$group <- 'Vehicle'
s3_final$celltype <- 'Microglia'

s1 <- Read10X_h5("C://Bioinf/Penny_project/filtered_feature_bc_matrix_s4.h5", use.names = TRUE, unique.features = TRUE)
s1 <- CreateSeuratObject(counts = s1, project = "s1", min.cells = 3, min.features = 200)
s1 [["percent.rb"]] <- PercentageFeatureSet(s1 , pattern = "^Rps|^Rpl|^Mrps|^Mrpl", assay = 'RNA')
s1 <- CellCycleScoring(s1, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^mt-")
VlnPlot(s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
s1 <- subset(s1, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20 & percent.rb < 30)
s1 <- ScaleData(s1, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))
s1 <- ProcessSeu(s1)

s1 <- RDoublet(s1)

s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][12] == 'Singlet')])
s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][13] == 'Singlet')])

s1 <- ProcessSeu(s1)

s4_final <- s1
s4_final$condition <- 's4'
s4_final$group <- 'Vehicle'
s4_final$celltype <- 'RGC'

s1 <- Read10X_h5("C://Bioinf/Penny_project/filtered_feature_bc_matrix_s5.h5", use.names = TRUE, unique.features = TRUE)
s1 <- CreateSeuratObject(counts = s1, project = "s1", min.cells = 3, min.features = 200)
s1 [["percent.rb"]] <- PercentageFeatureSet(s1 , pattern = "^Rps|^Rpl|^Mrps|^Mrpl", assay = 'RNA')
s1 <- CellCycleScoring(s1, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^mt-")
VlnPlot(s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
s1 <- subset(s1, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20 & percent.rb < 30)
s1 <- ScaleData(s1, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))
s1 <- ProcessSeu(s1)

s1 <- RDoublet(s1)

s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][12] == 'Singlet')])
s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][13] == 'Singlet')])

s1 <- ProcessSeu(s1)

s5_final <- s1
s5_final$condition <- 's5'
s5_final$group <- 'BA'
s5_final$celltype <- 'Microglia'

s1 <- Read10X_h5("C://Bioinf/Penny_project/filtered_feature_bc_matrix_s6.h5", use.names = TRUE, unique.features = TRUE)
s1 <- CreateSeuratObject(counts = s1, project = "s1", min.cells = 3, min.features = 200)
s1 [["percent.rb"]] <- PercentageFeatureSet(s1 , pattern = "^Rps|^Rpl|^Mrps|^Mrpl", assay = 'RNA')
s1 <- CellCycleScoring(s1, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^mt-")
VlnPlot(s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
s1 <- subset(s1, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20 & percent.rb < 30)
s1 <- ScaleData(s1, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))
s1 <- ProcessSeu(s1)

s1 <- RDoublet(s1)

s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][12] == 'Singlet')])
s1 <- subset(s1, cells = colnames(s1)[which(s1[[]][13] == 'Singlet')])

s1 <- ProcessSeu(s1)

s6_final <- s1
s6_final$condition <- 's6'
s6_final$group <- 'BA'
s6_final$celltype <- 'RGC'

SaveH5Seurat(s1_final, 'C://Bioinf/Penny_project/s1.h5Seurat')
SaveH5Seurat(s2_final, 'C://Bioinf/Penny_project/s2.h5Seurat')
SaveH5Seurat(s3_final, 'C://Bioinf/Penny_project/s3.h5Seurat')
SaveH5Seurat(s4_final, 'C://Bioinf/Penny_project/s4.h5Seurat')
SaveH5Seurat(s5_final, 'C://Bioinf/Penny_project/s5.h5Seurat')
SaveH5Seurat(s6_final, 'C://Bioinf/Penny_project/s6.h5Seurat')
