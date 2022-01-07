
# Load the ccRCC data -----------------------------------------------------
load('Data/Single/KIDNEY/rennes2grenoble.rda')
tmp <- ls()



# Load the packages and helper functions ----------------------------------
source("Scripts/Packages.R")
source("Scripts/Helper.R")




# For multithreading ------------------------------------------------------
plan("multiprocess", workers = 4)




# Creates Seurat Object for ccRCC data ------------------------------------
gc()
ccRCC <- CreateSeuratObject(counts = matrix.data, project = "ccRCC", 
                                  min.cells = 3, min.features = 200, 
                                  meta.data = meta.data)


# Split the dataset into Tumor and Normal  --------------------------------
ccRCC <- SplitObject(ccRCC, split.by = c("sample"))


# Subset the dataset to Tumor Only ---------------------------------------
ccRCC <- ccRCC[names(ccRCC) %in% c('RCC5 (Kidney, Tumor)')]
ccRCC <- ccRCC$`RCC5 (Kidney, Tumor)`


plots <- list()

ccRCC[["percent.mt"]] <- PercentageFeatureSet(ccRCC, pattern = "^MT-")


plots[["Violin1"]] <- VlnPlot(ccRCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plots[["Scatter1"]] <- FeatureScatter(ccRCC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plots[["Scatter2"]] <- FeatureScatter(ccRCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plots$Violin1 + plots$Scatter1 + plots$Scatter2

ccRCC <- subset(ccRCC, subset = nFeature_RNA > 0 
                & nFeature_RNA < 5000 & percent.mt < 5)
plots[["Violin2"]] <- VlnPlot(ccRCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plots[["Scatter3"]] <- FeatureScatter(ccRCC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plots[["Scatter4"]] <- FeatureScatter(ccRCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plots$Violin2 + plots$Scatter3 + plots$Scatter4


options(future.globals.maxSize= 10000*1024^2)
ccRCC <- NormalizeData(ccRCC, normalization.method = "LogNormalize",
                       scale.factor = 10000)

ccRCC <- FindVariableFeatures(ccRCC, selection.method = "vst", nfeatures = 5000)

top20 <- head(VariableFeatures(ccRCC), 20)
plot1 <- VariableFeaturePlot(ccRCC)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2


all.genes <- rownames(ccRCC)
ccRCC <- ScaleData(ccRCC, features = all.genes)
ccRCC <- RunPCA(ccRCC, features = VariableFeatures(object = ccRCC))

ElbowPlot(ccRCC)
ccRCC <- FindNeighbors(ccRCC, reduction = "pca", dims = 1:15)
ccRCC <- FindClusters(ccRCC, resolution = 0.5)



ccRCC <- RunUMAP(object = ccRCC, dims = 1:15)

DimPlot(ccRCC, reduction = "umap", group.by = "cell.association", label = TRUE)




library(Rmagic)


var.genes <- VariableFeatures(ccRCC)
ccRCC <- magic(ccRCC,genes='all_genes',t='auto')


DefaultAssay(ccRCC) <- "MAGIC_RNA"
ccRCC.counts.all <- ccRCC@assays$MAGIC_RNA@data

ccRCC.counts <- ccRCC.counts[var.genes,]

sc.cor.matrix <- cor(t(ccRCC.counts),t(ccRCC.counts),use="pairwise.complete.obs")
sc.cor.all <- cor(t(ccRCC.counts.all),t(ccRCC.counts.all),use="pairwise.complete.obs")


sc.cor.df <- sc.cor.matrix %>% as.table %>% as.data.frame

load("0630-ccRCC-5.RData")


sc.cor.df <- sc.cor.df %>% mutate(Abs = abs(Freq))
sc.cor.df %>% filter(Abs < 1) %>% arrange(Freq)

sc_bulk <- inner_join(cor.rna$cor.pairs,sc.cor.df, by = 'Var1')

write.csv(sc.cor.df,"sc.cor.csv")

source("BulkSeq.R")
