setwd('/bettik/lbundalian-ext/Practice/RNASeq.Half.Cell')

library(edgeR)
library(tidyverse)
library(vroom)
library(ggplot2)
library(fgsea)
library(DOSE)
library(clusterProfiler)
library(future)
library(data.table)
library(enrichplot)
plan("multiprocess", workers = 4)

files <- list.files(path = "Data/supp",pattern = ".xlsx", full.names = TRUE)
dataset <- lapply(files, readxl::read_xlsx) 
labels <- list.files(path = "Data/supp",
                        pattern = ".xlsx") %>% gsub(".xlsx.*", "", .)
names(dataset) <- labels

for (e in 1:length(dataset)) {

  
  if (e %in% c(1:3)) {
   colnames(dataset[[e]]) <- dataset[[e]][2,]
   dataset[[e]] <- dataset[[e]][-c(1,2),]
 } else {
   
   colnames(dataset[[e]]) <- dataset[[e]][1,]
   dataset[[e]] <- dataset[[e]][-c(1),]
 }
  
}

mrna <- vroom("Data/genes/mrna.gct", skip=2) %>% as.data.frame
mirna <- vroom("Data/genes/mirna.gct") %>% as.data.frame


assay <- list()
meta <- list()

assay[['MIRNA']] <- mirna[,-c(1,2)] %>% as.matrix
rownames(assay[['MIRNA']]) <- mirna$Name
meta[['MIRNA']] <- colnames(mirna[,-c(1,2)])

assay[['MRNA']] <- mrna[,-c(1,2)] %>% as.matrix
rownames(assay[['MRNA']]) <- mrna$NAME
meta[['MRNA']] <- colnames(mrna[,-c(1,2)])

mirna.idx <- c(1:19)
mrna.idx <- c(1:5,7:20)
assay$MIRNA <- assay$MIRNA[,mirna.idx]
assay$MRNA <- assay$MRNA[,mrna.idx]
colnames(assay$MIRNA)
colnames(assay$MRNA)
idx <- NULL
for(i in colnames(assay$MRNA)){
    idx <- append(idx,paste0('HC_',gsub(".*?([0-9]+).*", "\\1", i)))
}
colnames(assay$MIRNA) <- idx
colnames(assay$MRNA) <- idx

tables <- list()
# Converts the assay to table
tables[['MIRNA']] <- as.table(assay$MIRNA) %>% as.data.frame 
colnames(tables[['MIRNA']]) <- c("GENES","HC","COUNTS")
tables[['MRNA']] <- as.table(assay$MRNA) %>% as.data.frame
colnames(tables[['MRNA']]) <- c("GENES","HC","COUNTS")

# Combine duplicates by values_fn = mean ( Averaging duplicate genes)
tables[['MIRNA_WIDE']] <- tables$MIRNA %>% 
  pivot_wider(names_from=HC, values_from=COUNTS, values_fn = mean) 
tables[['MRNA_WIDE']] <- tables$MRNA %>% 
  pivot_wider(names_from=HC, values_from=COUNTS, values_fn = mean) 

# reassign to assay
assay$MIRNA <- tables$MIRNA_WIDE %>% column_to_rownames('GENES') %>% as.matrix
assay$MRNA <- tables$MRNA_WIDE %>% column_to_rownames('GENES') %>% as.matrix

# assign assay with log base 2 for mRNA since it is not log normalized
assay[['MRNA.LOG']] <- log(base=2,assay$MRNA)

# assign default assay to log since miRNA is already log normalized
assay[['MIRNA.LOG']] <- assay$MIRNA

# calculate the base 10 for original miRNA assay
assay[['MIRNA']] <- 2^assay$MIRNA

hist(rowSums(assay$MIRNA),breaks=150, main = "miRNA Counts", xlab = "miRNA")
hist(rowSums(assay$MIRNA.LOG),breaks=150, main = "miRNA log 2 Counts", xlab = "miRNA")
hist(rowSums(assay$MRNA),breaks=150, main = "mRNA Counts", xlab = "mRNA")
hist(rowSums(assay$MRNA.LOG),breaks=150, main = "mRNA log2 Counts", xlab = "mRNA")

cor.matrix <- list()
cor.matrix[['RAW']] <- cor(t(assay$MIRNA),t(assay$MRNA),use='pairwise.complete.obs')
cor.matrix[['LOG2']] <- cor(t(assay$MIRNA.LOG),t(assay$MRNA.LOG),use='pairwise.complete.obs')

cor.matrix$RAW  <- cor.matrix$RAW %>% as.table %>% as.data.frame
cor.matrix$LOG2  <- cor.matrix$LOG2 %>% as.table %>% as.data.frame

selection <- list()
selection[['log2_miR92']] <- cor.matrix$LOG2 %>% filter(Var1 == 'hsa-miR-92a-3p') %>% 
column_to_rownames("Var2") %>% dplyr::select(Freq) %>% as.matrix %>%
.[,1]
selection[['raw_miR92']] <- cor.matrix$RAW %>% filter(Var1 == 'hsa-miR-92a-3p') %>% 
column_to_rownames("Var2") %>% dplyr::select(Freq) %>% as.matrix %>%
.[,1]

targets <- list()
targets[['ALL.MIR92']] <- vroom('Data/targets/all_target_mir92a.txt') %>% dplyr::select(c(1,15))
colnames(targets[['ALL.MIR92']]) <- c('GENE','SCORE')
targets[['CON.MIR92']] <- vroom('Data/targets/conserved_target_mir92a.txt') %>% dplyr::select(c(1,15))
colnames(targets[['CON.MIR92']]) <- c('GENE','SCORE')

plots <- list()

options(repr.plot.width=10, repr.plot.height=10)
for (thr in c(-0.7, -0.5, -0.3, -0.1)){
  target.list <- targets$ALL.MIR92 %>% as.data.frame %>% filter(SCORE < -0.1) %>%
    dplyr::select(GENE) %>% mutate(TERM="GENE") %>% dplyr::select(TERM,GENE)
 
  result <- GSEA(
            sort(selection$raw_miR92, decreasing = T),
            exponent = 1,
            minGSSize = 100,
            maxGSSize = nrow(target.list),
            eps = 1e-10,
            pvalueCutoff = 10,
            pAdjustMethod = "fdr",
            target.list,
            verbose = TRUE,
            seed = FALSE,
            by = "fgsea"
            )
    
  # broad institute 
    plots[[paste0("THR",thr)]] <- gseaplot2(result, geneSetID = 1, title = "Enrichment for hsa-miR-92a-3p" , 
          subplots=1:3,  pvalue_table = TRUE,  rel_heights = c(1.5, 0.25, .25))     
}




for(p in plots){
    print(p)
}
