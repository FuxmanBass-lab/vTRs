library(ggplot2)
library(ArchR)
library(hexbin)
library(Biostrings)
library(stringr)
library(stringi)
library(Biostrings)
library(parallel)

#set.seed(1)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 4)

tsv_files <- c("/media/jecorban/Expansion/FullvTR_fragments/fullvTR_fragments-1.tsv.gz", "/media/jecorban/Expansion/FullvTR_fragments/fullvTR_fragments-2.tsv.gz", "/media/jecorban/Expansion/FullvTR_fragments/fullvTR_fragments-3.tsv.gz", "/media/jecorban/Expansion/FullvTR_fragments/fullvTR_fragments-4.tsv.gz", "/media/jecorban/Expansion/FullvTR_fragments/fullvTR_fragments-5.tsv.gz", "/media/jecorban/Expansion/FullvTR_fragments/fullvTR_fragments-6.tsv.gz", "/media/jecorban/Expansion/FullvTR_fragments/fullvTR_fragments-7.tsv.gz", "/media/jecorban/Expansion/FullvTR_fragments/fullvTR_fragments-8.tsv.gz")
               
names(tsv_files) <- c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8")

FF_ArrowFiles <- createArrowFiles(
  inputFiles = tsv_files,
  sampleNames = names(tsv_files),
  minTSS = 7, 
  minFrags = 40000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = FALSE 
)

#Note I need to create the doublet scores otherwise subsequent ArchR calls will complain
#but I will _not_ ever filter on them or use them
FF_doubScores <- addDoubletScores(
  input = FF_ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  force = FALSE
)

FF_proj_AllSamples <- ArchRProject(
  ArrowFiles = FF_ArrowFiles,
  copyArrows = FALSE
)

FF_proj_AllSamples_Filt <- FF_proj_AllSamples[which(FF_proj_AllSamples$TSSEnrichment >= 7 & FF_proj_AllSamples$nFrags >= 40000),]


FF_proj_AllSamples_Filt <- addIterativeLSI(ArchRProj = FF_proj_AllSamples_Filt, useMatrix = "TileMatrix", name = "IterativeLSI")

FF_proj_AllSamples_Filt <- addClusters(input = FF_proj_AllSamples_Filt, reducedDims = "IterativeLSI")

FF_proj_AllSamples_Filt <- addUMAP(ArchRProj = FF_proj_AllSamples_Filt, reducedDims = "IterativeLSI")

FF_AllSample_plt <- plotEmbedding(ArchRProj = FF_proj_AllSamples_Filt, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", labelMeans = FALSE)

plotPDF(FF_AllSample_plt, name = "FullvTR-UMAP-Clusters.pdf",
        ArchRProj = FF_proj_AllSamples_Filt, addDOC = FALSE, width = 5, height = 5)


sample1_ass <- read.csv("/media/jecorban/Expansion/FullvTR_MergedAssociations/Sample1_merged_associations.csv", header = TRUE)
sample2_ass <- read.csv("/media/jecorban/Expansion/FullvTR_MergedAssociations/Sample2_merged_associations.csv", header = TRUE)
sample3_ass <- read.csv("/media/jecorban/Expansion/FullvTR_MergedAssociations/Sample3_merged_associations.csv", header = TRUE)
sample4_ass <- read.csv("/media/jecorban/Expansion/FullvTR_MergedAssociations/Sample4_merged_associations.csv", header = TRUE)
sample5_ass <- read.csv("/media/jecorban/Expansion/FullvTR_MergedAssociations/Sample5_merged_associations.csv", header = TRUE)
sample6_ass <- read.csv("/media/jecorban/Expansion/FullvTR_MergedAssociations/Sample6_merged_associations.csv", header = TRUE)
sample7_ass <- read.csv("/media/jecorban/Expansion/FullvTR_MergedAssociations/Sample7_merged_associations.csv", header = TRUE)
sample8_ass <- read.csv("/media/jecorban/Expansion/FullvTR_MergedAssociations/Sample8_merged_associations.csv", header = TRUE)

FF_proj_AllSamples_Filt@cellColData$cellBC_id = stringr::str_split_fixed(rownames(FF_proj_AllSamples_Filt@cellColData),"-",n=2)[,1]

FF_proj_AllSamples_Filt@cellColData$cellBC = stringr::str_split_fixed(FF_proj_AllSamples_Filt@cellColData$cellBC_id,"#",n=2)[,2]

FF_proj_AllSamples_Filt@cellColData$condition = stringr::str_split_fixed(FF_proj_AllSamples_Filt@cellColData$cellBC_id,"#",n=2)[,1]

newDF = data.frame(cellBC = FF_proj_AllSamples_Filt@cellColData$cellBC, condition = FF_proj_AllSamples_Filt@cellColData$condition, order= 1:length(FF_proj_AllSamples_Filt@cellColData$cellBC))

sample1_ass$revComp <- sapply(sample1_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample2_ass$revComp <- sapply(sample2_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample3_ass$revComp <- sapply(sample3_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample4_ass$revComp <- sapply(sample4_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample5_ass$revComp <- sapply(sample5_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample6_ass$revComp <- sapply(sample6_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample7_ass$revComp <- sapply(sample7_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample8_ass$revComp <- sapply(sample8_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))


sample1_ass_SUB = sample1_ass[sample1_ass $revComp %in% newDF$cellBC,]

sample2_ass_SUB = sample2_ass[sample2_ass $revComp %in% newDF$cellBC,]

sample3_ass_SUB = sample3_ass[sample3_ass $revComp %in% newDF$cellBC,]

sample4_ass_SUB = sample4_ass[sample4_ass $revComp %in% newDF$cellBC,]

sample5_ass_SUB = sample5_ass[sample5_ass $revComp %in% newDF$cellBC,]

sample6_ass_SUB = sample6_ass[sample6_ass $revComp %in% newDF$cellBC,]

sample7_ass_SUB = sample7_ass[sample7_ass $revComp %in% newDF$cellBC,]

sample8_ass_SUB = sample8_ass[sample8_ass $revComp %in% newDF$cellBC,]


sample1_ass_SUB_test = sample1_ass[sample1_ass $revComp %in% newDF$cellBC[newDF$condition == "sample1"],]

sample2_ass_SUB_test = sample2_ass[sample2_ass $revComp %in% newDF$cellBC[newDF$condition == "sample2"],]

sample3_ass_SUB_test = sample3_ass[sample3_ass $revComp %in% newDF$cellBC[newDF$condition == "sample3"],]

sample4_ass_SUB_test = sample4_ass[sample4_ass $revComp %in% newDF$cellBC[newDF$condition == "sample4"],]

sample5_ass_SUB_test = sample5_ass[sample5_ass $revComp %in% newDF$cellBC[newDF$condition == "sample5"],]

sample6_ass_SUB_test = sample6_ass[sample6_ass $revComp %in% newDF$cellBC[newDF$condition == "sample6"],]

sample7_ass_SUB_test = sample7_ass[sample7_ass $revComp %in% newDF$cellBC[newDF$condition == "sample7"],]

sample8_ass_SUB_test = sample8_ass[sample8_ass $revComp %in% newDF$cellBC[newDF$condition == "sample8"],]

newDF$vTR = "None"

for(i in 1:nrow(newDF)){
  if(newDF$condition[i] == "sample1"){
    sub = sample1_ass_SUB$vTR[sample1_ass_SUB$revComp == newDF$cellBC[i]]
    if(length(sub) == 1){ newDF$vTR[i] = sub }
  }
  if(newDF$condition[i] == "sample2"){
    sub = sample2_ass_SUB$vTR[sample2_ass_SUB$revComp == newDF$cellBC[i]]
    if(length(sub) == 1){ newDF$vTR[i] = sub }
  }
  if(newDF$condition[i] == "sample3"){
    sub = sample3_ass_SUB$vTR[sample3_ass_SUB$revComp == newDF$cellBC[i]]
    if(length(sub) == 1){ newDF$vTR[i] = sub }
  }
  if(newDF$condition[i] == "sample4"){
    sub = sample4_ass_SUB$vTR[sample4_ass_SUB$revComp == newDF$cellBC[i]]
    if(length(sub) == 1){ newDF$vTR[i] = sub }
  }
  if(newDF$condition[i] == "sample5"){
    sub = sample5_ass_SUB$vTR[sample5_ass_SUB$revComp == newDF$cellBC[i]]
    if(length(sub) == 1){ newDF$vTR[i] = sub }
  }
  if(newDF$condition[i] == "sample6"){
    sub = sample6_ass_SUB$vTR[sample6_ass_SUB$revComp == newDF$cellBC[i]]
    if(length(sub) == 1){ newDF$vTR[i] = sub }
  }
  if(newDF$condition[i] == "sample7"){
    sub = sample7_ass_SUB$vTR[sample7_ass_SUB$revComp == newDF$cellBC[i]]
    if(length(sub) == 1){ newDF$vTR[i] = sub }
  }
  if(newDF$condition[i] == "sample8"){
    sub = sample8_ass_SUB$vTR[sample8_ass_SUB$revComp == newDF$cellBC[i]]
    if(length(sub) == 1){ newDF$vTR[i] = sub }
  }
}

newDF$vTR_modified = newDF$vTR
newDF$vTR_modified[is.na(newDF$vTR_modified) | newDF$vTR_modified == "None"] = "Unknown"
FF_proj_AllSamples_Filt@cellColData$vTR = newDF$vTR_modified


plt_vTR_color <- plotEmbedding(ArchRProj = FF_proj_AllSamples_Filt, colorBy = "cellColData", name = "vTR", embedding = "UMAP", labelMeans = FALSE)

plotPDF(plt_vTR_color, name = "FullvTR-UMAP-Colors-MyCellRanger.pdf",
        ArchRProj = FF_proj_AllSamples_Filt, addDOC = FALSE, width = 5, height = 5)

FF_proj_AllSamples_Filt_known <- FF_proj_AllSamples_Filt[which(FF_proj_AllSamples_Filt@cellColData$vTR != "Unknown"),]

plt_vTR_color_known <- plotEmbedding(ArchRProj = FF_proj_AllSamples_Filt_known, colorBy = "cellColData", name = "vTR", embedding = "UMAP", labelMeans = FALSE)

plotPDF(plt_vTR_color_known, name = "FullvTR-UMAP-Colors-KnownOnly-MyCellRanger.pdf",
        ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE, width = 5, height = 5)

df_forUMAP <- FF_proj_AllSamples_Filt@embeddings$UMAP$df

write.csv(df_forUMAP, file ="FullvTR_df_forUMAP_MyCellRanger.csv")


new_df <- data.frame(FF_proj_AllSamples_Filt_known@cellColData$vTR, FF_proj_AllSamples_Filt_known@cellColData$Clusters)
write.csv(new_df, file = "vTR_cluster_identity_df.csv")

plt_tst <- plotEmbedding(ArchRProj = FF_proj_AllSamples_Filt_known, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", labelMeans = FALSE)

plotPDF(plt_tst, name = "MAP-ByClusters-KnownOnly-MyCellRanger.pdf",
        ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE, width = 5, height = 5)



pathToMacs2 <- "/home/jecorban/miniconda3/bin/macs3"

FF_proj_AllSamples_Filt_known <- addGroupCoverages(ArchRProj = FF_proj_AllSamples_Filt_known, groupBy = "vTR")

FF_proj_AllSamples_Filt_known <- addReproduciblePeakSet(
  ArchRProj = FF_proj_AllSamples_Filt_known, 
  groupBy = "vTR", 
  pathToMacs2 = pathToMacs2,
)

#PercentinClusters csv making

p1 <- plotEmbedding(ArchRProj = FF_proj_AllSamples_Filt_known, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

plotPDF(p1, name = "Plot-UMAP-vTRClusters.pdf", ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE, width = 5, height = 5)

test_df <- as.data.frame(FF_proj_AllSamples_Filt_known@cellColData)

test_df2 <- test_df[,c(16,20)]

library(dplyr)

total <- test_df2 %>% group_by(Clusters) %>% summarize(clust_total_cells = n())

grouped <- test_df2 %>% group_by(Clusters,vTR) %>% summarize(number = n())

merge_df <- merge(total, grouped, by.x = "Clusters", by.y = "Clusters", all.y = TRUE)

merge_df$PercentvTRPerCluster <- merge_df$number / merge_df$clust_total_cells

write.csv(merge_df, file = "PercentvTRInEachCluster.csv")

total2 <- test_df2 %>% group_by(vTR) %>% summarize(vTR_total_cells = n())

merge_df2 <- merge(total2, grouped, by.x = "vTR", by.y = "vTR", all.y = TRUE)

merge_df2$PercentofClusterInvTR<- merge_df2$number / merge_df2$vTR_total_cells

write.csv(merge_df2, file = "PercentClusterInEachvTR.csv")



##ClusterPercentagePlot

# using simple calculation
ggplot(merge_df2, aes(x = , y = value / sum(value) * 100, fill = factor(item))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Sector", y = "Percentage", title = "Stacked 100% Bar Plot by Sector") + coord_flip()


##


FF_proj_AllSamples_Filt_known <- addPeakMatrix(ArchRProj = FF_proj_AllSamples_Filt_known)

FF_proj_AllSamples_Filt_known <- FF_proj_AllSamples_Filt_known[which(FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR12" & FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR55" & FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR81" & FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR82" & FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR83" & FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR84" & FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR85" & FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR92" & FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR101" & FF_proj_AllSamples_Filt_known@cellColData$vTR != "J_vTR102"),]


markersPeaks <- getMarkerFeatures(
  ArchRProj = FF_proj_AllSamples_Filt_known, 
  useMatrix = "PeakMatrix", 
  groupBy = "vTR",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

##Stop here if just using ArchRProj to make comparison volcano plots

test_heatmap_mat <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,
  nLabel = 0,
  binaryClusterRows = TRUE,
  returnMatrix = TRUE
)

write.csv(test_heatmap_mat, "FullvTR_Peak_Heatmap_Matrix_MyCellRanger.csv", row.names=T)

test_heatmap <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,
  nLabel = 0,
  binaryClusterRows = TRUE
)


# Making all pairwise comparisons 

devtools::install_github("GreenleafLab/chromVARmotifs")

FF_proj_AllSamples_Filt_known <- addMotifAnnotations(ArchRProj = FF_proj_AllSamples_Filt_known, motifSet = "homer", name = "Motif")

vTR_list <- c(unique(FF_proj_AllSamples_Filt_known@cellColData$vTR))



comp_list <- vector(mode='list', length=length(vTR_list))
pma_list <- vector(mode='list', length=length(vTR_list))
pv_list <- vector(mode='list', length=length(vTR_list))
motif_up_list <- vector(mode='list', length=length(vTR_list))
motif_down_list <- vector(mode='list', length=length(vTR_list))
df_down_combo <- NULL
df_up_combo <- NULL

for(i in 1:length(vTR_list)){
  comp_list[[i]] <- getMarkerFeatures(
    ArchRProj = FF_proj_AllSamples_Filt_known, 
    useMatrix = "PeakMatrix",
    groupBy = "vTR",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = vTR_list[[i]],
    bgdGroups = "Empty"
  )
  
pma_list[[i]] <- plotMarkers(seMarker = comp_list[[i]], name = vTR_list[i], cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
  
pv_list[[i]] <- plotMarkers(seMarker = comp_list[[i]], name = vTR_list[i], cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
  
motifsUp <- peakAnnoEnrichment(
  seMarker = comp_list[[i]],
  ArchRProj = FF_proj_AllSamples_Filt_known,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
  
  temp_df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
  temp_df <- temp_df[order(temp_df$mlog10Padj, decreasing = TRUE),]
  temp_df$rank <- seq_len(nrow(temp_df))
  temp_df$vTR <- vTR_list[[i]]
  df_up_combo <- rbind(df_up_combo, temp_df)
  
  rm(temp_df)
  
  motifsDo <- peakAnnoEnrichment(
    seMarker = comp_list[[i]],
    ArchRProj = FF_proj_AllSamples_Filt_known,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )
  temp_df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
  temp_df <- temp_df[order(temp_df$mlog10Padj, decreasing = TRUE),]
  temp_df$rank <- seq_len(nrow(temp_df))
  temp_df$vTR <- vTR_list[[i]]
  
  df_down_combo <- rbind(df_down_combo, temp_df)
  rm(temp_df)
  
}

write.csv(df_down_combo, file = "FullvTR_AllComparisons_DownMotifs_Dataframe.csv")

write.csv(df_up_combo, file = "FullvTR_AllComparisons_UpMotifs_Dataframe.csv")

pma_df_combo <- NULL
pma_df_combo <- as.data.frame(pma_df_combo)
pma_df_temp <- NULL

for(i in 1:length(pma_list)){
  pma_df_temp$Log2FC <- comp_list[[i]]@assays@data$Log2FC
  pma_df_temp$FDR <- comp_list[[i]]@assays@data$FDR
  pma_df_temp$Mean <- comp_list[[i]]@assays@data$Mean
  pma_df_temp$vTR <- vTR_list[[i]]
  pma_df_temp <- as.data.frame(pma_df_temp)
  colnames(pma_df_temp) <- c("Log2FC", "FDR", "Mean", "vTR")
  pma_df_combo <- rbind(pma_df_combo, pma_df_temp)
  rm(pma_df_temp)
  pma_df_temp <- NULL
}

write.csv(pma_df_combo, file = "FullvTR_AllComparisons_PMAandPV_Dataframe.csv")


##QC plots

df <- getCellColData(FF_proj_AllSamples_Filt_known, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(4.5, 5.05),
  ylim = c(6, 20)
) + geom_hline(yintercept = 7, lty = "dashed") + geom_vline(xintercept = 4.6, lty = "dashed")

plotPDF(p, name = "260313_TSS-vs-Frags-FullvTR-OnlyKnownGenotypeCells.pdf", ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE)


p2 <- plotGroups(
  ArchRProj = FF_proj_AllSamples_Filt_known, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(p2, name = "260313_ViolinPlot-TSS-FullvTR-NoFiltering.pdf", ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE)

p3 <- plotGroups(
  ArchRProj = FF_proj_AllSamples_Filt_known, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

plotPDF(p3, name = "260313_RidgePlot-Frags-FullvTR-NoFiltering.pdf", ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE)

p4 <- plotGroups(
  ArchRProj = FF_proj_AllSamples_Filt_known, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(p4, name = "260313_ViolinPlot-Frags-FullvTR-NoFiltering.pdf", ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE)

p1 <- plotFragmentSizes(ArchRProj = FF_proj_AllSamples_Filt_known)
p2 <- plotTSSEnrichment(ArchRProj = FF_proj_AllSamples_Filt_known)
plotPDF(p1,p2, name = "260313_C-Sample-FragSizes-TSSProfile-FullvTR.pdf", ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE, width = 5, height = 5)

#With vTRs
p4 <- plotGroups(
  ArchRProj = FF_proj_AllSamples_Filt_known, 
  groupBy = "vTR", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(p4, name = "ViolinPlot-Frags-FullvTR-NoFiltering-ByvTR.pdf", ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE)

p2 <- plotGroups(
  ArchRProj = FF_proj_AllSamples_Filt_known, 
  groupBy = "vTR", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(p2, name = "ViolinPlot-TSS-FullvTR-NoFiltering-ByvTR.pdf", ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE)

qc_df <- data.frame(FF_proj_AllSamples_Filt_known@cellColData$vTR, FF_proj_AllSamples_Filt_known@cellColData$TSSEnrichment, FF_proj_AllSamples_Filt_known@cellColData$nFrags, FF_proj_AllSamples_Filt_known@cellColData$ReadsInPeaks)

write.csv(qc_df, file ="260313_FF_FullvTR_QC_Dataframe.csv")

##MakingQCFigures

qc_df <- fread("260313_FullvTR_QC_Dataframe.csv")
qc_df <- qc_df[,-1]
colnames(qc_df)<- c("vTR", "TSS_Enrichment", "nFrags", "Reads_in_Peaks")
TSS_plot <- ggplot(qc_df, aes(x=vTR, y=TSS_Enrichment, fill = vTR)) + 
  geom_violin() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"),
                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                        legend.position = "none") + xlab("")
ggsave(TSS_plot, filename = "260313_FullvTR_TSS_Violin_Plot.pdf", 
       width = 20, height = 8, units ="in")


nFrag_plot <- ggplot(qc_df, aes(x=vTR, y=nFrags, fill = vTR)) + 
  geom_violin() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"),
                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                        legend.position = "none") + xlab("")
ggsave(nFrag_plot, filename = "260313_FullvTR_nFrag_Violin_Plot.pdf", 
       width = 20, height = 8, units ="in")

ReadsinPeaks_plt <- ggplot(qc_df, aes(x=vTR, y=Reads_in_Peaks, fill = vTR)) + 
  geom_violin() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"),
                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                        legend.position = "none") + xlab("") + ylab("Reads in Peaks")
ggsave(ReadsinPeaks_plt, filename = "260313_FullvTR_ReadsInPeaks_Violin_Plot.pdf", 
       width = 20, height = 8, units ="in")


# Check that IterativeLSI is available
FF_proj_AllSamples_Filt_known <- addCoAccessibility(
  ArchRProj = FF_proj_AllSamples_Filt_known,
  reducedDims = "IterativeLSI"
)


##Here we set correlation cutoffs and resolution for loops. Higher (eg. 10000) reduced overplotting, 100-1000 is good depending on range selected
cA <- getCoAccessibility(
  ArchRProj = FF_proj_AllSamples_Filt_known,
  corCutOff = 0.5,
  resolution = 2500,
  returnLoops = TRUE
)


markerGenes <- c("RET", "DUSP6", "CCND1")

markerGenes <- c("MET")

markerGenes <- c("TP73")


p <- plotBrowserTrack(
  ArchRProj = FF_proj_AllSamples_Filt_known, 
  groupBy = "vTR", 
  geneSymbol = markerGenes_setfull, 
  features = getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
  upstream = 70000,
  downstream = 70000,
  loops = getCoAccessibility(FF_proj_AllSamples_Filt_known)
)


