library(BiocManager)
library(ArchR)
library(httpgd)
library(pheatmap)
library(devtools)
set.seed(1)
 
addArchRThreads(threads = 16)
 
inputFiles <- getTutorialData("Hematopoiesis")
 
addArchRGenome("hg38")
 

 
names(inputFiles1) <-c("w21_dc1r3_r1","w21_dc2r2_r1","w21_dc2r2_r2")
 
#create arrow file
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles1,
  sampleNames = names(inputFiles1),
  minTSS  = 5, #Dont set this too high because you can always increase later
  minFrags = 700,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles
 


#creating an ArchR project
project2 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "data",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

#add doublet score
doubScores <- addDoubletScores(
  input = project2,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

#### Collect all samples into a joint data structure
project2

### Save the project in Raw folder after adding the Doublet score
saveArchRProject(ArchRProj= project2 , outputDirectory = "Raw")

getCellColData(project2)

##### 1.5 Quality control
quality <- getCellColData(project2, select = c("log10(nFrags)", "TSSEnrichment"))
quality
quality_plot <- ggPoint(
    x = quality[,1], 
    y = quality[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(quality[,1], probs = 0.99)),
    ylim = c(0, quantile(quality[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
quality_plot

ggsave(
  filename = "path.../quality_plot.png", 
  plot = quality_plot,           # The plot object
  width = 10,              # Width in inches
  height = 8,              # Height in inches
  dpi = 300                # Resolution (DPI)
)

#cells number of sample 1
idxSample1 <- BiocGenerics::which(project2$Sample %in% "w21_dc1r3_r1")
cellsSample <- project2$cellNames[idxSample1]
project2[cellsSample, ]

#cells number of sample 2
idxSample2 <- BiocGenerics::which(project2$Sample %in% "w21_dc2r2_r1")
cellsSample2 <- project2$cellNames[idxSample2]
project2[cellsSample2, ]

#cells number of sample 3
idxSample3 <- BiocGenerics::which(project2$Sample %in% "w21_dc2r2_r2")
cellsSample3 <- project2$cellNames[idxSample3]
project2[cellsSample3, ]


#plotting sample statistics 
Plot_1 <- plotGroups(
    ArchRProj = project2, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
png("path..../plot1.png", 
    width = 1024,    # Set width (in inches)
    height = 800,    # Set height (in inches)
    res = 300)
plot(Plot_1)  # Replace with your plotting code
dev.off()   


plot_2 <- plotGroups(
    ArchRProj = project2, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
ggsave(
  filename = "path.../plot2.png", 
  plot = plot_2,           # The plot object
  width = 10,              # Width in inches
  height = 8,              # Height in inches
  dpi = 300                # Resolution (DPI)
)


# plotting Fragment Size Distribution and TSS Enrichment
Fragemnt_plot <- plotFragmentSizes(ArchRProj = project2)
ggsave(
  filename = "path.../plot3.png", 
  plot = Fragemnt_plot,           # The plot object
  width = 10,              # Width in inches
  height = 8,              # Height in inches
  dpi = 300                # Resolution (DPI)
)

TSS_plot <- plotTSSEnrichment(ArchRProj = project2)
ggsave(
  filename = "path.../plot4.png", 
  plot = TSS_plot,           # The plot object
  width = 10,              # Width in inches
  height = 8,              # Height in inches
  dpi = 300                # Resolution (DPI)
)


##### 1.6 Filter the dataset
project2 <- filterDoublets(doubScores)
project2 <- project2[which(project2$TSSEnrichment >7 & project2$nFrags > 1000 & project2$nFrags < 50000)]
project2
getCellColData(project2)

##### Week 2
####2 Dimensionality Reduction (5P)
###2.1 Iterative LSI

project2 <- addIterativeLSI(
    ArchRProj = project2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)
#add the Umap and plot it for sample/Tss and nfragment
project2 <- addUMAP(
    ArchRProj = project2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
 
#plot UMAP before correction
sample_plot_Umap <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
TSS_plot_Umap <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
nFrag_plot_Umap <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")

ggsave(
  filename = "path.../nFrag_plot_Umap.png", 
  plot = nFrag_plot_Umap,           # The plot object
  width = 10,              # Width in inches
  height = 8,              # Height in inches
  dpi = 300                # Resolution (DPI)
)

###2.3 Dealing with batch effects

project2 <- addHarmony(
    ArchRProj = project2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)
project2 <- addUMAP(
    ArchRProj = project2, 
    reducedDims = "Harmony", 
    name = "HarmonyUMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

#plot UMAP after correction
sample_plot_UmapAF <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
TSS_plot_UmapAF <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
nFrag_plot_UmapAF <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")

ggsave(
  filename = "path.../sample_plot_UmapAF.png", 
  plot = sample_plot_UmapAF,           # The plot object
  width = 10,              # Width in inches
  height = 8,              # Height in inches
  dpi = 300                # Resolution (DPI)
)

#####3 Clustering 
#add clustring using Seurats
project2 <- addClusters(
    input = project2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

table(project2$Clusters)

#plot UMAP clustring 
clustring_Umap <- plotEmbedding(ArchRProj = project2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggsave(
  filename = "path.../clustring_Umap.png", 
  plot = clustring_Umap,           # The plot object
  width = 10,              # Width in inches
  height = 8,              # Height in inches
  dpi = 300                # Resolution (DPI)
)
clustring_Umap
#confusionMatrix for heat map and plot it 
Heat_matrix <- confusionMatrix(paste0(project2$Clusters), paste0(project2$Sample))
Heat_matrix

##### 4.1 Peaks calling

project3 <- addGroupCoverages(ArchRProj = project2, groupBy = "Clusters")
PeakCallspath <- ("path.../data/PeakCalls")

projHemeTmp <- addReproduciblePeakSet(
    ArchRProj = project3, 
    groupBy = "Clusters",
    peakMethod = "Tiles",
    pathToMacs2 = PeakCallspath,
    threads=1
)

#Adding matrix
project4 <- addPeakMatrix(projHemeTmp)
getAvailableMatrices(project4)


markersPeaks  <- getMarkerFeatures(
    ArchRProj = project4, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  threads=1
)

markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

heat_map_marker <- draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

png("heatmap_output.png", width = 800, height = 600)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

#Plot genes
plot_TOP2A<- plotBrowserTrack(
    ArchRProj = project4, 
    groupBy = "Clusters", 
    geneSymbol = c("TOP2A"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)

png("plot_TOP2A.png", width=800, height=400)
 
grid::grid.newpage()
grid::grid.draw(plot_TOP2A$TOP2A)
dev.off()


plot_MKI67<- plotBrowserTrack(
    ArchRProj = project4, 
    groupBy = "Clusters", 
    geneSymbol = c("MKI67"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)

png("plot_MKI67.png", width=800, height=400)
 
grid::grid.newpage()
grid::grid.draw(plot_MKI67$MKI67)
dev.off()


plot_AURKA<- plotBrowserTrack(
    ArchRProj = project4, 
    groupBy = "Clusters", 
    geneSymbol = c("AURKA"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)

png("plot_AURKA.png", width=800, height=400)
 
grid::grid.newpage()
grid::grid.draw(plot_AURKA$AURKA)
dev.off()


plot_SATB2<- plotBrowserTrack(
    ArchRProj = project4, 
    groupBy = "Clusters", 
    geneSymbol = c("SATB2"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)

png("plot_SATB2.png", width=800, height=400)
 
grid::grid.newpage()
grid::grid.draw(plot_SATB2$SATB2)
dev.off()

plot_SLC12A7<- plotBrowserTrack(
    ArchRProj = project4, 
    groupBy = "Clusters", 
    geneSymbol = c("SLC12A7"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)

png("plot_SLC12A7.png", width=800, height=400)
 
grid::grid.newpage()
grid::grid.draw(plot_SLC12A7$SLC12A7)
dev.off()

#Week 3
#Identify Marker genes
markers_genes <- getMarkerFeatures(
    ArchRProj = project4, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon" , threads=1
)

markers_for_genes <- getMarkers(markers_genes, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")#we used the parameter for the gene markers we got before for each cluster
markers_for_genes

heatmap_MG <- markerHeatmap(
  seMarker = markers_genes, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  transpose = TRUE
)

First5_MG  <- c("PERM1",  "MIR6726", "VWA1", "LINC01770",  "PANK4")

plot_First5_MG  <- plotEmbedding(
    ArchRProj = project4, 
    colorBy = "GeneScoreMatrix", 
    name = First5_MG, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)


library(patchwork)

combined_plot <- wrap_plots(plot_First5_MG, ncol = 2) # Combine plots into a grid
ggsave(
  filename = "path.../plot_First5_MG_combined.png", 
  plot = combined_plot,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)

#Using MAGIC
project5 <- addImputeWeights(project4)

MAGIC_plot <- plotEmbedding(
    ArchRProj = project5, 
    colorBy = "GeneScoreMatrix", 
    name = First5_MG, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(project5)
)

MAGIC_plot <- wrap_plots(MAGIC_plot, ncol = 3) # Combine plots into a grid
ggsave(
  filename = "path.../MAGIC_plot.png", 
  plot = MAGIC_plot,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)

#Transcription Factor motif activity
library(chromVARmotifs)
project6 <- addMotifAnnotations(ArchRProj = project5, motifSet = "cisbp",genome = "hg38")

project6 <- addBgdPeaks(project6)

project6 <- addDeviationsMatrix(
  ArchRProj = project6, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev_proj6 <- getVarDeviations(project6, name = "MotifMatrix", plot = TRUE)

getAvailableMatrices(project6)

motifs <- c("TAL1_62","TAL2_822")
markerMotifs <- getFeatures(project6, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% motifs]
markerMotifs

# 6.2 Plot UMAP embeddings for marker TFs
motifs_genes <- getMarkerFeatures(
    ArchRProj = project6, 
    useMatrix = "MotifMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon" , threads=1
)

motifs_plot <- plotEmbedding(
    ArchRProj = project6, 
    colorBy = "MotifMatrix", 
    name = markerMotifs, 
    embedding = "UMAP",
    imputeWeights =  getImputeWeights(project6)
)

motifs_plot <- wrap_plots(motifs_plot, ncol = 2) # Combine plots into a grid
ggsave(
  filename = "path.../motifs_plot.png", 
  plot = motifs_plot,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)
# 6.3 Motif activity
p_distribution <- plotGroups(ArchRProj = project6, 
  groupBy = "Clusters", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(project6)
)

p_distribution_plot <- wrap_plots(p_distribution, ncol = 2) # Combine plots into a grid
ggsave(
  filename = "path.../p_distribution_plot.png", 
  plot = p_distribution_plot,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)

#7.1 Data integration

seRNA <- readRDS("new_pbmc.rds")
seRNA
project6 <- addGeneIntegrationMatrix(
    ArchRProj = project6, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force = TRUE,
    groupRNA = "Cluster.Name",  
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)
getAvailableMatrices(project6)

project6 <- addImputeWeights(project6)

markerGenes  <- c(
    "ID4",
    "NFIA", 
    "ASCL1"
  )

p1 <- plotEmbedding(
    ArchRProj = project6, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    continuousSet = "horizonExtra",
    embedding = "UMAP",
    imputeWeights = getImputeWeights(project6)
)

p2 <- plotEmbedding(
    ArchRProj = project6,
    colorBy = "GeneScoreMatrix",
    continuousSet = "horizonExtra",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(project6)
)
p1c <- wrap_plots(p1, ncol = 3) # Combine plots into a grid
ggsave(
  filename = "path.../p1.png", 
  plot = p1c,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)

p2c <- wrap_plots(p2, ncol = 3) # Combine plots into a grid
ggsave(
  filename = "path.../p2.png", 
  plot = p2c,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)
#Correlation Coefficients
corGIM_MM <- correlateMatrices(
    ArchRProj = project6,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    threads=1
)

#7.2 Correlation Coefficients
cM <- confusionMatrix(project6$Clusters, project6$predictedGroup)
labelOld <- rownames(cM)
labelOld
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
project6$Clusters <- mapLabels(project6$Clusters, newLabels = labelNew, oldLabels = labelOld)
p1 <- plotEmbedding(project6, colorBy = "cellColData", name = "Clusters")
p1
order_vale <- sort(corGIM_MM)
head(order_vale,4)
tail(order_vale,4)

#peak-gene linkage 
project6 <- addPeak2GeneLinks(
    ArchRProj = project6,
    reducedDims = "IterativeLSI",
    threads=1
)

peak_to_gene <- getPeak2GeneLinks(
    ArchRProj = project6,
    corCutOff = 0.45,
    resolution = 1000,
    returnLoops = TRUE
)

plot_peak_to_gene <- plotBrowserTrack(
    ArchRProj = project6, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(project6)
)

plotPDF(plotList = plot_peak_to_gene, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
    ArchRProj = project6, 
    addDOC = FALSE, width = 10, height = 10)

heatmap_peak_to_gene <- plotPeak2GeneHeatmap(ArchRProj = project6, groupBy = "Clusters")

png("heatmap_peak_to_gene.png", width = 800, height = 600)
draw(heatmap_peak_to_gene, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()
#----------------------------------------------------------------------------------------------------------------------------
#Week 4
#9.1 Differential peak accessibility
Mono_T_mrkertest <- getMarkerFeatures(
  ArchRProj = project6, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "GluN2",
  bgdGroups = "IN2"
)

MA <- markerPlot(seMarker = Mono_T_mrkertest, name = "GluN2", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
ggsave(
  filename = "path.../MA.png", 
  plot = MA,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)

Vol <- markerPlot(seMarker = Mono_T_mrkertest, name = "GluN2", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
ggsave(
  filename = "path.../Vol.png", 
  plot = Vol,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)

#9.2 TF motif enrichment
motifsUp <- peakAnnoEnrichment(
    seMarker = Mono_T_mrkertest,
    ArchRProj = project6,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

head(df)
gg_motifsup <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggsave(
  filename = "path.../gg_motifsup.png", 
  plot = gg_motifsup,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)


motifsDo <- peakAnnoEnrichment(
    seMarker = Mono_T_mrkertest,
    ArchRProj = project6,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

motifsDo

dfDo <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
dfDo <- dfDo[order(dfDo$mlog10Padj, decreasing = TRUE),]
dfDo$rank <- seq_len(nrow(dfDo))

head(dfDo)

gg_motifsDo <- ggplot(dfDo, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dfDo[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggsave(
  filename = "path.../gg_motifsDo.png", 
  plot = gg_motifsDo,
  width = 20,  # Adjust dimensions as needed
  height = 16,
  dpi = 300
)


Mono_TM <- getMarkerFeatures(

  ArchRProj = project6,

  useMatrix = "PeakMatrix",

  groupBy = "Clusters",

  testMethod = "wilcoxon",

  bias = c("TSSEnrichment", "log10(nFrags)"),threads=1)


enrichMotifs <- peakAnnoEnrichment(

    seMarker = markersPeaks,

    ArchRProj = project6,

    peakAnnotation = "Motif",

    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",

    background = "all"

)

enrichMotifs

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

png("heatmapEM.png", width = 800, height = 600)
draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

motifPositions <- getPositions(project6)

