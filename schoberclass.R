Idents(obj)<- "snn_0.4"
DimPlot(obj)
View(obj@meta.data)
Idents(obj)<- "snn_res.0.4"
DimPlot(obj)
library(celldex)
library(SingleR)
ref <- MouseRNAseqData()
pred <- SingleR(test= new@assays$RNA@data, ref = ref, assay.type.test = 1, labels = ref$label.fine )
View(obj@meta.data)
obj[["pred"]] <- pred$labels
View(obj@meta.data)
DimPlot(obj, group.by = "pred")

DimPlot(obj, group.by = "pred", split.by = "orig.ident")
pred.fine <- SingleR(test= obj@assays$RNA@data, ref = ref, assay.type.test = 1, labels = ref$label.fine )
obj[["pred.fine"]] <- pred.fine$labels
View(obj@meta.data)
DimPlot(obj, group.by = "pred.fine")
Idents(obj)<- "pred"
DimPlot(obj)
fb <- subset(obj, ident= "Fibroblasts"
fb
obj
DimPlot(fb, group.by = pred.fine)
DimPlot(fb, group.by = 'pred.fine')
obj[["pred.fine"]] <- pred.fine$labels
fb <- subset(obj, ident= "Fibroblasts"
fb
obj
DimPlot(fb, group.by = pred.fine)
DimPlot(fb, group.by = 'pred.fine')
obj[["pred.fine"]] <- pred.fine$labels
fbp<- SingleR(fb, ref = ref, assay.type.test = 1, labels = ref$label.fine)
fbp<- SingleR(fb@assays$RNA@data, ref = ref, assay.type.test = 1, labels = ref$label.fine)
fb[["pred_fine"]]<- fbp$labels
View(fb@meta.data)
Idents(fb) <- "pred_fine"
DimPlot(fb)
fb <- FindVariableFeatures(fb, selection.method = "vst", nfeatures = 2000)
fb <- ScaleData(fb, features = rownames(fb))
fb<- RunPCA(fb, features = VariableFeatures(fb))
fb<- FindNeighbors(fb, dims = 1:20)
fb<- FindClusters(fb, resolution = .4)
fb<- RunUMAP(fb, dims = 1:20)
DimPlot(fb)
predfb<- SingleR(fb, ref = ref, assay.type.test = 1, labels = ref$label.main)
predfb<- SingleR(fb@assays$RNA@data, ref = ref, assay.type.test = 1, labels = ref$label.main)
DimPlot(fb)
predfb<- SingleR(fb, ref = ref, assay.type.test = 1, labels = ref$label.main)
predfb<- SingleR(fb@assays$RNA@data, ref = ref, assay.type.test = 1, labels = ref$label.main)
fb[["pred_main"]] <- predfb$labels

#next - > Find DEGs List
#     - > Plot Dotplot/Heatmaps of DEGS by cluster/orig.ident
#     - > Rerun PredList


#Load in the Reference 
library(celldex)
ref <- celldex::MouseRNAseqData()
#reference <- 


#Prepare the reference using CHETAH's Protocol
celltypes_hn <- ref$label.main
counts_hn <- assay(ref)
class(counts_hn)
counts_hn[1:5, 1:5]

#Prepare Counts Matrix from Seurat Object
counts<- GetAssayData(object = obj)

#Prepare Dimensions Matrix from SCE Object
sce<- as.SingleCellExperiment(obj)
umap<- reducedDim(sce, type = 'UMAP')
input_mel1 <- SingleCellExperiment(assays = list(counts = counts), reducedDims = SimpleList(UMAP = umap))

#Follow Protocol
class(counts)
umap[1:5,]
all.equal(rownames(umap), colnames(counts))
headneck_ref <- SingleCellExperiment(assays = list(counts = counts_hn), colData = DataFrame(celltypes = celltypes_hn))
input_mel2 <- CHETAHclassifier(input = input_mel1, ref_cells = headneck_ref)


###Attempting GSEA###
library(escape)
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Seurat))
library(dittoSeq)
Idents(obj)<- "snn_res.0.4"
sce<- as.SingleCellExperiment(new)
GS <- getGeneSets(species = 'Mus musculus', library = "C7")

###Prepare the SC Object###
counts<- GetAssayData(object = new)
#View(counts)
umap<- reducedDim(sce, type = 'UMAP')
input_mel1 <- SingleCellExperiment(assays = list(counts = counts), reducedDims = SimpleList(UMAP = umap))
View(input_mel1)

###Enrich the Object###
ES <- enrichIt(obj = input_mel1, gene.sets = GS, groups = 1000, cores = 12)
View(ES)
obj<- AddMetaData(obj, ES)

###Set up Heatmap
colors <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
dittoHeatmap(object = obj, genes = NULL, metas = names(ES), annot.by = "groups", fontsize = 7, cluster_cols = TRUE, heatmap.colors = rev(colors(50)))






