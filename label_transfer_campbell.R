# Label transfer to call celltypes using campbell hypothalamus dataset
library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)

############### Data edit for Seurat obj ####################
dat <- fread("./dat/merged_starmap_Nef_cellsrenamed_total.coutt.tsv")
meta <- fread("./dat/tomcells_metadata.csv")
rownames(meta) <- meta$V1
# Edit gene column
dat %>% separate(V1, into = c('all', 'pro_code'), sep = '__') %>% separate(all, into = c('gene', 'protein'), sep = '_') -> dat

# Restrict to only ensemble ids
dat <- dat[grepl('ENSMUS', gene),]
dat[1:10,1:10]

dat %>% select(matches('tom')) -> dat_edit 

# Set genecolumn to rownames
dat_edit <- as.sparse(dat_edit)
rownames(dat_edit) <- dat$gene

dat_edit[1:10,1:10]

############### Seurat workflow ####################
seur_obj <- CreateSeuratObject(counts = dat_edit)
seur_obj <- AddMetaData(seur_obj, meta)
seur_obj@meta.data

# run sctransform
seur_obj <- SCTransform(seur_obj)

# save object
save(seur_obj, file = "./dat/tomcells_r_adan_seurat_obj.RData")

############### Campbell edit ####################
# loading campbell dataset
camp.dat <- fread('/scratch/data-for_fast_access/pub-others/campbell2017/campbell.umi.csv.gz')
#### Edit to ensembl ids
# Change var.featues for campbell into ensembl ids
id2ensembl <- fread('/home/cbmr/qwn903/ygg-projects/amj/tune/timshel-bmicelltypes/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz')

# Translate to mouse ensembl
colnames(camp.dat)[1] <- 'gene_name_optimal'
id_join <- left_join(camp.dat, id2ensembl)
camp.dat <- id_join[, c(20923, 2:20922)] %>% na.omit %>% column_to_rownames('ensembl_gene_id')


cam.meta <- fread('/scratch/data-for_fast_access/pub-others/campbell2017/campbell.cell_metadata.csv')  %>% 
  select(-nGene,-nUMI, -orig.ident) %>% 
  column_to_rownames('cell_id')


cam_obj <- CreateSeuratObject(counts = camp.dat, meta.data = cam.meta)
cam_obj <- SCTransform(cam_obj)

# Save object
save(cam_obj, file = "./dat/campbell_seurat_obj.RData")

############### Label transfer - Campbell ####################
load('./dat/campbell_seurat_obj.RData')
load('./dat/tomcells_r_adan_seurat_obj.RData')
# Remove cells with no cell-calling from R. Adan
seur_obj <- subset(seur_obj, cells = na.omit(seur_obj@meta.data$V1))

# Find labels from campbell
hyp.anchors <- FindTransferAnchors(reference = cam_obj, normalization.method = "SCT",
                                   query = seur_obj, dims = 1:30, reference.assay = "SCT", 
                                   query.assay = "SCT", reduction = "cca")
predictions_neurons <- TransferData(anchorset = hyp.anchors, 
                                    refdata = cam_obj$cell_type_all_lvl2, dims = 1:30, weight.reduction = "cca")
predictions_neurons2 <- TransferData(anchorset = hyp.anchors, 
                                     refdata = cam_obj$cell_type_all_lvl1, dims = 1:30, weight.reduction = "cca")
predictions_taxo <- TransferData(anchorset = hyp.anchors, 
                                 refdata = cam_obj$taxonomy_lvl2, dims = 1:30, weight.reduction = "cca")

# Editing predictions together
predictions_neurons %>% rownames_to_column("V1") %>% filter(prediction.score.max > 0.5) %>% 
  rename(max_pred_neuron1 = predicted.id,
         pred_score_neuron1 = prediction.score.max) %>%
  select(V1, max_pred_neuron1, pred_score_neuron1) -> pred1

predictions_neurons %>% rownames_to_column("V1") %>% filter(prediction.score.max > 0.85) %>% 
  rename(max_pred_neuron1.5 = predicted.id,
         pred_score_neuron1.5 = prediction.score.max) %>%
  select(V1, max_pred_neuron1.5, pred_score_neuron1.5) -> pred1.5

predictions_neurons2 %>% rownames_to_column("V1") %>% filter(prediction.score.max > 0.5) %>% 
  rename(max_pred_neuron2 = predicted.id,
         pred_score_neuron2 = prediction.score.max) %>%
  select(V1, max_pred_neuron2, pred_score_neuron2) -> pred2

predictions_taxo %>% rownames_to_column("V1") %>% filter(prediction.score.max > 0.5) %>% 
  rename(max_pred_taxo = predicted.id,
         pred_score_taxo = prediction.score.max) %>%
  select(V1, max_pred_taxo, pred_score_taxo) -> pred3

predictions <- left_join(pred1, pred2) %>% left_join(., pred3) %>% left_join(., pred1.5) %>% column_to_rownames('V1')


seur_obj <- AddMetaData(seur_obj, metadata = predictions)
head(seur_obj@meta.data)

# Save object
save(seur_obj, file = "./dat/tomcells_w.labels_r_adan_seurat_obj.RData")


############### PCA and visualization ####################
load('./dat/tomcells_w.labels_r_adan_seurat_obj.RData')

# These are now standard steps in the Seurat workflow for visualization and clustering
seur_obj <- RunPCA(seur_obj, verbose = F)
seur_obj <- RunUMAP(seur_obj, dims = 1:30)

seur_obj <- FindNeighbors(seur_obj, dims = 1:30)
seur_obj <- FindClusters(seur_obj, resolution = 2.5)

p1 <- DimPlot(seur_obj, label = TRUE, group.by = 'SCT_snn_res.2.5', label.size = 6, pt.size = 3, repel = T) + 
  theme(legend.position = "bottom") + ggtitle("Clustering res. 2.5")

p2 <- DimPlot(seur_obj, label = TRUE, group.by = 'cell_types', label.size = 6, pt.size = 3, repel = T) + 
  theme(legend.position = "bottom") + ggtitle("R Adan Celltypes")

p1 + p2

###### Adan vs campbell taxonomy
head(seur_obj@meta.data)
p3 <- DimPlot(seur_obj, label = T, group.by = 'max_pred_taxo', label.size = 6, pt.size = 3, repel = T) + 
  theme(legend.position = "bottom") + ggtitle("Campbell Taxonomy")

p2 + p3

######## Subset neuronal cluster 
# Neuronal clusters only
seur_obj_neuro <- subset(seur_obj, cells = seur_obj@meta.data[grepl('Neuron', seur_obj@meta.data$max_pred_taxo), 'V1'])
# These are now standard steps in the Seurat workflow for visualization and clustering
seur_obj_neuro <- RunPCA(seur_obj_neuro, verbose = F)
seur_obj_neuro <- RunUMAP(seur_obj_neuro, dims = 1:30)
seur_obj_neuro <- FindNeighbors(seur_obj_neuro, dims = 1:30)
seur_obj_neuro <- FindClusters(seur_obj_neuro, resolution = 1.7)

p_n1 <- DimPlot(seur_obj_neuro, group.by = "cell_types", label = T, label.size = 4, pt.size = 3, repel = T) + 
  theme(legend.position = "bottom")

p_n2 <- DimPlot(seur_obj_neuro, group.by = "max_pred_neuron1", label = T, label.size = 4, pt.size = 3, repel = T) + 
  theme(legend.position = "bottom")

###### Adan vs campbell neuronal 1
p4 <- DimPlot(seur_obj, label = T, group.by = 'max_pred_neuron2', label.size = 6, pt.size = 3, repel = T) + 
  theme(legend.position = "bottom") + ggtitle("Campbell Celltypes 1")

###### Adan vs campbell neuronal 2
p5 <- DimPlot(seur_obj, label = T, group.by = 'max_pred_neuron1', label.size = 6, pt.size = 3, repel = T) + 
  theme(legend.position = "bottom") + ggtitle("Campbell Celltypes 2")

###### Adan vs campbell neuronal 2 - max pred 85%
p6 <- DimPlot(seur_obj, label = T, group.by = 'max_pred_neuron1.5', label.size = 6, pt.size = 3, repel = T) + 
  theme(legend.position = "bottom") + ggtitle("Campbell Celltypes 2, max.pred 85%")

############### POMC visualization ####################
ex_1 <- DimPlot(seur_obj_neuro, group.by = "cell_types", pt.size = 3, label = T, repel = T) + NoLegend() + 
  (DimPlot(seur_obj_neuro, group.by = "max_pred_neuron1", pt.size = 3, label = T, repel = T) + NoLegend()) + 
  FeaturePlot(seur_obj_neuro, "ENSMUSG00000020660", pt.size = 3) 

ex_2 <- VlnPlot(seur_obj_neuro, features = c("ENSMUSG00000020660", "ENSMUSG00000115958", "ENSMUSG00000025400"), group.by = "max_pred_neuron1", assay = "RNA")

pdf("extra_plots.pdf", width = 15)
ex_1
ex_2
dev.off()


############### PDF output ####################

pdf('./Label_transfer_campbell.pdf', width = 17, height = 10)
# R_adan celltypes
p2 + p1 + plot_annotation(
  title = 'Seurat clustering comparison',
  subtitle = 'Comparing initial clustering from seurat - resolution 2.5 - with cell-calling from R. Adan')
# R_adan vs campbell taxonomy
p2 + p3 + plot_annotation(
  title = 'Cell-calling comparison 1',
  subtitle = 'Comparing cell-calling from R. Adan with Campbell taxonomy classification - max.pred 60%')
# R_adan vs campbell neuron 1
p2 + p4 + plot_annotation(
  title = 'Cell-calling comparison 2',
  subtitle = 'Comparing cell-calling from R. Adan with Campbell celltype 1 classification - max.pred 60%')
# R_adan vs campbell neuron 2
p2 + p5 + plot_annotation(
  title = 'Cell-calling comparison 3',
  subtitle = 'Comparing cell-calling from R. Adan with Campbell celltype 2 classification - max.pred 60%')
# R_adan vs campbell neuron 2, max predition 85%
p2 + p6 + plot_annotation(
  title = 'Cell-calling comparison 4',
  subtitle = 'Comparing cell-calling from R. Adan with Campbell celltype 2 classification - max.pred 85%')
# R_adan vs campbell neuronal sub-cluster
p_n1 + p_n2 + plot_annotation(
  title = 'Cell-calling comparison, Neuronal clusters only',
  subtitle = 'Comparing cell-calling from R. Adan with Campbell celltype 2 classification - max.pred 60%, neuronal subcluster')
dev.off()


pdf('./short_label.transfer_campbell.pdf', width = 17, height = 10)
# R_adan vs campbell taxonomy
p2 + p3 + plot_annotation(
  title = 'Cell-calling comparison 1',
  subtitle = 'Comparing cell-calling from R. Adan with Campbell taxonomy classification - max.pred 60%')
# R_adan vs campbell neuronal sub-cluster
p_n1 + p_n2 + plot_annotation(
  title = 'Cell-calling comparison, Neuronal clusters only',
  subtitle = 'Comparing cell-calling from R. Adan with Campbell celltype 2 classification - max.pred 60%, neuronal subcluster')
dev.off()









