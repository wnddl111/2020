
library(dplyr)

library(Seurat)

library(patchwork)

#�ʿ��� �͵��� Ȱ��ȭ ��Ų��

# Load the PBMC dataset

#data�� wd�� �ٿ� �޴� �� ���ϰ� �Ʒ� �ڵ忡 �ּҰ��� �װ� �°� �����Ѵ�.

pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")



'�� filtered_gene_bc_matrices/hg19/ ���͸� �ȿ��� �� �� ���� ������ ����� ����

�� barcodes.tsv : �� single cell���� ���ڵ� ����

�� genes.tsv : �ΰ� �������� ���ĸ�Ī�� ���� 

�� matrix.mtx : �� ���ڵ�� �����ڿ��� �� �������� �������� sparse matrix�� ��Ÿ�� ��



32738 1047 2572919



 32,738 : �ΰ� �������� ������ �������� �ʾƵ� ��

 1,047�� ���ڵ��� ������ ���迡 ���� �����ؾ� �ϴ� ����

 2,572,919�� 0�� �ƴ� �������� ������ ���迡 ���� �����ؾ� �ϴ� ����

'

#read10x �Լ��� cell ranger�� ����� �޾� umi ī��Ʈ ��Ʈ������ ��ȯ��

#�� ��Ʈ������ ���� ��-���� ��-�� �������� ����� Ư¡�� ���� ������ ���� ��Ÿ�� 

#cell ranger - �м� ����, ��� �� �� �ϳ��� ���ڵ� ��Ʈ������ ����



#����� umi count matrix ���� ����



# Initialize the Seurat object with the raw (non-normalized data).

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc

#���� �״�� �ʱ�ȭ�ϴ� �����ε� ������ count ������ �ʿ��ϰ� ������Ʈ���� ���÷� �����ִ� �� ����

# min.features���� �� ���ڴ� �۾ƾ� �����Ǵ� ������ ����/

# ��Ȯ���� ������ min.cells �� ����� ������ �� ���� 



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#���Ⱑ ���� qc�� �ϴ� �ܰ��ε� ppt�� ������ �� ó�� �����ܵ帮�Ƹ� �� ������ (scRNA seq������)

#�׷��� �����ܵ帮�� ������� ����ϴ� �Լ��� percentagefeatureset ���





# Visualize QC metrics as a violin plot

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Violin plot�� knerdensity plot�� ���� �߽����� �� ���̵忡 ǥ�� �� ��  

#knerdensity plot�� non-parametric density estimation �� �ϳ�

#����) �� ��� ��� ���� ������ ���� ���� �����ϰ� ������ �����͸����� Ȯ���е��Լ��� �����ؾ� �ϴµ� �̸� non-parametric density estimation�� �θ���.





# FeatureScatter is typically used to visualize feature-feature relationships, but can be used

# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.



plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot1 + plot2

#���踦 ��Ÿ���� ������? ���� �� ���� 





pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#�̰Ŵ� 200�� �̸� 2500�� �ʰ��ϴ� ������ ��� ������ ���� ���� ���͸� �ϴ� ��

#�̶� ���� ���� 5%�Ѵ� �� ������ / �� ���� ���� ������ ������Ű�� ���� ����

pbmc <- NormalizeData(pbmc)

'

������ �ʴ� ���� ������ �� ���� �ܰ� -> ����ȭ

scale factor ��⵵ ppt�� �־��ݾ� �̰� ����� ������?�� �ٷ�� ���ϱ� 

scaling �۾��� �ʿ��� 

���⼭�� lognormalize ��� ��� - ���� ���ϰ� �� ����� �α׷� ��ȯ�ϴ� ���� �����ϸ� ����ȭ ����̷�

 

 

�̰Ÿ��� �����ϰ� ���� �� �ִµ� �װ� default ���̶� �̷��� �ƹ��͵� ���� �ص� ����ȭ��

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

 

'



'

�� ������ �ϴ� �� �� �� ������ ū Ư¡ - ��𼭴� ���� ��𼭴� ����

 

'







pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#�ٸ� ��� 2000�� ǥ���ش޶�� �ǰ�?

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(pbmc), 10)



# plot variable features with and without labels

plot1 <- VariableFeaturePlot(pbmc)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# plot1 �����ٰ� ������ ǥ���� �� top 10�̰� repel�� ���� ���̴� ��

# ���⼭ ��� repel ������ ������ ���µ� ��Ȯ�� ������ �𸣰���

plot1 + plot2





# ���⼭ �ϴ� �Ŵ� ����ȭ ����� 0 �л��� 1�� �ǵ��� �� �������� ������ ������

# scaling �ϴ� ����/ ����ϰ� �������� �ϴϱ� 

all.genes <- rownames(pbmc)

pbmc <- ScaleData(pbmc, features = all.genes)

#�̰� ���ư��� ��Ȯ�� ������ �𸣰��� 







#����� pca�� �����ϴ� ��/ �̰� �ð�ȭ �ϴ� ���� 3���� ���� 

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#���⼭ dims�� dataframe�� ��� ���� ���̸� �ǹ� 

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")









# �ؿ��Ŵ� ������ ���� 3���� ��� �� �������ε� ������ �ҽ��� ���� Ž���� �� ���� 

# pca ������ ���� ���ĵȴ� 

DimPlot(pbmc, reduction = "pca")



DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)



DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)





#pca�� 12���� ������ ������ 12������ �ʿ��ѵ� �̶� �� ���� ��Ÿ���� 2���� ������ ��Ÿ���� �ǵ� �״ϱ� ������� 

#�� ����� �غ��ϱ� ���� pca������ ������� ���� Ŭ�����͸� ��

# NOTE: This process can take a long time for big datasets, comment out for expediency. More

# approximate techniques such as those implemented in ElbowPlot() can be used to reduce

# computation time

pbmc <- JackStraw(pbmc, num.replicate = 100)

#p value ������ ���� ����(�׸��� �ļ�) �� ���� �� �ִ� �ð�ȭ ������ ������ 

#�ٵ� �̰� ��Ȯ�� ���� �� �𸣰ھ� 

pbmc <- ScoreJackStraw(pbmc, dims = 1:20)



JackStrawPlot(pbmc, dims = 1:15)



ElbowPlot(pbmc)







#clustering�� �� �� ���� �����ϴٰ� �ư� �Ȱǵ�

#���� single cell rna seq�� ������ ���� ��Ȳ���� clustering�� �ؾ��ϴϱ� ������ �ּҰŸ��� �� �����  



pbmc <- FindNeighbors(pbmc, dims = 1:10)

pbmc <- FindClusters(pbmc, resolution = 0.5)



# Look at cluster IDs of the first 5 cells

head(Idents(pbmc), 5)





#���� ���� ���

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =

# 'umap-learn')

pbmc <- RunUMAP(pbmc, dims = 1:10)



#�̰͵� �ڼ��� �� �� �𸣰ڵ� 







# note that you can set `label = TRUE` or use the LabelClusters function to help label

# individual clusters

DimPlot(pbmc, reduction = "umap")

saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")









# find all markers of cluster 1

cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)

head(cluster1.markers, n = 5)



# find all markers distinguishing cluster 5 from clusters 0 and 3

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)

head(cluster5.markers, n = 5)



# find markers for every cluster compared to all remaining cells, report only the positive ones

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



VlnPlot(pbmc, features = c("MS4A1", "CD79A"))



# you can plot raw counts as well

VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)



FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               
                               "CD8A"))

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(pbmc, features = top10$gene) + NoLegend()







#��� ������ ��Ʈ�� ��� ǥ�� ��Ŀ�� ����Ͽ� ������� �ʴ� Ŭ�����͸��� �����Ͽ� �� ������ ���� ��ġ ����

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     
                     "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(pbmc) # ���� �� �̸��� �������� ���� 

#level�� factor���� ��ȯ�ϴ°� 

pbmc <- RenameIdents(pbmc, new.cluster.ids)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



saveRDS(pbmc, file = "../output/pbmc3k_final.rds")





































































































������

�ְ����� ���� 2 ������ü����