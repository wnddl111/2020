
library(dplyr)

library(Seurat)

library(patchwork)

#필요한 것들을 활성화 시킨다

# Load the PBMC dataset

#data는 wd에 다운 받는 게 편하고 아래 코드에 주소값은 그게 맞게 변형한다.

pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")



'⑴ filtered_gene_bc_matrices/hg19/ 디렉터리 안에는 총 세 개의 파일이 저장돼 있음

① barcodes.tsv : 각 single cell들의 바코드 네임

② genes.tsv : 인간 유전자의 정식명칭과 별명 

③ matrix.mtx : 각 바코드와 유전자에서 그 유전자의 발현량을 sparse matrix로 나타낸 것



32738 1047 2572919



 32,738 : 인간 유전자의 개수로 변경하지 않아도 됨

 1,047은 바코드의 개수로 실험에 따라 변경해야 하는 값임

 2,572,919는 0이 아닌 데이터의 개수로 실험에 따라 변경해야 하는 값임

'

#read10x 함수는 cell ranger의 출력을 받아 umi 카운트 매트릭스를 반환함

#이 매트릭스의 값는 열-세포 행-각 세포에서 검출된 특징에 대한 분자의 수를 나타냄 

#cell ranger - 분석 도구, 결과 값 중 하나로 바코드 매트릭스를 만듬



#결론은 umi count matrix 같은 거지



# Initialize the Seurat object with the raw (non-normalized data).

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc

#설명 그대로 초기화하는 과정인데 느낌상 count 파일이 필요하고 프로젝트명은 음시로 정해주는 것 같고

# min.features에서 이 숫자는 작아야 누락되는 세포가 적대/

# 정확하지 않지만 min.cells 도 비슷한 느낌인 것 같음 



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#여기가 이제 qc를 하는 단계인데 ppt에 정리한 것 처럼 미토콘드리아를 잘 봐야해 (scRNA seq에서는)

#그래서 미토콘드리아 백분율을 계산하는 함수인 percentagefeatureset 사용





# Visualize QC metrics as a violin plot

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Violin plot은 knerdensity plot을 축을 중심으로 각 사이드에 표시 한 것  

#knerdensity plot은 non-parametric density estimation 중 하나

#참고) 이 경우 어떠한 사전 정보나 지식 없이 순수하게 관측된 데이터만으로 확률밀도함수를 추정해야 하는데 이를 non-parametric density estimation라 부른다.





# FeatureScatter is typically used to visualize feature-feature relationships, but can be used

# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.



plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot1 + plot2

#관계를 나타내는 산점도? 같은 것 같음 





pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#이거는 200개 미만 2500개 초과하는 고유한 기능 개수를 가진 셀을 필터링 하는 것

#이때 미콘 수가 5%넘는 건 여과함 / 헷갈 ㄴㄴ 저기 기준이 여과시키고 싶은 기준

pbmc <- NormalizeData(pbmc)

'

원하지 않는 셀을 제거한 후 다음 단계 -> 정규화

scale factor 얘기도 ppt에 있었잖아 이게 희박한 데이터?를 다루다 보니까 

scaling 작업이 필요해 

여기서는 lognormalize 방법 사용 - 만을 곱하고 그 결과를 로그로 변환하는 전역 스케일링 정규화 방법이래

 

 

이거말고 복잡하게 쓰는 법 있는데 그게 default 값이라 이렇게 아무것도 없이 해도 정규화됨

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

 

'



'

이 다음에 하는 게 셀 간 변동이 큰 특징 - 어디서는 높고 어디서는 낮고

 

'







pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#다른 기능 2000개 표시해달라는 건가?

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(pbmc), 10)



# plot variable features with and without labels

plot1 <- VariableFeaturePlot(pbmc)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# plot1 위에다가 점으로 표시할 건 top 10이고 repel은 라벨을 붙이는 것

# 여기서 계속 repel 때문에 오류가 나는데 명확한 이유는 모르겠음

plot1 + plot2





# 여기서 하는 거는 정규화 평균을 0 분산을 1이 되도록 각 유전자의 발현을 조정함

# scaling 하는 거지/ 비슷하게 만들어줘야 하니까 

all.genes <- rownames(pbmc)

pbmc <- ScaleData(pbmc, features = all.genes)

#이게 돌아가는 정확한 원리는 모르겠음 







#여기는 pca를 수행하는 것/ 이걸 시각화 하는 법이 3가지 있음 

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#여기서 dims는 dataframe의 행과 열의 길이를 의미 

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")









# 밑에거는 위에서 말한 3가지 방법 중 마지막인데 이질성 소스를 쉽게 탐색할 수 있음 

# pca 점수에 따라서 정렬된대 

DimPlot(pbmc, reduction = "pca")



DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)



DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)





#pca는 12개의 변수가 있으면 12차원이 필요한데 이때 음 제일 나타내는 2개의 축으로 나타내는 건데 그니까 녹음들어 

#음 노이즈를 극복하기 위해 pca점수를 기반으로 셀을 클러스터링 함

# NOTE: This process can take a long time for big datasets, comment out for expediency. More

# approximate techniques such as those implemented in ElbowPlot() can be used to reduce

# computation time

pbmc <- JackStraw(pbmc, num.replicate = 100)

#p value 분포를 균일 분포(그림의 파선) 와 비교할 수 있는 시각화 도구를 제공함 

#근데 이게 정확히 뭔지 잘 모르겠어 

pbmc <- ScoreJackStraw(pbmc, dims = 1:20)



JackStrawPlot(pbmc, dims = 1:15)



ElbowPlot(pbmc)







#clustering을 할 때 논문 공부하다가 아게 된건데

#보톤 single cell rna seq은 정보가 없는 상황에서 clustering을 해야하니까 군집내 최소거리법 을 사용해  



pbmc <- FindNeighbors(pbmc, dims = 1:10)

pbmc <- FindClusters(pbmc, resolution = 0.5)



# Look at cluster IDs of the first 5 cells

head(Idents(pbmc), 5)





#비선형 차원 축소

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =

# 'umap-learn')

pbmc <- RunUMAP(pbmc, dims = 1:10)



#이것도 자세한 건 잘 모르겠따 







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







#몇몇 데이터 세트의 경우 표준 마커를 사용하여 편향되지 않는 클러스터링을 적용하여 셀 유형에 쉽게 일치 가능

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     
                     "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(pbmc) # 오른 쪽 이름을 왼쪽으로 변경 

#level은 factor들을 반환하는것 

pbmc <- RenameIdents(pbmc, new.cluster.ids)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



saveRDS(pbmc, file = "../output/pbmc3k_final.rds")





































































































맨위로

주고받은 메일 2 도움말전체보기
